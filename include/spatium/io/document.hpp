#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/error.hpp>
#  include <spatium/io/field.hpp>
#  include <spatium/io/field_mathjson.hpp>
#  include <spatium/io/field_parse.hpp>
#  include <spatium/io/field_text.hpp>
#  include <spatium/io/json.hpp>
#  include <algorithm>
#  include <cctype>
#  include <format>
#  include <optional>
#  include <string>
#  include <string_view>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::io::build {

// ── A document: named definitions, one a line ────────────────────
//
// What a person writes to say what is wanted, before anything is computed: numbers and
// expressions with names, parameters that range over an interval, roots of equations,
// integrals and shapes, each a line, each reading only the names above it.
//
//   r = 0.5                              a value: a number or an expression in the names above
//   p in [0.5, 4]                        a parameter: any number in the interval
//   x: x^2 - p = 0 on [0, 10]            a root: the x on [0, 10] for which the equation holds
//   y: 2y^2 - 3y - 5 = 0                 a root with no bracket: every y it holds for, in closed form
//   I = integral exp(-p * t^2) dt on [0, 1]  an integral over t, named by its dt
//   S0 = torus(2, r)                     a shape: an operation on numbers and the shapes above
//   S1 = offset(S0, h)
//   # a comment, to the end of the line
//
// A parameter is chosen (`in`); a root and an integral are found on an interval (`on`), and a
// root written without one is solved in closed form or not at all. A value or a root that reads a
// parameter is a function of it: the calculator answers such a definition for the whole interval,
// a scene draws it at the parameter's current value.
//
// The document is a description and computes nothing: what each definition is worth is the
// next layer's (spatium/answer/evaluate.hpp gives a definition its value, by the route it takes;
// io/document_trace.hpp builds its shapes as a scene). Each field keeps the text it was
// read from, so a document is printed back as it was written, and its JSON carries the same
// texts, read through the same reader.

struct Written {
    Field<double> field;
    std::string text;   // as written; for a field built in code, its `to_text`
};

enum class StatementKind { Value, Parameter, Root, Integral, Shape };

inline constexpr StatementKind statement_kinds[] = {StatementKind::Value, StatementKind::Parameter, StatementKind::Root,
                                                    StatementKind::Integral, StatementKind::Shape};

inline const char* statement_kind_name(StatementKind k) {
    switch (k) {
        case StatementKind::Value: return "value";
        case StatementKind::Parameter: return "parameter";
        case StatementKind::Root: return "root";
        case StatementKind::Integral: return "integral";
        case StatementKind::Shape: return "shape";
    }
    return "?";
}

// The operations a shape is made by, and what each argument is: `n` a number, `f` a field (both
// expressions in the names above), `s` the name of a shape above. Arguments past `required` may
// be left out; a `variadic` operation takes its last kind any number of times.
struct ShapeOperation {
    std::string_view name, arguments;
    std::size_t required;
    bool variadic;
    std::string_view caption;
};

inline constexpr ShapeOperation shape_operations[] = {
    {"torus", "nn", 2, false, "a torus: the radius to the tube's centre, and the tube's"},
    {"cylinder", "nn", 2, false, "a closed cylinder on the z axis: its radius and its height"},
    {"sphere", "n", 1, false, "a sphere at the origin: its radius"},
    {"flake", "nnnn", 3, false, "a shallow cap in a slab: the slab's half extents, and how much wider its sphere is"},
    {"offset", "sf", 2, false, "a closed surface pushed along its normal by a thickness"},
    {"scatter", "ssnn", 3, false, "copies of a shape on a surface: the shape, the surface, how many, and a seed"},
    {"union", "s", 1, true, "the shapes together"},
};

// The relations a root's sides may stand in, as written.
inline constexpr std::string_view relations[] = {"=", "<", "<=", ">", ">="};

namespace detail {

inline const ShapeOperation* shape_operation(std::string_view name) {
    for (const auto& op : shape_operations)
        if (op.name == name) return &op;
    return nullptr;
}

// what argument k of an operation is
inline char argument_kind(const ShapeOperation& op, std::size_t k) {
    return k < op.arguments.size() ? op.arguments[k] : op.arguments.back();
}

} // namespace detail

struct Statement {
    StatementKind kind = StatementKind::Value;
    std::string name;
    Written value;            // a value's expression, an integral's integrand, a root's left side
    Written rhs;              // a root's right side
    std::string variable;     // an integral's variable (a root's unknown is its own name)
    Written lo, hi;           // the interval of a parameter, a root or an integral; a root's may be unwritten
    std::string operation;    // a shape's operation, one of `shape_operations`
    std::vector<Written> args;   // and its arguments, a shape's by its name
    // a root's relation: "=" for an equation, or an inequality's "<", "<=", ">", ">=", whose
    // solution is the set where it holds (`x: x^2 - 4 > 0`)
    std::string relation = "=";
    // a system: the other unknowns of a root, each with its equation, solved together with it
    // (`x, y: x + 2y = 5, 3x - y = 1`)
    std::vector<Statement> also;
};

struct Document {
    std::vector<Statement> statements;

    // the definition of a name: its own, or the system it is an unknown of
    const Statement* find(std::string_view name) const {
        for (const auto& d : statements)
            if (d.name == name || std::ranges::any_of(d.also, [&](const Statement& e) { return e.name == name; })) return &d;
        return nullptr;
    }
};

// ── The forms a statement is written in ──────────────────────────
//
// One row a form: the line as typed, the row of a formula editor, and the tail a root on a
// bracket adds to each. A form is a template: `<slot>` a part of the statement, the rest words
// and marks as written, spaces as the line is printed. Reading a line, printing it back, its LaTeX
// and its JSON all walk these rows, so no two of them can say a statement differently. A form
// with no `<name>` is a shorthand, what a calculator's line says (`read_shorthand`): the name
// is found from what it reads. Where two forms start alike, the reader takes the first that reads
// a word or a mark the others have not (`in`, `integral`, `:`, a system's `,`, a shape's operation).
//
// A row of the editor comes back as MathJSON, which `mathjson` reads (`mathjson_row_text`): the
// trees a form is written as, a JSON array of patterns -- `"<slot>"` a part of the statement, a
// string at the head of a node the head it is (`"Limits|Tuple"` either), a node exactly that many
// items. A row is the line of the first form one of whose patterns it fits.
struct StatementForm {
    StatementKind kind;
    std::string_view line, latex = {};          // latex empty: a form no row of an editor writes
    std::string_view on = {}, on_latex = {};    // a root's bracket, where it has one
    std::string_view mathjson = {};
};

inline constexpr StatementForm statement_forms[] = {
    {.kind = StatementKind::Integral,
     .line = "<name> = integral <integrand> d<variable> on [<lo>, <hi>]",
     .latex = "<named>\\int_{<lo>}^{<hi>} <integrand>\\,d<variable>",
     .mathjson = R"([["Equal", "<name>", ["Integrate", "<integrand>", ["Limits|Tuple", "<variable>", "<lo>", "<hi>"]]],
                     ["Integrate", "<integrand>", ["Limits|Tuple", "<variable>", "<lo>", "<hi>"]]])"},
    {.kind = StatementKind::Shape, .line = "<name> = <operation>(<args>)"},
    {.kind = StatementKind::Value, .line = "<name> = <value>", .latex = "<named><value>"},
    {.kind = StatementKind::Root,
     .line = "<name>: <lhs> <relation> <rhs>",
     .latex = "<lhs> <relation> <rhs>",
     .on = " on [<lo>, <hi>]",
     .on_latex = ", <name> \\in \\left[<lo>, <hi>\\right]",
     .mathjson = R"([["Sequence", ["<relation>", "<lhs>", "<rhs>"], ["Element", "<name>", ["Interval", "<lo>", "<hi>"]]]])"},
    {.kind = StatementKind::Parameter,
     .line = "<name> in [<lo>, <hi>]",
     .latex = "<name> \\in \\left[<lo>, <hi>\\right]",
     .mathjson = R"([["Element", "<name>", ["Interval", "<lo>", "<hi>"]]])"},
    {.kind = StatementKind::Root, .line = "<names>: <equations>"},   // a system, an equation a row (`to_latex_rows`)
    // shorthands: an equation or inequality, on a bracket or not; an expression, `2 + 2` or `2 + 2 =`
    {.kind = StatementKind::Root, .line = "<lhs> <relation> <rhs>", .on = " on [<lo>, <hi>]",
     .mathjson = R"([["<relation>", "<lhs>", "<rhs>"]])"},
    {.kind = StatementKind::Value, .line = "<value>", .mathjson = R"(["<value>"])"},
    {.kind = StatementKind::Value, .line = "<value> ="},
};

// What a node of a row is, said where a row that holds to a pattern departs from it.
inline constexpr std::pair<std::string_view, std::string_view> mathjson_nodes[] = {
    {"Element", "a name and its interval are written x \\in [a, b]"},
    {"Interval", "an interval is written [a, b]"},
    {"Integrate", "an integral is written with its ends and its d, \\int_a^b f\\,dt"},
    {"Limits|Tuple", "an integral here has both ends, \\int_a^b"},
    {"Sequence", "a row holds one statement, or an equation or inequality and its unknown's interval: x^2 = 2, x \\in [0, 2]"},
};

namespace detail {

// the slots that hold an expression, and the field of the statement each is
struct WrittenSlot {
    std::string_view slot;
    Written Statement::* field;
};
inline constexpr WrittenSlot written_slots[] = {{"value", &Statement::value}, {"lhs", &Statement::value}, {"integrand", &Statement::value},
                                                {"rhs", &Statement::rhs},     {"lo", &Statement::lo},     {"hi", &Statement::hi}};

inline Written Statement::* written_slot(std::string_view slot) {
    for (const auto& w : written_slots)
        if (w.slot == slot) return w.field;
    return nullptr;
}

inline bool is_named(const StatementForm& f) { return f.line.starts_with("<name"); }
inline bool is_system(const StatementForm& f) { return f.line.starts_with("<names>"); }

// the form a statement is written in
inline const StatementForm& form_of(StatementKind kind, bool system) {
    for (const auto& f : statement_forms)
        if (f.kind == kind && is_named(f) && is_system(f) == system) return f;
    return statement_forms[0];   // every kind has a named form; held by the round-trip tests
}

// The parts of a form in order: a slot (`<lo>`), a word (`in`), a mark (`[`); spaces only part them.
struct FormPart {
    std::string_view text;
    bool slot;
};

inline std::vector<FormPart> form_parts(std::string_view form) {
    std::vector<FormPart> out;
    for (std::size_t i = 0; i < form.size();) {
        if (form[i] == ' ') {
            ++i;
        } else if (form[i] == '<') {
            const std::size_t end = form.find('>', i);
            out.push_back({form.substr(i + 1, end - i - 1), true});
            i = end + 1;
        } else {
            std::size_t end = i + 1;
            while (std::isalpha(static_cast<unsigned char>(form[i])) && end < form.size() && std::isalpha(static_cast<unsigned char>(form[end]))) ++end;
            out.push_back({form.substr(i, end - i), false});
            i = end;
        }
    }
    return out;
}

// A form written out, each slot replaced by what `part` gives for it.
template <class Part>
inline Result<std::string> fill(std::string_view form, Part part) {
    std::string out;
    for (std::size_t i = 0; i < form.size();) {
        const std::size_t open = form.find('<', i);
        out += form.substr(i, open == std::string_view::npos ? std::string_view::npos : open - i);
        if (open == std::string_view::npos) break;
        const std::size_t close = form.find('>', open);
        Result<std::string> p = part(form.substr(open + 1, close - open - 1));
        if (!p) return std::unexpected(p.error());
        out += *p;
        i = close + 1;
    }
    return out;
}

inline std::string_view trim(std::string_view s) {
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) s.remove_prefix(1);
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) s.remove_suffix(1);
    return s;
}

inline bool reserved_name(std::string_view name) {
    if (name == "in" || name == "on" || name == "integral") return true;
    if (detail::shape_operation(name)) return true;
    if (std::ranges::any_of(functions, [&](const Function& f) { return f.name == name; })) return true;
    return std::ranges::any_of(constants, [&](const Constant& c) { return c.name == name; });
}

// an error in a line; line 0 is a document that came from no text, whose error names no line
inline std::unexpected<Error> doc_error(std::size_t line, std::string message) {
    return std::unexpected(Error(ErrorCode::ParseError, line ? std::format("line {}: {}", line, message) : message));
}

// an expression read where the reader stands, with the text it was read from
inline Result<Written> written(FieldReader& r, std::string_view text) {
    r.word();   // past the spaces, so the text starts where the expression does
    const std::size_t from = r.position();
    auto f = r.expression();
    if (!f) return std::unexpected(f.error());
    return Written{std::move(*f), std::string(trim(text.substr(from, r.position() - from)))};
}

// The relation between a root's sides: `=`, or `<`, `<=`, `>`, `>=`; none when the next is none.
inline std::string relation(FieldReader& r) {
    if (r.take('=')) return "=";
    for (const char c : {'<', '>'})
        if (r.take(c)) return r.take('=') ? std::string{c, '='} : std::string{c};
    return {};
}

// A relation as LaTeX: \le and \ge for the two that hold at their ends.
inline std::string relation_latex(std::string_view rel) {
    return rel == "<=" ? "\\le" : rel == ">=" ? "\\ge" : std::string(rel);
}

// the variable an integral names by its d: the word d<name> standing before `on [`
inline std::string integral_variable(std::string_view text) {
    for (std::size_t at = text.find(" d"); at != std::string_view::npos; at = text.find(" d", at + 1)) {
        std::size_t end = at + 2;
        while (end < text.size() && (std::isalnum(static_cast<unsigned char>(text[end])) || text[end] == '_')) ++end;
        if (end == at + 2) continue;
        if (trim(text.substr(end)).starts_with("on")) return std::string(text.substr(at + 2, end - at - 2));
    }
    return {};
}

} // namespace detail

// Whether a root is sought on a bracket; one written without is solved in closed form.
inline bool has_bracket(const Statement& d) { return !d.lo.text.empty(); }

// The names a definition reads that it does not bind itself.
inline std::vector<std::string> reads(const Statement& d) {
    std::vector<std::string> out;
    const auto add = [&](const Written& w, std::string_view bound = {}) {
        for (auto& n : w.field.free_symbols())
            if (n != bound && std::ranges::find(out, n) == out.end()) out.push_back(std::move(n));
    };
    // every expression it holds, an integral's integrand less its variable, a root's equations less
    // the unknowns they solve for; then its bracket, and a shape's arguments
    add(d.value, d.variable), add(d.rhs);
    for (const auto& e : d.also) add(e.value), add(e.rhs);
    if (d.kind == StatementKind::Root)
        std::erase_if(out, [&](const std::string& n) { return n == d.name || std::ranges::any_of(d.also, [&](const Statement& e) { return e.name == n; }); });
    add(d.lo), add(d.hi);
    for (const auto& a : d.args) add(a);
    return out;
}

// That a document reads well: every name defined once and not a word of the language, every
// definition reading only the names above it, a parameter's ends numbers, a root reading its
// unknown, a shape's arguments as many and of the kinds its operation takes, and a shape read
// only where a shape is taken. `line_of` says which line each definition came from, for the
// message.
inline Result<void> check(const Document& doc, const std::vector<std::size_t>& line_of = {}) {
    const auto line = [&](std::size_t i) { return i < line_of.size() ? line_of[i] : std::size_t{0}; };
    for (std::size_t i = 0; i < doc.statements.size(); ++i) {
        const Statement& d = doc.statements[i];
        const auto defined_above = [&](std::string_view n) {
            return std::ranges::any_of(doc.statements.begin(), doc.statements.begin() + static_cast<std::ptrdiff_t>(i), [&](const Statement& e) {
                return e.name == n || std::ranges::any_of(e.also, [&](const Statement& a) { return a.name == n; });
            });
        };
        for (std::size_t k = 0; k <= d.also.size(); ++k) {
            const std::string& n = k ? d.also[k - 1].name : d.name;
            if (detail::reserved_name(n)) return detail::doc_error(line(i), std::format("'{}' is a word of the language, not a name", n));
            if (defined_above(n) || (k && (n == d.name || std::ranges::any_of(d.also.begin(), d.also.begin() + static_cast<std::ptrdiff_t>(k - 1),
                                                                                [&](const Statement& a) { return a.name == n; }))))
                return detail::doc_error(line(i), std::format("'{}' is defined twice", n));
        }
        if (d.kind == StatementKind::Integral) {
            if (d.variable.empty()) return detail::doc_error(line(i), "an integral names its variable by d, as in exp(-t^2) dt");
            if (d.variable == d.name || defined_above(d.variable))
                return detail::doc_error(line(i), std::format("the variable '{}' of the integral is a name already", d.variable));
        }
        if (d.kind == StatementKind::Parameter)
            for (const Written* w : {&d.lo, &d.hi})
                if (const auto open = w->field.free_symbols(); !open.empty())
                    return detail::doc_error(line(i), std::format("an end of the parameter {} reads '{}'; it is a number", d.name, open.front()));
        if (d.kind == StatementKind::Root && d.also.empty()) {
            const auto both = (d.value.field - d.rhs.field).free_symbols();
            if (std::ranges::find(both, d.name) == both.end())
                return detail::doc_error(line(i), std::format("the equation of {} does not read {}", d.name, d.name));
        }
        // a system's unknowns are each read by one of its equations at least
        for (std::size_t k = 0; k <= d.also.size() && !d.also.empty(); ++k) {
            const std::string& n = k ? d.also[k - 1].name : d.name;
            const auto reads_it = [&](const Statement& e) {
                const auto both = (e.value.field - e.rhs.field).free_symbols();
                return std::ranges::find(both, n) != both.end();
            };
            if (!reads_it(d) && std::ranges::none_of(d.also, reads_it))
                return detail::doc_error(line(i), std::format("no equation of the system reads {}", n));
        }
        for (const auto& n : reads(d))
            if (!defined_above(n))
                return detail::doc_error(line(i), std::format("{} reads '{}', which no line above defines", d.name, n));
        const auto is_shape = [&](std::string_view n) {
            const Statement* s = doc.find(n);
            return s && s->kind == StatementKind::Shape;
        };
        // a shape is taken by its name alone, where an operation takes one, and nowhere else
        const auto numbers_only = [&](const Field<double>& f, std::string_view what) -> Result<void> {
            for (const auto& n : f.free_symbols())
                if (is_shape(n)) return detail::doc_error(line(i), std::format("{} reads the shape {} as a number", what, n));
            return {};
        };
        if (d.kind != StatementKind::Shape) {
            for (const Written* w : {&d.value, &d.rhs, &d.lo, &d.hi})
                if (auto ok = numbers_only(w->field, d.name); !ok) return ok;
            continue;
        }
        const ShapeOperation* op = detail::shape_operation(d.operation);
        if (!op) return detail::doc_error(line(i), std::format("'{}' is not an operation a shape is made by", d.operation));
        if (d.args.size() < op->required || (!op->variadic && d.args.size() > op->arguments.size()))
            return detail::doc_error(line(i), std::format("{} takes {} arguments, given {}", op->name,
                                                          op->required == op->arguments.size() || op->variadic
                                                              ? std::format("{}{}", op->required, op->variadic ? " or more" : "")
                                                              : std::format("{} to {}", op->required, op->arguments.size()),
                                                          d.args.size()));
        for (std::size_t k = 0; k < d.args.size(); ++k) {
            const Written& a = d.args[k];
            const std::string what = std::format("argument {} of {}", k + 1, op->name);
            if (detail::argument_kind(*op, k) == 's') {
                if (!is_shape(a.text)) return detail::doc_error(line(i), std::format("{} is a shape, given '{}'", what, a.text));
            } else if (auto ok = numbers_only(a.field, what); !ok) {
                return ok;
            }
        }
    }
    return {};
}

namespace detail {

// One line read by one form: the statement, or how far the form went (`reader`, standing where it
// stopped, and the mark it `expected` there when that was one of the form's own words or marks),
// and whether it went past a part the forms that start alike do not share (`committed`).
struct FormRead {
    Result<Statement> statement;
    FieldReader reader;
    std::string expected;
    bool committed = false;
};

inline FormRead read_form(const StatementForm& form, std::string_view text) {
    const std::string variable = integral_variable(text);
    std::vector<std::string> reserved{"in", "on", "integral"};
    if (!variable.empty()) reserved.push_back("d" + variable);
    FormRead out{Statement{}, FieldReader(text, "", "", std::move(reserved)), {}, false};
    FieldReader& r = out.reader;
    Statement& d = *out.statement;
    d.kind = form.kind;
    std::optional<Error> error;
    const auto stop = [&](std::unexpected<Error> e, std::string expected = {}) {
        error = std::move(e.error()), out.expected = std::move(expected);
        return false;
    };
    const auto name = [&](std::string& into) {
        into = std::string(r.word());
        if (into.empty() || std::isdigit(static_cast<unsigned char>(into.front()))) return stop(r.fail("a name"));
        return r.take_word(into), true;
    };
    const auto expression = [&](Written& into) {
        auto w = written(r, text);
        if (!w) return stop(std::unexpected(w.error()));
        return into = std::move(*w), true;
    };
    const auto slot = [&](std::string_view p) {
        if (Written Statement::* w = written_slot(p)) return expression(d.*w);
        if (p == "name") return name(d.name);
        if (p == "relation") {
            d.relation = relation(r);
            return !d.relation.empty() || stop(r.fail("'=', '<' or '>'"));
        }
        if (p == "operation") {
            const ShapeOperation* op = shape_operation(r.word());
            if (!op) return stop(r.fail("an operation a shape is made by"));
            r.take_word(op->name), d.operation = std::string(op->name), out.committed = true;
            return true;
        }
        if (p == "args") {
            if (r.sees(')')) return true;
            for (Written a;; d.args.push_back(std::move(a))) {
                if (!expression(a)) return false;
                if (r.sees(')')) return d.args.push_back(std::move(a)), true;
                if (!r.take(',')) return stop(r.fail("',' or ')'"));
            }
        }
        if (p == "names") {   // a system's unknowns, two or more
            if (!name(d.name)) return false;
            while (r.take(',')) {
                out.committed = true;
                Statement& e = d.also.emplace_back();
                e.kind = StatementKind::Root;
                if (!name(e.name)) return false;
            }
            return out.committed || stop(r.fail("','"), "','");
        }
        if (p == "equations") {   // as many as the unknowns, each an equation
            for (std::size_t k = 0; k <= d.also.size(); ++k) {
                Statement& e = k ? d.also[k - 1] : d;
                if (k && !r.take(','))
                    return stop(r.fail(std::format("',' and the equation of {}: a system has an equation for each unknown", e.name)));
                if (!expression(e.value)) return false;
                if (!r.take('=')) return stop(r.fail("'='"));
                if (!expression(e.rhs)) return false;
            }
            return r.done() || stop(r.fail("the end of the line: a system has as many equations as unknowns, and no bracket"));
        }
        return stop(r.fail(std::format("<{}>", p)));   // a slot no reader knows: a fault in the table
    };
    const auto walk = [&](std::string_view part_of_form) {
        const auto parts = form_parts(part_of_form);
        std::string_view word_before;   // the word just read, which a missing mark is said to follow
        for (std::size_t k = 0; k < parts.size(); ++k) {
            const auto [p, is_slot] = parts[k];
            if (is_slot) {
                if (!slot(p)) return false;
                word_before = p == "operation" ? std::string_view(d.operation) : std::string_view{};
            } else if (p == "d" && k + 1 < parts.size() && parts[k + 1].text == "variable") {
                // an integral's d and its variable, one word
                if (variable.empty() || !r.take_word("d" + variable)) return stop(r.fail("the variable, as in dt"));
                d.variable = variable, ++k;
            } else {
                const bool word = std::isalpha(static_cast<unsigned char>(p.front()));
                if (word ? !r.take_word(p) : !r.take(p.front())) {
                    std::string want = std::format("'{}'{}", p, word_before.empty() ? "" : std::format(" after {}", word_before));
                    return stop(r.fail(want), want);
                }
                out.committed = out.committed || p != "=";
                word_before = word ? p : std::string_view{};
            }
        }
        return true;
    };
    if (walk(form.line) && (form.on.empty() || r.done() || walk(form.on)) && (r.done() || stop(r.fail("the end of the line"))))
        return out;
    out.statement = std::unexpected(std::move(*error));
    return out;
}

// One line, which is not empty and not a comment, read by the forms with a name (`named`) or the
// shorthands: the first form that reads it whole, or the error of the first that committed to it;
// where none did, the error of the one that read furthest, its marks and those of the forms that
// stopped at the same place named together ("expected '=', ':' or 'in'").
inline Result<Statement> read_statement(std::string_view text, bool named = true) {
    std::optional<FormRead> furthest;
    std::vector<std::string> expected;
    for (const auto& form : statement_forms) {
        if (is_named(form) != named) continue;
        FormRead f = read_form(form, text);
        if (f.statement || f.committed) return std::move(f.statement);
        const std::size_t at = f.reader.position();
        if (furthest && at < furthest->reader.position()) continue;
        if (!furthest || at > furthest->reader.position()) furthest = f, expected.clear();
        if (!f.expected.empty() && std::ranges::find(expected, f.expected) == expected.end()) expected.push_back(f.expected);
    }
    if (expected.size() < 2) return std::move(furthest->statement);
    std::string marks = expected.front();
    for (std::size_t k = 1; k < expected.size(); ++k) marks += (k + 1 < expected.size() ? ", " : " or ") + expected[k];
    return furthest->reader.fail(marks);
}

} // namespace detail

// A document from its text, a definition a line; blank lines and `#` comments are skipped. An
// error names the line, and a parse error its column in that line.
inline Result<Document> read_document(std::string_view text) {
    Document doc;
    std::vector<std::size_t> line_of;
    std::size_t number = 0;
    while (!text.empty()) {
        const std::size_t end = text.find('\n');
        std::string_view line = text.substr(0, end);
        text = end == std::string_view::npos ? std::string_view{} : text.substr(end + 1);
        ++number;
        if (const auto hash = line.find('#'); hash != std::string_view::npos) line = line.substr(0, hash);
        if (detail::trim(line).empty()) continue;
        auto d = detail::read_statement(line);
        if (!d) return detail::doc_error(number, d.error().message);
        doc.statements.push_back(std::move(*d));
        line_of.push_back(number);
    }
    if (auto ok = check(doc, line_of); !ok) return std::unexpected(ok.error());
    return doc;
}

namespace detail {

// a statement's form, and whether it is written with its bracket
inline const StatementForm& form_of(const Statement& d) { return form_of(d.kind, !d.also.empty()); }
inline bool written_on(const Statement& d) { return !d.lo.text.empty() || !d.hi.text.empty(); }

// what a slot of a statement's form holds, as the line writes it
inline Result<std::string> text_part(const Statement& d, std::string_view slot) {
    if (Written Statement::* w = written_slot(slot)) return (d.*w).text;
    std::string out;
    if (slot == "name") out = d.name;
    else if (slot == "relation") out = d.relation;
    else if (slot == "variable") out = d.variable;
    else if (slot == "operation") out = d.operation;
    else if (slot == "args")
        for (const auto& a : d.args) out += (out.empty() ? "" : ", ") + a.text;
    else if (slot == "names")
        for (std::size_t k = 0; k <= d.also.size(); ++k) out += (k ? ", " + d.also[k - 1].name : d.name);
    else if (slot == "equations")
        for (std::size_t k = 0; k <= d.also.size(); ++k) {
            const Statement& e = k ? d.also[k - 1] : d;
            out += (k ? ", " : "") + e.value.text + " = " + e.rhs.text;
        }
    return out;
}

// and as a row of a formula editor writes it; a row the editor wrote with no name is read as
// `answer`, and written back with none (`named`)
inline Result<std::string> latex_part(const Statement& d, std::string_view slot) {
    if (Written Statement::* w = written_slot(slot)) return to_latex((d.*w).field);
    if (slot == "name") return name_latex(d.name);
    if (slot == "named") return d.name == "answer" ? std::string() : name_latex(d.name) + " = ";
    if (slot == "relation") return relation_latex(d.relation);
    if (slot == "variable") return name_latex(d.variable);
    return std::unexpected(Error(ErrorCode::InvalidArgument, std::format("{}: no row of an editor writes its {}", d.name, slot)));
}

} // namespace detail

// A definition as the line it was written as.
inline std::string to_text(const Statement& d) {
    const StatementForm& form = detail::form_of(d);
    const auto part = [&](std::string_view slot) { return detail::text_part(d, slot); };
    return *detail::fill(form.line, part) + (detail::written_on(d) ? *detail::fill(form.on, part) : "");
}

// A definition as the row of a formula editor that reads back as it (`mathjson_row_text`): the
// question a page shows is the library's own, `p \in [0.5, 4]`, `x^{2} - p = 0, x \in [0, 10]`,
// `I = \int_{0}^{1} f\,dt`. A shape is a scene's, and has no row; a system is several rows
// (`to_latex_rows`), and refused here rather than written as its first equation.
inline Result<std::string> to_latex(const Statement& d) {
    if (!d.also.empty())
        return std::unexpected(Error(ErrorCode::InvalidArgument, std::format("{} is a system, an equation a row: to_latex_rows", d.name)));
    const StatementForm& form = detail::form_of(d);
    if (form.latex.empty())
        return std::unexpected(Error(ErrorCode::InvalidArgument, std::format("{} is a shape, which a scene shows, not a row", d.name)));
    const auto part = [&](std::string_view slot) { return detail::latex_part(d, slot); };
    auto row = detail::fill(form.latex, part);
    if (!row || !has_bracket(d)) return row;
    auto on = detail::fill(form.on_latex, part);
    if (!on) return std::unexpected(on.error());
    return *row + *on;
}

// A definition as the rows of a formula editor that read back as it: one row, or a system's
// equations one a row, as several rows of equations read (`respond`); a shape has none.
inline Result<std::vector<std::string>> to_latex_rows(const Statement& d) {
    if (d.kind == StatementKind::Shape) return std::vector<std::string>{};
    if (d.also.empty()) {
        auto row = to_latex(d);
        if (!row) return std::unexpected(row.error());
        return std::vector<std::string>{std::move(*row)};
    }
    std::vector<std::string> rows;
    for (std::size_t k = 0; k <= d.also.size(); ++k) {
        const Statement& e = k ? d.also[k - 1] : d;
        auto l = to_latex(e.value.field);
        if (!l) return std::unexpected(l.error());
        auto r = to_latex(e.rhs.field);
        if (!r) return std::unexpected(r.error());
        rows.push_back(*l + " = " + *r);
    }
    return rows;
}

inline std::string to_text(const Document& doc) {
    std::string out;
    for (const auto& d : doc.statements) out += to_text(d) + "\n";
    return out;
}

// The document as JSON: each definition its name, its kind, and a key for each slot of its form
// with the text written there (`relation` only where it is not '='; a system's other unknowns
// and their equations in `also`). Read back, the keys fill the same form and the line it writes
// is read by the same reader as a typed line, so the two forms cannot mean different things. A
// field may also be MathJSON, what a LaTeX editor on a page writes; it is read as the text of the
// language first (field_mathjson.hpp), so it means the same again.
inline JsonValue to_json(const Document& doc) {
    JsonValue::ArrayType defs;
    for (const auto& d : doc.statements) {
        JsonValue::ObjectType o{{"name", d.name}, {"kind", statement_kind_name(d.kind)}};
        const StatementForm& form = detail::form_of(d);
        for (const std::string_view part_of_form : {form.line, detail::written_on(d) ? form.on : std::string_view{}})
            for (const auto& [slot, is_slot] : detail::form_parts(part_of_form)) {
                if (!is_slot || slot == "name" || slot == "names" || (slot == "relation" && d.relation == "=")) continue;
                if (slot == "args") {
                    JsonValue::ArrayType args;
                    for (const auto& a : d.args) args.emplace_back(a.text);
                    o.emplace_back("args", std::move(args));
                } else if (slot == "equations") {
                    JsonValue::ArrayType also;
                    for (const auto& e : d.also) also.emplace_back(JsonValue::ObjectType{{"name", e.name}, {"lhs", e.value.text}, {"rhs", e.rhs.text}});
                    o.emplace_back("lhs", d.value.text), o.emplace_back("rhs", d.rhs.text), o.emplace_back("also", std::move(also));
                } else {
                    o.emplace_back(std::string(slot), *detail::text_part(d, slot));
                }
            }
        defs.push_back(std::move(o));
    }
    return JsonValue::ObjectType{{"statements", std::move(defs)}};
}

inline Result<Document> document_from_json(const JsonValue& json) {
    const JsonValue* defs = json.find("statements");
    if (!defs || !defs->is_array()) return std::unexpected(Error(ErrorCode::ParseError, "document: no 'statements' array"));
    std::string text;
    for (const auto& d : defs->as_array()) {
        const std::string kind = d.string_or("kind", "");
        Statement st;
        st.name = d.string_or("name", "");
        const auto known = std::ranges::find_if(statement_kinds, [&](StatementKind k) { return kind == statement_kind_name(k); });
        if (known == std::end(statement_kinds))
            return std::unexpected(Error(ErrorCode::ParseError, std::format("document: '{}' is not a kind of definition", kind)));
        st.kind = *known;
        // a field is its text, or MathJSON (field_mathjson.hpp) read as that text; a root's whole
        // equation may stand as its left side, ["Equal", a, b], with no right side
        std::optional<Error> unread;
        const auto read = [&](std::string_view key, const JsonValue& v) -> std::string {
            if (v.is_string()) return v.as_string();
            auto t = mathjson_text(v);
            if (!t && !unread) unread = Error(ErrorCode::ParseError, std::format("document: {}, {}: {}", st.name, key, t.error().message));
            return t ? *t : std::string();
        };
        std::vector<std::pair<std::string, std::string>> sides;
        if (const JsonValue* lhs = d.find("lhs"); st.kind == StatementKind::Root && lhs && !d.find("rhs") && lhs->is_array()) {
            auto eq = mathjson_equation(*lhs);
            if (!eq) return std::unexpected(Error(ErrorCode::ParseError, std::format("document: {}, lhs: {}", st.name, eq.error().message)));
            sides = {{"lhs", eq->first}, {"rhs", eq->second}};
        }
        const auto s = [&](std::string_view key) -> std::string {
            for (const auto& [k, v] : sides)
                if (k == key) return v;
            const JsonValue* v = d.find(key);
            return v ? read(key, *v) : std::string();
        };
        const JsonValue* also = d.find("also");
        const StatementForm& form = detail::form_of(st.kind, also && also->is_array());
        for (const std::string_view part_of_form : {form.line, d.find("lo") || d.find("hi") ? form.on : std::string_view{}})
            for (const auto& [slot, is_slot] : detail::form_parts(part_of_form)) {
                if (!is_slot || slot == "name" || slot == "names") continue;
                if (Written Statement::* w = detail::written_slot(slot)) st.*w = {{}, s(slot)};
                else if (slot == "relation") st.relation = d.string_or("relation", "=");
                else if (slot == "variable") st.variable = s(slot);
                else if (slot == "operation") st.operation = s(slot);
                else if (slot == "args") {
                    if (const JsonValue* a = d.find("args"); a && a->is_array())
                        for (const auto& x : a->as_array()) st.args.push_back({{}, x.is_string() ? x.as_string() : std::string("?")});
                } else if (slot == "equations") {
                    st.value = {{}, s("lhs")}, st.rhs = {{}, s("rhs")};
                    for (const auto& e : also->as_array())
                        st.also.push_back({.kind = StatementKind::Root, .name = e.string_or("name", "?"),
                                           .value = {{}, e.string_or("lhs", "?")}, .rhs = {{}, e.string_or("rhs", "?")}});
                }
            }
        if (unread) return std::unexpected(*unread);
        text += to_text(st) + "\n";
    }
    return read_document(text);
}

// A line written as no statement, read as the one it stands for, with the names `above`
// defined: an expression, `2 + 2` (or `2 + 2 =`), the value `answer`; an equation or inequality
// the root of it in the one name it reads that nothing above defines, on a bracket or not; and a
// new name alone on the left of '=' that the right does not read the value of that name, `f =
// a^2 + 1`. The statement alone, which the caller puts below `above` and checks.
inline Result<Statement> read_shorthand(std::string_view line, const Document& above) {
    auto read = detail::read_statement(line, false);
    if (!read) return read;
    Statement d = std::move(*read);
    const auto refuse = [](std::string why) { return std::unexpected(Error(ErrorCode::InvalidArgument, std::move(why))); };
    const auto open = [&](const Field<double>& f) {
        std::vector<std::string> out;
        for (auto& n : f.free_symbols())
            if (!above.find(n)) out.push_back(std::move(n));
        return out;
    };
    if (d.kind == StatementKind::Value) {
        if (const auto names = open(d.value.field); !names.empty())
            return refuse(std::format("the expression reads '{}', which has no value; to find it, write an equation ({} = 0)",
                                      names.front(), d.value.text));
        d.name = "answer";
        return d;
    }
    const std::string left(detail::trim(d.value.text));
    if (const auto right_reads = d.rhs.field.free_symbols(); d.relation == "=" && !has_bracket(d) && !above.find(left) &&
        d.value.field.free_symbols() == std::vector<std::string>{left} && !std::ranges::contains(right_reads, left))
        return Statement{.kind = StatementKind::Value, .name = left, .value = std::move(d.rhs)};
    const auto names = open(d.value.field - d.rhs.field);
    if (names.size() != 1)
        return refuse(std::format("{} is solved for one unknown; this one reads {}", d.relation == "=" ? "an equation" : "an inequality", names.size()));
    d.name = names.front();
    return d;
}

// Lines that are each an equation, two or more, reading as many unknowns as there are lines: the
// system they are, `x, y: x + 2y = 5, 3x - y = 1`, its unknowns in the order the lines read them;
// none where they are not, or where each line defines a new name by the names above it (`a = 2`
// above `f = a^2 + 1`), as each alone would (`read_shorthand`).
inline std::optional<Statement> read_system(std::string_view text) {
    std::vector<std::string> names;
    std::vector<Statement> equations;
    Document defined;   // the lines so far, while each is a definition
    while (!text.empty()) {
        const auto end = text.find('\n');
        const std::string_view line = detail::trim(text.substr(0, end));
        text = end == std::string_view::npos ? std::string_view{} : text.substr(end + 1);
        if (line.empty() || line.starts_with('#')) continue;
        auto e = detail::read_statement(line, false);
        if (!e || e->kind != StatementKind::Root || e->relation != "=" || detail::written_on(*e)) return std::nullopt;
        if (auto s = read_shorthand(line, defined); s && s->kind == StatementKind::Value && defined.statements.size() == equations.size()) {
            defined.statements.push_back(std::move(*s));
            if (!check(defined)) defined.statements.pop_back();
        }
        for (auto& n : (e->value.field - e->rhs.field).free_symbols())
            if (std::ranges::find(names, n) == names.end()) names.push_back(std::move(n));
        equations.push_back(std::move(*e));
    }
    if (equations.size() < 2 || names.size() != equations.size() || defined.statements.size() == equations.size()) return std::nullopt;
    Statement d = std::move(equations.front());
    d.name = names.front();
    for (std::size_t k = 1; k < equations.size(); ++k) d.also.push_back(std::move(equations[k])), d.also.back().name = names[k];
    if (!check(Document{{d}})) return std::nullopt;
    return d;
}

namespace detail {

// the relation a MathJSON head is, as the line writes it, or none
inline std::string_view mathjson_relation_of(const JsonValue& head) {
    static constexpr std::pair<std::string_view, std::string_view> heads[]{
        {"Equal", "="}, {"Less", "<"}, {"LessEqual", "<="}, {"Greater", ">"}, {"GreaterEqual", ">="}};
    for (const auto& [name, written] : heads)
        if (head.is_string() && head.as_string() == name) return written;
    return {};
}

// Where a row departs from a pattern it holds to: the innermost node with a head about the place.
struct Departure {
    std::string path;
    std::string head;
};

// Whether a row fits a pattern, its slots aside: each node's head and number of items. A row
// holds to the pattern (`held`) once it matched a head other than `Equal`, which an equation and
// an integral both start with.
inline bool mathjson_fits(const JsonValue& p, const JsonValue& v, const std::string& path, bool& held, std::optional<Departure>& off) {
    if (!p.is_array()) return true;   // a slot
    const auto &pa = p.as_array(), &va = v.as_array();   // empty when v is not an array
    const std::string& head = pa[0].as_string();
    const bool literal = !head.starts_with('<');
    const auto depart = [&] {
        if (!off && literal) off = Departure{path, head};
        return false;
    };
    // a head of the pattern is one name, or several "A|B", or a relation's slot
    const auto one_of = [&](std::string_view name) {
        for (std::size_t from = 0;; from = head.find('|', from) + 1) {
            if (std::string_view(head).substr(from, head.find('|', from) - from) == name) return true;
            if (head.find('|', from) == std::string::npos) return false;
        }
    };
    const bool head_fits = !va.empty() && (literal ? va[0].is_string() && one_of(va[0].as_string()) : !mathjson_relation_of(va[0]).empty());
    if (!head_fits || va.size() != pa.size()) return depart();
    held = held || (literal && head != "Equal");
    for (std::size_t k = 1; k < pa.size(); ++k)
        if (!mathjson_fits(pa[k], va[k], std::format("{}[{}]", path, k), held, off)) return depart();
    return true;
}

// The slots of a pattern a row fits, read into the statement: an expression as the text the
// language's printer writes for it, a name as a name, a relation by its head.
inline Result<void> mathjson_bind(const JsonValue& p, const JsonValue& v, const std::string& path, Statement& d) {
    if (p.is_string()) {
        const std::string_view slot = std::string_view(p.as_string()).substr(1, p.as_string().size() - 2);
        if (Written Statement::* w = written_slot(slot)) {
            auto t = mathjson_side(v, path);
            if (!t) return std::unexpected(t.error());
            (d.*w).text = std::move(*t);
            return {};
        }
        auto n = mathjson_name(v, path);
        if (!n) return std::unexpected(n.error());
        (slot == "variable" ? d.variable : d.name) = std::move(*n);
        return {};
    }
    const auto &pa = p.as_array(), &va = v.as_array();
    if (pa[0].as_string() == "<relation>") d.relation = std::string(mathjson_relation_of(va[0]));
    for (std::size_t k = 1; k < pa.size(); ++k)
        if (auto ok = mathjson_bind(pa[k], va[k], std::format("{}[{}]", path, k), d); !ok) return ok;
    return {};
}

} // namespace detail

// A row of a LaTeX editor, as MathJSON, read as the line of a document it is (or the question a
// page asks): the line of the first form one of whose patterns it fits (`StatementForm`), written
// as that form writes it -- `p \in [0.5, 4]` the parameter `p in [0.5, 4]`; an equation and its
// unknown's interval, `x^2 - p = 0, x \in [0, 10]`, the root `x: x^2 - p = 0 on [0, 10]`; `I =
// \int_0^1 f\,dt` the integral `I = integral f dt on [0, 1]` (with no name, `answer`); an equation
// or inequality as written; anything else an expression. Which of these a line is, and what it
// reads, the document's reader decides, so a row means what the same line typed means.
inline Result<std::string> mathjson_row_text(const JsonValue& row) {
    // the editor writes a list as ["Delimiter", ["Sequence", ...], ","]: the list is the row
    const auto& r = row.as_array();
    const bool listed = r.size() >= 2 && r[0].is_string() && r[0].as_string() == "Delimiter" && r[1].is_array() &&
                        !r[1].as_array().empty() && r[1].as_array()[0].is_string() && r[1].as_array()[0].as_string() == "Sequence";
    const JsonValue& v = listed ? r[1] : row;
    const std::string at = listed ? "[1]" : "";
    for (const auto& form : statement_forms) {
        if (form.mathjson.empty()) continue;
        const auto patterns = spatium::io::parse(form.mathjson);
        for (const auto& p : patterns->as_array()) {
            bool held = false;
            std::optional<detail::Departure> off;
            if (!detail::mathjson_fits(p, v, at, held, off)) {
                if (!held) continue;
                const auto said = std::ranges::find(mathjson_nodes, off ? off->head : "", &std::pair<std::string_view, std::string_view>::first);
                return detail::mathjson_error(off ? off->path : at, said != std::end(mathjson_nodes) ? std::string(said->second)
                                                                                                      : "a row holds one statement");
            }
            Statement d;
            d.kind = form.kind;
            if (auto ok = detail::mathjson_bind(p, v, at, d); !ok) return std::unexpected(ok.error());
            if (detail::is_named(form) && d.name.empty()) d.name = "answer";
            const auto part = [&](std::string_view slot) { return detail::text_part(d, slot); };
            return *detail::fill(form.line, part) + (detail::written_on(d) ? *detail::fill(form.on, part) : "");
        }
    }
    return detail::mathjson_error(at, "no form of a statement reads the row");   // the last form reads any expression
}

// Whether two documents say the same: the same names, kinds, variables and operations in the
// same order, and every field the same function.
inline bool same_document(const Document& a, const Document& b) {
    if (a.statements.size() != b.statements.size()) return false;
    for (std::size_t i = 0; i < a.statements.size(); ++i) {
        const Statement &x = a.statements[i], &y = b.statements[i];
        if (!same_document(Document{x.also}, Document{y.also})) return false;
        if (x.kind != y.kind || x.name != y.name || x.variable != y.variable || x.operation != y.operation || x.relation != y.relation) return false;
        const auto same = [](const Written& p, const Written& q) {
            return p.text.empty() == q.text.empty() && (p.text.empty() || same_value(p.field, q.field));
        };
        if (!same(x.value, y.value) || !same(x.rhs, y.rhs) || !same(x.lo, y.lo) || !same(x.hi, y.hi)) return false;
        if (x.args.size() != y.args.size()) return false;
        for (std::size_t k = 0; k < x.args.size(); ++k)
            if (!same(x.args[k], y.args[k])) return false;
    }
    return true;
}

} // namespace spatium::io::build
