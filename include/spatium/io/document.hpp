#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/error.hpp>
#  include <spatium/io/field.hpp>
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
// expressions with names, parameters that range over an interval, roots of equations and
// integrals, each a line, each reading only the names above it.
//
//   r = 0.5                              a value: a number or an expression in the names above
//   p in [0.5, 4]                        a parameter: any number in the interval
//   x: x^2 - p = 0 on [0, 10]            a root: the x on [0, 10] for which the equation holds
//   I = integral exp(-p * t^2) dt on [0, 1]  an integral over t, named by its dt
//   # a comment, to the end of the line
//
// A parameter is chosen (`in`); a root and an integral are found on an interval (`on`). A value
// or a root that reads a parameter is a function of it: the calculator answers such a definition
// for the whole interval, a scene draws it at the parameter's current value.
//
// The document is a description and computes nothing: what each definition is worth is the
// next layer's (rsc/spatium/query.hpp reads a root or an integral of one parameter into a
// question). Each field keeps the text it was read from, so a document is printed back as it
// was written, and its JSON carries the same texts, read through the same reader.

struct Written {
    Field<double> field;
    std::string text;   // as written; for a field built in code, its `to_text`
};

enum class StatementKind { Value, Parameter, Root, Integral };

inline const char* statement_kind_name(StatementKind k) {
    switch (k) {
        case StatementKind::Value: return "value";
        case StatementKind::Parameter: return "parameter";
        case StatementKind::Root: return "root";
        case StatementKind::Integral: return "integral";
    }
    return "?";
}

struct Statement {
    StatementKind kind = StatementKind::Value;
    std::string name;
    Written value;            // a value's expression, an integral's integrand, a root's left side
    Written rhs;              // a root's right side
    std::string variable;     // an integral's variable (a root's unknown is its own name)
    Written lo, hi;           // the interval of a parameter, a root or an integral
};

struct Document {
    std::vector<Statement> statements;

    const Statement* find(std::string_view name) const {
        for (const auto& d : statements)
            if (d.name == name) return &d;
        return nullptr;
    }
};

namespace detail {

inline std::string_view trim(std::string_view s) {
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) s.remove_prefix(1);
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) s.remove_suffix(1);
    return s;
}

inline bool reserved_name(std::string_view name) {
    if (name == "in" || name == "on" || name == "integral") return true;
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

// `on [a, b]` or `in [a, b]`
inline Result<std::pair<Written, Written>> interval(FieldReader& r, std::string_view text, std::string_view word) {
    if (!r.take_word(word)) return r.fail(std::format("'{}'", word));
    if (!r.take('[')) return r.fail("'['");
    auto a = written(r, text);
    if (!a) return std::unexpected(a.error());
    if (!r.take(',')) return r.fail("','");
    auto b = written(r, text);
    if (!b) return std::unexpected(b.error());
    if (!r.take(']')) return r.fail("']'");
    return std::pair{std::move(*a), std::move(*b)};
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

// The names a definition reads that it does not bind itself.
inline std::vector<std::string> reads(const Statement& d) {
    std::vector<std::string> out;
    const auto add = [&](const Field<double>& f, std::string_view bound) {
        for (auto& n : f.free_symbols())
            if (n != bound && std::ranges::find(out, n) == out.end()) out.push_back(std::move(n));
    };
    switch (d.kind) {
        case StatementKind::Value: add(d.value.field, {}); break;
        case StatementKind::Parameter: add(d.lo.field, {}), add(d.hi.field, {}); break;
        case StatementKind::Root: add(d.value.field, d.name), add(d.rhs.field, d.name), add(d.lo.field, {}), add(d.hi.field, {}); break;
        case StatementKind::Integral: add(d.value.field, d.variable), add(d.lo.field, {}), add(d.hi.field, {}); break;
    }
    return out;
}

// That a document reads well: every name defined once and not a word of the language, every
// definition reading only the names above it, a parameter's ends numbers, a root reading its
// unknown. `line_of` says which line each definition came from, for the message.
inline Result<void> check(const Document& doc, const std::vector<std::size_t>& line_of = {}) {
    const auto line = [&](std::size_t i) { return i < line_of.size() ? line_of[i] : std::size_t{0}; };
    for (std::size_t i = 0; i < doc.statements.size(); ++i) {
        const Statement& d = doc.statements[i];
        const auto defined_above = [&](std::string_view n) {
            return std::ranges::any_of(doc.statements.begin(), doc.statements.begin() + static_cast<std::ptrdiff_t>(i),
                                       [&](const Statement& e) { return e.name == n; });
        };
        if (detail::reserved_name(d.name)) return detail::doc_error(line(i), std::format("'{}' is a word of the language, not a name", d.name));
        if (defined_above(d.name)) return detail::doc_error(line(i), std::format("'{}' is defined twice", d.name));
        if (d.kind == StatementKind::Integral) {
            if (d.variable.empty()) return detail::doc_error(line(i), "an integral names its variable by d, as in exp(-t^2) dt");
            if (d.variable == d.name || defined_above(d.variable))
                return detail::doc_error(line(i), std::format("the variable '{}' of the integral is a name already", d.variable));
        }
        if (d.kind == StatementKind::Parameter)
            for (const Written* w : {&d.lo, &d.hi})
                if (const auto open = w->field.free_symbols(); !open.empty())
                    return detail::doc_error(line(i), std::format("an end of the parameter {} reads '{}'; it is a number", d.name, open.front()));
        if (d.kind == StatementKind::Root) {
            const auto both = (d.value.field - d.rhs.field).free_symbols();
            if (std::ranges::find(both, d.name) == both.end())
                return detail::doc_error(line(i), std::format("the equation of {} does not read {}", d.name, d.name));
        }
        for (const auto& n : reads(d))
            if (!defined_above(n))
                return detail::doc_error(line(i), std::format("{} reads '{}', which no line above defines", d.name, n));
    }
    return {};
}

namespace detail {

// One line, which is not empty and not a comment.
inline Result<Statement> read_statement(std::string_view text) {
    Statement d;
    const std::string variable = detail::integral_variable(text);
    std::vector<std::string> reserved{"in", "on", "integral"};
    if (!variable.empty()) reserved.push_back("d" + variable);
    FieldReader r(text, "", "", reserved);
    d.name = std::string(r.word());
    if (d.name.empty() || std::isdigit(static_cast<unsigned char>(d.name.front()))) return r.fail("a name");
    r.take_word(d.name);
    if (r.word() == "in") {
        d.kind = StatementKind::Parameter;
        auto ends = detail::interval(r, text, "in");
        if (!ends) return std::unexpected(ends.error());
        std::tie(d.lo, d.hi) = std::move(*ends);
    } else if (r.take(':')) {
        d.kind = StatementKind::Root;
        auto lhs = detail::written(r, text);
        if (!lhs) return std::unexpected(lhs.error());
        if (!r.take('=')) return r.fail("'='");
        auto rhs = detail::written(r, text);
        if (!rhs) return std::unexpected(rhs.error());
        d.value = std::move(*lhs), d.rhs = std::move(*rhs);
        auto ends = detail::interval(r, text, "on");
        if (!ends) return std::unexpected(ends.error());
        std::tie(d.lo, d.hi) = std::move(*ends);
    } else if (r.take('=')) {
        if (r.take_word("integral")) {
            d.kind = StatementKind::Integral;
            auto f = detail::written(r, text);
            if (!f) return std::unexpected(f.error());
            d.value = std::move(*f);
            if (variable.empty() || !r.take_word("d" + variable)) return r.fail("the variable, as in dt");
            d.variable = variable;
            auto ends = detail::interval(r, text, "on");
            if (!ends) return std::unexpected(ends.error());
            std::tie(d.lo, d.hi) = std::move(*ends);
        } else {
            d.kind = StatementKind::Value;
            auto f = detail::written(r, text);
            if (!f) return std::unexpected(f.error());
            d.value = std::move(*f);
        }
    } else {
        return r.fail("'=', ':' or 'in'");
    }
    if (!r.done()) return r.fail("the end of the line");
    return d;
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

// A definition as the line it was written as.
inline std::string to_text(const Statement& d) {
    switch (d.kind) {
        case StatementKind::Value: return std::format("{} = {}", d.name, d.value.text);
        case StatementKind::Parameter: return std::format("{} in [{}, {}]", d.name, d.lo.text, d.hi.text);
        case StatementKind::Root:
            return std::format("{}: {} = {} on [{}, {}]", d.name, d.value.text, d.rhs.text, d.lo.text, d.hi.text);
        case StatementKind::Integral:
            return std::format("{} = integral {} d{} on [{}, {}]", d.name, d.value.text, d.variable, d.lo.text, d.hi.text);
    }
    return {};
}

inline std::string to_text(const Document& doc) {
    std::string out;
    for (const auto& d : doc.statements) out += to_text(d) + "\n";
    return out;
}

// The document as JSON: each definition its name, its kind and the texts of its fields, which
// are read back through the same reader as a line is, so the two forms cannot mean different
// things.
inline JsonValue to_json(const Document& doc) {
    JsonValue::ArrayType defs;
    for (const auto& d : doc.statements) {
        JsonValue::ObjectType o{{"name", d.name}, {"kind", statement_kind_name(d.kind)}};
        switch (d.kind) {
            case StatementKind::Value: o.emplace_back("value", d.value.text); break;
            case StatementKind::Parameter: o.emplace_back("lo", d.lo.text), o.emplace_back("hi", d.hi.text); break;
            case StatementKind::Root:
                o.emplace_back("lhs", d.value.text), o.emplace_back("rhs", d.rhs.text);
                o.emplace_back("lo", d.lo.text), o.emplace_back("hi", d.hi.text);
                break;
            case StatementKind::Integral:
                o.emplace_back("integrand", d.value.text), o.emplace_back("variable", d.variable);
                o.emplace_back("lo", d.lo.text), o.emplace_back("hi", d.hi.text);
                break;
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
        const std::string kind = d.string_or("kind", ""), name = d.string_or("name", "");
        const auto s = [&](const char* key) { return d.string_or(key, ""); };
        if (kind == "value") text += std::format("{} = {}\n", name, s("value"));
        else if (kind == "parameter") text += std::format("{} in [{}, {}]\n", name, s("lo"), s("hi"));
        else if (kind == "root") text += std::format("{}: {} = {} on [{}, {}]\n", name, s("lhs"), s("rhs"), s("lo"), s("hi"));
        else if (kind == "integral")
            text += std::format("{} = integral {} d{} on [{}, {}]\n", name, s("integrand"), s("variable"), s("lo"), s("hi"));
        else return std::unexpected(Error(ErrorCode::ParseError, std::format("document: '{}' is not a kind of definition", kind)));
    }
    return read_document(text);
}

// Whether two documents say the same: the same names, kinds and variables in the same order,
// and every field the same function.
inline bool same_document(const Document& a, const Document& b) {
    if (a.statements.size() != b.statements.size()) return false;
    for (std::size_t i = 0; i < a.statements.size(); ++i) {
        const Statement &x = a.statements[i], &y = b.statements[i];
        if (x.kind != y.kind || x.name != y.name || x.variable != y.variable) return false;
        const auto same = [](const Written& p, const Written& q) {
            return p.text.empty() == q.text.empty() && (p.text.empty() || same_value(p.field, q.field));
        };
        if (!same(x.value, y.value) || !same(x.rhs, y.rhs) || !same(x.lo, y.lo) || !same(x.hi, y.hi)) return false;
    }
    return true;
}

} // namespace spatium::io::build
