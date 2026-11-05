#pragma once

// A question typed as text, and the answer as the JSON a page reads.
//
// A question is a document: the rows of a formula editor, each a line of it (`respond_rows`;
// `p \in [0.5, 4]` the parameter `p in [0.5, 4]`, `x^2 - p = 0, x \in [0, 10]` the root
// `x: x^2 - p = 0 on [0, 10]`), or its lines typed. One line that is not a statement is a
// shorthand `respond` writes out and evaluates (evaluate.hpp): `2 + 2` is the value
// `answer = 2 + 2`, `2x^2 - 3x - 5 = 0` the root `x: 2x^2 - 3x - 5 = 0`. A row reading a
// parameter above is answered as a function of it, synthesised with its certificate. A text that
// is none of these is refused with the reason, which is the library's own (a parse error says its
// column). The answer gives the question back as rows the library writes (`question_latex`),
// carries the document it was read as and every route that was tried, so it is a whole record of
// what was asked and how it was answered. The page composes nothing: what it shows of an answer is in this JSON, and the
// gallery is written by the same function, so a card and a live answer cannot differ.
//
// `evaluate_request` is the same door for a document a page or a builder writes: the document as
// text or JSON, the name it asks for and the parameters it chooses, in one JSON object.
//
// An error is written as the power of ten above it (`error_below`): a bound a person reads, and
// one that does not move with the last digits of a libm.

#include <rsc/spatium/closed_form.hpp>
#include <rsc/spatium/evaluate.hpp>
#include <rsc/spatium/query.hpp>
#include <rsc/spatium/synthesis.hpp>
#include <spatium/io/document.hpp>
#include <spatium/io/field_mathjson.hpp>
#include <spatium/io/plot_glsl.hpp>
#include <spatium/io/field_parse.hpp>
#include <spatium/io/field_text.hpp>
#include <spatium/io/json.hpp>

#include <algorithm>
#include <cctype>
#include <functional>
#include <optional>
#include <cmath>
#include <expected>
#include <format>
#include <string>
#include <string_view>

namespace rsc {

namespace respond_detail {

using spatium::io::JsonValue;
using Object = JsonValue::ObjectType;
using Array = JsonValue::ArrayType;

inline double decade_above(double e) { return e > 0 ? std::pow(10.0, std::ceil(std::log10(e))) : 0.0; }

inline std::string_view trimmed(std::string_view s) {
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) s.remove_prefix(1);
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) s.remove_suffix(1);
    return s;
}

} // namespace respond_detail

// An answer by synthesis as JSON: its pieces written out when there are at most `written` of them,
// else counted, and the certificate they share.
inline spatium::io::JsonValue::ObjectType answer_json(const Answer& a, const std::string& parameter, std::size_t written = 2) {
    using namespace respond_detail;
    namespace bd = spatium::io::build;
    Object out{{"answered", a.answered}};
    if (!a.refused.empty()) out.emplace_back("refused", a.refused);
    double worst = 0;
    for (const auto& s : a.pieces) worst = std::max(worst, s.error);
    out.emplace_back("pieces", a.pieces.size());
    if (!a.pieces.empty()) {
        const Certificate& c = a.pieces.front().certificate;
        out.emplace_back("certificate", Object{{"strength", c.strength},
                                               {"norm", c.relative ? "largest relative error" : "largest absolute error"},
                                               {"points", c.points},
                                               {"scalar", c.scalar},
                                               {"computed", c.computed},
                                               {"error_below", decade_above(worst)}});
    }
    if (!a.pieces.empty() && a.pieces.size() <= written) {
        Array pieces;
        for (const auto& s : a.pieces) {
            if (s.form == Form::Chebyshev) {
                // the polynomial is a long field; the page shows what it is and its coefficients
                const std::string p = bd::detail::name_latex(parameter);
                Array coefficients;
                for (double c : s.coefficients) coefficients.push_back(c);
                pieces.push_back(Object{{"lo", s.lo}, {"hi", s.hi}, {"form", "chebyshev"}, {"degree", s.terms - 1},
                                        {"latex", std::format("\\sum_{{k=0}}^{{{}}} c_k\\,T_k\\!\\left(\\frac{{2{} - {} - {}}}{{{}}}\\right)",
                                                              s.terms - 1, p, bd::detail::number_text(s.lo, bd::detail::TextStyle::Latex),
                                                              bd::detail::number_text(s.hi, bd::detail::TextStyle::Latex),
                                                              bd::detail::number_text(s.hi - s.lo, bd::detail::TextStyle::Latex))},
                                        {"coefficients", coefficients}, {"accepted", s.accepted}});
                continue;
            }
            pieces.push_back(Object{{"lo", s.lo}, {"hi", s.hi}, {"form", "expression"},
                                    {"latex", *bd::to_latex(s.expr, parameter)},
                                    {"text", *bd::to_text(s.expr, parameter)}, {"terms", s.terms},
                                    {"accepted", s.accepted}});
        }
        out.emplace_back("answer", pieces);
    }
    return out;
}

// An evaluation as JSON: the route its value came by and the value as that route gives it, or
// why there is none, and every candidate with its cost.
// A computed number on a page: a fraction as one, otherwise the arithmetic's own digits, and a
// double as the printer writes it.
inline std::string number_latex(const Candidate& c, double value) {
    namespace bd = spatium::io::build;
    if (c.text.empty()) return bd::detail::number_text(value, bd::detail::TextStyle::Latex);
    if (const auto slash = c.text.find('/'); slash != std::string::npos) {
        const bool negative = c.text.front() == '-';
        return std::format("{}\\frac{{{}}}{{{}}}", negative ? "-" : "", c.text.substr(negative, slash - negative),
                           c.text.substr(slash + 1));
    }
    return c.text;
}

// A number's proved digits: how far from zero only where none is proved, since elsewhere that is
// the value's last digits.
inline spatium::io::JsonValue::ObjectType proved_json(const spatium::Proved& p) {
    using namespace respond_detail;
    if (p.exact) return {{"digits", p.digits}, {"text", p.text}, {"exact", true}};
    if (p.digits == 0) return {{"digits", 0}, {"bound", decade_above(p.bound)}};
    return {{"digits", p.digits}, {"text", p.text}, {"exact", p.exact}};
}

namespace respond_detail {

// An inequality's axis as the page draws it: each part's kind, its ends, and a root's value.
inline Object axis_json(const Axis& a) {
    static constexpr const char* kinds[]{"holds", "fails", "root", "break", "undecided"};
    Array parts;
    for (const auto& p : a.parts) {
        Object o{{"kind", kinds[static_cast<int>(p.kind)]}, {"lo", p.lo}, {"hi", p.hi}};
        if (!p.text.empty()) o.emplace_back("text", p.text);
        parts.push_back(std::move(o));
    }
    return Object{{"relation", a.relation}, {"line", a.line}, {"parts", parts}};
}

// The set where an inequality holds, as LaTeX: its intervals joined, an end in a bracket where it
// belongs to the set -- a root of a relation that holds at equality, or a bracket's end -- and in a
// parenthesis where it does not; the empty set, or every real number.
inline std::string set_latex(const Axis& a, const std::string& name) {
    namespace bd = spatium::io::build;
    const auto number = [](double v) { return bd::detail::number_text(v, bd::detail::TextStyle::Latex); };
    const auto value = [](const std::string& t) {
        const auto slash = t.find('/');
        if (slash == std::string::npos) return t;
        const bool minus = t.starts_with('-');
        return std::string(minus ? "-" : "") + "\\frac{" + t.substr(minus, slash - minus) + "}{" + t.substr(slash + 1) + "}";
    };
    const bool at_equality = a.relation.size() == 2;
    using K = Axis::Kind;
    std::vector<std::string> spans;
    std::string lo;
    bool open = false, whole = true;
    for (std::size_t i = 0; i < a.parts.size(); ++i) {
        const auto& p = a.parts[i];
        const K prev = i ? a.parts[i - 1].kind : K::Fails, next = i + 1 < a.parts.size() ? a.parts[i + 1].kind : K::Fails;
        if (p.kind == K::Holds) {
            if (open) continue;
            open = true;
            if (i == 0) lo = a.line ? "(-\\infty" : "[" + number(p.lo);
            else if (prev == K::Root) lo = (at_equality ? "[" : "(") + value(a.parts[i - 1].text);
            else if (prev == K::Break && !a.parts[i - 1].text.empty()) lo = "(" + value(a.parts[i - 1].text);
            else lo = "(" + number(p.lo);
            continue;
        }
        whole = false;
        if (p.kind == K::Root && at_equality) {
            if (open && next == K::Holds) continue;   // a root inside the set, as x^2 >= 0 has at 0
            if (open) spans.push_back(lo + ", " + value(p.text) + "]"), open = false;
            else if (next != K::Holds) spans.push_back("\\{" + value(p.text) + "\\}");
            continue;
        }
        if (open) spans.push_back(lo + ", " + (!p.text.empty() ? value(p.text) : number(p.lo)) + ")"), open = false;
    }
    if (open) spans.push_back(lo + ", " + (a.line ? std::string("\\infty)") : number(a.parts.back().hi) + "]"));
    const std::string x = bd::detail::name_latex(name) + " \\in ";
    if (spans.empty()) return x + "\\varnothing";
    if (whole && a.line) return x + "\\mathbb{R}";
    std::string out;
    for (const auto& sp : spans) out += (out.empty() ? "" : " \\cup ") + sp;
    return x + out;
}

} // namespace respond_detail


inline spatium::io::JsonValue::ObjectType evaluation_json(const Evaluation& e, std::size_t written = 2) {
    using namespace respond_detail;
    namespace bd = spatium::io::build;
    Object out;
    if (const auto r = e.route()) out.emplace_back("route", route_name(*r));
    // a synthesis that found no value still says how far it got, its refusal among it
    if (!e.taken && e.answer)
        for (auto& [k, v] : answer_json(*e.answer, e.parameter, written)) out.emplace_back(std::move(k), std::move(v));
    else if (!e.refused.empty())
        out.emplace_back("refused", e.refused);
    if (e.taken) switch (e.candidates[*e.taken].route) {
        case Route::Evaluate: {
            const Candidate& c = e.candidates[*e.taken];
            out.emplace_back("value", e.values.front());
            if (!c.text.empty()) out.emplace_back("text", c.text);
            // an exact value as its fraction and, where it ends, its decimal: = \frac{391}{4} = 97.75
            std::string value = number_latex(c, e.values.front());
            if (e.exact && c.arithmetic == Arithmetic::Double) {
                const auto d = e.exact->decimal();
                value = e.exact->den() == 1 ? *d
                        : spatium::io::build::detail::fraction_latex(*e.exact) + (d ? " = " + *d : "");
            }
            out.emplace_back("latex", *bd::to_latex(e.written) + " = " + value);
            if (!e.steps.empty()) {
                Array steps;
                for (const auto& st : e.steps) steps.push_back(Object{{"label", st.label}, {"latex", st.latex}});
                out.emplace_back("steps", steps);
            }
            if (e.proved) out.emplace_back("proved", proved_json(*e.proved));
            else if (!e.unproved.empty())
                out.emplace_back("unproved", e.unproved);
            break;
        }
        case Route::Expression:
            out.emplace_back("reads", e.free.front());
            out.emplace_back("latex", *bd::to_latex(e.written));
            out.emplace_back("text", *bd::to_text(e.written));
            break;
        case Route::ClosedForm: {
            const auto& c = *e.closed;
            Array steps;
            for (const auto& s : c.steps) steps.push_back(Object{{"label", s.label}, {"latex", s.latex}});
            out.emplace_back("unknown", c.unknown);
            out.emplace_back("steps", steps);
            if (c.check)
                out.emplace_back("check", Object{{"by", "Vieta's formulas"},
                                                 {"holds", c.check->error <= 1e-12},
                                                 {"error_below", decade_above(c.check->error)}});
            // what the steps were computed through, as numbers
            Array roots;
            for (const auto& z : c.solution.roots) roots.push_back(Array{z.re, z.im});
            out.emplace_back("solution", Object{{"a", c.solution.a}, {"b", c.solution.b}, {"c", c.solution.c},
                                                {"discriminant", c.solution.discriminant}, {"roots", roots}});
            if (c.exact) out.emplace_back("exact", Object{{"text", c.exact->text}, {"latex", c.exact->latex}});
            Array values;
            for (double v : e.values) values.push_back(v);
            out.emplace_back("values", values);
            break;
        }
        case Route::Intervals: {
            out.emplace_back("axis", axis_json(*e.axis));
            out.emplace_back("set", set_latex(*e.axis, e.name));
            break;
        }
        case Route::Linear: {
            const auto& l = *e.linear;
            Array steps;
            for (const auto& s : l.steps) steps.push_back(Object{{"label", s.label}, {"latex", s.latex}});
            Object solution;
            for (std::size_t j = 0; j < l.unknowns.size(); ++j) solution.emplace_back(l.unknowns[j], l.values[j]);
            out.emplace_back("steps", steps);
            out.emplace_back("solution", solution);
            if (!l.exact.empty()) {
                Object exact;
                for (std::size_t j = 0; j < l.unknowns.size(); ++j) exact.emplace_back(l.unknowns[j], l.exact[j].text());
                out.emplace_back("exact", exact);
            }
            break;
        }
        case Route::Enclosed: {
            if (e.system) {
                // a system's roots, each as its unknowns' proved digits
                Array values;
                std::string latex;
                for (const auto& r : e.system->roots) {
                    std::string names, digits;
                    for (std::size_t j = 0; j < r.size(); ++j) {
                        values.push_back(r[j].mid());
                        names += (j ? ", " : "") + bd::detail::name_latex(e.system->unknowns[j]);
                        digits += (j ? ", " : "") + spatium::proved(r[j]).text;
                    }
                    latex += (latex.empty() ? "" : ";\\quad ") + std::format("({}) = ({})", names, digits);
                }
                out.emplace_back("values", values);
                out.emplace_back("latex", latex);
                break;
            }
            // each root as its proved digits, or within a bound of zero where none is proved
            Array values;
            for (double v : e.values) values.push_back(v);
            std::string latex;
            const auto& roots = e.enclosed->roots;
            for (std::size_t k = 0; k < roots.size(); ++k) {
                const auto p = spatium::proved(roots[k]);
                const std::string name = roots.size() == 1 ? e.name : std::format("{}_{{{}}}", e.name, k + 1);
                latex += (k ? ",\\quad " : "") +
                         (p.digits > 0 ? name + " = " + p.text
                                       : std::format("\\left|{}\\right| \\le {}", name,
                                                     bd::detail::number_text(decade_above(p.bound), bd::detail::TextStyle::Latex)));
            }
            out.emplace_back("unknown", e.name);
            out.emplace_back("values", values);
            out.emplace_back("latex", latex);
            break;
        }
        case Route::Bracket: out.emplace_back("value", e.values.front()); break;
        case Route::Integral:
            out.emplace_back("value", e.values.front());
            out.emplace_back("error_below", decade_above(*e.error));
            break;
        case Route::Synthesis:
            for (auto& [k, v] : answer_json(*e.answer, e.parameter, written)) out.emplace_back(std::move(k), std::move(v));
            break;
    }
    // the roots on the box searched, each proved alone in its interval, whichever route was taken
    if (e.enclosed) {
        const auto& f = *e.enclosed;
        Array roots, undecided;
        for (const auto& r : f.roots) {
            auto p = proved_json(spatium::proved(r));
            p.emplace_back("lo", r.lo);
            p.emplace_back("hi", r.hi);
            roots.push_back(std::move(p));
        }
        for (const auto& u : f.undecided) undecided.push_back(Array{u.lo, u.hi});
        out.emplace_back("enclosed", Object{{"box", Array{f.box.lo, f.box.hi}}, {"roots", roots},
                                            {"undecided", undecided}, {"complete", f.complete()}});
    }
    Array candidates;
    for (std::size_t i = 0; i < e.candidates.size(); ++i) {
        const Candidate& c = e.candidates[i];
        Object o{{"route", route_name(c.route)}, {"arithmetic", arithmetic_name(c.arithmetic)}, {"cost", c.cost},
                 {"taken", e.taken == i}};
        if (!c.refused.empty()) o.emplace_back("refused", c.refused);
        if (!c.values.empty()) {
            Array values;
            for (double v : c.values) values.push_back(v);
            o.emplace_back("values", values);
        }
        if (!c.text.empty()) o.emplace_back("text", c.text);
        candidates.push_back(std::move(o));
    }
    out.emplace_back("candidates", candidates);
    return out;
}

namespace respond_detail {

// The shorthand of one line read as the document it stands for, with the names `above` defined
// before it: the definition asked for, or the route it asked by and why it does not read.
struct Shorthand {
    std::expected<AskedDocument, std::pair<const char*, std::string>> asked;
};

inline Shorthand shorthand(std::string_view q, const spatium::io::build::Document& above) {
    namespace bd = spatium::io::build;
    const auto refuse = [](const char* route, std::string why) {
        return Shorthand{std::unexpected(std::pair{route, std::move(why)})};
    };
    // the statements above, then the one asked; checked as one document
    const auto with_above = [&](bd::Document asked, std::string wanted, const char* route) {
        asked.statements.insert(asked.statements.begin(), above.statements.begin(), above.statements.end());
        if (auto ok = bd::check(asked); !ok) return refuse(route, ok.error().message);
        return Shorthand{AskedDocument{std::move(asked), std::move(wanted)}};
    };
    const auto unknown = [&](const auto& f) {
        std::vector<std::string> open;
        for (auto& n : f.free_symbols())
            if (!above.find(n)) open.push_back(std::move(n));
        return open;
    };
    // an equation, or an expression
    const std::string_view equation = q;
    // no relation but '=' at the end (2 + 2 =): an expression to evaluate, as a calculator does
    if (const auto eq = equation.find_first_of("=<>"); eq == std::string_view::npos || (equation[eq] == '=' && trimmed(equation.substr(eq + 1)).empty())) {
        const std::string_view expression = trimmed(equation.substr(0, eq));
        const auto f = bd::parse(expression, "", "");
        if (!f) return refuse("evaluate", f.error().message);
        if (const auto open = unknown(*f); !open.empty())
            return refuse("evaluate", "evaluate: the expression reads '" + open.front() +
                                          "', which has no value; to find it, write an equation (" + std::string(expression) + " = 0)");
        bd::Statement value;
        value.name = "answer", value.value = {*f, std::string(expression)};
        return with_above({{std::move(value)}}, "answer", "evaluate");
    }
    // the roots of the equation in its one unknown: in closed form, and proved on its bracket or the whole line
    bd::FieldReader r(equation, "", "", {"on"});
    bd::Statement root;
    root.kind = bd::StatementKind::Root;
    auto lhs = bd::detail::written(r, equation);
    if (!lhs) return refuse("closed form", lhs.error().message);
    root.relation = bd::detail::relation(r);
    if (root.relation.empty()) return refuse("closed form", r.fail("'=', '<' or '>'").error().message);
    auto rhs = bd::detail::written(r, equation);
    if (!rhs) return refuse("closed form", rhs.error().message);
    // on [a, b]: the box its roots are sought on, as in a document
    if (!r.done()) {
        auto ends = bd::detail::interval(r, equation, "on");
        if (!ends) return refuse("closed form", ends.error().message);
        std::tie(root.lo, root.hi) = std::move(*ends);
        if (!r.done()) return refuse("closed form", r.fail("the end").error().message);
    }
    // a new name on the left that the right does not read is a definition, f = a^2 + 1, not an equation in f
    const std::string left(trimmed(lhs->text));
    if (const auto right_reads = rhs->field.free_symbols(); root.relation == "=" && !left.empty() && !above.find(left) &&
        lhs->field.free_symbols() == std::vector<std::string>{left} && !std::ranges::contains(right_reads, left)) {
        bd::Statement value;
        value.name = left, value.value = std::move(*rhs);
        return with_above({{std::move(value)}}, left, "evaluate");
    }
    const auto names = unknown(lhs->field - rhs->field);
    const char* route = root.relation == "=" ? "closed form" : "interval method";
    if (names.size() != 1 && root.relation != "=")
        return refuse(route, std::format("an inequality is solved for one unknown; this one reads {}", names.size()));
    if (names.size() != 1) return refuse(route, solve_closed_form(lhs->field - rhs->field, "").error().message);
    root.name = names.front(), root.value = std::move(*lhs), root.rhs = std::move(*rhs);
    return with_above({{std::move(root)}}, names.front(), route);
}

// What a page draws beside an answer, as the request `plot_request` takes: a value of one
// parameter, the document it is defined in and the parameter along the x axis. A root in closed
// form is drawn as its equation's left side less its right along the unknown, over its real
// roots and a margin on each side (around the vertex when it has none); a synthesised answer of
// one expression as that expression along its parameter; a value reading one parameter as
// itself. Anything else has no graph.
inline std::optional<Object> graph_of(const AskedDocument& asked, const Object& out) {
    namespace bd = spatium::io::build;
    const auto field = [&](std::string_view key) -> const JsonValue* {
        for (const auto& [k, v] : out)
            if (k == key) return &v;
        return nullptr;
    };
    if (field("refused")) return std::nullopt;
    const bd::Statement* wanted = asked.document.find(asked.wanted);
    if (!wanted) return std::nullopt;
    std::vector<const bd::Statement*> parameters;
    for (const auto& d : asked.document.statements)
        if (d.kind == bd::StatementKind::Parameter) parameters.push_back(&d);
    const JsonValue* route = field("route");
    const std::string how = route && route->is_string() ? route->as_string() : "";
    const auto request = [](std::string document, std::string name, std::string along) {
        return Object{{"document", std::move(document)}, {"name", std::move(name)}, {"along", std::move(along)}};
    };
    if ((how == "closed form" || how == "every root, proved") && wanted->kind == bd::StatementKind::Root) {
        std::vector<double> roots;
        if (const JsonValue* v = field("values"); v && v->is_array())
            for (const auto& x : v->as_array())
                if (x.is_number() && std::isfinite(x.as_number())) roots.push_back(x.as_number());
        double lo = 0, hi = 0;
        if (!roots.empty()) {
            const auto [a, b] = std::ranges::minmax(roots);
            const double margin = std::max(1.0, 0.5 * (b - a));
            lo = a - margin, hi = b + margin;
        } else if (const JsonValue* s = field("solution"); s && s->is_object() && s->number_or("a", 0) != 0) {
            const double a = s->number_or("a", 0), b = s->number_or("b", 0), D = s->number_or("discriminant", 0);
            const double centre = -b / (2 * a), half = std::max(2.0, std::sqrt(std::abs(D)) / std::abs(a));
            lo = centre - half, hi = centre + half;
        } else {
            return std::nullopt;
        }
        std::string text;
        for (const auto& d : asked.document.statements) {
            if (&d == wanted) break;
            text += bd::to_text(d) + "\n";
        }
        const std::string name = asked.document.find("graph") ? "graph_1" : "graph";
        text += std::format("{} in [{}, {}]\n{} = ({}) - ({})", wanted->name, lo, hi, name, wanted->value.text, wanted->rhs.text);
        return request(std::move(text), name, wanted->name);
    }
    if (parameters.size() != 1) return std::nullopt;
    const std::string& p = parameters.front()->name;
    if (how == "synthesis") {
        const JsonValue* pieces = field("answer");
        if (!pieces || !pieces->is_array() || pieces->as_array().size() != 1) return std::nullopt;
        const JsonValue& piece = pieces->as_array().front();
        const std::string expression = piece.string_or("text", "");
        if (piece.string_or("form", "") != "expression" || expression.empty()) return std::nullopt;
        return request(std::format("{} in [{}, {}]\n{} = {}", p, piece.number_or("lo", 0), piece.number_or("hi", 0), wanted->name, expression),
                       wanted->name, p);
    }
    if (wanted->kind == bd::StatementKind::Value) return request(bd::to_text(asked.document), wanted->name, p);
    return std::nullopt;
}


// ── The answer as one sheet ───────────────────────────────────────
//
// Every answer, whatever route gave it, as the same sheet a page draws (`sheet`): what it is
// (`outcome`: a number, roots, a function, none, a parameter, or refused), how strongly it holds
// (`guarantee`: exact, proved digits, proved complete, estimated, sampled, computed, none) in
// words (`verdict`: its `strength`, "proved to 16 digits", and `notes` on how, "closed form",
// "no other real root"; `label` the two as one line), the result in one line of LaTeX, and its
// sections in one order -- Input, Steps, Roots, Check, Properties, How it was found -- each a list
// of steps, of facts or a table, with whether it starts open and a `gist` where a few words say
// what is inside. A page draws these shapes and never asks which route answered; a route's own
// fields stay beside it as the record.

// a step as a sheet writes it: its label, its formula, and the words a page sets in bold with the
// rest beside them, where the step has them
inline JsonValue step_json(const spatium::io::build::Step& st) {
    Object o{{"label", st.label}, {"latex", st.latex}};
    if (!st.lead.empty()) o.emplace_back("lead", st.lead);
    if (!st.note.empty()) o.emplace_back("note", st.note);
    return o;
}

inline const JsonValue* field_of(const Object& o, std::string_view key) {
    for (const auto& [k, v] : o)
        if (k == key) return &v;
    return nullptr;
}

inline Object sheet_of(const Evaluation& e, spatium::io::build::StatementKind kind, const Object& out) {
    namespace bd = spatium::io::build;
    Array sections;
    const auto section = [&](std::string title, bool open, std::string shape, JsonValue body, Object extra = {}) {
        Object o{{"title", std::move(title)}, {"open", open}, {std::move(shape), std::move(body)}};
        for (auto& kv : extra) o.push_back(std::move(kv));
        sections.push_back(std::move(o));
    };
    const auto fact = [](std::string label, std::string value, bool latex = false, std::optional<bool> ok = std::nullopt) {
        Object f{{"label", std::move(label)}, {latex ? "latex" : "text", std::move(value)}};
        if (ok) f.emplace_back("ok", *ok);
        return JsonValue(f);
    };
    const auto text_of = [](double v) { return bd::detail::number_text(v, bd::detail::TextStyle::Plain); };
    const auto interval_text = [&](const spatium::Interval<double>& v) { return std::format("[{}, {}]", text_of(v.lo), text_of(v.hi)); };
    const auto box_text = [&](const std::vector<spatium::Interval<double>>& b) {
        std::string t;
        for (const auto& v : b) t += (t.empty() ? "" : " × ") + interval_text(v);
        return t;
    };
    // a piece nothing was proved of, and why where it is past what a double holds
    const auto undecided_text = [](double lo, double hi) {
        const bool overflows = std::min(std::abs(lo), std::abs(hi)) >= 1e300 && (lo > 0) == (hi > 0);
        return std::format("[{:.4g}, {:.4g}]{}", lo, hi, overflows ? ", where the equation overflows a double" : "");
    };
    // a power of ten as LaTeX, for a fact, and as text with raised digits, 10⁻¹⁵, for a label
    const auto decade_text = [](double v) {
        return v > 0 ? std::format("10^{{{}}}", static_cast<int>(std::lround(std::log10(v)))) : std::string("0");
    };
    const auto decade_label = [](double v) {
        if (!(v > 0)) return std::string("0");
        static constexpr const char* raised[] = {"⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹"};
        std::string out = "10";
        for (char c : std::to_string(static_cast<int>(std::lround(std::log10(v))))) out += c == '-' ? "⁻" : raised[c - '0'];
        return out;
    };
    std::string outcome = "number", guarantee = "computed", label, strength, result;
    std::vector<std::string> notes;
    bool ok = true;
    // Input: the question as the library read it, and what it asks for
    if (const JsonValue* rows = field_of(out, "question_latex"); rows && rows->is_array() && !rows->as_array().empty()) {
        std::string read;
        for (const auto& r : rows->as_array()) read += (read.empty() ? "" : ",\\quad ") + r.as_string();
        Array facts{fact("Read as", read, true)};
        if (e.name != "answer") facts.push_back(fact("Asks for", bd::detail::name_latex(e.name), true));
        section("Input", false, "facts", facts);
    }
    const JsonValue* latex = field_of(out, "latex");
    if (latex && latex->is_string()) result = latex->as_string();
    const std::optional<Route> route = e.route();

    if (!e.taken) {
        outcome = "refused", guarantee = "none", ok = false;
        strength = "not answered";
        if (route && e.candidates.size() <= 1) notes.push_back(std::format("by {}", route_name(*route)));
        result.clear();
        // a synthesis that answered on part of the interval still shows its pieces below
    }
    if (kind == spatium::io::build::StatementKind::Parameter && e.taken) {
        outcome = "parameter", guarantee = "exact", strength = "a parameter", ok = false;
        if (const JsonValue* rows = field_of(out, "question_latex"); rows && rows->is_array() && !rows->as_array().empty())
            result = rows->as_array().back().as_string();
    } else if (e.taken) switch (*route) {
        case Route::Evaluate: {
            // the value alone, or the name and its value: the question stands above it, and a long
            // one would push the answer off a phone's screen
            const Candidate& c = e.candidates[*e.taken];
            std::string value = number_latex(c, e.values.front());
            if (e.exact && c.arithmetic == Arithmetic::Double) {
                const auto d = e.exact->decimal();
                value = e.exact->den() == 1 ? *d : bd::detail::fraction_latex(*e.exact) + (d ? " = " + *d : "");
            } else if (e.proved && !e.proved->exact && e.proved->digits > 0) {
                value = "\\approx " + e.proved->text;
            }
            result = (e.name == "answer" ? std::string() : bd::detail::name_latex(e.name) + (value.starts_with("\\approx") ? " " : " = ")) + value;
            if (e.proved && e.proved->exact) guarantee = "exact", strength = "exact";
            else if (e.proved && e.proved->digits > 0)
                guarantee = "proved", strength = std::format("proved to {} digit{}", e.proved->digits, e.proved->digits == 1 ? "" : "s");
            else if (e.proved)
                guarantee = "none", ok = false, strength = "no digit proved",
                notes.push_back(std::format("within {} of zero", decade_label(e.proved->bound)));
            else guarantee = "computed", ok = false, strength = "computed in double";
            if (e.steps.size() > 1) {
                Array steps;
                for (const auto& st : e.steps) steps.push_back(step_json(st));
                section("Steps", true, "steps", steps, Object{{"gist", std::format("{} steps", e.steps.size())}});
            }
            Array facts{fact("In double", text_of(e.values.front()))};
            if (e.proved && !e.proved->exact && e.proved->digits > 0) facts.push_back(fact("Proved", e.proved->text, false, true));
            if (!e.proved && !e.unproved.empty()) facts.push_back(fact("Not proved", e.unproved, false, false));
            section("Check", false, "facts", facts);
            break;
        }
        case Route::Expression: {
            outcome = "function", guarantee = "exact", ok = false;
            strength = "exact", notes.push_back(std::format("a function of {}, as written", e.free.empty() ? std::string("its parameter") : e.free.front()));
            break;
        }
        case Route::ClosedForm: {
            outcome = "roots";
            const auto& c = *e.closed;
            std::vector<bd::Step> chain, checks;
            for (const auto& st : c.steps) (st.label.starts_with("check") ? checks : chain).push_back(st);
            const auto named = std::ranges::find_if(chain, [](const bd::Step& st) { return st.label == "roots" || st.label == "root"; });
            result = named != chain.end() ? named->latex : (chain.empty() ? result : chain.back().latex);
            Array steps;
            for (const auto& st : chain) steps.push_back(step_json(st));
            section("Steps", true, "steps", steps, Object{{"gist", std::format("closed form, {} steps", chain.size())}});
            Array facts;
            const bool vieta = c.check && c.check->error <= 1e-12;
            for (const auto& st : checks) facts.push_back(fact("Vieta's formulas", st.latex, true, vieta));
            // the second route: each root of the formula inside an interval proved to hold one
            bool agree = false, complete = false;
            if (e.enclosed) {
                std::string roots;
                for (const auto& r : e.enclosed->roots) {
                    const auto p = spatium::proved(r);
                    roots += (roots.empty() ? "" : ", ") + (p.digits > 0 ? std::format("{} ({} digits)", p.text, p.digits) : std::string("0"));
                }
                agree = e.enclosed->roots.size() == e.values.size() &&
                        std::ranges::all_of(e.values, [&](double v) {
                            return std::ranges::any_of(e.enclosed->roots, [&](const auto& r) {
                                return std::nextafter(r.lo, -INFINITY) <= v && v <= std::nextafter(r.hi, INFINITY);
                            });
                        });
                complete = agree && e.enclosed->complete();
                facts.push_back(fact("Every real root, proved", roots.empty() ? std::string("none") : roots, false, agree));
                for (const auto& u : e.enclosed->undecided) facts.push_back(fact("Undecided", undecided_text(u.lo, u.hi), false, false));
            }
            const bool proved_all = agree;
            section("Check", true, "facts", facts);
            Array props{fact("Discriminant", text_of(c.solution.discriminant))};
            section("Properties", false, "facts", props);
            // proved digits are of real roots; complex ones are computed, and checked by Vieta's formulas
            guarantee = c.exact ? "exact" : proved_all && !e.values.empty() ? "proved" : "computed";
            strength = guarantee, notes.push_back("closed form");
            if (vieta) notes.push_back("checked by Vieta's formulas");
            if (complete) notes.push_back(e.values.empty() ? "no real root" : "no other real root");
            break;
        }
        case Route::Intervals: {
            const Axis& a = *e.axis;
            outcome = "set", result = set_latex(a, e.name);
            const auto undecided = std::ranges::count_if(a.parts, [](const Axis::Part& p) { return p.kind == Axis::Kind::Undecided; });
            guarantee = "proved", ok = undecided == 0;
            strength = "proved", notes.push_back("by the interval method");
            if (undecided) notes.push_back(std::format("except on {} piece{} undecided", undecided, undecided == 1 ? "" : "s"));
            section("Axis", true, "axis", axis_json(a));
            Array facts;
            for (const auto& p : a.parts)
                if (p.kind == Axis::Kind::Root) facts.push_back(fact("End, a root", std::format("{} in [{}, {}]", p.text, text_of(p.lo), text_of(p.hi)), false, true));
                else if (p.kind == Axis::Kind::Break) facts.push_back(fact("Break, no value", std::format("in [{}, {}]", text_of(p.lo), text_of(p.hi))));
                else if (p.kind == Axis::Kind::Undecided) facts.push_back(fact("Undecided", std::format("[{}, {}]", text_of(p.lo), text_of(p.hi)), false, false));
            if (!facts.empty()) section("Ends", false, "facts", facts);
            break;
        }
        case Route::Linear: {
            const auto& l = *e.linear;
            outcome = "system";
            const bool exact = !l.exact.empty();
            guarantee = exact ? "exact" : "computed", ok = exact;
            strength = guarantee, notes.push_back(exact ? "by elimination in fractions" : "by elimination in doubles");
            result = l.steps.back().latex;
            Array steps;
            for (const auto& st : l.steps) steps.push_back(step_json(st));
            section("Steps", true, "steps", steps, Object{{"gist", std::format("elimination, {} steps", l.steps.size())}});
            std::string residuals;
            for (double r : l.residuals) residuals += (residuals.empty() ? "" : ", ") + text_of(r);
            const bool holds = std::ranges::all_of(l.residuals, [](double r) { return std::abs(r) <= 1e-12; });
            section("Check", !holds, "facts", Array{fact("Each equation at the solution, left less right", residuals, false, holds)});
            break;
        }
        case Route::Enclosed: {
            outcome = "roots";
            if (e.system) {
                const auto& sys = *e.system;
                const std::size_t n = sys.roots.size();
                guarantee = sys.undecided.empty() ? "complete" : "proved";
                strength = "proved", notes.push_back(std::format("{} root{}", n, n == 1 ? "" : "s"));
                notes.push_back(sys.undecided.empty() ? "no other on the box" : "elsewhere undecided");
                Array columns, rows;
                for (const auto& u : sys.unknowns) columns.push_back(u);
                columns.push_back("digits"), columns.push_back("box");
                for (const auto& r : sys.roots) {
                    Array row;
                    int digits = 17;
                    for (const auto& v : r) {
                        const auto p = spatium::proved(v);
                        row.push_back(p.text), digits = std::min(digits, p.digits);
                    }
                    row.push_back(std::to_string(digits)), row.push_back(box_text(r));
                    rows.push_back(std::move(row));
                }
                section("Roots", true, "table", Object{{"columns", columns}, {"rows", rows}},
                        Object{{"gist", std::format("{} root{}", n, n == 1 ? "" : "s")}});
                // the second route: each equation on its root's box, an interval that holds zero
                Array checks;
                bool hold = true;
                for (std::size_t k = 0; k < sys.at_roots.size(); ++k) {
                    const bool zero = std::ranges::all_of(sys.at_roots[k], [](const auto& v) { return v.lo <= 0 && 0 <= v.hi; });
                    hold = hold && zero;
                    std::string t;
                    for (const auto& v : sys.at_roots[k]) t += (t.empty() ? "" : ", ") + std::format("[{:.3g}, {:.3g}]", v.lo, v.hi);
                    checks.push_back(fact(std::format("Root {}, each equation on its box", k + 1), t, false, zero));
                }
                if (!checks.empty()) section("Check", !hold, "facts", checks);
                Array facts{fact("Searched", "every double of each unknown"),
                            fact("Cut into", std::format("{} pieces, {} proved to hold no root", sys.pieces, sys.excluded))};
                for (std::size_t k = 0; k < sys.undecided.size() && k < 8; ++k) facts.push_back(fact("Undecided", box_text(sys.undecided[k]), false, false));
                if (sys.undecided.size() > 8) facts.push_back(fact("Undecided", std::format("and {} more", sys.undecided.size() - 8), false, false));
                section("Properties", false, "facts", facts);
                break;
            }
            const auto& f = *e.enclosed;
            const std::size_t n = f.roots.size();
            guarantee = f.complete() ? "complete" : "proved";
            strength = "proved", notes.push_back(std::format("{} root{}", n, n == 1 ? "" : "s"));
            notes.push_back(f.complete() ? "no other on the box" : "elsewhere undecided");
            Array rows;
            for (const auto& r : f.roots) {
                const auto p = spatium::proved(r);
                rows.push_back(Array{p.digits > 0 ? p.text : std::string("0"), p.digits > 0 ? std::to_string(p.digits) : std::string("-"),
                                     std::format("[{}, {}]", text_of(r.lo), text_of(r.hi))});
            }
            section("Roots", true, "table", Object{{"columns", Array{"root", "digits", "interval"}}, {"rows", rows}},
                    Object{{"gist", std::format("{} root{}", n, n == 1 ? "" : "s")}});
            const bool line = f.box.hi == std::numeric_limits<double>::max() && f.box.lo == -std::numeric_limits<double>::max();
            Array facts{fact("Searched", line ? std::string("every double") : std::format("[{}, {}]", text_of(f.box.lo), text_of(f.box.hi))),
                        fact("Cut into", std::format("{} pieces, {} proved to hold no root", f.pieces, f.excluded))};
            for (const auto& u : f.undecided) facts.push_back(fact("Undecided", undecided_text(u.lo, u.hi), false, false));
            section("Properties", false, "facts", facts);
            break;
        }
        case Route::Bracket: {
            result = std::format("{} = {}", bd::detail::name_latex(e.name), bd::detail::number_text(e.values.front(), bd::detail::TextStyle::Latex));
            guarantee = "computed", strength = "computed", notes.push_back("a root on its bracket, by bisection to the last bit");
            break;
        }
        case Route::Integral: {
            result = std::format("{} = {}", bd::detail::name_latex(e.name), bd::detail::number_text(e.values.front(), bd::detail::TextStyle::Latex));
            guarantee = "estimated", strength = std::format("error below {}, estimated", decade_label(e.error ? *e.error : 0));
            break;
        }
        case Route::Synthesis: break;
    }
    // a synthesis, answered or in part: its pieces and their certificate
    if (route == Route::Synthesis && e.answer) {
        outcome = e.taken ? "function" : outcome;
        const JsonValue* pieces = field_of(out, "answer");
        const JsonValue* cert = field_of(out, "certificate");
        // one piece as its formula, an expression or a Chebyshev sum; more as how many there are
        if (e.taken && pieces && pieces->is_array() && pieces->as_array().size() == 1)
            result = std::format("{} = {}", bd::detail::name_latex(e.name), pieces->as_array().front().string_or("latex", ""));
        else if (e.taken) result = std::format("\\text{{piecewise, {} pieces}}", e.answer->pieces.size());
        if (cert && cert->is_object()) {
            const double below = cert->number_or("error_below", 0);
            const bool relative = cert->string_or("norm", "") == "largest relative error";
            if (e.taken) {
                guarantee = "sampled";
                strength = below > 0 ? std::format("{}error below {}", relative ? "relative " : "", decade_label(below)) : "exact on the grid";
                notes.push_back(std::format("certified on {} points", cert->number_or("points", 0)));
            }
            section("Certificate", true, "facts",
                    Array{fact("Strength", std::format("{}, on {} points", cert->string_or("strength", ""), cert->number_or("points", 0))),
                          fact(relative ? "Relative error below" : "Absolute error below", decade_text(below), true),
                          fact("Scalar", cert->string_or("scalar", ""))});
        }
        if (pieces && pieces->is_array() && pieces->as_array().size() > 1) {
            Array rows;
            for (const auto& pc : pieces->as_array())
                rows.push_back(Array{std::format("[{}, {}]", text_of(pc.number_or("lo", 0)), text_of(pc.number_or("hi", 0))),
                                     pc.string_or("form", ""), pc.string_or("text", pc.string_or("latex", ""))});
            section("Pieces", false, "table", Object{{"columns", Array{"on", "form", "expression"}}, {"rows", rows}});
        }
    }
    // a refusal no route was tried for says why here; one with routes says it route by route below
    if (outcome == "refused" && e.candidates.empty()) section("Check", true, "facts", Array{fact("Why", e.refused, false, false)});
    // How it was found: every route tried, its cost and what it gave or why not
    Array rows;
    for (std::size_t i = 0; i < e.candidates.size(); ++i) {
        const Candidate& c = e.candidates[i];
        std::string gave = c.refused;
        if (gave.empty() && c.route == Route::Synthesis && e.answer)
            gave = e.answer->pieces.size() == 1 && e.answer->pieces.front().form == Form::Expression
                       ? *bd::to_text(e.answer->pieces.front().expr, e.parameter) : std::format("{} pieces", e.answer->pieces.size());
        else if (gave.empty()) {
            for (double v : c.values) gave += (gave.empty() ? "" : ", ") + text_of(v);
            if (!c.text.empty()) gave = c.text;
            if (gave.empty() && c.route == Route::ClosedForm) gave = "no real root";
        }
        rows.push_back(Array{std::string(route_name(c.route)) + (e.taken == i ? ", taken" : ""), std::to_string(c.cost), gave});
    }
    if (!rows.empty()) {
        Object extra{{"gist", std::format("{} route{}", rows.size(), rows.size() == 1 ? "" : "s")}};
        if (route) extra.emplace_back("note", route_caption(*route));
        section("How it was found", false, "table", Object{{"columns", Array{"route", "cost", "gave"}}, {"rows", rows}}, std::move(extra));
    }
    // Every route that answered, each with its own chain, in the order they were tried: the taken
    // one's is the Steps above; a closed form's check by Vieta's formulas is a chain of its own
    Array routes;
    const auto chain_of = [](const std::vector<bd::Step>& steps) {
        Array out;
        for (const auto& st : steps) out.push_back(step_json(st));
        return out;
    };
    const auto add_route = [&](std::string name, std::size_t cost, bool taken, Array steps) {
        if (steps.empty()) return;
        const std::size_t n = steps.size();
        routes.push_back(Object{{"route", std::move(name)}, {"cost", cost}, {"taken", taken},
                                {"gist", std::format("{} step{}", n, n == 1 ? "" : "s")}, {"steps", std::move(steps)}});
    };
    const auto x_of = [&](const std::string& name) { return bd::detail::name_latex(name); };
    for (std::size_t i = 0; i < e.candidates.size(); ++i) {
        const Candidate& c = e.candidates[i];
        if (!c.refused.empty()) continue;
        const bool taken = e.taken == i;
        const std::string name = route_name(c.route);
        const auto one = [&](std::string label, std::string latex) { return Array{Object{{"label", std::move(label)}, {"latex", std::move(latex)}}}; };
        switch (c.route) {
            case Route::Evaluate:
                if (e.steps.size() > 1) add_route(name, c.cost, taken, chain_of(e.steps));
                else if (taken && !result.empty()) add_route(name, c.cost, taken, one("value", result));
                break;
            case Route::ClosedForm: {
                if (!e.closed) break;
                std::vector<bd::Step> chain, checks;
                for (const auto& st : e.closed->steps) (st.label.starts_with("check") ? checks : chain).push_back(st);
                add_route(name, c.cost, taken, chain_of(chain));
                // Vieta's formulas: the equation and its coefficients, then the sum and the product of the roots
                if (!checks.empty() && chain.size() >= 2) {
                    std::vector<bd::Step> vieta{chain[0], chain[1]};
                    vieta.insert(vieta.end(), checks.begin(), checks.end());
                    add_route("Vieta's formulas", 0, false, chain_of(vieta));
                }
                break;
            }
            case Route::Linear:
                if (e.linear) add_route(name, c.cost, taken, chain_of(e.linear->steps));
                break;
            case Route::Intervals: {
                if (!e.axis) break;
                Array steps;
                for (const auto& part : e.axis->parts)
                    if (part.kind == Axis::Kind::Root)
                        steps.push_back(Object{{"label", "an end: a root, proved alone by Newton's step"}, {"latex", std::format("{} \\approx {}", x_of(e.name), part.text)}});
                    else if (part.kind == Axis::Kind::Break)
                        steps.push_back(Object{{"label", "a break: no value"}, {"latex", std::format("{} \\in \\left[{}, {}\\right]", x_of(e.name), text_of(part.lo), text_of(part.hi))}});
                steps.push_back(Object{{"label", "where it holds"}, {"latex", set_latex(*e.axis, e.name)}});
                add_route(name, c.cost, taken, std::move(steps));
                break;
            }
            case Route::Enclosed: {
                Array steps;
                const auto searched = [&](std::string latex) { steps.push_back(Object{{"label", "searched"}, {"latex", std::move(latex)}}); };
                const auto cut = [&](std::size_t pieces, std::size_t excluded, std::size_t undecided) {
                    steps.push_back(Object{{"label", "pieces cut where undecided: with no root + with one + undecided"},
                                           {"latex", std::format("{} = {} + {} + {}", pieces, excluded, pieces - excluded - undecided, undecided)}});
                };
                if (e.system) {
                    const auto& sys = *e.system;
                    std::string names;
                    for (const auto& u : sys.unknowns) names += (names.empty() ? "" : ", ") + x_of(u);
                    searched(std::format("{} \\in \\mathbb{{R}}", names));
                    cut(sys.pieces, sys.excluded, sys.undecided.size());
                    for (const auto& r : sys.roots) {
                        std::string at;
                        for (std::size_t j = 0; j < r.size(); ++j) at += (j ? ",\\ " : "") + x_of(sys.unknowns[j]) + " \\approx " + spatium::proved(r[j]).text;
                        steps.push_back(Object{{"label", "Krawczyk's step lands inside: one root"}, {"latex", at}});
                    }
                } else if (e.enclosed) {
                    const auto& f = *e.enclosed;
                    const bool line = f.box.hi == std::numeric_limits<double>::max() && f.box.lo == -std::numeric_limits<double>::max();
                    searched(line ? x_of(e.name) + " \\in \\mathbb{R}"
                                  : std::format("{} \\in \\left[{}, {}\\right]", x_of(e.name), text_of(f.box.lo), text_of(f.box.hi)));
                    cut(f.pieces, f.excluded, f.undecided.size());
                    for (const auto& r : f.roots) {
                        const auto p = spatium::proved(r);
                        steps.push_back(Object{{"label", "Newton's step lands inside: one root"},
                                               {"latex", std::format("{} \\approx {}", x_of(e.name), p.digits > 0 ? p.text : std::string("0"))}});
                    }
                }
                add_route(name, c.cost, taken, std::move(steps));
                break;
            }
            case Route::Bracket:
                if (!c.values.empty())
                    add_route(name, c.cost, taken, one("bisection to the last bit", std::format("{} = {}", x_of(e.name), text_of(c.values.front()))));
                break;
            case Route::Expression: case Route::Integral: case Route::Synthesis:
                if (taken && !result.empty()) add_route(name, c.cost, taken, one(name, result));
                break;
        }
    }
    // Next: what else this question can become, each written as a question of its own -- the same
    // equation as an inequality, on an interval about its roots, its integral between two roots;
    // an inequality as its equation, or the other side of it. Only for a root alone on its page.
    Array next;
    const JsonValue* doc = field_of(out, "document");
    const JsonValue* statements = doc ? doc->find("statements") : nullptr;
    if (kind == bd::StatementKind::Root && statements && statements->is_array() && statements->as_array().size() == 1) {
        const JsonValue& st = statements->as_array().front();
        const std::string name = st.string_or("name", ""), lhs = st.string_or("lhs", ""), rhs = st.string_or("rhs", ""),
                          rel = st.string_or("relation", "=");
        const auto offer = [&](std::string caption, std::string text) {
            const auto asked = bd::read_document(text);
            if (!asked || asked->statements.size() != 1) return;
            if (const auto latex = bd::to_latex(asked->statements.front()))
                next.push_back(Object{{"caption", std::move(caption)}, {"text", std::move(text)}, {"latex", *latex}});
        };
        if (!st.find("also") && !name.empty()) {
            const std::string bracket = st.find("lo") ? std::format(" on [{}, {}]", st.string_or("lo", ""), st.string_or("hi", "")) : "";
            if (rel == "=") {
                offer("Where it is greater", std::format("{}: {} > {}{}", name, lhs, rhs, bracket));
                std::vector<double> real = e.values;
                std::ranges::sort(real);
                if (bracket.empty() && !real.empty() && std::isfinite(real.front()) && std::isfinite(real.back()))
                    offer("Only on an interval about its roots", std::format("{}: {} = {} on [{}, {}]", name, lhs, rhs,
                                                                             std::floor(real.front()) - 1, std::ceil(real.back()) + 1));
                if (real.size() >= 2 && real.front() < real.back())
                    offer("Its integral between the outer roots", std::format("I = integral {} d{} on [{}, {}]",
                                                                              rhs == "0" ? lhs : std::format("({}) - ({})", lhs, rhs), name,
                                                                              detail::shortest(real.front()), detail::shortest(real.back())));
            } else {
                const std::string other = rel == "<" ? ">=" : rel == "<=" ? ">" : rel == ">" ? "<=" : "<";
                offer("Where it is equal", std::format("{}: {} = {}{}", name, lhs, rhs, bracket));
                offer("Where it does not hold", std::format("{}: {} {} {}{}", name, lhs, other, rhs, bracket));
            }
        }
    }
    if (!next.empty()) section("Next", true, "next", next, Object{{"gist", std::format("{} question{}", next.size(), next.size() == 1 ? "" : "s")}});
    // the label is the verdict in one line: "Proved, 4 roots, no other on the box"
    label = strength;
    if (!label.empty()) label[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(label[0])));
    Array note_list;
    for (const auto& n : notes) label += ", " + n, note_list.push_back(n);
    return Object{{"outcome", outcome}, {"guarantee", guarantee}, {"verdict", Object{{"strength", strength}, {"notes", note_list}}},
                  {"label", label}, {"ok", ok}, {"result", result}, {"sections", sections}, {"routes", routes}};
}

// How a request asks for its values to be computed: `engine`, {"arithmetic": "rational",
// "strict": false, "tolerance": 1e-8}, each optional, a `tolerance` beside it the synthesis's
// too; or why it names an arithmetic the engine does not have.
inline std::expected<Engine, std::string> engine_of(const JsonValue& request) {
    Engine engine{.tolerance = request.number_or("tolerance", 1e-8)};
    const JsonValue* e = request.find("engine");
    if (!e || !e->is_object()) return engine;
    const std::string a = e->string_or("arithmetic", "double");
    const auto known = std::ranges::find_if(arithmetics, [&](Arithmetic x) { return a == arithmetic_name(x); });
    if (known == std::end(arithmetics)) {
        std::string names;
        for (Arithmetic x : arithmetics) names += std::string(names.empty() ? "" : ", ") + arithmetic_name(x);
        return std::unexpected(std::format("the engine has no arithmetic '{}'; it has {}", a, names));
    }
    engine.arithmetic = *known;
    if (const JsonValue* strict = e->find("strict"); strict && strict->is_bool()) engine.strict = strict->as_bool();
    engine.tolerance = e->number_or("tolerance", engine.tolerance);
    return engine;
}

// Lines that are each an equation, two or more, reading as many unknowns as there are lines,
// written as the system they are: `x, y: x + 2y = 5, 3x - y = 1`, asked for its first unknown.
inline std::optional<AskedDocument> system_of(std::string_view text) {
    namespace bd = spatium::io::build;
    std::vector<std::string> names, equations;
    while (!text.empty()) {
        const auto end = text.find('\n');
        const std::string_view line = trimmed(text.substr(0, end));
        text = end == std::string_view::npos ? std::string_view{} : text.substr(end + 1);
        if (line.empty() || line.starts_with('#')) continue;
        bd::FieldReader r(line, "", "");
        auto lhs = r.expression();
        if (!lhs || !r.take('=')) return std::nullopt;
        auto rhs = r.expression();
        if (!rhs || !r.done()) return std::nullopt;
        for (auto& n : (*lhs - *rhs).free_symbols())
            if (std::ranges::find(names, n) == names.end()) names.push_back(std::move(n));
        equations.emplace_back(line);
    }
    if (equations.size() < 2 || names.size() != equations.size()) return std::nullopt;
    std::string line;
    for (std::size_t k = 0; k < names.size(); ++k) line += (k ? ", " : "") + names[k];
    for (std::size_t k = 0; k < equations.size(); ++k) line += (k ? ", " : ": ") + equations[k];
    auto doc = bd::read_document(line);
    if (!doc) return std::nullopt;
    return AskedDocument{std::move(*doc), names.front()};
}

} // namespace respond_detail

// What a person types, answered. One line is a shorthand (`2 + 2`, `2x^2 - 3x - 5 = 0`, `solve
// ... for ...`), or else a line of a document (`x: x^2 - 2 = 0 on [0, 2]`). Several lines are a
// document and its question: the lines above the last define names, and the last is the
// shorthand read with them (`p in [0.5, 4]` above `solve`-less `x^2 - p = 0` asks for x as a
// function of p), or else the document's last definition.
inline spatium::io::JsonValue respond(std::string_view text, Engine engine = {}) {
    using namespace respond_detail;
    namespace bd = spatium::io::build;
    const std::string_view all = trimmed(text);
    Object out{{"question", std::string(all)}};
    const auto answered = [&](const AskedDocument& asked) {
        const Evaluation e = evaluate(asked.document, asked.wanted, {}, engine);
        for (auto& [k, v] : evaluation_json(e)) out.emplace_back(std::move(k), std::move(v));
        out.emplace_back("document", bd::to_json(asked.document));
        // the question as the rows of a formula editor, written by the library
        Array rows;
        for (const auto& d : asked.document.statements) {
            if (auto row = bd::to_latex(d)) rows.push_back(*row);
            for (const auto& e : d.also)   // a system: an equation a row, as a card writes them
                if (auto row = bd::to_latex(e)) rows.push_back(*row);
        }
        out.emplace_back("question_latex", rows);
        if (auto g = graph_of(asked, out)) out.emplace_back("graph", std::move(*g));
        const bd::Statement* wanted = asked.document.find(asked.wanted);
        out.emplace_back("sheet", sheet_of(e, wanted ? wanted->kind : bd::StatementKind::Value, out));
        return JsonValue(std::move(out));
    };
    // the rows above the last, each a line of the document: a statement as written, or a row a
    // page writes as no statement -- an expression, an equation -- as the same shorthand the last
    // row is, an expression named answer_k by its row so a row below can read it
    bd::Document above;
    // rows that are all equations, as many as the unknowns they read: one system, every unknown
    // solved together (`x + 2y = 5` above `3x - y = 1`)
    if (auto system = system_of(all)) return answered(*system);
    std::string_view q = all;
    if (const auto nl = all.rfind('\n'); nl != std::string_view::npos) {
        std::string_view rest = all.substr(0, nl);
        for (std::size_t k = 1; !rest.empty(); ++k) {
            const auto end = rest.find('\n');
            const std::string_view line = trimmed(rest.substr(0, end));
            rest = end == std::string_view::npos ? std::string_view{} : rest.substr(end + 1);
            if (line.empty() || line.starts_with('#')) continue;
            const auto refuse = [&](std::string why) {
                out.emplace_back("route", "document");
                out.emplace_back("refused", std::format("line {}: {}", k, why));
                return JsonValue(std::move(out));
            };
            if (auto doc = bd::read_document(bd::to_text(above) + std::string(line))) {
                above = std::move(*doc);
                continue;
            }
            auto row = shorthand(line, above);
            if (!row.asked) return refuse(row.asked.error().second);
            bd::Statement st = row.asked->document.statements.back();
            if (st.name == "answer") st.name = std::format("answer_{}", k);
            above.statements.push_back(std::move(st));
        }
        q = trimmed(all.substr(nl + 1));
    }
    auto s = shorthand(q, above);
    if (s.asked) return answered(*s.asked);
    // a shorthand that does not read, but a document that does
    if (auto whole = bd::read_document(bd::to_text(above) + std::string(q)); whole && !whole->statements.empty()) {
        std::string wanted = whole->statements.back().name;
        return answered({std::move(*whole), std::move(wanted)});
    }
    // a line in a document's own words, a name and ':' or 'in', is refused as the document refuses it
    if (q.find(':') != std::string_view::npos || q.find(" in ") != std::string_view::npos) {
        if (auto whole = bd::read_document(bd::to_text(above) + std::string(q)); !whole) {
            out.emplace_back("route", "document");
            out.emplace_back("refused", whole.error().message);
            return JsonValue(std::move(out));
        }
    }
    // refused under the route the shorthand asked for
    out.emplace_back("route", s.asked.error().first);
    out.emplace_back("refused", std::move(s.asked.error().second));
    return JsonValue(std::move(out));
}

// The rows of a LaTeX editor, each a formula as MathJSON, answered as the lines they read as
// (`mathjson_row_text`): `{"rows": [...]}`, computed as its `engine` asks (`engine_of`). The
// answer is the one `respond` gives the same lines typed, and a row that does not read is
// refused with its number, `row`, counted from 1.
inline spatium::io::JsonValue respond_rows(std::string_view request) {
    using namespace respond_detail;
    namespace bd = spatium::io::build;
    const auto json = spatium::io::parse(request);
    const JsonValue* rows = json && json->is_object() ? json->find("rows") : nullptr;
    if (!rows || !rows->is_array() || rows->as_array().empty())
        return JsonValue(Object{{"question", ""}, {"route", "document"}, {"refused", "the request has no rows"}});
    std::string text;
    for (std::size_t k = 0; k < rows->as_array().size(); ++k) {
        // a row is MathJSON, or a line already written in the language, {"text": "..."}
        const JsonValue& row = rows->as_array()[k];
        const JsonValue* typed = row.is_object() ? row.find("text") : nullptr;
        auto line = typed && typed->is_string() ? spatium::Result<std::string>(typed->as_string()) : bd::mathjson_row_text(row);
        if (!line)
            return JsonValue(Object{{"question", text}, {"route", "document"}, {"row", static_cast<double>(k + 1)},
                                    {"refused", std::format("row {}: {}", k + 1, line.error().message)}});
        text += (k ? "\n" : "") + *line;
    }
    const auto engine = engine_of(*json);
    if (!engine) return JsonValue(Object{{"question", text}, {"route", "document"}, {"refused", engine.error()}});
    return respond(text, *engine);
}

// A document given to the door as one JSON object: `document` its text or its JSON, `name` the
// definition wanted, `chosen` the parameters given a value ({"p": 2}), and `engine`, how it is
// computed (`engine_of`). The answer is `name` and the evaluation, or what could not be read.
inline spatium::io::JsonValue evaluate_request(std::string_view request) {
    using namespace respond_detail;
    namespace bd = spatium::io::build;
    const auto refuse = [](std::string why) { return JsonValue(Object{{"refused", std::move(why)}}); };
    const auto json = spatium::io::parse(request);
    if (!json || !json->is_object()) return refuse("the request is not a JSON object");
    const JsonValue* given = json->find("document");
    spatium::Result<bd::Document> doc = std::unexpected(spatium::Error(spatium::ErrorCode::ParseError, "the request has no document"));
    if (given && given->is_string()) doc = bd::read_document(given->as_string());
    else if (given && given->is_object()) doc = bd::document_from_json(*given);
    if (!doc) return refuse(doc.error().message);
    Chosen chosen;
    if (const JsonValue* c = json->find("chosen"); c && c->is_object())
        for (const auto& [k, v] : c->as_object()) {
            if (!v.is_number()) return refuse(std::format("the value chosen for {} is not a number", k));
            chosen.emplace_back(k, v.as_number());
        }
    const auto engine = engine_of(*json);
    if (!engine) return refuse(engine.error());
    const std::string name = json->string_or("name", "");
    Object out{{"name", name}};
    for (auto& [k, v] : evaluation_json(evaluate(*doc, name, chosen, *engine)))
        out.emplace_back(std::move(k), std::move(v));
    return out;
}

// A value of a document drawn as the graph of a function of one of its parameters, as one JSON
// object: `document` (text or JSON), `name`, and `along`, the parameter on the x axis. The other
// parameters it reads stay symbols, read from the shader's slots, and come back with their
// intervals so a page can give each a slider. With `grid` -- {"width", "height", "view": [x0, y0,
// x1, y1], "slots": [...]} -- the answer also holds F at every pixel's centre as the shader's
// check mode computes it, run in fp32 on the host: what a browser holds its WebGL2 frame against.
// With `curve` -- {"x0", "x1", "n", "slots": [...]} -- it holds F at the centres of n columns
// from x0 to x1 in doubles: what a page draws where it has zoomed past the steps of fp32.
inline spatium::io::JsonValue plot_request(std::string_view request) {
    using namespace respond_detail;
    namespace bd = spatium::io::build;
    const auto refuse = [](std::string why) { return JsonValue(Object{{"refused", std::move(why)}}); };
    const auto json = spatium::io::parse(request);
    if (!json || !json->is_object()) return refuse("the request is not a JSON object");
    std::string name = json->string_or("name", ""), along = json->string_or("along", "");
    Expr f;
    std::function<Array(const std::string&)> interval;
    spatium::Result<bd::Document> doc = std::unexpected(spatium::Error(spatium::ErrorCode::ParseError, "the request has no document"));
    if (const JsonValue* n = json->find("net"); n && n->is_object()) {
        // a net of one input, y = net(x): its width and its weights in `TanhNet`'s order
        std::vector<double> weights;
        if (const JsonValue* w = n->find("weights"); w && w->is_array())
            for (const auto& x : w->as_array()) weights.push_back(x.as_number());
        auto made = Expr::net_of(std::move(weights), static_cast<std::size_t>(n->number_or("hidden", 0)), {Expr::u()});
        if (!made) return refuse(made.error().message);
        f = std::move(*made);
        if (name.empty()) name = "net";
        if (along.empty()) along = "x";
        Array range{-1.0, 1.0};
        if (const JsonValue* r = n->find("interval"); r && r->is_array() && r->as_array().size() == 2)
            range = Array{r->as_array()[0].as_number(), r->as_array()[1].as_number()};
        interval = [range](const std::string&) { return range; };
    } else {
        const JsonValue* given = json->find("document");
        if (given && given->is_string()) doc = bd::read_document(given->as_string());
        else if (given && given->is_object()) doc = bd::document_from_json(*given);
        if (!doc) return refuse(doc.error().message);
        const bd::Statement* d = doc->find(name);
        if (!d) return refuse(std::format("the document defines no '{}'", name));
        if (d->kind != bd::StatementKind::Value)
            return refuse(std::format("{} is a {}; a graph is drawn of a value", d->name, bd::statement_kind_name(d->kind)));
        auto u = detail::unfolded(*doc, *d, {});
        if (!u) return refuse(u.error().message);
        if (std::ranges::find(u->free, along) == u->free.end())
            return refuse(std::format("{} does not read the parameter '{}'", d->name, along));
        f = u->f.bind(along, Expr::u());
        interval = [&doc](const std::string& p) {
            const bd::Statement* s = doc->find(p);
            return Array{s->lo.field(0.0, 0.0), s->hi.field(0.0, 0.0)};
        };
    }
    const spatium::Result<bd::PlotShader> shader = bd::plot_shader(f);
    if (!shader) return refuse(shader.error().message);
    Array symbols;
    for (const auto& s : shader->symbols) symbols.push_back(Object{{"name", s}, {"interval", interval(s)}});
    Array weights;
    for (float w : shader->weights) weights.push_back(static_cast<double>(w));
    Object out{{"name", name}, {"along", along}, {"interval", interval(along)}, {"vertex", shader->vertex},
               {"fragment", shader->fragment}, {"symbols", symbols}, {"weights", weights}};
    // the heights the graph reaches over its interval, the other parameters at their middles, in
    // the fp32 the shader computes in: where a page puts the y axis
    if (auto pod = bd::lower(f)) {
        const auto twin = bd::narrowed(*pod);
        const Array along_ends = interval(along);
        std::vector<float> middles;
        for (const auto& s : shader->symbols) {
            const Array ends = interval(s);
            middles.push_back(static_cast<float>(0.5 * (ends[0].as_number() + ends[1].as_number())));
        }
        constexpr int samples = 256;
        const float view[4]{static_cast<float>(along_ends[0].as_number()), 0.0f, static_cast<float>(along_ends[1].as_number()), 0.0f};
        const float size[2]{static_cast<float>(samples), 1.0f};
        double low = INFINITY, high = -INFINITY;
        for (int k = 0; k < samples; ++k) {
            const double F = bd::plot_value(twin, view, size, k, 0, middles).F;
            if (std::isfinite(F)) low = std::min(low, F), high = std::max(high, F);
        }
        if (low <= high) out.emplace_back("range", Array{low, high});
    }
    if (const JsonValue* g = json->find("grid"); g && g->is_object()) {
        const int w = static_cast<int>(g->number_or("width", 0)), h = static_cast<int>(g->number_or("height", 0));
        float view[4]{}, size[2]{static_cast<float>(w), static_cast<float>(h)};
        if (const JsonValue* v = g->find("view"); v && v->is_array() && v->as_array().size() == 4)
            for (std::size_t k = 0; k < 4; ++k) view[k] = static_cast<float>(v->as_array()[k].as_number());
        std::vector<float> slots;
        if (const JsonValue* s = g->find("slots"); s && s->is_array())
            for (const auto& x : s->as_array()) slots.push_back(static_cast<float>(x.as_number()));
        auto pod = bd::lower(f);
        if (!pod) return refuse(pod.error().message);
        const auto twin = bd::narrowed(*pod);
        Array values;
        for (int py = 0; py < h; ++py)
            for (int px = 0; px < w; ++px) {
                const bd::PlotValue at = bd::plot_value(twin, view, size, px, py, slots);
                values.push_back(static_cast<double>(at.F));
            }
        out.emplace_back("values", values);
    }
    if (const JsonValue* c = json->find("curve"); c && c->is_object()) {
        const double x0 = c->number_or("x0", 0), x1 = c->number_or("x1", 1);
        const int n = std::clamp(static_cast<int>(c->number_or("n", 0)), 0, 8192);
        auto pod = bd::lower(f);
        if (!pod) return refuse(pod.error().message);
        // the sliders' values, given in the shader's order, put in the program's slots by name
        std::vector<double> slots(pod->symbols.size(), 0.0);
        if (const JsonValue* s = c->find("slots"); s && s->is_array())
            for (std::size_t k = 0; k < s->as_array().size() && k < shader->symbols.size(); ++k)
                if (const auto at = std::ranges::find(pod->symbols, shader->symbols[k]); at != pod->symbols.end())
                    slots[static_cast<std::size_t>(at - pod->symbols.begin())] = s->as_array()[k].as_number();
        Array values;
        for (int k = 0; k < n; ++k) {
            const double F = bd::interpret(*pod, x0 + (k + 0.5) / n * (x1 - x0), 0.0, slots);
            values.push_back(std::isfinite(F) ? JsonValue(F) : JsonValue(nullptr));
        }
        out.emplace_back("curve", values);
    }
    return out;
}

} // namespace rsc
