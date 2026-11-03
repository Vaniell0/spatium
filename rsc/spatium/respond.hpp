#pragma once

// A question typed as text, and the answer as the JSON a page reads.
//
// What a person types into the calculator is a shorthand for a document of one definition, and
// `respond` writes it out and evaluates it (evaluate.hpp): `2 + 2` is the value `answer = 2 + 2`,
// `2x^2 - 3x - 5 = 0` the root `x: 2x^2 - 3x - 5 = 0`, solved in closed form with its steps, and
// `solve x^2 - p = 0 for x in [0, 10], p in [0.5, 4]` the parameter and the root it asks for,
// synthesised with its certificate. A text that is none of these is refused with the reason,
// which is the library's own (a parse error says its column). The answer carries the document it
// was read as and every route that was tried, so it is a whole record of what was asked and how
// it was answered. The page composes nothing: what it shows of an answer is in this JSON, and the
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

inline bool has_word(std::string_view s, std::string_view w) {
    for (std::size_t at = s.find(w); at != std::string_view::npos; at = s.find(w, at + 1)) {
        const bool before = at == 0 || std::isspace(static_cast<unsigned char>(s[at - 1]));
        const bool after = at + w.size() == s.size() || std::isspace(static_cast<unsigned char>(s[at + w.size()]));
        if (before && after) return true;
    }
    return false;
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
            out.emplace_back("latex", *bd::to_latex(e.written) + " = " + number_latex(c, e.values.front()));
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
        case Route::Bracket: out.emplace_back("value", e.values.front()); break;
        case Route::Integral:
            out.emplace_back("value", e.values.front());
            out.emplace_back("error_below", decade_above(*e.error));
            break;
        case Route::Synthesis:
            for (auto& [k, v] : answer_json(*e.answer, e.parameter, written)) out.emplace_back(std::move(k), std::move(v));
            break;
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
    if (has_word(q, "for")) {
        auto asked = question_document(q);
        if (!asked) return refuse("synthesis", asked.error().message);
        return with_above(std::move(asked->document), std::move(asked->wanted), "synthesis");
    }
    // an equation, with or without the word solve in front
    std::string_view equation = q;
    if (equation.starts_with("solve") && (equation.size() == 5 || std::isspace(static_cast<unsigned char>(equation[5]))))
        equation = trimmed(equation.substr(5));
    // no '=' but at the end (2 + 2 =): an expression to evaluate, as a calculator does
    if (const auto eq = equation.find('='); eq == std::string_view::npos || trimmed(equation.substr(eq + 1)).empty()) {
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
    // the root of the equation in its one unknown, with no bracket: solved in closed form
    bd::FieldReader r(equation, "", "");
    bd::Statement root;
    root.kind = bd::StatementKind::Root;
    auto lhs = bd::detail::written(r, equation);
    if (!lhs) return refuse("closed form", lhs.error().message);
    if (!r.take('=')) return refuse("closed form", r.fail("'='").error().message);
    auto rhs = bd::detail::written(r, equation);
    if (!rhs) return refuse("closed form", rhs.error().message);
    if (!r.done()) return refuse("closed form", r.fail("the end").error().message);
    // a new name on the left that the right does not read is a definition, f = a^2 + 1, not an equation in f
    const std::string left(trimmed(lhs->text));
    if (const auto right_reads = rhs->field.free_symbols(); !left.empty() && !above.find(left) &&
        lhs->field.free_symbols() == std::vector<std::string>{left} && !std::ranges::contains(right_reads, left)) {
        bd::Statement value;
        value.name = left, value.value = std::move(*rhs);
        return with_above({{std::move(value)}}, left, "evaluate");
    }
    const auto names = unknown(lhs->field - rhs->field);
    if (names.size() != 1) return refuse("closed form", solve_closed_form(lhs->field - rhs->field, "").error().message);
    root.name = names.front(), root.value = std::move(*lhs), root.rhs = std::move(*rhs);
    return with_above({{std::move(root)}}, names.front(), "closed form");
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
    if (how == "closed form" && wanted->kind == bd::StatementKind::Root) {
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

} // namespace respond_detail

// What a person types, answered. One line is a shorthand (`2 + 2`, `2x^2 - 3x - 5 = 0`, `solve
// ... for ...`), or else a line of a document (`x: x^2 - 2 = 0 on [0, 2]`). Several lines are a
// document and its question: the lines above the last define names, and the last is the
// shorthand read with them (`p in [0.5, 4]` above `solve`-less `x^2 - p = 0` asks for x as a
// function of p), or else the document's last definition.
inline spatium::io::JsonValue respond(std::string_view text, double tolerance = 1e-8) {
    using namespace respond_detail;
    namespace bd = spatium::io::build;
    const std::string_view all = trimmed(text);
    Object out{{"question", std::string(all)}};
    const auto answered = [&](const AskedDocument& asked) {
        for (auto& [k, v] : evaluation_json(evaluate(asked.document, asked.wanted, {}, Engine{.tolerance = tolerance})))
            out.emplace_back(std::move(k), std::move(v));
        out.emplace_back("document", bd::to_json(asked.document));
        if (auto g = graph_of(asked, out)) out.emplace_back("graph", std::move(*g));
        return JsonValue(std::move(out));
    };
    bd::Document above;
    std::string_view q = all;
    if (const auto nl = all.rfind('\n'); nl != std::string_view::npos) {
        auto doc = bd::read_document(all.substr(0, nl));
        if (!doc) {
            out.emplace_back("route", "document");
            out.emplace_back("refused", doc.error().message);
            return JsonValue(std::move(out));
        }
        above = std::move(*doc);
        q = trimmed(all.substr(nl + 1));
    }
    auto s = shorthand(q, above);
    if (s.asked) return answered(*s.asked);
    // a shorthand that does not read, but a document that does
    if (auto whole = bd::read_document(all); whole && !whole->statements.empty()) {
        std::string wanted = whole->statements.back().name;
        return answered({std::move(*whole), std::move(wanted)});
    }
    // refused under the route the shorthand asked for
    out.emplace_back("route", s.asked.error().first);
    out.emplace_back("refused", std::move(s.asked.error().second));
    return JsonValue(std::move(out));
}

// The rows of a LaTeX editor, each a formula as MathJSON, answered as the lines they read as
// (`mathjson_row_text`): `{"rows": [...], "tolerance": 1e-8}`. The answer is the one `respond`
// gives the same lines typed, and a row that does not read is refused with its number, `row`,
// counted from 1.
inline spatium::io::JsonValue respond_rows(std::string_view request) {
    using namespace respond_detail;
    namespace bd = spatium::io::build;
    const auto json = spatium::io::parse(request);
    const JsonValue* rows = json && json->is_object() ? json->find("rows") : nullptr;
    if (!rows || !rows->is_array() || rows->as_array().empty())
        return JsonValue(Object{{"question", ""}, {"route", "document"}, {"refused", "the request has no rows"}});
    std::string text;
    for (std::size_t k = 0; k < rows->as_array().size(); ++k) {
        auto line = bd::mathjson_row_text(rows->as_array()[k]);
        if (!line)
            return JsonValue(Object{{"question", text}, {"route", "document"}, {"row", static_cast<double>(k + 1)},
                                    {"refused", std::format("row {}: {}", k + 1, line.error().message)}});
        text += (k ? "\n" : "") + *line;
    }
    return respond(text, json->number_or("tolerance", 1e-8));
}

// A document given to the door as one JSON object: `document` its text or its JSON, `name` the
// definition wanted, `chosen` the parameters given a value ({"p": 2}), and `engine`, how it is
// computed ({"arithmetic": "rational", "strict": false, "tolerance": 1e-8}, each optional; a
// `tolerance` beside them is the synthesis's too). The answer is `name` and the evaluation, or
// what could not be read.
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
    Engine engine{.tolerance = json->number_or("tolerance", 1e-8)};
    if (const JsonValue* e = json->find("engine"); e && e->is_object()) {
        const std::string a = e->string_or("arithmetic", "double");
        const auto known = std::ranges::find_if(arithmetics, [&](Arithmetic x) { return a == arithmetic_name(x); });
        if (known == std::end(arithmetics)) {
            std::string names;
            for (Arithmetic x : arithmetics) names += std::string(names.empty() ? "" : ", ") + arithmetic_name(x);
            return refuse(std::format("the engine has no arithmetic '{}'; it has {}", a, names));
        }
        engine.arithmetic = *known;
        if (const JsonValue* strict = e->find("strict"); strict && strict->is_bool()) engine.strict = strict->as_bool();
        engine.tolerance = e->number_or("tolerance", engine.tolerance);
    }
    const std::string name = json->string_or("name", "");
    Object out{{"name", name}};
    for (auto& [k, v] : evaluation_json(evaluate(*doc, name, chosen, engine)))
        out.emplace_back(std::move(k), std::move(v));
    return out;
}

// A value of a document drawn as the graph of a function of one of its parameters, as one JSON
// object: `document` (text or JSON), `name`, and `along`, the parameter on the x axis. The other
// parameters it reads stay symbols, read from the shader's slots, and come back with their
// intervals so a page can give each a slider. With `grid` -- {"width", "height", "view": [x0, y0,
// x1, y1], "slots": [...]} -- the answer also holds F at every pixel's centre as the shader's
// check mode computes it, run in fp32 on the host: what a browser holds its WebGL2 frame against.
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
    return out;
}

} // namespace rsc
