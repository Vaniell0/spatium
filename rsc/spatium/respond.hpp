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
#include <spatium/io/field_parse.hpp>
#include <spatium/io/field_text.hpp>
#include <spatium/io/json.hpp>

#include <cctype>
#include <cmath>
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
        case Route::Evaluate:
            out.emplace_back("value", e.values.front());
            out.emplace_back("latex", *bd::to_latex(e.written) + " = " + bd::detail::number_text(e.values.front(), bd::detail::TextStyle::Latex));
            break;
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
        Object o{{"route", route_name(c.route)}, {"cost", c.cost}, {"taken", e.taken == i}};
        if (!c.refused.empty()) o.emplace_back("refused", c.refused);
        if (!c.values.empty()) {
            Array values;
            for (double v : c.values) values.push_back(v);
            o.emplace_back("values", values);
        }
        candidates.push_back(std::move(o));
    }
    out.emplace_back("candidates", candidates);
    return out;
}

inline spatium::io::JsonValue respond(std::string_view text, double tolerance = 1e-8) {
    using namespace respond_detail;
    namespace bd = spatium::io::build;
    const std::string_view q = trimmed(text);
    Object out{{"question", std::string(q)}};
    // a shorthand that does not read as a document is refused under the route it asked for
    const auto refuse = [&](const char* route, std::string why) {
        out.emplace_back("route", route);
        out.emplace_back("refused", std::move(why));
        return JsonValue(std::move(out));
    };
    const auto answered = [&](const AskedDocument& asked) {
        for (auto& [k, v] : evaluation_json(evaluate(asked.document, asked.wanted, {}, tolerance)))
            out.emplace_back(std::move(k), std::move(v));
        out.emplace_back("document", bd::to_json(asked.document));
        return JsonValue(std::move(out));
    };
    if (has_word(q, "for")) {
        auto asked = question_document(q);
        if (!asked) return refuse("synthesis", asked.error().message);
        return answered(*asked);
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
        if (const auto open = f->free_symbols(); !open.empty())
            return refuse("evaluate", "evaluate: the expression reads '" + open.front() +
                                          "', which has no value; to find it, write an equation (" + std::string(expression) + " = 0)");
        bd::Statement value;
        value.name = "answer", value.value = {*f, std::string(expression)};
        return answered({{{std::move(value)}}, "answer"});
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
    const auto names = (lhs->field - rhs->field).free_symbols();
    if (names.size() != 1) return refuse("closed form", solve_closed_form(lhs->field - rhs->field, "").error().message);
    root.name = names.front(), root.value = std::move(*lhs), root.rhs = std::move(*rhs);
    AskedDocument asked{{{std::move(root)}}, names.front()};
    if (auto ok = bd::check(asked.document); !ok) return refuse("closed form", ok.error().message);
    return answered(asked);
}

// A document given to the door as one JSON object: `document` its text or its JSON, `name` the
// definition wanted, `chosen` the parameters given a value ({"p": 2}), `tolerance` the
// synthesis's. The answer is `name` and the evaluation, or what could not be read.
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
    const std::string name = json->string_or("name", "");
    Object out{{"name", name}};
    for (auto& [k, v] : evaluation_json(evaluate(*doc, name, chosen, json->number_or("tolerance", 1e-8))))
        out.emplace_back(std::move(k), std::move(v));
    return out;
}

} // namespace rsc
