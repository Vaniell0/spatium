#pragma once

// One door for a question typed as text, and the answer as the JSON a page reads.
//
// `respond("2 + 2")` evaluates an expression with no unknown, as a calculator does;
// `respond("2x^2 - 3x - 5 = 0")` solves the equation in closed form and returns its steps;
// `respond("solve x^2 - p = 0 for x in [0, 10], p in [0.5, 4]")` reads a question with a
// parameter and returns the synthesised answer with its certificate; a text that is neither is
// refused with the reason, which is the library's own (a parse error says its column). The page
// composes nothing: what it shows of an answer is in this JSON, and the gallery is written by the
// same function, so a card and a live answer cannot differ.
//
// An error is written as the power of ten above it (`error_below`): a bound a person reads, and
// one that does not move with the last digits of a libm.

#include <rsc/spatium/closed_form.hpp>
#include <rsc/spatium/query.hpp>
#include <rsc/spatium/synthesis.hpp>
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

inline spatium::io::JsonValue respond(std::string_view text, double tolerance = 1e-8) {
    using namespace respond_detail;
    const std::string_view q = trimmed(text);
    Object out{{"question", std::string(q)}};
    if (has_word(q, "for")) {
        auto asked = ask(q);
        out.emplace_back("route", "synthesis");
        if (!asked) {
            out.emplace_back("refused", asked.error().message);
            return out;
        }
        asked->tolerance = tolerance;
        for (auto& [k, v] : answer_json(answer(*asked), asked->parameter))
            out.emplace_back(std::move(k), std::move(v));
        return out;
    }
    // an equation, with or without the word solve in front
    std::string_view equation = q;
    if (equation.starts_with("solve") && (equation.size() == 5 || std::isspace(static_cast<unsigned char>(equation[5]))))
        equation = trimmed(equation.substr(5));
    // no '=' but at the end (2 + 2 =): an expression to evaluate, as a calculator does
    if (const auto eq = equation.find('='); eq == std::string_view::npos || trimmed(equation.substr(eq + 1)).empty()) {
        const std::string_view expression = trimmed(equation.substr(0, eq));
        out.emplace_back("route", "evaluate");
        const auto f = spatium::io::build::parse(expression, "", "");
        if (!f) {
            out.emplace_back("refused", f.error().message);
            return out;
        }
        if (const auto open = f->free_symbols(); !open.empty()) {
            out.emplace_back("refused", "evaluate: the expression reads '" + open.front() +
                                            "', which has no value; to find it, write an equation (" +
                                            std::string(expression) + " = 0)");
            return out;
        }
        const double value = (*f)(0.0, 0.0);
        out.emplace_back("value", value);
        out.emplace_back("latex", *spatium::io::build::to_latex(*f) + " = " +
                                      spatium::io::build::detail::number_text(value,
                                          spatium::io::build::detail::TextStyle::Latex));
        return out;
    }
    out.emplace_back("route", "closed form");
    const auto solved = solve_closed_form(equation);
    if (!solved) {
        out.emplace_back("refused", solved.error().message);
        return out;
    }
    Array steps;
    for (const auto& s : solved->steps) steps.push_back(Object{{"label", s.label}, {"latex", s.latex}});
    out.emplace_back("unknown", solved->unknown);
    out.emplace_back("steps", steps);
    if (solved->check)
        out.emplace_back("check", Object{{"by", "Vieta's formulas"},
                                         {"holds", solved->check->error <= 1e-12},
                                         {"error_below", decade_above(solved->check->error)}});
    return out;
}

} // namespace rsc
