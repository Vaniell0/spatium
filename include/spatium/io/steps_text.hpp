#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/polynomial.hpp>
#  include <spatium/io/field_text.hpp>
#  include <cmath>
#  include <format>
#  include <string>
#  include <string_view>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::io::build {

// ── A solution as a chain of steps ────────────────────────────────
//
// What a page shows under an answer found in closed form: each step a label and a line of LaTeX,
// the formula first and then the same formula with the numbers in it, so a reader can follow
// from the equation to the roots. The steps are read off the solver's own record
// (`solve_quadratic_steps`), never recomputed here: what is written is what was computed.

struct Step {
    std::string label;   // what the step is, in words
    std::string latex;   // the step itself
};

namespace detail {

inline std::string num(double v) { return number_text(v, TextStyle::Latex); }

// a number as a factor or a term after an operator: in parentheses when negative
inline std::string grouped_num(double v) { return std::signbit(v) ? "\\left(" + num(v) + "\\right)" : num(v); }

// a x^2 + b x + c with the zero terms left out, a coefficient of 1 unwritten, and a negative one
// written as a subtraction
inline std::string polynomial_latex(double a, double b, double c, std::string_view x) {
    std::string out;
    const auto term = [&](double k, std::string_view power) {
        if (k == 0) return;
        const double m = std::abs(k);
        const std::string coef = (m == 1 && !power.empty()) ? "" : num(m);
        if (out.empty()) out = std::signbit(k) ? "-" : "";
        else out += std::signbit(k) ? " - " : " + ";
        out += coef + std::string(power);
    };
    const std::string name = name_latex(x);
    term(a, name + "^{2}");
    term(b, name);
    term(c, "");
    return out.empty() ? "0" : out;
}

inline std::string complex_latex(const Complex<double>& z) {
    if (z.im == 0) return num(z.re);
    const std::string im = std::abs(z.im) == 1 ? "i" : num(std::abs(z.im)) + "\\,i";
    if (z.re == 0) return (std::signbit(z.im) ? "-" : "") + im;
    return num(z.re) + (std::signbit(z.im) ? " - " : " + ") + im;
}

} // namespace detail

inline std::vector<Step> to_steps(const QuadraticSteps<double>& s, std::string_view unknown = "x") {
    using detail::num, detail::grouped_num;
    const std::string x = detail::name_latex(unknown);
    std::vector<Step> out;
    out.push_back({"equation", detail::polynomial_latex(s.a, s.b, s.c, unknown) + " = 0"});
    out.push_back({"coefficients", std::format("a = {},\\quad b = {},\\quad c = {}", num(s.a), num(s.b), num(s.c))});
    switch (s.kind) {
        case QuadraticCase::Identity:
            out.push_back({"a = b = c = 0", "0 = 0 \\text{ for every } " + x});
            return out;
        case QuadraticCase::None:
            out.push_back({"a = b = 0, c \\ne 0", std::format("{} = 0 \\text{{ holds for no }} {}", num(s.c), x)});
            return out;
        case QuadraticCase::Linear:
            out.push_back({"a = 0: the equation is linear",
                           std::format("{} = -\\frac{{c}}{{b}} = -\\frac{{{}}}{{{}}} = {}", x, num(s.c), num(s.b),
                                       num(s.roots[0].re))});
            return out;
        default: break;
    }
    out.push_back({"discriminant", std::format("D = b^{{2}} - 4ac = {}^{{2}} - 4 \\cdot {} \\cdot {} = {}", grouped_num(s.b),
                                               grouped_num(s.a), grouped_num(s.c), num(s.discriminant))});
    const std::string minus_b = num(-s.b), two_a = num(2 * s.a);
    if (s.kind == QuadraticCase::ComplexPair) {
        out.push_back({"D < 0: two complex conjugate roots",
                       std::format("{}_{{1,2}} = \\frac{{-b \\pm i\\sqrt{{-D}}}}{{2a}} = \\frac{{{} \\pm i\\sqrt{{{}}}}}{{{}}}", x,
                                   minus_b, num(-s.discriminant), two_a)});
    } else if (s.kind == QuadraticCase::Double) {
        out.push_back({"D = 0: one root, counted twice",
                       std::format("{} = -\\frac{{b}}{{2a}} = \\frac{{{}}}{{{}}}", x, minus_b, two_a)});
        out.push_back({"root", std::format("{} = {}", x, num(s.roots[0].re))});
        return out;
    } else {
        out.push_back({"D > 0: two real roots",
                       std::format("{}_{{1,2}} = \\frac{{-b \\pm \\sqrt{{D}}}}{{2a}} = \\frac{{{} \\pm \\sqrt{{{}}}}}{{{}}}", x,
                                   minus_b, num(s.discriminant), two_a)});
    }
    out.push_back({"roots", std::format("{}_{{1}} = {},\\quad {}_{{2}} = {}", x, detail::complex_latex(s.roots[0]), x,
                                        detail::complex_latex(s.roots[1]))});
    return out;
}

} // namespace spatium::io::build
