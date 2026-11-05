#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/polynomial.hpp>
#  include <spatium/core/interval.hpp>
#  include <spatium/core/rational.hpp>
#  include <spatium/io/field_parse.hpp>
#  include <spatium/io/field_pod.hpp>
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
// p / q on a page with its sign worked out: -\frac{4}{2}, never \frac{-4}{2} or -\frac{-4}{-2}
inline std::string signed_fraction(double p, double q) {
    return std::format("{}\\frac{{{}}}{{{}}}", std::signbit(p) != std::signbit(q) ? "-" : "", num(std::abs(p)), num(std::abs(q)));
}

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
            out.push_back({"a = b = c = 0: every value is a root", "0 = 0 \\text{ for every } " + x});
            return out;
        case QuadraticCase::None:
            out.push_back({"a = b = 0 and c is not 0: no root", std::format("{} = 0 \\text{{ holds for no }} {}", num(s.c), x)});
            return out;
        case QuadraticCase::Linear:
            out.push_back({"a = 0: the equation is linear",
                           std::format("{} = -\\frac{{c}}{{b}} = {} = {}", x, detail::signed_fraction(-s.c, s.b),
                                       num(s.roots[0].re))});
            return out;
        default: break;
    }
    out.push_back({"discriminant", std::format("D = b^{{2}} - 4ac = {}^{{2}} - 4 \\cdot {} \\cdot {} = {}", grouped_num(s.b),
                                               grouped_num(s.a), grouped_num(s.c), num(s.discriminant))});
    // (-b ± √D) / 2a with a negative 2a is (b ± √D) / |2a|, the same pair; a zero b is not written
    const double top_b = s.a < 0 ? s.b : -s.b;
    const std::string head = top_b == 0 ? "" : num(top_b) + " ", two_a = num(std::abs(2 * s.a));
    if (s.kind == QuadraticCase::ComplexPair) {
        out.push_back({"D < 0: two complex conjugate roots",
                       std::format("{}_{{1,2}} = \\frac{{-b \\pm i\\sqrt{{-D}}}}{{2a}} = \\frac{{{}\\pm i\\sqrt{{{}}}}}{{{}}}", x,
                                   head, num(-s.discriminant), two_a)});
    } else if (s.kind == QuadraticCase::Double) {
        out.push_back({"D = 0: one root, counted twice",
                       std::format("{} = -\\frac{{b}}{{2a}} = {}", x, detail::signed_fraction(-s.b, 2 * s.a))});
        out.push_back({"root", std::format("{} = {}", x, num(s.roots[0].re))});
        return out;
    } else {
        out.push_back({"D > 0: two real roots",
                       std::format("{}_{{1,2}} = \\frac{{-b \\pm \\sqrt{{D}}}}{{2a}} = \\frac{{{}\\pm \\sqrt{{{}}}}}{{{}}}", x,
                                   head, num(s.discriminant), two_a)});
    }
    out.push_back({"roots", std::format("{}_{{1}} = {},\\quad {}_{{2}} = {}", x, detail::complex_latex(s.roots[0]), x,
                                        detail::complex_latex(s.roots[1]))});
    return out;
}

// The check of the roots by Vieta's formulas as the last step: the sum and the product as
// computed, against -b/a and c/a; a check that does not hold says by how much.
inline Step to_step(const VietaCheck<double>& v, std::string_view unknown = "x") {
    const std::string x = detail::name_latex(unknown);
    const bool holds = v.error <= 1e-12;
    const std::string rel = holds ? "=" : "\\ne";
    // where it holds, a value off by the doubles' last bits is written as the one it is near, with
    // \approx: x_1 + x_2 \approx 0, not -2.220446049250313 \cdot 10^{-16}
    const auto shown = [&](const Complex<double>& got, double want) {
        const double scale = std::max({std::abs(want), std::abs(got.re), 1.0});
        if (!holds || got.im != 0 || got.re == want || std::abs(got.re - want) > 1e-12 * scale)
            return "= " + detail::complex_latex(got);
        return "\\approx " + detail::num(want);
    };
    std::string latex = std::format("{0}_{{1}} + {0}_{{2}} {1} {2} -\\frac{{b}}{{a}} = {3},\\quad {0}_{{1}} {0}_{{2}} {4} {2} "
                                    "\\frac{{c}}{{a}} = {5}",
                                    x, shown(v.sum, v.want_sum), rel, detail::num(v.want_sum),
                                    shown(v.product, v.want_product), detail::num(v.want_product));
    if (!holds) latex += std::format("\\quad\\text{{(relative error {})}}", detail::num(v.error));
    return {holds ? "check: Vieta's formulas hold" : "check: Vieta's formulas do not hold", latex};
}

// ── A computation as its steps ────────────────────────────────────
//
// What a value was computed through, a step an operation: the operation written on the values of
// what it reads, and the value it gave -- 3 + 4 = 7, then 7^{2} = 49 -- read off the lowered
// program as it ran, never worked out again from the answer. The values are fractions where every
// operation is one on fractions (`Rational`), and otherwise the digits an interval of the same
// program proves (core/interval.hpp), written with \approx. A leaf -- a number, a name -- is no
// step. The label is the operation's, from the language's own table (`functions`), the one the
// reference shows.

namespace detail {

inline std::string fraction_latex(const Rational& r) {
    if (r.den() == 1) return std::to_string(r.num());
    return std::format("{}\\frac{{{}}}{{{}}}", r.num() < 0 ? "-" : "", r.num() < 0 ? -r.num() : r.num(), r.den());
}

// what an operation is called: a function by the head of its caption in the language's table
inline std::string step_label(Op o) {
    switch (o) {
        case Op::Add: return "sum";
        case Op::Sub: return "difference";
        case Op::Mul: return "product";
        case Op::Div: return "quotient";
        case Op::Pow: return "power";
        default: break;
    }
    std::string name = op_name(o);
    for (auto& c : name) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    for (const auto& f : functions)
        if (f.name == name) return std::string(f.caption.substr(0, f.caption.find(',')));
    return name;
}

inline bool is_leaf(Op o) {
    switch (o) {
        case Op::Const: case Op::U: case Op::V: case Op::Symbol: case Op::Coord: case Op::Origin: case Op::Point:
        case Op::Id: case Op::Hash: case Op::Opaque:
            return true;
        default:
            return false;
    }
}

} // namespace detail

inline Result<std::vector<Step>> evaluation_steps(const Field<double>& whole) {
    using namespace detail;
    using I = Interval<double>;
    const Field<double> f = whole.reachable_from_root();
    const auto pod = lower(f);
    if (!pod) return std::unexpected(pod.error());
    if (const auto op = beyond<I>(*pod))
        return std::unexpected(Error(ErrorCode::InvalidArgument, std::format("an interval has no {} in it", op_name(*op))));
    const auto count = static_cast<std::uint32_t>(pod->ops.size());
    // each operation as written, and its value: a fraction where it is one, else proved digits
    std::vector<Spelled> written(count), value(count);
    for (std::uint32_t i = 0; i < count; ++i) written[i] = spell_op(f.op(i), written, "u", "v", TextStyle::Latex);
    std::vector<Rational> fraction(count);
    const bool fractions = !beyond<Rational>(*pod);
    if (fractions) {
        const auto in = converted<Rational>(*pod);
        interpret(in.ops.data(), count, in.tables.data(), in.points.data(), fraction.data(),
                  FieldInputs<Rational>{Rational(0), Rational(0)});
    }
    const auto in = converted<I>(*pod);
    std::vector<I> enclosed(count);
    interpret(in.ops.data(), count, in.tables.data(), in.points.data(), enclosed.data(), FieldInputs<I>{I{0}, I{0}});
    std::vector<bool> approximate(count, false);
    for (std::uint32_t i = 0; i < count; ++i) {
        if (is_leaf(f.op(i).op)) {
            value[i] = written[i];   // a number or a name as it was written
            continue;
        }
        if (fractions && fraction[i].ok()) {
            const Rational& r = fraction[i];
            value[i] = Spelled{fraction_latex(r), r.num() < 0 ? Binding::Sum : Binding::Atom, r.num() == 1 && r.den() == 1,
                               r.num() < 0, r.den() != 1};
        } else {
            // where no digit is proved, how far from zero the value may be: 10^20 + 1 - 10^20 is
            // 1, which a double's interval holds only within 2^14 of zero
            const Proved p = proved(enclosed[i]);
            const std::string text = p.digits > 0 ? p.text
                                                  : std::format("0 \\pm 10^{{{}}}", static_cast<int>(std::ceil(std::log10(p.bound))));
            const bool negative = text.front() == '-';
            value[i] = Spelled{text, negative ? Binding::Sum : Binding::Atom, p.exact && text == "1", negative, false};
            approximate[i] = !p.exact;
        }
        // a power written as a chain of products keeps its base, so the chain is one step
        value[i].base = written[i].base, value[i].power = written[i].power;
    }
    // a link of such a chain that the next one continues is no step of its own
    std::vector<bool> continued(count, false);
    for (std::uint32_t k = 0; k < count; ++k)
        if (const auto& n = f.op(k); n.op == Op::Mul && written[k].power > 1 && written[n.a].base == written[k].base)
            continued[n.a] = true;
    std::vector<Step> out;
    for (std::uint32_t i = 0; i < count; ++i) {
        const auto& n = f.op(i);
        if (is_leaf(n.op) || continued[i]) continue;
        const std::string on = spell_op(n, value, "u", "v", TextStyle::Latex).text;
        if (on == value[i].text) continue;   // a fraction as written, 1/4, is its own value
        // a chain of products of one base is a power, and named so
        const std::string label = written[i].power > 1 ? "power" : step_label(n.op);
        out.push_back({label, std::format("{} {} {}", on, approximate[i] ? "\\approx" : "=", value[i].text)});
    }
    // the answer, a fraction whose decimal ends, in that decimal too: \frac{391}{4} = 97.75
    if (const Rational& r = fraction[count - 1]; !out.empty() && fractions && r.ok() && r.den() != 1)
        if (const auto d = r.decimal()) out.back().latex += " = " + *d;
    return out;
}

} // namespace spatium::io::build
