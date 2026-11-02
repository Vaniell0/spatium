#pragma once

// An equation the library solves in closed form, typed as text, with the steps it was solved by.
//
// `2x^2 - 3x - 5 = 0` is read into a field in one unknown; its coefficients are the field's
// values at x = 0, 1 and -1 (c = f(0), b = (f(1) - f(-1)) / 2, a = (f(1) + f(-1)) / 2 - c), and
// that it is a polynomial of degree at most two is checked, not assumed: the quadratic with those
// coefficients must give the field's values at five more points, or the equation is refused by
// name. No expression is rewritten. The roots come from `solve_quadratic_steps`, the solver
// itself, so the steps a page shows are the ones the answer was computed through; the last step
// checks the two roots by Vieta's formulas, a second route that says when the first lost digits.
//
// Asked for exactly (`exact`), the coefficients are computed again as fractions -- the unknown a
// slot of the same lowered program, run in `Rational` -- and the degree checked by equality, not
// to a tolerance. The same solver, run on fractions, gives the roots as fractions where the
// discriminant is a square (5/2 and -1), and where it is not they are written (p ± q√m)/d with m
// free of square factors: (3 ± √17)/2, (-1 ± 2i)/1 as -1 ± 2i. Nothing is rounded on that route;
// what does not fit in 64 bits is refused by name, and the double's route stands.
//
// This is the closed-form side of the calculator; a question with a parameter whose answer has
// to be found goes to `ask` and synthesis (query.hpp).

#include <spatium/algebra/polynomial.hpp>
#include <spatium/core/rational.hpp>
#include <spatium/io/field_pod.hpp>
#include <spatium/io/field_parse.hpp>
#include <spatium/io/steps_text.hpp>

#include <array>
#include <cmath>
#include <format>
#include <numeric>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace rsc {

// The roots of a quadratic with fractional coefficients, as they are: fractions, or
// (p ± q√m)/d with i where the discriminant is negative.
struct ExactRoots {
    std::vector<spatium::Rational> fractions;   // when every root is one
    spatium::Rational p, q;                     // otherwise p ± q√m, imaginary when `i`
    long long m = 1;
    bool i = false;
    std::string text, latex;                    // `(3 ± √17)/2`, `\frac{3 \pm \sqrt{17}}{2}`
};

struct ClosedForm {
    std::string unknown;
    std::optional<ExactRoots> exact;                    // when asked for and the coefficients are fractions
    std::string inexact;                                // why not, when asked for and there are none
    spatium::QuadraticSteps<double> solution;
    std::optional<spatium::VietaCheck<double>> check;   // of the roots, when there are two
    std::vector<spatium::io::build::Step> steps;        // the solution, and the check last
    std::size_t evaluations = 0;                        // of the equation's field, the route's cost
};

namespace closed_form_detail {

using spatium::Rational;

// n = k^2 m with m free of square factors, by trial division to 2^20 (past it a square factor is
// left inside m: the root is still exact, only less simplified).
inline std::pair<long long, long long> square_free(long long n) {
    long long k = 1, m = 1;
    for (long long p = 2; p <= (1 << 20) && p * p <= n; ++p) {
        int e = 0;
        while (n % p == 0) n /= p, ++e;
        for (; e >= 2; e -= 2) k *= p;
        if (e) m *= p;
    }
    const auto r = static_cast<long long>(std::sqrt(static_cast<double>(n)));
    for (long long s : {r - 1, r, r + 1})
        if (s > 1 && s * s == n) return {k * s, m};
    return {k, m * n};
}

inline std::string fraction_latex(Rational r) {
    if (r.den() == 1) return std::to_string(r.num());
    return std::format("{}\\frac{{{}}}{{{}}}", r.num() < 0 ? "-" : "", r.num() < 0 ? -r.num() : r.num(), r.den());
}

// The roots of a x^2 + b x + c = 0 exactly, through the library's own solver run on fractions.
inline spatium::Result<ExactRoots> exact_roots(Rational a, Rational b, Rational c) {
    using spatium::Error, spatium::ErrorCode;
    const auto s = spatium::solve_quadratic_steps(a, b, c);
    const auto too_big = [] { return std::unexpected(Error(ErrorCode::InvalidArgument, "the exact roots do not fit in 64 bits")); };
    ExactRoots out;
    bool exact = true;
    for (const auto& z : s.roots) exact = exact && z.re.ok() && z.im.ok() && z.im == Rational(0);
    if (exact) {
        for (const auto& z : s.roots) out.fractions.push_back(z.re);
        for (const auto& r : out.fractions) {
            out.text += (out.text.empty() ? "" : ", ") + r.text();
            out.latex += (out.latex.empty() ? "" : ",\\quad ") + fraction_latex(r);
        }
        return out;
    }
    // roots = -b/2a ± sqrt(D)/2a: sqrt(n/d) = sqrt(n d)/d, and n d = k^2 m
    const Rational D = s.discriminant, two_a = Rational(2) * a;
    out.i = D < Rational(0);
    const Rational size = out.i ? -D : D;
    const __int128 nd = __int128(size.num()) * size.den();
    if (!D.ok() || nd > INT64_MAX) return too_big();
    const auto [k, m] = square_free(static_cast<long long>(nd));
    out.m = m;
    out.p = -b / two_a;
    out.q = abs(Rational(k) / (Rational(size.den()) * two_a));
    if (!out.p.ok() || !out.q.ok()) return too_big();
    // over one denominator: (P ± Q√m)/d
    const long long d = std::lcm(out.p.den(), out.q.den());
    const __int128 P = __int128(out.p.num()) * (d / out.p.den()), Q = __int128(out.q.num()) * (d / out.q.den());
    if (P > INT64_MAX || P < -INT64_MAX || Q > INT64_MAX) return too_big();
    const std::string q_text = (Q == 1 && (out.i || m > 1) ? "" : std::to_string(static_cast<long long>(Q)));
    const std::string root_text = q_text + (out.i ? "i" : "") + (m > 1 ? "√" + std::to_string(m) : "");
    const std::string root_latex = q_text + (out.i ? "i" : "") + (m > 1 ? std::format("\\sqrt{{{}}}", m) : "");
    const std::string head_text = P != 0 ? std::to_string(static_cast<long long>(P)) + " ± " : "±";
    const std::string head_latex = P != 0 ? std::to_string(static_cast<long long>(P)) + " \\pm " : "\\pm ";
    out.text = d == 1 ? head_text + root_text : "(" + head_text + root_text + ")/" + std::to_string(d);
    out.latex = d == 1 ? head_latex + root_latex : std::format("\\frac{{{}}}{{{}}}", head_latex + root_latex, d);
    return out;
}

// The coefficients of f, a polynomial of degree at most two in its slot, as fractions: from its
// values at 0, 1 and -1, checked equal at five more points.
inline spatium::Result<std::array<Rational, 3>> exact_coefficients(const spatium::io::build::Field<double>& f, const std::string& x) {
    using spatium::Error, spatium::ErrorCode;
    namespace bd = spatium::io::build;
    const auto pod = bd::lower(f);
    if (!pod) return std::unexpected(pod.error());
    if (const auto op = bd::beyond<Rational>(*pod))
        return std::unexpected(Error(ErrorCode::InvalidArgument, std::format("the equation reads {}, which has no fraction", bd::op_name(*op))));
    const auto program = bd::converted<Rational>(*pod);
    const auto at = [&](Rational v) { return bd::interpret(program, Rational(0), Rational(0), {v}); };
    const Rational c = at(0), up = at(1), down = at(-1);
    const Rational a = (up + down) / Rational(2) - c, b = (up - down) / Rational(2);
    for (Rational v : {Rational(2), Rational(-2), Rational(1, 2), Rational(3), Rational(-7)}) {
        const Rational want = at(v);
        if (!want.ok() || want != (a * v + b) * v + c)
            return std::unexpected(Error(ErrorCode::InvalidArgument, std::format(
                "the equation is not a quadratic in {} with fractions for coefficients (at {} = {} it is {})", x, x, v.text(), want.text())));
    }
    if (!a.ok() || !b.ok() || !c.ok()) return std::unexpected(Error(ErrorCode::InvalidArgument, "the coefficients do not fit in 64 bits"));
    return std::array<Rational, 3>{a, b, c};
}

} // namespace closed_form_detail

// The equation f = 0 in its one unknown, solved in closed form. Given `exact` -- the same
// equation as written, its constants not folded into doubles -- its roots as they are too.
inline spatium::Result<ClosedForm> solve_closed_form(const spatium::io::build::Field<double>& f, const std::string& x,
                                                     const spatium::io::build::Field<double>* exact = nullptr) {
    using spatium::Error, spatium::ErrorCode;
    namespace bd = spatium::io::build;
    if (const auto names = f.free_symbols(); names.size() != 1 || names.front() != x)
        return std::unexpected(Error(ErrorCode::InvalidArgument,
            std::format("closed form: the equation must read exactly one unknown, it reads {}", names.size())));
    std::size_t evaluations = 0;
    const auto at = [&](double v) { return ++evaluations, f.bind(x, bd::Field<double>(v))(0.0, 0.0); };
    const double c = at(0), up = at(1), down = at(-1);
    const double a = (up + down) / 2 - c, b = (up - down) / 2;
    for (double v : {2.0, -2.0, 0.5, 3.0, -7.0}) {
        const double want = at(v), got = (a * v + b) * v + c;
        const double scale = std::max({std::abs(want), std::abs(a * v * v), std::abs(b * v), std::abs(c), 1.0});
        if (!(std::abs(want - got) <= 1e-12 * scale))
            return std::unexpected(Error(ErrorCode::InvalidArgument,
                std::format("closed form: the equation is not a polynomial of degree at most two in {} "
                            "(at {} = {} it is {}, a quadratic through x = 0, 1, -1 gives {})", x, x, v, want, got)));
    }
    ClosedForm out{x, std::nullopt, {}, spatium::solve_quadratic_steps(a, b, c), std::nullopt, {}, evaluations};
    out.check = spatium::check_vieta(out.solution);
    out.steps = bd::to_steps(out.solution, x);
    if (exact) {
        auto k = closed_form_detail::exact_coefficients(*exact, x);
        auto r = k ? closed_form_detail::exact_roots((*k)[0], (*k)[1], (*k)[2]) : std::unexpected(k.error());
        if (r) out.exact = std::move(*r);
        else out.inexact = r.error().message;
        // the roots' step says them as they are: one root, two fractions, or a pair (p ± q√m)/d
        if (out.exact && !out.steps.empty()) {
            using closed_form_detail::fraction_latex;
            const auto& e = *out.exact;
            const std::string name = bd::detail::name_latex(x);
            auto& last = out.steps.back();
            if (last.label == "root" && !e.fractions.empty())
                last.latex = name + " = " + fraction_latex(e.fractions[0]);
            else if (last.label == "roots" && e.fractions.size() == 2)
                last.latex = std::format("{0}_{{1}} = {1},\\quad {0}_{{2}} = {2}", name, fraction_latex(e.fractions[0]),
                                         fraction_latex(e.fractions[1]));
            else if (last.label == "roots")
                last.latex = name + "_{1,2} = " + e.latex;
        }
    }
    if (out.check) out.steps.push_back(bd::to_step(*out.check, x));
    return out;
}

inline spatium::Result<ClosedForm> solve_closed_form(std::string_view text) {
    using spatium::Error, spatium::ErrorCode;
    namespace bd = spatium::io::build;
    bd::FieldReader r(text, "", "");
    auto lhs = r.expression();
    if (!lhs) return std::unexpected(lhs.error());
    if (!r.take('=')) return r.fail("'='");
    auto rhs = r.expression();
    if (!rhs) return std::unexpected(rhs.error());
    if (!r.done()) return r.fail("the end");
    const auto f = std::move(*lhs) - *rhs;
    const auto names = f.free_symbols();
    return solve_closed_form(f, names.size() == 1 ? names.front() : std::string());
}

} // namespace rsc
