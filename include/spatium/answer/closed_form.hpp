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

#include <spatium/algebra/linear_solve.hpp>
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

namespace spatium::answer {

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
    std::vector<spatium::io::build::Step> steps;        // the solution; `check` is a step of its own (to_step)
    std::size_t evaluations = 0;                        // of the equation's field, the route's cost
};

// A system of equations, each linear in the unknowns, solved: as fractions where every
// coefficient is one, else in doubles.
struct LinearSystem {
    std::vector<std::string> unknowns;
    std::vector<double> values;                         // the solution, in the unknowns' order
    std::vector<spatium::Rational> exact;               // the same as fractions, when the coefficients are
    std::vector<spatium::io::build::Step> steps;        // the system as a matrix, then its solution
    std::vector<double> residuals;                      // each equation at the solution, in doubles: the check
    std::size_t evaluations = 0;                        // of the equations' fields, the route's cost
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

// A x = b for A an n x n matrix of rows, by the library's elimination at the size it has
// (solve_direct, n from 1 to 4); in fractions it is exact, and a singular A is refused.
template<class T, std::size_t N = 1>
spatium::Result<std::vector<T>> solved(const std::vector<std::vector<T>>& A, const std::vector<T>& b) {
    if constexpr (N > 4) {
        return std::unexpected(spatium::Error(spatium::ErrorCode::InvalidArgument, "a system of more than four equations"));
    } else {
        if (A.size() != N) return solved<T, N + 1>(A, b);
        spatium::Matrix<T, N, N> m{};
        spatium::Vec<T, N> v{};
        for (std::size_t i = 0; i < N; ++i) {
            v[i] = b[i];
            for (std::size_t j = 0; j < N; ++j) m(i, j) = A[i][j];
        }
        auto x = spatium::solve_direct(m, v);
        if (!x) return std::unexpected(spatium::Error(spatium::ErrorCode::SingularMatrix,
                                                      "the system has no single solution: its equations are not independent"));
        return std::vector<T>(x->data.begin(), x->data.end());
    }
}

// The coefficients of each equation f = 0 in the unknowns, a row of A and its b = -f(0): from
// its values at 0 and at each unit vector, checked at twice each and at each pair together (a
// square or a product of two unknowns shows there) and at one point more. In T: Rational checks
// by equality, a double to a tolerance.
template<class T>
spatium::Result<std::pair<std::vector<std::vector<T>>, std::vector<T>>> affine(const std::vector<spatium::io::build::Field<double>>& equations,
                                                                                const std::vector<std::string>& unknowns, std::size_t& evaluations) {
    using spatium::Error, spatium::ErrorCode;
    namespace bd = spatium::io::build;
    const std::size_t n = unknowns.size();
    std::vector<std::vector<T>> A;
    std::vector<T> b;
    for (std::size_t k = 0; k < equations.size(); ++k) {
        const auto pod = bd::lower(equations[k]);
        if (!pod) return std::unexpected(pod.error());
        if constexpr (std::is_same_v<T, Rational>)
            if (const auto op = bd::beyond<Rational>(*pod))
                return std::unexpected(Error(ErrorCode::InvalidArgument, std::format("equation {} reads {}, which has no fraction", k + 1, bd::op_name(*op))));
        const auto program = bd::converted<T>(*pod);
        const auto at = [&](const std::vector<T>& x) {
            std::vector<T> slots(program.symbols.size(), T(0));
            for (std::size_t j = 0; j < n; ++j)
                if (const auto s = std::ranges::find(program.symbols, unknowns[j]); s != program.symbols.end()) slots[static_cast<std::size_t>(s - program.symbols.begin())] = x[j];
            return ++evaluations, bd::interpret(program, T(0), T(0), slots);
        };
        const auto unit = [&](std::size_t j, T by) { std::vector<T> x(n, T(0)); x[j] = by; return x; };
        const T c = at(std::vector<T>(n, T(0)));
        std::vector<T> row;
        for (std::size_t j = 0; j < n; ++j) row.push_back(at(unit(j, T(1))) - c);
        std::vector<std::vector<T>> checks;
        for (std::size_t i = 0; i < n; ++i) {
            checks.push_back(unit(i, T(2)));
            for (std::size_t j = i + 1; j < n; ++j) checks.push_back(unit(i, T(1))), checks.back()[j] = T(1);
        }
        checks.emplace_back();
        for (std::size_t j = 0; j < n; ++j) checks.back().push_back(T(static_cast<int>(3 * j + 2)) / T(7) - T(1));
        for (const auto& x : checks) {
            T want = c, got = at(x);
            for (std::size_t j = 0; j < n; ++j) want += row[j] * x[j];
            bool same;
            if constexpr (std::is_same_v<T, Rational>) same = got.ok() && got == want;
            else same = std::abs(got - want) <= 1e-12 * std::max({std::abs(got), std::abs(want), 1.0});
            if (!same) return std::unexpected(Error(ErrorCode::InvalidArgument, std::format("equation {} is not linear in the unknowns", k + 1)));
        }
        A.push_back(std::move(row)), b.push_back(T(0) - c);
    }
    return std::pair{std::move(A), std::move(b)};
}

} // namespace closed_form_detail

// A system of equations f_k = 0, as many as unknowns, each linear in them, solved by elimination:
// in fractions where every coefficient is one, so x + 2y = 5, 3x - y = 1 is x = 1, y = 2 exactly,
// else in doubles. Given `exact`, the equations as written, their constants not folded.
inline spatium::Result<LinearSystem> solve_linear_system(const std::vector<spatium::io::build::Field<double>>& equations,
                                                         const std::vector<std::string>& unknowns,
                                                         const std::vector<spatium::io::build::Field<double>>* exact = nullptr) {
    using closed_form_detail::fraction_latex;
    namespace bd = spatium::io::build;
    if (equations.size() != unknowns.size())
        return std::unexpected(spatium::Error(spatium::ErrorCode::InvalidArgument,
            std::format("linear system: {} equations in {} unknowns", equations.size(), unknowns.size())));
    LinearSystem out{unknowns, {}, {}, {}, {}, 0};
    // the check, by a second route: each equation, as a field, at the solution
    const auto checked = [&] {
        for (const auto& f : equations) {
            auto at = f;
            for (std::size_t j = 0; j < unknowns.size(); ++j) at = at.bind(unknowns[j], bd::Field<double>(out.values[j]));
            out.residuals.push_back(at(0.0, 0.0));
        }
        return out;
    };
    const auto matrix = [&](const auto& A, const auto& b, const auto& tex) {
        std::string a, x, r;
        for (std::size_t i = 0; i < A.size(); ++i) {
            for (std::size_t j = 0; j < A[i].size(); ++j) a += (j ? " & " : "") + tex(A[i][j]);
            const std::string end = i + 1 < A.size() ? " \\\\ " : "";
            a += end, x += bd::detail::name_latex(unknowns[i]) + end, r += tex(b[i]) + end;
        }
        return std::format("\\begin{{pmatrix}} {} \\end{{pmatrix}} \\begin{{pmatrix}} {} \\end{{pmatrix}} = \\begin{{pmatrix}} {} \\end{{pmatrix}}", a, x, r);
    };
    const auto solution = [&](const auto& x, const auto& tex) {
        std::string s;
        for (std::size_t j = 0; j < x.size(); ++j) s += (j ? ",\\quad " : "") + bd::detail::name_latex(unknowns[j]) + " = " + tex(x[j]);
        return s;
    };
    // in fractions first, and in doubles where they are not
    if (auto k = closed_form_detail::affine<spatium::Rational>(exact ? *exact : equations, unknowns, out.evaluations)) {
        auto x = closed_form_detail::solved(k->first, k->second);
        if (!x) return std::unexpected(x.error());
        if (std::ranges::all_of(*x, [](const spatium::Rational& r) { return r.ok(); })) {
            out.exact = *x;
            for (const auto& r : *x) out.values.push_back(r.to_double());
            out.steps = {{"coefficients", matrix(k->first, k->second, fraction_latex)}, {"solution", solution(*x, fraction_latex)}};
            return checked();
        }
    }
    auto k = closed_form_detail::affine<double>(equations, unknowns, out.evaluations);
    if (!k) return std::unexpected(k.error());
    auto x = closed_form_detail::solved(k->first, k->second);
    if (!x) return std::unexpected(x.error());
    const auto number = [](double v) { return bd::detail::number_text(v, bd::detail::TextStyle::Plain); };
    out.values = *x;
    out.steps = {{"coefficients", matrix(k->first, k->second, number)}, {"solution", solution(*x, number)}};
    return checked();
}

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
        // the roots' step says them as they are, then as decimals: one root, two fractions, or a pair
        // (p ± q√m)/d -- 5/2 = 2.5, and (3 ± √17)/2 ≈ 3.56, -0.56
        if (out.exact && !out.steps.empty()) {
            using closed_form_detail::fraction_latex;
            const auto& e = *out.exact;
            const std::string name = bd::detail::name_latex(x);
            const auto shown = [](spatium::Rational r) {
                if (r.den() == 1) return fraction_latex(r);
                const auto d = r.decimal();
                return fraction_latex(r) + (d ? " = " + *d : std::format(" \\approx {:.12g}", r.to_double()));
            };
            auto& last = out.steps.back();
            if (last.label == "root" && !e.fractions.empty())
                last.latex = name + " = " + shown(e.fractions[0]);
            else if (last.label == "roots" && e.fractions.size() == 2)
                last.latex = std::format("{0}_{{1}} = {1},\\quad {0}_{{2}} = {2}", name, shown(e.fractions[0]), shown(e.fractions[1]));
            else if (last.label == "roots" && !e.i) {
                // twelve digits, and a pair about zero as one: ±√2 ≈ ±1.41421356237
                const auto twelve = [](double v) { return std::format("{:.12g}", v); };
                const double r0 = out.solution.roots[0].re, r1 = out.solution.roots[1].re;
                const std::string near = e.p == spatium::Rational(0) ? "\\pm " + twelve(std::abs(r0)) : twelve(r0) + ",\\ " + twelve(r1);
                last.latex = std::format("{}_{{1,2}} = {} \\approx {}", name, e.latex, near);
            }
            else if (last.label == "roots")
                last.latex = name + "_{1,2} = " + e.latex;
        }
    }
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

} // namespace spatium::answer
