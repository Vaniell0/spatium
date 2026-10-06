#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/calculus.hpp>
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <algorithm>
#  include <array>
#  include <cmath>
#  include <cstddef>
#  include <limits>
#  include <numbers>
#  include <utility>
#endif

SPATIUM_EXPORT namespace spatium {
inline namespace algebra {

// Quadrature over a domain that is a type.
//
// `integrate` (calculus.hpp) takes an interval as a finite pair, samples both
// ends, and so cannot do [0, inf) or an integrable singularity at an end: the
// NaN it returned there was a property of the model, not of the rule. The cure
// is a change of variable, and which change depends on what the domain is, so
// the domain is a type:
//
//   Finite<T>{a, b}       tanh-sinh: x = tanh(pi/2 sinh t). Its nodes crowd
//                         doubly exponentially toward the ends, which are never
//                         sampled, so an integrable singularity there costs
//                         nothing.
//   HalfLine<T>{from}     exp-sinh: x = exp(pi/2 sinh t), onto [from, inf) -- or
//                         (-inf, from] with `toward_infinity = false`.
//   WholeLine             sinh-sinh: x = sinh(pi/2 sinh t).
//
// All three are formulas in T, so the rule is right on float, double, Real50
// and Dual alike with a tolerance that follows the scalar (eps^(3/4) by
// default): no table of nodes to copy to fifty digits.
//
// The answer is an `IntegralResult`, as `integrate`'s. A doubly exponential
// rule gains digits by doubling per level, so its estimate is the difference of
// the last two levels and its status reads how that difference fell:
// `Converged` when it met the tolerance and fell by at least half, `Suspicious`
// when it met it without falling (two levels agreeing by coincidence),
// `DepthCapped` when the levels ran out, `Failed` when a bound of a finite
// domain, or every value, was not finite.
//
// One limit worth knowing: a node close to an end is stored as a point of the
// scalar, so near an end b != 0 it cannot be closer than an ulp of b, and an
// integrable singularity AT such an end (1/sqrt(b - x)) is truncated there --
// worse than the digits asked for. At an end that is 0, or for a weak
// singularity (log), nothing is lost.
//
// A second witness for a finite interval that samples completely differently:
// `gauss_kronrod`, G7K15 adaptive, up to long double. Two rules that agree are
// two witnesses; where they disagree by more than their estimates, one of them
// has been fooled, and which is a question for a third.

template<Scalar T> struct Finite { T a, b; };
template<Scalar T> struct HalfLine { T from; bool toward_infinity = true; };
struct WholeLine {};

namespace quadrature_detail {

template<Scalar T>
bool usable(const T& x) { return std::isfinite(primal_double(x)); }

// The abscissa and weight at the parameter t, and whether the node is still
// usable: false once it has overflowed, underflowed or reached an end.
template<Scalar T>
struct Node { T x{}; T w{}; bool ok = false; };

// The doubly exponential engine, for any node function. Levels halve the step
// h from 1; level 0 walks t = k h for every k (k = 0 once), later levels only
// the new odd multiples, and the sum is the previous level's halved plus them.
// Each side (t > 0, t < 0) is walked until its nodes stop being usable or its
// values stop being finite -- an overflow where the weights are already nothing.
template<Scalar T, typename F, typename NodeFn>
IntegralResult<T> doubly_exponential(F& f, NodeFn node, T scale, T tol, int max_levels) {
    using std::abs;
    const T nan = T(std::numeric_limits<double>::quiet_NaN());
    const T eps = T(machine_epsilon<T>());
    long evaluations = 0;
    T sum{0}, l1{0};                    // h * sum w f, h * sum w |f|
    T total{0}, last_total{0};
    T delta{0}, previous_delta{0};
    bool any = false;
    T h{1};

    for (int level = 0; level <= max_levels; ++level) {
        T add{0}, add_l1{0};
        const auto visit = [&](T t, bool centre) -> bool {      // false: this side is done
            const Node<T> nd = node(t);
            if (!nd.ok) return false;
            const T fx = f(nd.x);
            ++evaluations;
            const T term = nd.w * fx;
            if (!usable(fx) || !usable(term)) return centre;
            add = add + term;
            add_l1 = add_l1 + nd.w * abs(fx);
            any = true;
            return true;
        };
        if (level == 0) visit(T{0}, true);
        const int stride = (level == 0) ? 1 : 2;
        const long k_max = 10L << level;                        // t up to 10: past every usable node of a double, and of Real50 too
        for (const int sign : {+1, -1})
            for (long k = 1; k <= k_max; k += stride)
                if (!visit(T(static_cast<double>(k)) * h * T(sign), false)) break;

        if (level == 0) { sum = h * add; l1 = h * add_l1; }
        else            { sum = sum / T{2} + h * add; l1 = l1 / T{2} + h * add_l1; }
        total = scale * sum;
        if (!any) return {nan, nan, evaluations, IntegralStatus::Failed};   // nothing finite anywhere: no 0 == 0 'convergence'
        if (level > 0) {
            previous_delta = delta;
            delta = abs(T(total - last_total));
            const T bar = tol * abs(scale) * l1;
            if (level >= 2 && usable(delta) && delta <= bar) {
                const bool fell = delta <= previous_delta / T{2} || delta <= T{8} * eps * abs(scale) * l1;
                return {total, delta, evaluations, fell ? IntegralStatus::Converged : IntegralStatus::Suspicious};
            }
        }
        last_total = total;
        h = h / T{2};
    }
    return {total, delta, evaluations, IntegralStatus::DepthCapped};
}

template<Scalar T>
T default_tolerance(double requested) {
    return T(requested > 0 ? requested : std::pow(machine_epsilon<T>(), 0.75));
}

}  // namespace quadrature_detail

// ── Finite: tanh-sinh ─────────────────────────────────────────

template<Scalar T, typename F>
    requires Function<F, T, T>
IntegralResult<T> quadrature(F&& f, Finite<T> d, double tolerance = 0.0, int max_levels = 9) {
    using std::exp; using std::sinh; using std::cosh;
    using quadrature_detail::Node;
    const T nan = T(std::numeric_limits<double>::quiet_NaN());
    if (!quadrature_detail::usable(d.a) || !quadrature_detail::usable(d.b))
        return {nan, nan, 0, IntegralStatus::Failed};
    if (d.a == d.b) return {T{0}, T{0}, 0, IntegralStatus::Converged};
    T orientation{1};
    if (d.b < d.a) { std::swap(d.a, d.b); orientation = T{-1}; }
    const T half = (d.b - d.a) / T{2};
    const T pi = T(std::numbers::pi);
    const T eps = T(machine_epsilon<T>());
    const T tiny = eps * eps;

    // x = a +- half * dist, with dist = 1 - |tanh(u)| = 2q/(1+q), q = exp(-2u),
    // 2u = pi sinh|t|: written through the distance to the near end so that a
    // node close to it is not lost to the rounding of 1 - tanh. t < 0 is the
    // left end, t > 0 the right; t = 0 is the middle.
    const auto node = [&](T t) {
        if (t == T{0}) return Node<T>{T(d.a + half), pi / T{2}, true};
        const T s = sinh(t < T{0} ? -t : t);
        const T q = exp(-pi * s);
        if (q < tiny) return Node<T>{};
        const T dist = T{2} * q / (T{1} + q);
        const T w = (pi / T{2}) * cosh(t) * T{4} * q / ((T{1} + q) * (T{1} + q));
        const T x = (t > T{0}) ? T(d.b - half * dist) : T(d.a + half * dist);
        if (!(x > d.a) || !(x < d.b)) return Node<T>{};          // reached an end
        return Node<T>{x, w, true};
    };
    auto r = quadrature_detail::doubly_exponential<T>(f, node, half, quadrature_detail::default_tolerance<T>(tolerance), max_levels);
    r.value = r.value * orientation;
    return r;
}

// ── Half line: exp-sinh ───────────────────────────────────────

template<Scalar T, typename F>
    requires Function<F, T, T>
IntegralResult<T> quadrature(F&& f, HalfLine<T> d, double tolerance = 0.0, int max_levels = 9) {
    using std::exp; using std::sinh; using std::cosh;
    using quadrature_detail::Node;
    const T nan = T(std::numeric_limits<double>::quiet_NaN());
    if (!quadrature_detail::usable(d.from)) return {nan, nan, 0, IntegralStatus::Failed};
    const T pi = T(std::numbers::pi);
    const T direction = d.toward_infinity ? T{1} : T{-1};
    // x = from + direction * exp(pi/2 sinh t): t -> -inf is the end, t -> +inf infinity.
    const auto node = [&](T t) {
        const T y = exp((pi / T{2}) * sinh(t));
        const T w = (pi / T{2}) * cosh(t) * y;
        if (!quadrature_detail::usable(y) || !quadrature_detail::usable(w) || y == T{0}) return Node<T>{};
        const T x = T(d.from + direction * y);
        if (x == d.from) return Node<T>{};                       // reached the end
        return Node<T>{x, w, true};
    };
    // The integral over the set [from, inf) or (-inf, from], positive for a positive f either way.
    return quadrature_detail::doubly_exponential<T>(f, node, T{1}, quadrature_detail::default_tolerance<T>(tolerance), max_levels);
}

// ── Whole line: sinh-sinh ─────────────────────────────────────

// The line has no endpoint to carry the scalar, so it is named: `quadrature<Real50>(f, WholeLine{})`;
// double where it is not.
template<Scalar T = double, typename F>
    requires Function<F, T, T>
IntegralResult<T> quadrature(F&& f, WholeLine, double tolerance = 0.0, int max_levels = 9) {
    using std::sinh; using std::cosh;
    using quadrature_detail::Node;
    const T pi = T(std::numbers::pi);
    const auto node = [&](T t) {
        const T u = (pi / T{2}) * sinh(t);
        const T x = sinh(u);
        const T w = (pi / T{2}) * cosh(t) * cosh(u);
        if (!quadrature_detail::usable(x) || !quadrature_detail::usable(w)) return Node<T>{};
        return Node<T>{x, w, true};
    };
    return quadrature_detail::doubly_exponential<T>(f, node, T{1}, quadrature_detail::default_tolerance<T>(tolerance), max_levels);
}

// ── A second witness for a finite interval: Gauss-Kronrod G7K15 ──
//
// Its nodes are fixed and its error estimate is |K15 - G7|, a different thing
// from tanh-sinh's difference of levels: the two are independent. Constants to
// 33 digits, so the scalar may be float, double or long double (or a Dual of
// one); beyond that the rule would be exact only to long double's digits and is
// refused at compile time -- tanh-sinh is the rule for those.

namespace quadrature_detail {

struct GK15 {
    // abscissae of K15 in decreasing order; x[1], x[3], x[5], x[7] are the G7 nodes too
    static constexpr std::array<long double, 8> x = {
        0.991455371120812639206854697526329L, 0.949107912342758524526189684047851L,
        0.864864423359769072789712788640926L, 0.741531185599394439863864773280788L,
        0.586087235467691130294144838258730L, 0.405845151377397166906606412076961L,
        0.207784955007898467600689403773245L, 0.0L};
    static constexpr std::array<long double, 8> wk = {
        0.022935322010529224963732008058970L, 0.063092092629978553290700663189204L,
        0.104790010322250183839876322541518L, 0.140653259715525918745189590510238L,
        0.169004726639267902826583426598550L, 0.190350578064785409913256402421014L,
        0.204432940075298892414161999234649L, 0.209482141084727828012999174891714L};
    // Gauss 7-point weights at x[1], x[3], x[5], x[7]
    static constexpr std::array<long double, 4> wg = {
        0.129484966168869693270611432679082L, 0.279705391489276667901467771423780L,
        0.381830050505118944950369775488975L, 0.417959183673469387755102040816327L};
};

template<Scalar T, typename F>
void gk15_panel(F& f, T a, T b, T& kronrod, T& gauss, T& l1, long& evaluations) {
    using std::abs;
    const T centre = (a + b) / T{2}, half = (b - a) / T{2};
    const T fc = f(centre);
    T k = T(GK15::wk[7]) * fc;
    T g = T(GK15::wg[3]) * fc;
    T a1 = T(GK15::wk[7]) * abs(fc);
    for (std::size_t j = 0; j < 7; ++j) {
        const T dx = half * T(GK15::x[j]);
        const T fl = f(T(centre - dx)), fr = f(T(centre + dx));
        const T s = fl + fr;
        a1 = a1 + T(GK15::wk[j]) * (abs(fl) + abs(fr));
        k = k + T(GK15::wk[j]) * s;
        if (j % 2 == 1) g = g + T(GK15::wg[j / 2]) * s;
    }
    evaluations += 15;
    kronrod = k * half;
    gauss = g * half;
    l1 = a1 * abs(half);
}

template<Scalar T, typename F>
T gk15_adaptive(F& f, T a, T b, T tol, int depth, T& error, long& evaluations, bool& capped, bool& finite_values) {
    using std::abs;
    T k{0}, g{0}, l1{0};
    gk15_panel(f, a, b, k, g, l1, evaluations);
    if (!usable(k) || !usable(g)) { finite_values = false; return T{0}; }
    const T err = abs(T(k - g));
    if (depth <= 0 || err <= tol) {
        if (err > tol) capped = true;
        error = error + err;
        return k;
    }
    const T m = (a + b) / T{2};
    return gk15_adaptive(f, a, m, tol / T{2}, depth - 1, error, evaluations, capped, finite_values) +
           gk15_adaptive(f, m, b, tol / T{2}, depth - 1, error, evaluations, capped, finite_values);
}

}  // namespace quadrature_detail

template<Scalar T, typename F>
    requires Function<F, T, T>
IntegralResult<T> gauss_kronrod(F&& f, T a, T b, double tolerance = 0.0) {
    static_assert(machine_epsilon<T>() >= 1e-19,
                  "Gauss-Kronrod's constants go to long double's digits; use quadrature(f, Finite<T>{a, b}) beyond");
    const T nan = T(std::numeric_limits<double>::quiet_NaN());
    if (!quadrature_detail::usable(a) || !quadrature_detail::usable(b)) return {nan, nan, 0, IntegralStatus::Failed};
    if (a == b) return {T{0}, T{0}, 0, IntegralStatus::Converged};
    // The tolerance is relative to what the integrand weighs over the whole
    // interval (the L1 norm), as tanh-sinh's is, so the two are asked the same.
    long evaluations = 0;
    T k0{0}, g0{0}, l1{0};
    quadrature_detail::gk15_panel(f, a, b, k0, g0, l1, evaluations);
    if (!quadrature_detail::usable(k0) || !quadrature_detail::usable(l1))
        return {nan, nan, evaluations, IntegralStatus::Failed};
    const T tol = quadrature_detail::default_tolerance<T>(tolerance) * (l1 > T{0} ? l1 : T{1});
    T error{0};
    bool capped = false, finite_values = true;
    const T value = quadrature_detail::gk15_adaptive(f, a, b, tol, 20, error, evaluations, capped, finite_values);
    if (!finite_values) return {nan, nan, evaluations, IntegralStatus::Failed};
    return {value, error, evaluations, capped ? IntegralStatus::DepthCapped : IntegralStatus::Converged};
}

}  // namespace algebra
}  // namespace spatium
