#pragma once

// What the pairs written by rsc/spatium/rsc_find.py share, and what the finder's own runs use.
//
// A route found in the typed graph (rsc/spatium/graph.json) is turned into C++ by rsc/spatium/emit.py:
// one function per edge, the last one the route. Nothing here is per route. This header is the part
// that is the same for every route:
//
//   - draw<Cell>(i)      the i-th input: points and tangents from the connectivity matrix's sampler for
//                        the space (tests/connectivity/probes.hpp), seeded from the Halton sequence so a
//                        run is the same on every machine;
//   - up_* / down        a sample into another scalar and an answer back to double;
//   - dist, mag, show    how far apart two answers are (the largest coordinate), how large an input is,
//                        and how to print one;
//   - the bound          what each path promises, as pair.hpp asks: a sum of two allowances, each a
//                        multiple of the arithmetic's epsilon in the input's size. A path that goes
//                        through an inverse (log, distance and what is built on them) is held to the
//                        matrix's own tolerance, sqrt(eps) * 64, because acos and acosh lose half the
//                        digits near coincident points; the rest to 1000 eps;
//   - compare            a pair run by the finder: the same as symmetry::check, and it also keeps the
//                        worst input, the calls a route makes (its cost) and its time.
//
// The cost of a route is the number of operation calls it makes for one input: every emitted edge
// counts one, except the conversions of a scalar lift (calls()).

#include "pair.hpp"

#include "../connectivity/probes.hpp"

#include <spatium/algebra/groups/se3.hpp>
#include <spatium/algebra/groups/so3.hpp>
#include <spatium/core/access.hpp>
#include <spatium/core/epsilon.hpp>
#include <spatium/spaces/lie_group.hpp>
#if SPATIUM_HAS_BOOST_MULTIPRECISION
#  include <spatium/core/precision.hpp>
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

namespace rsc_pair {

using spatium::Matrix;
using spatium::Vec;

// ── Cost ───────────────────────────────────────────────────────

inline std::size_t& calls() {
    static thread_local std::size_t n = 0;
    return n;
}

// ── Values, not expressions ────────────────────────────────────
// Over Real50 a sum or a square root is an expression template that refers to the locals it was made
// from; returned out of a function it dangles. Every emitted function returns val(...), the number.

template<class X>
auto val(X&& x) {
    using D = std::remove_cvref_t<X>;
    if constexpr (!std::is_arithmetic_v<D> && requires { typename D::result_type; }) return typename D::result_type(std::forward<X>(x));
    else return D(std::forward<X>(x));
}

// ── Scalars: a sample into another scalar, an answer back ──────

template<class R, class X> auto up(const X& x) { return connectivity::cast<R>(x); }

template<class R, class S, class P>
auto up_pt(const S& s, const P& p) {   // as the matrix does it: put back on the space, or every cell would measure double's rounding
    using X = spatium::spaces::point_t<S>;
    auto q = connectivity::cast<R>(p);
    if constexpr (requires { spatium::spaces::project(s, q); }) return X{spatium::spaces::project(s, q)};
    else return X{q};
}
template<class R, class V> auto up_tan(const V& v) { return connectivity::cast<R>(v); }
template<class R> R up_sc(double t) { return connectivity::from_double<R>(t); }

// To double: a scalar of any depth (Real50, Dual) by its primal, a Vec, a Matrix or a product point entry by entry.
template<class X> double down(const X& x) { return spatium::primal_double(x); }
template<class T, std::size_t N> Vec<double, N> down(const Vec<T, N>& a) {
    Vec<double, N> r;
    for (std::size_t i = 0; i < N; ++i) r.data[i] = down(a.data[i]);
    return r;
}
template<class T, std::size_t R, std::size_t C> Matrix<double, R, C> down(const Matrix<T, R, C>& a) {
    Matrix<double, R, C> r;
    for (std::size_t i = 0; i < R * C; ++i) r.data[i] = down(a.data[i]);
    return r;
}
template<class A, class B> auto down(const spatium::product_detail::Pair<A, B>& a) {
    return spatium::product_detail::Pair<decltype(down(a.first)), decltype(down(a.second))>{down(a.first), down(a.second)};
}

// How far apart two answers are: the largest coordinate of the difference.
inline double dist_d(double a, double b) { return std::abs(a - b); }
template<std::size_t N> double dist_d(const Vec<double, N>& a, const Vec<double, N>& b) {
    double m = 0;
    for (std::size_t i = 0; i < N; ++i) m = std::max(m, std::abs(a.data[i] - b.data[i]));
    return m;
}
template<std::size_t R, std::size_t C> double dist_d(const Matrix<double, R, C>& a, const Matrix<double, R, C>& b) {
    double m = 0;
    for (std::size_t i = 0; i < R * C; ++i) m = std::max(m, std::abs(a.data[i] - b.data[i]));
    return m;
}
template<class A, class B> double dist_d(const spatium::product_detail::Pair<A, B>& a, const spatium::product_detail::Pair<A, B>& b) {
    return std::max(dist_d(a.first, b.first), dist_d(a.second, b.second));
}
template<class A, class B> double dist(const A& a, const B& b) { return dist_d(down(a), down(b)); }

// How large an input is.
inline double mag(double x) { return std::abs(x); }
template<class T, std::size_t N> double mag(const Vec<T, N>& a) {
    double m = 0;
    for (std::size_t i = 0; i < N; ++i) m = std::max(m, std::abs(spatium::primal_double(a.data[i])));
    return m;
}
template<class T, std::size_t R, std::size_t C> double mag(const Matrix<T, R, C>& a) {
    double m = 0;
    for (std::size_t i = 0; i < R * C; ++i) m = std::max(m, std::abs(spatium::primal_double(a.data[i])));
    return m;
}
template<class A, class B> double mag(const spatium::product_detail::Pair<A, B>& a) { return std::max(mag(a.first), mag(a.second)); }

inline std::string show(double x) { char b[40]; std::snprintf(b, sizeof b, "%.17g", x); return b; }
template<class T, std::size_t N> std::string show(const Vec<T, N>& a) {
    std::string s = "[";
    for (std::size_t i = 0; i < N; ++i) s += (i ? " " : "") + show(spatium::primal_double(a.data[i]));
    return s + "]";
}
template<class T, std::size_t R, std::size_t C> std::string show(const Matrix<T, R, C>& a) {
    std::string s = "[";
    for (std::size_t i = 0; i < R * C; ++i) s += (i ? " " : "") + show(spatium::primal_double(a.data[i]));
    return s + "]";
}
template<class A, class B> std::string show(const spatium::product_detail::Pair<A, B>& a) { return "(" + show(a.first) + ", " + show(a.second) + ")"; }

// ── Inputs: the matrix's samplers, one draw per index ──────────

template<class Cell>
struct Draw {
    using Pt = typename decltype(Cell::points(std::declval<connectivity::Rand&>()))::value_type;
    using Tv = decltype(Cell::tangent(std::declval<connectivity::Rand&>(), std::declval<const Pt&>()));
    Pt p, q;   // two points of the space
    Tv u, v;   // two tangent vectors at p
};

template<class Cell>
Draw<Cell> draw(std::size_t i) {
    connectivity::Rand r;
    // the Halton sequence picks the seed, so the inputs are the same on every machine and well spread in index
    r.rng.seed(20260928ULL + static_cast<std::uint64_t>(symmetry::halton(i, 0) * 4294967296.0));
    const auto pts = Cell::points(r);
    Draw<Cell> d;
    d.p = pts[0];
    d.q = pts[1];
    d.u = Cell::tangent(r, d.p);
    d.v = Cell::tangent(r, d.p);
    return d;
}

template<class Cell>
std::string show_draw(const Draw<Cell>& d) {
    return "p=" + show(d.p) + " q=" + show(d.q) + " u=" + show(d.u) + " v=" + show(d.v);
}

template<class Cell>
double scale_of(const Draw<Cell>& d) { return 1.0 + std::max({mag(d.p), mag(d.q), mag(d.u), mag(d.v)}); }

// ── What each path promises ────────────────────────────────────
// smooth: 1000 eps in the size of the input. inverse: the matrix's tolerance, sqrt(eps) * 64.
// The pair's bound is the two paths' allowances added, as pair.hpp has it.

enum class Tier { Smooth, Inverse };

inline constexpr double allowance(Tier t) {
    constexpr double eps = std::numeric_limits<double>::epsilon();
    return t == Tier::Smooth ? 1000.0 * eps : 64.0 * 1.4901161193847656e-08;   // sqrt(eps) = 2^-26
}

template<class Cell>
double bound(Tier t, const Draw<Cell>& d) { return 2.0 * allowance(t) * scale_of(d); }

// ── A pair, run by the finder ──────────────────────────────────

struct Row {
    std::size_t samples = 0, failures = 0, first_failure = 0, worst_index = 0;
    double worst_ratio = 0, worst_dist = 0, worst_bound = 0;
    std::size_t calls = 0, ref_calls = 0;
    double ns = 0, ref_ns = 0;   // mean time of one evaluation
    std::string worst_input;
};

template<class Cell, class Ref, class Route, class Bound>
Row compare(std::size_t n, Ref ref, Route route, Bound bnd) {
    using Clock = std::chrono::steady_clock;
    Row row;
    row.samples = n;
    double t_route = 0, t_ref = 0;
    for (std::size_t i = 1; i <= n; ++i) {
        const Draw<Cell> d = draw<Cell>(i);
        calls() = 0;
        const auto r0 = Clock::now();
        const auto a = ref(d);
        const auto r1 = Clock::now();
        if (i == 1) row.ref_calls = calls();
        calls() = 0;
        const auto b = route(d);
        const auto r2 = Clock::now();
        if (i == 1) row.calls = calls();
        t_ref += std::chrono::duration<double, std::nano>(r1 - r0).count();
        t_route += std::chrono::duration<double, std::nano>(r2 - r1).count();

        const double dd = dist(a, b);
        const double tol = bnd(d);
        const bool ok = std::isfinite(dd) && dd <= tol;
        if (!ok) {
            ++row.failures;
            if (row.first_failure == 0) row.first_failure = i;
        }
        const double ratio = std::isfinite(dd) ? dd / tol : std::numeric_limits<double>::infinity();
        if (ratio > row.worst_ratio || row.worst_index == 0) {
            row.worst_ratio = ratio;
            row.worst_dist = dd;
            row.worst_bound = tol;
            row.worst_index = i;
            row.worst_input = show_draw(d);
        }
    }
    row.ns = n ? t_route / static_cast<double>(n) : 0;
    row.ref_ns = n ? t_ref / static_cast<double>(n) : 0;
    return row;
}

// ── Samplers for the spaces that are not cells of the matrix ───
// A Lie group as a space: rotations and rigid motions from small exponentials, so that the relative
// motion between two stays inside the injectivity radius of log.

struct LieSO3 {
    static constexpr const char* name = "LieSO3";
    template<class T> static spatium::LieGroupManifold<spatium::algebra::SO3<T>> make() { return {}; }
    static std::vector<Matrix<double, 3, 3>> points(connectivity::Rand& r) {
        const spatium::algebra::SO3<double> g;
        std::vector<Matrix<double, 3, 3>> p;
        for (int i = 0; i < 6; ++i) p.push_back(g.exp(Vec<double, 3>{r(-0.7, 0.7), r(-0.7, 0.7), r(-0.7, 0.7)}));
        return p;
    }
    static Vec<double, 3> tangent(connectivity::Rand& r, const Matrix<double, 3, 3>&) { return {r(), r(), r()}; }
};

struct LieSE3 {
    static constexpr const char* name = "LieSE3";
    template<class T> static spatium::LieGroupManifold<spatium::algebra::SE3<T>> make() { return {}; }
    static std::vector<Matrix<double, 4, 4>> points(connectivity::Rand& r) {
        const spatium::algebra::SE3<double> g;
        std::vector<Matrix<double, 4, 4>> p;
        for (int i = 0; i < 6; ++i)
            p.push_back(g.exp(Vec<double, 6>{r(-0.7, 0.7), r(-0.7, 0.7), r(-0.7, 0.7), r(-1, 1), r(-1, 1), r(-1, 1)}));
        return p;
    }
    static Vec<double, 6> tangent(connectivity::Rand& r, const Matrix<double, 4, 4>&) { return {r(), r(), r(), r(), r(), r()}; }
};

}  // namespace rsc_pair
