#pragma once
// A ball on a manifold -- every point within geodesic distance r of a
// centre -- and a geodesic ray against it, in closed form where the space
// has one. The bound a tree over a curved space is built from, and the
// test a ray or a swept body is walked through it by: an axis-aligned box
// needs coordinates, and a chart does not preserve distances, so a node
// that misses in a chart can be hit on the manifold. A geodesic ball is
// the bound that means the same thing in every chart.
//
// A ray here is a unit-speed geodesic: origin p, unit tangent v in the
// space's own metric, gamma(t) = exp_map(p, v, t). Three closed forms:
//
//   Euclidean<N>   |p + t v - c| <= r, a quadratic;
//   Sphere<N>      cos(d(t)/R) = A cos(t/R - phi), with A and phi from
//                  <p, c> and <v, c>; the ray is periodic in 2 pi R, so it
//                  passes the ball once a turn and the first entry is
//                  returned;
//   Hyperbolic<N>  cosh d(t) = a cosh t + b sinh t, a quadratic in e^t,
//                  solved in the form that does not cancel.
//
// Also here, per space: the lower distance from a point to a ball (exact in
// any metric space, by the triangle inequality), the merge of two balls
// (along the geodesic between their centres while that geodesic is the
// shortest one, and by the triangle inequality alone past that), the
// measure a surface-area heuristic weighs a ball by, and the geometry a
// tree builder needs to know about the space.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/spaces/euclidean.hpp>
#  include <spatium/spaces/hyperbolic.hpp>
#  include <spatium/spaces/sphere.hpp>
#  include <algorithm>
#  include <cmath>
#  include <limits>
#  include <numbers>
#  include <optional>
#  include <utility>
#endif

SPATIUM_EXPORT namespace spatium::spatial {

template<typename Space>
struct GeodesicBall {
    using ScalarType = typename Space::ScalarType;
    typename Space::PointType c{};
    ScalarType r{};
};

// What a tree builder has to know about a space and cannot ask it: how far
// a geodesic stays the unique shortest path (a merge along one is valid
// only inside that), and the range of sectional curvature. Infinity for
// "never ends"; unknown spaces have no overload and get the triangle
// inequality alone.
template<Scalar T>
struct GeometryBounds {
    T injectivity_radius;
    T curvature_min, curvature_max;
};

template<std::size_t N, Scalar T>
GeometryBounds<T> geometry_bounds(const Euclidean<N, T>&) {
    return {std::numeric_limits<T>::infinity(), T{0}, T{0}};
}

template<std::size_t N, Scalar T>
GeometryBounds<T> geometry_bounds(const Sphere<N, T>& s) {
    const T k = T{1} / (s.radius * s.radius);
    return {std::numbers::pi_v<T> * s.radius, k, k};
}

template<std::size_t N, Scalar T>
GeometryBounds<T> geometry_bounds(const Hyperbolic<N, T>&) {
    return {std::numeric_limits<T>::infinity(), T{-1}, T{-1}};
}

// No point of the ball is closer to p than this, in any metric space.
template<typename Space>
typename Space::ScalarType lower_distance(const Space& space, const GeodesicBall<Space>& b,
                                          const typename Space::PointType& p) {
    const auto d = space.distance(p, b.c) - b.r;
    return d > decltype(d){0} ? d : decltype(d){0};
}

// A ball containing both. Along the geodesic from one centre to the other
// the smallest such ball has its centre where the two far edges are
// equidistant; that holds while the geodesic is the shortest path, which
// is what the injectivity radius bounds. Past it, the first ball grown by
// the distance and the second radius, which the triangle inequality alone
// guarantees. A pad of a few ulps keeps containment true in floating point.
template<typename Space>
GeodesicBall<Space> merge(const Space& space, const GeodesicBall<Space>& a, const GeodesicBall<Space>& b) {
    using T = typename Space::ScalarType;
    const T d = space.distance(a.c, b.c);
    if (d + b.r <= a.r) return a;
    if (d + a.r <= b.r) return b;
    const T pad = std::numeric_limits<T>::epsilon() * T{16} * (T{1} + d + a.r + b.r);
    if (d < geometry_bounds(space).injectivity_radius) {
        const T r = (d + a.r + b.r) / T{2};
        return {space.exp_map(a.c, space.log_map(a.c, b.c), (r - a.r) / d), r + pad};
    }
    return {a.c, d + b.r + pad};
}

// What a random geodesic's chance of crossing the ball scales with: the
// measure of its boundary, by Cauchy-Crofton (Santalo in constant
// curvature). Only ratios within one tree are used.
template<std::size_t N, Scalar T>
T sah_measure(const Euclidean<N, T>&, const GeodesicBall<Euclidean<N, T>>& b) {
    constexpr T pi = std::numbers::pi_v<T>;
    if constexpr (N == 2) return T{2} * pi * b.r;
    else return T{4} * pi * b.r * b.r;
}

template<std::size_t N, Scalar T>
T sah_measure(const Sphere<N, T>& s, const GeodesicBall<Sphere<N, T>>& b) {
    using std::sin;
    constexpr T pi = std::numbers::pi_v<T>;
    const T rs = s.radius * sin(std::min(b.r / s.radius, pi / T{2}));
    if constexpr (N == 2) return T{2} * pi * rs;
    else return T{4} * pi * rs * rs;
}

template<std::size_t N, Scalar T>
T sah_measure(const Hyperbolic<N, T>&, const GeodesicBall<Hyperbolic<N, T>>& b) {
    using std::sinh;
    constexpr T pi = std::numbers::pi_v<T>;
    if constexpr (N == 2) return T{2} * pi * sinh(b.r);
    else return T{4} * pi * sinh(b.r) * sinh(b.r);
}

// ── A geodesic ray against a ball ──────────────────────────────
//
// The parameter interval the unit-speed ray from `p` along `v` spends in
// the ball, starting no earlier than 0, if it enters before `tmax`. The
// end of the interval may lie past `tmax`.

template<std::size_t N, Scalar T>
std::optional<std::pair<T, T>> ray_interval(const Euclidean<N, T>&, const Vec<T, N>& p, const Vec<T, N>& v,
                                            const GeodesicBall<Euclidean<N, T>>& b, T tmax) {
    using std::sqrt;
    const Vec<T, N> oc{p - b.c};
    const T half_b = v.dot(oc);
    const T cc = oc.dot(oc) - b.r * b.r;
    const T disc = half_b * half_b - cc;
    if (disc < T{0}) return std::nullopt;
    const T s = sqrt(disc);
    const T t1 = -half_b + s;
    if (t1 < T{0}) return std::nullopt;
    const T t0 = std::max(-half_b - s, T{0});
    if (t0 > tmax) return std::nullopt;
    return std::pair{t0, t1};
}

// gamma(t) = cos(t/R) p + R sin(t/R) v, so <gamma, c> = a cos + R b sin with
// a = <p, c>, b = <v, c>, and cos(d/R) = <gamma, c>/R^2 = A cos(theta - phi).
// The ball is where cos(theta - phi) >= cos(r/R)/A: an arc of theta of half
// width alpha around phi, once every turn.
template<std::size_t N, Scalar T>
std::optional<std::pair<T, T>> ray_interval(const Sphere<N, T>& s, const Vec<T, N + 1>& p,
                                            const Vec<T, N + 1>& v, const GeodesicBall<Sphere<N, T>>& b,
                                            T tmax) {
    using std::acos; using std::atan2; using std::cos; using std::floor; using std::sqrt;
    constexpr T two_pi = T{2} * std::numbers::pi_v<T>;
    const T R = s.radius;
    const T a = p.dot(b.c), bb = v.dot(b.c) * R;
    const T A = sqrt(a * a + bb * bb) / (R * R);
    const T cr = cos(std::min(b.r / R, std::numbers::pi_v<T>));
    if (A <= T{0}) {   // the whole great circle is at a quarter turn from c
        return cr <= T{0} ? std::optional{std::pair{T{0}, two_pi * R}} : std::nullopt;
    }
    const T k = cr / A;
    if (k > T{1}) return std::nullopt;
    if (k <= T{-1}) return std::pair{T{0}, std::numeric_limits<T>::infinity()};   // always inside
    const T phi = atan2(bb, a);
    const T alpha = acos(k);
    // The arc [phi - alpha, phi + alpha] repeats every 2 pi; take the first
    // copy that has not ended by theta = 0.
    T lo = phi - alpha, hi = phi + alpha;
    const T shift = floor(-hi / two_pi) + T{1};
    lo += shift * two_pi;
    hi += shift * two_pi;
    if (hi - two_pi >= T{0}) { lo -= two_pi; hi -= two_pi; }
    const T t0 = std::max(lo, T{0}) * R, t1 = hi * R;
    if (t0 > tmax) return std::nullopt;
    return std::pair{t0, t1};
}

// gamma(t) = cosh t p + sinh t v, so cosh d(t) = a cosh t + b sinh t with
// a = -<p, c>_L >= 1 and b = -<v, c>_L, |b| < a. Inside the ball is
// a cosh t + b sinh t <= cosh r; in u = e^t that is
// (a + b) u^2 - 2 cosh(r) u + (a - b) <= 0. Its smaller root is taken as
// (a - b) / (C + sqrt(disc)) rather than (C - sqrt(disc)) / (a + b), which
// loses every digit when the ball is small and far.
template<std::size_t N, Scalar T>
std::optional<std::pair<T, T>> ray_interval(const Hyperbolic<N, T>&, const Vec<T, N + 1>& p,
                                            const Vec<T, N + 1>& v, const GeodesicBall<Hyperbolic<N, T>>& b,
                                            T tmax) {
    using std::cosh; using std::log; using std::sqrt;
    using H = Hyperbolic<N, T>;
    const T a = -H::minkowski(p, b.c), bb = -H::minkowski(v, b.c);
    const T C = cosh(b.r);
    const T disc = C * C - (a * a - bb * bb);
    if (disc < T{0}) return std::nullopt;
    const T sq = sqrt(disc);
    const T u_hi = (C + sq) / (a + bb);
    if (u_hi < T{1}) return std::nullopt;   // left the ball before t = 0
    const T u_lo = (a - bb) / (C + sq);
    const T t0 = u_lo > T{1} ? log(u_lo) : T{0};
    if (t0 > tmax) return std::nullopt;
    return std::pair{t0, log(u_hi)};
}

} // namespace spatium::spatial
