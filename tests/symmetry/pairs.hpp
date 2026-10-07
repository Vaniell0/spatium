#pragma once

// The registry of pairs of paths. One place to read which answers the library
// computes two ways, how far apart the two may be, and how many inputs it was
// checked on. A pair is added here when it is written, not discovered later.
//
// Each pair says where its tolerance comes from, and what the first run measured.

#include "pair.hpp"

#include <spatium/algebra/calculus.hpp>
#include <spatium/algebra/dual.hpp>
#include <spatium/algebra/vector.hpp>
#include <spatium/core/access.hpp>
#include <spatium/core/precision.hpp>
#include <spatium/geometry/geometry.hpp>
#include <spatium/geometry/make.hpp>
#include <spatium/geometry/ray_parametric.hpp>
#include <spatium/geometry/ray_surface.hpp>
#include <spatium/physics/mechanics/narrow_phase.hpp>
#include <spatium/spaces/parametric.hpp>
#include <spatium/spaces/spd.hpp>
#include <spatium/spaces/sphere.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <type_traits>
#include <vector>

namespace symmetry {

using spatium::Vec;

// ── ray_torus: the closed form on double against the same code on Real50 ──
//
// A ray from a point of the box around the torus (R = 2, r = 0.5), aimed at a
// point inside the torus's own box, so most rays hit and some graze. The
// reference is the same function on fifty digits; the route is `double`. The
// bound is 1e-9 of the scale of the ray. A hit with two coinciding roots could
// only promise sqrt(eps) = 1.5e-8, but none of the 6 000 rays is that close, the
// worst distance is 8e-13, and a bound nearer the observation is what lets the
// pair see a fault of 1e-6.
//
// First run (20 000 rays, 15 846 with hits): `ray_torus<Real50>` did not compile
// (a qualified std::abs in torus_basis), and once it did, one ray (the 5 037th)
// had both hits off by 4.3e-4 while the other roots were 3.8 away. The cause was
// the closed forms' absolute branch tests; the roots are now polished by Newton
// steps on the polynomial itself, and the worst error over the 20 000 is 2e-11.
struct TorusRay { Vec<double, 3> o, d; };

inline TorusRay torus_ray(std::size_t i) {
    const auto u = [&](std::size_t k) { return halton(i, k); };
    Vec<double, 3> o{-6 + 12 * u(0), -6 + 12 * u(1), -2 + 4 * u(2)};
    const Vec<double, 3> target{-2.6 + 5.2 * u(3), -2.6 + 5.2 * u(4), -0.6 + 1.2 * u(5)};
    Vec<double, 3> d = target - o;
    d = d / d.norm();
    return {o, d};
}

inline std::vector<double> torus_hits_double(const TorusRay& s) {
    using namespace spatium::geometry;
    Torus<double> torus{.major_radius = 2.0, .minor_radius = 0.5};
    std::vector<double> t;
    for (const auto& h : ray_torus(Ray<3, double>{s.o, s.d}, torus)) t.push_back(h.t);
    return t;
}

inline std::vector<double> torus_hits_real50(const TorusRay& s) {
    using namespace spatium::geometry;
    using R = spatium::Real50;
    Torus<R> torus{.major_radius = R(2.0), .minor_radius = R(0.5)};
    const Vec<R, 3> o{R(s.o[0]), R(s.o[1]), R(s.o[2])};
    const Vec<R, 3> d{R(s.d[0]), R(s.d[1]), R(s.d[2])};
    std::vector<double> t;
    for (const auto& h : ray_torus(Ray<3, R>{o, d}, torus)) t.push_back(spatium::primal_double(h.t));
    return t;
}

inline Entry ray_torus_precision() {
    return {"ray_torus: double against Real50", 6000, [] {
        return check<TorusRay>(
            "ray_torus: double against Real50", Kind::Function, 6000, torus_ray,
            torus_hits_real50, torus_hits_double,
            [](const std::vector<double>& a, const std::vector<double>& b) {
                if (a.size() != b.size()) return std::numeric_limits<double>::infinity();          // a hit gained or lost is not a rounding
                double worst = 0.0;
                for (std::size_t k = 0; k < a.size(); ++k) worst = std::max(worst, std::abs(a[k] - b[k]));
                return worst;
            },
            [](const TorusRay& s) { return 1e-9 * (1.0 + s.o.norm()); });
    }};
}

// ── gradient: the exact derivative through Dual against central differences ──
//
// The bound is the central difference's own: a step of cbrt(eps) leaves a truncation of
// order h^2 f''' and a rounding of order eps / h, both near 1e-11 for this f; 1e-9 of the
// size of the gradient is the promise, and the worst distance measured is 6e-11 of it.
template<class T>
T gradient_test_function(const Vec<T, 3>& x) {
    using std::sin; using std::exp;
    return sin(x[0] * x[1]) + exp(x[2]) * x[0];
}

inline Entry gradient_dual_against_differences() {
    return {"gradient: Dual against central differences", 2000, [] {
        return check<Vec<double, 3>>(
            "gradient: Dual against central differences", Kind::Function, 2000,
            [](std::size_t i) { return Vec<double, 3>{-2 + 4 * halton(i, 0), -2 + 4 * halton(i, 1), -2 + 4 * halton(i, 2)}; },
            [](const Vec<double, 3>& x) {
                return spatium::gradient([](const Vec<spatium::Dual<double>, 3>& v) { return gradient_test_function(v); }, x);
            },
            [](const Vec<double, 3>& x) {
                const double h = std::cbrt(std::numeric_limits<double>::epsilon());
                Vec<double, 3> g{};
                for (std::size_t k = 0; k < 3; ++k) {
                    Vec<double, 3> hi = x, lo = x;
                    const double step = h * std::max(1.0, std::abs(x[k]));
                    hi[k] += step; lo[k] -= step;
                    g[k] = (gradient_test_function(hi) - gradient_test_function(lo)) / (2 * step);
                }
                return g;
            },
            [](const Vec<double, 3>& a, const Vec<double, 3>& b) { return (a - b).norm(); },
            [](const Vec<double, 3>& x) { return 1e-9 * (1.0 + spatium::gradient([](const Vec<spatium::Dual<double>, 3>& v) { return gradient_test_function(v); }, x).norm()); });
    }};
}

// ── exp then log on the sphere: a relation, not a second path ──
//
// log_p(exp_p(v)) = v for |v| inside the injectivity radius (pi r). Checked on tangent
// vectors of length up to 2.75 r, away from the cut locus where acos loses its footing.
// The bound is 1e-11 r: acos near +-1 is what limits it and the vectors are kept off both,
// and the worst distance measured over 2 000 draws is 2e-14.
struct SphereVector { Vec<double, 3> p, v; };

inline SphereVector sphere_vector(std::size_t i, double radius) {
    const double w = 2 * halton(i, 0) - 1, phi = 2 * std::numbers::pi * halton(i, 1);
    const double rho = std::sqrt(std::max(0.0, 1 - w * w));
    const Vec<double, 3> n{rho * std::cos(phi), rho * std::sin(phi), w};
    const Vec<double, 3> a{2 * halton(i, 2) - 1, 2 * halton(i, 3) - 1, 2 * halton(i, 4) - 1};
    Vec<double, 3> t = a - n * a.dot(n);
    const double length = radius * (0.05 + 2.7 * halton(i, 5));
    t = t * (length / t.norm());
    return {n * radius, t};
}

inline Entry sphere_exp_log() {
    return {"Sphere: log after exp is the identity", 2000, [] {
        const spatium::Sphere<2> s{.radius = 2.0};
        return check<SphereVector>(
            "Sphere: log after exp is the identity", Kind::Relation, 2000,
            [](std::size_t i) { return sphere_vector(i, 2.0); },
            [](const SphereVector& x) { return x.v; },
            [&](const SphereVector& x) {
                const auto q = spatium::spaces::exp_map(s, x.p, x.v, 1.0);
                return Vec<double, 3>{spatium::spaces::log_map(s, x.p, q)};
            },
            [](const Vec<double, 3>& a, const Vec<double, 3>& b) { return (a - b).norm(); },
            [](const SphereVector&) { return 1e-11 * 2.0; });
    }};
}

// ── ray_torus: the quartic against the general path through the chart ──
//
// The closed form against the same torus as a ParametricSurface, hit by ray_parametric, a
// Newton search over the chart that knows nothing about quartics. Two paths that share no
// code, on the rays of the first pair. The route's own promise is its Newton tolerance, 1e-6;
// over the 3 000 rays the worst distance measured is 3.5e-7 and no hit is gained or lost.
inline Entry ray_torus_against_chart() {
    return {"ray_torus: the quartic against the parametric search", 3000, [] {
        const auto surface = spatium::make_torus<double>(2.0, 0.5);
        return check<TorusRay>(
            "ray_torus: the quartic against the parametric search", Kind::Function, 3000, torus_ray,
            torus_hits_double,
            [&](const TorusRay& s) {
                std::vector<double> t;
                for (const auto& h : spatium::geometry::ray_parametric(spatium::geometry::Ray<3, double>{s.o, s.d}, surface))
                    t.push_back(h.t);
                return t;
            },
            [](const std::vector<double>& a, const std::vector<double>& b) {
                if (a.size() != b.size()) return std::numeric_limits<double>::infinity();
                double worst = 0.0;
                for (std::size_t k = 0; k < a.size(); ++k) worst = std::max(worst, std::abs(a[k] - b[k]));
                return worst;
            },
            [](const TorusRay& s) { return 1e-6 * (1.0 + s.o.norm()); });
    }};
}

// ── SPD(3): the eigenvalues through the characteristic cubic, double against Real50 ──
//
// S = Q diag(1, 1 + gap, 3) Q^T with Q from three Euler angles and the gap log-uniform in
// [1e-10, 1], rounded once to double and handed to both scalars, so the two paths see one
// matrix. The bound is the conditioning of a cubic root: the coefficients carry a rounding of
// eps S^3, a root with gap g to its neighbour moves by that over its slope there, about g, and
// the slope cannot be smaller than sqrt(eps) S before the root is a double one, where the error
// saturates at sqrt(eps) S. Thirty times that is the promise (about four times the worst measured, 2.2e-7 at gaps
// near 1e-7, against a bound of 4.3e-6 there).
//
// First run: `eigen_sym<Real50>` did not compile (a unary minus on a Boost number is a lazy
// expression and solve_cubic<T> cannot deduce T from four different types).
struct SpdCase { std::array<double, 9> m; double gap; };

inline SpdCase spd_case(std::size_t i) {
    const double a = 2 * std::numbers::pi * halton(i, 0), b = std::numbers::pi * halton(i, 1),
                 c = 2 * std::numbers::pi * halton(i, 2);
    const double gap = std::pow(10.0, -10.0 * halton(i, 3));
    const double D[3] = {1.0, 1.0 + gap, 3.0};
    const double ca = std::cos(a), sa = std::sin(a), cb = std::cos(b), sb = std::sin(b), cc = std::cos(c), sc = std::sin(c);
    const double Rz1[9] = {ca, -sa, 0, sa, ca, 0, 0, 0, 1}, Rx[9] = {1, 0, 0, 0, cb, -sb, 0, sb, cb},
                 Rz2[9] = {cc, -sc, 0, sc, cc, 0, 0, 0, 1};
    const auto mul = [](const double* A, const double* B, double* C) {
        for (int r = 0; r < 3; ++r)
            for (int k = 0; k < 3; ++k) {
                double sum = 0;
                for (int j = 0; j < 3; ++j) sum += A[r * 3 + j] * B[j * 3 + k];
                C[r * 3 + k] = sum;
            }
    };
    double T1[9], Q[9];
    mul(Rz1, Rx, T1);
    mul(T1, Rz2, Q);
    SpdCase out{};
    out.gap = gap;
    for (int r = 0; r < 3; ++r)
        for (int k = 0; k < 3; ++k) {
            double sum = 0;
            for (int j = 0; j < 3; ++j) sum += Q[r * 3 + j] * D[j] * Q[k * 3 + j];
            out.m[r * 3 + k] = sum;
        }
    for (int r = 0; r < 3; ++r)
        for (int k = r + 1; k < 3; ++k) out.m[r * 3 + k] = out.m[k * 3 + r] = 0.5 * (out.m[r * 3 + k] + out.m[k * 3 + r]);
    return out;
}

template<class T>
std::array<double, 3> spd_eigenvalues(const SpdCase& s) {
    spatium::Matrix<T, 3, 3> M;
    for (int r = 0; r < 3; ++r)
        for (int k = 0; k < 3; ++k) M(r, k) = T(s.m[r * 3 + k]);
    const auto e = spatium::detail::eigen_sym(M);
    std::array<double, 3> v{spatium::primal_double(e.values[0]), spatium::primal_double(e.values[1]),
                            spatium::primal_double(e.values[2])};
    std::sort(v.begin(), v.end());
    return v;
}

inline Entry spd3_eigenvalues_precision() {
    return {"SPD(3) eigenvalues: double against Real50", 4000, [] {
        return check<SpdCase>(
            "SPD(3) eigenvalues: double against Real50", Kind::Function, 4000, spd_case,
            spd_eigenvalues<spatium::Real50>, spd_eigenvalues<double>,
            [](const std::array<double, 3>& a, const std::array<double, 3>& b) {
                double worst = 0.0;
                for (int k = 0; k < 3; ++k) worst = std::max(worst, std::abs(a[k] - b[k]));
                return worst;
            },
            [](const SpdCase& s) {
                const double eps = std::numeric_limits<double>::epsilon(), S = 4.0;
                return 30.0 * eps * S * S * S / std::max(s.gap, std::sqrt(eps) * S);
            });
    }};
}

// ── Distances and hits between primitives: double against Real50 ──
//
// Ten answers from one draw of 40 numbers in [-2, 2]: a point against a segment, a triangle
// and a box; two segments, a segment and a triangle, two triangles, two lines; a ray against a
// triangle and against a box. The same function templates on fifty digits are the reference.
// A hit gained or lost is a failure (the flag is part of the answer); the bound is 1e-9, the
// worst measured distance over the draws is below it by orders of magnitude.
//
// First run: none of this compiled on Real50. Qualified std::abs / std::clamp / std::min /
// std::max / std::swap in the segment-segment and line-line distances, in the ray-triangle and
// ray-box intersections and in the barycentric coordinates of a triangle blocked the scalar's
// own functions, and expression templates (`auto t1 = (a - b) * c`) broke std::swap and
// std::max. A pair that exercises a function is what makes it worth fixing.
template<class T>
Vec<T, 3> geometry_point(const double* a) { return Vec<T, 3>{T(a[0]), T(a[1]), T(a[2])}; }

template<class T>
std::vector<double> geometry_answers(const std::vector<double>& u) {
    using namespace spatium::geometry;
    const double* d = u.data();
    std::vector<double> out;
    const auto push = [&](const T& x) { out.push_back(spatium::primal_double(x)); };
    const Vec<T, 3> p = geometry_point<T>(d);
    const Segment<3, T> s1{geometry_point<T>(d + 3), geometry_point<T>(d + 6)};
    const Segment<3, T> s2{geometry_point<T>(d + 9), geometry_point<T>(d + 12)};
    const Triangle<3, T> t1(geometry_point<T>(d + 15), geometry_point<T>(d + 18), geometry_point<T>(d + 21));
    const Triangle<3, T> t2(geometry_point<T>(d + 24), geometry_point<T>(d + 27), geometry_point<T>(d + 30));
    const Vec<T, 3> lo = geometry_point<T>(d + 33);
    const Box<3, T> box{lo, lo + Vec<T, 3>{T(1), T(1), T(1)}};
    push(distance(p, s1));
    push(distance(p, t1));
    push(distance(p, box));
    push(distance(s1, s2));
    push(distance(s1, t1));
    push(distance(t1, t2));
    const auto l1 = Line<3, T>::from(geometry_point<T>(d + 3), geometry_point<T>(d + 6));
    const auto l2 = Line<3, T>::from(geometry_point<T>(d + 9), geometry_point<T>(d + 12));
    if (l1 && l2) push(distance(*l1, *l2)); else out.push_back(-1.0);
    const auto ray = Ray<3, T>::from(p, geometry_point<T>(d + 36) - p);
    if (ray) {
        const auto h = intersect(*ray, t1);
        out.push_back(h ? 1.0 : 0.0);
        if (h) push((*h)[0]);
        const auto hb = intersect(*ray, box);
        out.push_back(hb ? 1.0 : 0.0);
        if (hb) push((*hb)[0]);
    }
    return out;
}

inline Entry geometry_precision() {
    return {"geometry: distances and hits, double against Real50", 3000, [] {
        return check<std::vector<double>>(
            "geometry: distances and hits, double against Real50", Kind::Function, 3000,
            [](std::size_t i) {
                std::vector<double> u(40);
                for (std::size_t k = 0; k < 40; ++k) u[k] = -2 + 4 * halton(k < 32 ? i : i + 104729, k % 32);
                return u;
            },
            geometry_answers<spatium::Real50>, geometry_answers<double>,
            [](const std::vector<double>& a, const std::vector<double>& b) {
                if (a.size() != b.size()) return std::numeric_limits<double>::infinity();
                double worst = 0.0;
                for (std::size_t k = 0; k < a.size(); ++k) worst = std::max(worst, std::abs(a[k] - b[k]));
                return worst;
            },
            [](const std::vector<double>&) { return 1e-9; });
    }};
}

// ── The derivative of a distance to a primitive: Dual against central differences ──
//
// The gradient of the distance from a point to a segment, a triangle and a box with respect to
// the point, taken through `Dual` and by differences. Away from the places where the nearest
// feature switches (a measure-zero set, and a step of 1e-6 does not reach it on these draws) the
// distance is smooth, and the two must agree to what a difference allows, 1e-8 here (the worst measured is 7e-10).
template<class T>
std::vector<double> distance_gradients(const std::vector<double>& u, const Vec<double, 3>& p0, int axis, double shift) {
    using namespace spatium::geometry;
    const double* d = u.data();
    Vec<T, 3> p{T(p0[0]), T(p0[1]), T(p0[2])};
    if constexpr (std::is_same_v<T, spatium::Dual<double>>) {
        p[axis] = T::variable(p0[axis]);
    } else {
        p[axis] = T(p0[axis] + shift);
    }
    const Segment<3, T> s{geometry_point<T>(d + 3), geometry_point<T>(d + 6)};
    const Triangle<3, T> t(geometry_point<T>(d + 15), geometry_point<T>(d + 18), geometry_point<T>(d + 21));
    const Vec<T, 3> lo = geometry_point<T>(d + 33);
    const Box<3, T> box{lo, lo + Vec<T, 3>{T(1), T(1), T(1)}};
    const T results[3] = {distance(p, s), distance(p, t), distance(p, box)};
    std::vector<double> out;
    for (const T& r : results) {
        if constexpr (std::is_same_v<T, spatium::Dual<double>>) out.push_back(r.deriv);
        else out.push_back(r);
    }
    return out;
}

inline Entry distance_derivatives() {
    return {"geometry: the derivative of a distance, Dual against differences", 2000, [] {
        return check<std::vector<double>>(
            "geometry: the derivative of a distance, Dual against differences", Kind::Function, 2000,
            [](std::size_t i) {
                std::vector<double> u(40);
                for (std::size_t k = 0; k < 40; ++k) u[k] = -2 + 4 * halton(k < 32 ? i : i + 104729, k % 32);
                return u;
            },
            [](const std::vector<double>& u) {
                std::vector<double> g;
                const Vec<double, 3> p{u[0], u[1], u[2]};
                for (int axis = 0; axis < 3; ++axis) {
                    const auto v = distance_gradients<spatium::Dual<double>>(u, p, axis, 0.0);
                    g.insert(g.end(), v.begin(), v.end());
                }
                return g;
            },
            [](const std::vector<double>& u) {
                std::vector<double> g;
                const Vec<double, 3> p{u[0], u[1], u[2]};
                const double h = 1e-6;
                for (int axis = 0; axis < 3; ++axis) {
                    const auto hi = distance_gradients<double>(u, p, axis, +h);
                    const auto lo = distance_gradients<double>(u, p, axis, -h);
                    for (std::size_t k = 0; k < hi.size(); ++k) g.push_back((hi[k] - lo[k]) / (2 * h));
                }
                return g;
            },
            [](const std::vector<double>& a, const std::vector<double>& b) {
                double worst = 0.0;
                for (std::size_t k = 0; k < a.size(); ++k) worst = std::max(worst, std::abs(a[k] - b[k]));
                return worst;
            },
            [](const std::vector<double>&) { return 1e-8; });
    }};
}

// ── Polygons: hull, areas of a boolean, clipping and distances, double against Real50 ──
//
// Two convex hulls of six points each in [-1, 1]^2, their areas, the area of their intersection,
// of their difference and of their symmetric difference, a point clipped to a hull and to a box,
// and the distance from the point to a hull and between the hulls. A result that is an error in
// one scalar and a value in the other changes the length of the answer and fails the pair.
//
// First run: `Polygon::measure` called `std::abs` qualified, so none of this compiled on Real50.
template<class T>
std::vector<double> polygon_answers(const std::vector<double>& u) {
    using namespace spatium::geometry;
    std::vector<double> out;
    const auto V = [&](std::size_t k) { return Vec<T, 2>{T(u[2 * k]), T(u[2 * k + 1])}; };
    std::vector<Vec<T, 2>> pa, pb;
    for (std::size_t k = 0; k < 6; ++k) pa.push_back(V(k));
    for (std::size_t k = 6; k < 12; ++k) pb.push_back(V(k));
    const auto ha = convex_hull(pa);
    const auto hb = convex_hull(pb);
    if (!ha || !hb) return out;
    const auto push = [&](const T& x) { out.push_back(spatium::primal_double(x)); };
    push(ha->area());
    push(hb->area());
    const auto ir = intersection_region(*ha, *hb);
    if (ir) push(ir->area()); else out.push_back(-1.0);
    const auto da = difference_area<2, T>(*ha, *hb);
    if (da) push(*da); else out.push_back(-1.0);
    const auto sa = symmetric_difference_area<2, T>(*ha, *hb);
    if (sa) push(*sa); else out.push_back(-1.0);
    const Vec<T, 2> p = V(12);
    const auto c = clip(p, *ha);
    if (c) { push((*c)[0]); push((*c)[1]); } else out.push_back(-1.0);
    const Box<2, T> box{Vec<T, 2>{T(-0.5), T(-0.5)}, Vec<T, 2>{T(0.5), T(0.5)}};
    const auto cb = clip(p, box);
    if (cb) { push((*cb)[0]); push((*cb)[1]); } else out.push_back(-1.0);
    push(distance(p, *ha));
    push(distance(*ha, *hb));
    return out;
}

inline Entry polygon_precision() {
    return {"geometry: polygons, double against Real50", 3000, [] {
        return check<std::vector<double>>(
            "geometry: polygons, double against Real50", Kind::Function, 3000,
            [](std::size_t i) {
                std::vector<double> u(26);
                for (std::size_t k = 0; k < 26; ++k) u[k] = -1 + 2 * halton(i, k);
                return u;
            },
            polygon_answers<spatium::Real50>, polygon_answers<double>,
            [](const std::vector<double>& a, const std::vector<double>& b) {
                if (a.size() != b.size()) return std::numeric_limits<double>::infinity();
                double worst = 0.0;
                for (std::size_t k = 0; k < a.size(); ++k) worst = std::max(worst, std::abs(a[k] - b[k]));
                return worst;
            },
            [](const std::vector<double>&) { return 1e-9; });
    }};
}

// ── The force of the contact barrier is minus the gradient of its energy ──
//
// A point a distance d in [0.01, 0.3) off a unit sphere and off a torus (R = 2, r = 0.5), inside
// the barrier's band (d_hat = 0.3). One path is `ipc_contact_force`, the closed-form force along the
// normal of the contact query; the other is the energy `ipc_contact_energy` differentiated through
// `Dual` with respect to the point, negated. Two paths that share the barrier's formula and nothing
// else: the query's normal must be the gradient of its distance, and `ipc_barrier_grad` the
// derivative of `ipc_barrier`. The bound is 1e-12 of the force; the worst measured is 4e-16.
//
// First run: `point_to_torus<Dual>` did not compile (a qualified std::sqrt).
struct ContactDraw { Vec<double, 3> on_sphere, on_torus; };

inline ContactDraw contact_draw(std::size_t i) {
    const double pi = std::numbers::pi;
    const double d = 0.01 + 0.29 * halton(i, 2);
    const double u = 2 * pi * halton(i, 0), v = 2 * pi * halton(i, 1);
    const double rho = 2.0 + (0.5 + d) * std::cos(v);
    const double w = 2 * halton(i, 3) - 1, ph = 2 * pi * halton(i, 4), q = std::sqrt(1 - w * w);
    return {Vec<double, 3>{(1 + d) * q * std::cos(ph), (1 + d) * q * std::sin(ph), (1 + d) * w},
            Vec<double, 3>{rho * std::cos(u), rho * std::sin(u), (0.5 + d) * std::sin(v)}};
}

inline std::vector<double> contact_forces(const ContactDraw& s) {
    namespace m = spatium::physics::mechanics;
    const Vec<double, 3> origin{};
    const auto fs = m::ipc_contact_force(m::point_to_sphere(s.on_sphere, origin, 1.0), 0.3, 1.0);
    const auto ft = m::ipc_contact_force(m::point_to_torus(s.on_torus, spatium::geometry::Torus<double>{.major_radius = 2.0, .minor_radius = 0.5}), 0.3, 1.0);
    return {fs[0], fs[1], fs[2], ft[0], ft[1], ft[2]};
}

inline std::vector<double> contact_forces_from_energy(const ContactDraw& s) {
    namespace m = spatium::physics::mechanics;
    using D = spatium::Dual<double>;
    const Vec<D, 3> origin{D(0.0), D(0.0), D(0.0)};
    const spatium::geometry::Torus<D> torus{.major_radius = D(2.0), .minor_radius = D(0.5)};
    const auto gs = spatium::gradient([&](const Vec<D, 3>& p) {
        return m::ipc_contact_energy(m::point_to_sphere(p, origin, D(1.0)), D(0.3), D(1.0)); }, s.on_sphere);
    const auto gt = spatium::gradient([&](const Vec<D, 3>& p) {
        return m::ipc_contact_energy(m::point_to_torus(p, torus), D(0.3), D(1.0)); }, s.on_torus);
    return {-gs[0], -gs[1], -gs[2], -gt[0], -gt[1], -gt[2]};
}

inline Entry contact_force_is_minus_energy_gradient() {
    return {"contact: the barrier's force is minus the gradient of its energy", 3000, [] {
        return check<ContactDraw>(
            "contact: the barrier's force is minus the gradient of its energy", Kind::Function, 3000, contact_draw,
            contact_forces, contact_forces_from_energy,
            [](const std::vector<double>& a, const std::vector<double>& b) {
                double worst = 0.0, scale = 1.0;
                for (std::size_t k = 0; k < a.size(); ++k) { worst = std::max(worst, std::abs(a[k] - b[k])); scale = std::max(scale, std::abs(a[k])); }
                return worst / scale;
            },
            [](const ContactDraw&) { return 1e-12; });
    }};
}

inline std::vector<Entry> all_pairs() {
    return {ray_torus_precision(),         ray_torus_against_chart(), spd3_eigenvalues_precision(),
            geometry_precision(),          polygon_precision(),  distance_derivatives(), contact_force_is_minus_energy_gradient(),    gradient_dual_against_differences(),
            sphere_exp_log()};
}

}  // namespace symmetry
