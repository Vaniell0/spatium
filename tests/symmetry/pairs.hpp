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
#include <spatium/geometry/make.hpp>
#include <spatium/geometry/ray_parametric.hpp>
#include <spatium/geometry/ray_surface.hpp>
#include <spatium/spaces/parametric.hpp>
#include <spatium/spaces/spd.hpp>
#include <spatium/spaces/sphere.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
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

inline std::vector<Entry> all_pairs() {
    return {ray_torus_precision(),         ray_torus_against_chart(), spd3_eigenvalues_precision(),
            gradient_dual_against_differences(), sphere_exp_log()};
}

}  // namespace symmetry
