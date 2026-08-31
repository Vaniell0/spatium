#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/matrix.hpp>
#  include <spatium/algebra/polynomial.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/error.hpp>
#  include <spatium/geometry/line.hpp>
#  include <algorithm>
#  include <variant>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::geometry {

// ── Quadric surface ───────────────────────────────────────────
// General quadric in 3D: p^T Q p = 0 where p = (x,y,z,1)^T
// and Q is a symmetric 4x4 matrix.

template<Scalar T = double>
struct Quadric {
    Matrix<T, 4, 4> Q;

    // Evaluate: positive outside, negative inside, zero on surface
    T operator()(const Vec<T, 3>& p) const {
        Vec<T, 4> h{p[0], p[1], p[2], T{1}};
        auto Qh = Q * h;
        return h.dot(Qh);
    }

    // Surface normal at point (gradient of quadric form)
    Vec<T, 3> normal(const Vec<T, 3>& p) const {
        Vec<T, 4> h{p[0], p[1], p[2], T{1}};
        auto Qh = Q * h;
        // gradient = 2 * (Q + Q^T) * h, but Q symmetric → 2*Q*h
        Vec<T, 3> grad{Qh[0], Qh[1], Qh[2]};
        auto n = grad.norm();
        return n > epsilon<T>() ? grad / n : Vec<T, 3>{};
    }

    // ── Factories ─────────────────────────────────────────────

    // x² + y² + z² - r² = 0
    static Quadric sphere(T radius) {
        auto Q = Matrix<T, 4, 4>{};
        Q(0, 0) = T{1}; Q(1, 1) = T{1}; Q(2, 2) = T{1};
        Q(3, 3) = -(radius * radius);
        return {Q};
    }

    // x² + y² - r² = 0 (infinite cylinder along Z)
    static Quadric cylinder_z(T radius) {
        auto Q = Matrix<T, 4, 4>{};
        Q(0, 0) = T{1}; Q(1, 1) = T{1};
        Q(3, 3) = -(radius * radius);
        return {Q};
    }

    // x² + y² - z² = 0 (cone along Z, half-angle 45°)
    static Quadric cone_z() {
        auto Q = Matrix<T, 4, 4>{};
        Q(0, 0) = T{1}; Q(1, 1) = T{1}; Q(2, 2) = T{-1};
        return {Q};
    }

    // x²/a² + y²/b² + z²/c² - 1 = 0
    static Quadric ellipsoid(T a, T b, T c) {
        auto Q = Matrix<T, 4, 4>{};
        Q(0, 0) = T{1} / (a * a);
        Q(1, 1) = T{1} / (b * b);
        Q(2, 2) = T{1} / (c * c);
        Q(3, 3) = T{-1};
        return {Q};
    }
};

// ── Ray-quadric intersection result ──────────────────────────

template<Scalar T>
struct RayHit {
    T t;
    Vec<T, 3> point;
    Vec<T, 3> normal;
};

// Proximity info when ray misses: how close it came
template<Scalar T>
struct RayProximity {
    T closest_t;        // ray parameter of closest approach
    T miss_distance;    // |imaginary part| — measures how far ray misses
    Vec<T, 3> closest_point;
};

// ── Ray-quadric intersection ─────────────────────────────────
// Substitutes p(t) = o + t*d into p^T Q p = 0.
// Gets at² + bt + c = 0 where:
//   a = d^T Q3 d       (Q3 = upper-left 3x3 of Q)
//   b = 2*(d^T Q3 o + q3^T d)   where q3 = column 3 of Q, rows 0-2
//   c = o^T Q3 o + 2*q3^T o + Q(3,3)

// internal — do not use, no API stability. The QuadricCoeffs POD
// and quadric_coeffs() helper are an implementation detail shared
// between ray_quadric() and ray_quadric_proximity(); call those
// instead and treat anything in this namespace as private.
namespace detail {

template<Scalar T>
struct QuadricCoeffs { T a, b, c; };

// Build the (a, b, c) of the at² + bt + c = 0 substitution.
// Shared by `ray_quadric` and `ray_quadric_proximity` so a single
// inlining covers the hit-or-miss decision.
template<Scalar T>
inline QuadricCoeffs<T> quadric_coeffs(const Ray<3, T>& ray, const Quadric<T>& q) {
    auto& o = ray.origin;
    auto& d = ray.direction;
    auto& Q = q.Q;

    Vec<T, 3> Qd{Q(0,0)*d[0] + Q(0,1)*d[1] + Q(0,2)*d[2],
                  Q(1,0)*d[0] + Q(1,1)*d[1] + Q(1,2)*d[2],
                  Q(2,0)*d[0] + Q(2,1)*d[1] + Q(2,2)*d[2]};

    Vec<T, 4> o_h{o[0], o[1], o[2], T{1}};
    auto Qo_h = Q * o_h;
    Vec<T, 3> Qo3{Qo_h[0], Qo_h[1], Qo_h[2]};

    return { d.dot(Qd), T{2} * d.dot(Qo3), o_h.dot(Qo_h) };
}

} // namespace detail

template<Scalar T>
std::vector<RayHit<T>> ray_quadric(const Ray<3, T>& ray, const Quadric<T>& q) {
    auto [a_coeff, b_coeff, c_coeff] = detail::quadric_coeffs(ray, q);
    auto roots = solve_quadratic(a_coeff, b_coeff, c_coeff);

    std::vector<RayHit<T>> hits;
    for (auto& root : roots) {
        if (root.is_real() && root.re >= T{0}) {
            auto pt = ray.origin + ray.direction * root.re;
            hits.push_back({root.re, pt, q.normal(pt)});
        }
    }

    std::sort(hits.begin(), hits.end(), [](auto& a, auto& b) { return a.t < b.t; });
    return hits;
}

// Proximity: when ray misses, complex roots tell how close it came.
// |imaginary part| of the roots ∝ miss distance.
template<Scalar T>
Result<RayProximity<T>> ray_quadric_proximity(const Ray<3, T>& ray, const Quadric<T>& q) {
    auto [a_coeff, b_coeff, c_coeff] = detail::quadric_coeffs(ray, q);
    auto roots = solve_quadratic(a_coeff, b_coeff, c_coeff);

    // If real roots exist, ray actually hits — no proximity needed
    for (auto& r : roots) {
        if (r.is_real() && r.re >= T{0})
            return std::unexpected(Error{ErrorCode::DegenerateInput, "ray intersects surface"});
    }

    // Complex roots: real part = closest approach t, |imag| = miss metric
    auto closest_t = std::max(roots[0].re, T{0});
    auto miss = std::abs(roots[0].im);
    auto closest_pt = ray.origin + ray.direction * closest_t;

    return RayProximity<T>{closest_t, miss, closest_pt};
}

// Unified: returns either hits or proximity
template<Scalar T>
std::variant<std::vector<RayHit<T>>, RayProximity<T>>
ray_quadric_full(const Ray<3, T>& ray, const Quadric<T>& q) {
    auto hits = ray_quadric(ray, q);
    if (!hits.empty()) return hits;

    auto prox = ray_quadric_proximity(ray, q);
    if (prox) return *prox;

    // Fallback: no hits, no proximity (degenerate case)
    return std::vector<RayHit<T>>{};
}

// ── Torus ─────────────────────────────────────────────────────
// Quartic surface: (|p|² + R² − r²)² − 4R²(p_x² + p_y²) = 0
// with axis along +Z in the local frame, center at origin.
// General torus is given by (center, axis); ray is transformed to local frame.

template<Scalar T = double>
struct Torus {
    Vec<T, 3> center{};
    Vec<T, 3> axis{T{0}, T{0}, T{1}};   // unit vector along tube axis
    T major_radius{T{1}};               // R: centerline circle radius
    T minor_radius{T{0.25}};            // r: tube radius
};

// Orthonormal basis (u, v, w=axis) with w given (unit length).
template<Scalar T>
inline void torus_basis(const Vec<T, 3>& w, Vec<T, 3>& u, Vec<T, 3>& v) {
    Vec<T, 3> helper = std::abs(w[2]) < T{0.9}
        ? Vec<T, 3>{T{0}, T{0}, T{1}}
        : Vec<T, 3>{T{1}, T{0}, T{0}};
    u = w.cross(helper);
    u = u / u.norm();
    v = w.cross(u);
}

// internal — do not use, no API stability. TorusLocalRay and the
// torus_local_ray()/torus_quartic_coeffs() helpers back the
// quartic torus-ray solver in ray_torus(); they may be replaced
// or moved at any time.
namespace detail {

template<Scalar T>
struct TorusLocalRay {
    Vec<T, 3> u, v;          // basis vectors (axis = w stored in torus.axis)
    Vec<T, 3> o, d;          // ray re-expressed in the local frame
};

template<Scalar T>
inline TorusLocalRay<T> torus_local_ray(const Ray<3, T>& ray, const Torus<T>& torus) {
    Vec<T, 3> w = torus.axis;
    Vec<T, 3> u, v;
    torus_basis(w, u, v);
    Vec<T, 3> delta = ray.origin - torus.center;
    return { u, v,
             { delta.dot(u),         delta.dot(v),         delta.dot(w) },
             { ray.direction.dot(u), ray.direction.dot(v), ray.direction.dot(w) } };
}

template<Scalar T>
struct TorusQuartic { T c4, c3, c2, c1, c0; };

// Coefficients of the quartic (|p|² + R² − r²)² − 4R²(x² + y²) = 0
// after substituting p(t) = o + t·d in the torus' local frame.
template<Scalar T>
inline TorusQuartic<T> torus_quartic_coeffs(const Vec<T, 3>& o, const Vec<T, 3>& d,
                                            T R, T r) {
    T alpha  = d.dot(d);
    T beta   = T{2} * o.dot(d);
    T gamma  = o.dot(o);
    T alpha2 = d[0]*d[0] + d[1]*d[1];
    T beta2  = T{2} * (o[0]*d[0] + o[1]*d[1]);
    T gamma2 = o[0]*o[0] + o[1]*o[1];
    T mu     = gamma + R*R - r*r;
    return { alpha * alpha,
             T{2} * alpha * beta,
             T{2} * alpha * mu + beta * beta - T{4} * R * R * alpha2,
             T{2} * beta * mu - T{4} * R * R * beta2,
             mu * mu - T{4} * R * R * gamma2 };
}

} // namespace detail

// Ray-torus hits, sorted by t ascending. Non-unit ray.direction supported.
template<Scalar T>
std::vector<RayHit<T>> ray_torus(const Ray<3, T>& ray, const Torus<T>& torus) {
    auto loc = detail::torus_local_ray(ray, torus);
    Vec<T, 3> w = torus.axis;
    Vec<T, 3> u = loc.u, v = loc.v, o = loc.o, d = loc.d;
    T R = torus.major_radius;
    T r = torus.minor_radius;

    auto k = detail::torus_quartic_coeffs(o, d, R, r);
    auto roots = solve_quartic(k.c4, k.c3, k.c2, k.c1, k.c0);

    std::vector<RayHit<T>> hits;
    hits.reserve(4);
    for (auto& root : roots) {
        if (!root.is_real()) continue;
        T t = root.re;
        if (t < T{0}) continue;
        Vec<T, 3> p_local = o + d * t;
        Vec<T, 3> p_world = ray.origin + ray.direction * t;

        // Normal: gradient in local frame, then rotate to world
        T lp2 = p_local.dot(p_local);
        T s = lp2 + R * R - r * r;
        Vec<T, 3> grad_local{
            p_local[0] * (lp2 - R*R - r*r),
            p_local[1] * (lp2 - R*R - r*r),
            p_local[2] * s
        };
        // world_normal = u*gx + v*gy + w*gz
        Vec<T, 3> n = u * grad_local[0] + v * grad_local[1] + w * grad_local[2];
        T nlen = n.norm();
        if (nlen > epsilon<T>()) n = n / nlen;

        hits.push_back({t, p_world, n});
    }

    std::sort(hits.begin(), hits.end(), [](auto& a, auto& b) { return a.t < b.t; });
    return hits;
}

// Proximity for misses via minimum |imaginary part| of quartic roots.
template<Scalar T>
Result<RayProximity<T>> ray_torus_proximity(const Ray<3, T>& ray, const Torus<T>& torus) {
    auto loc = detail::torus_local_ray(ray, torus);
    auto k = detail::torus_quartic_coeffs(loc.o, loc.d,
                                          torus.major_radius, torus.minor_radius);
    auto roots = solve_quartic(k.c4, k.c3, k.c2, k.c1, k.c0);

    for (auto& root : roots) {
        if (root.is_real() && root.re >= T{0})
            return std::unexpected(Error{ErrorCode::DegenerateInput, "ray intersects torus"});
    }

    // Pick the complex pair with smallest |imag| — closest approach
    T best_im = std::numeric_limits<T>::max();
    T best_t  = T{0};
    for (auto& root : roots) {
        T im = std::abs(root.im);
        if (im > T{0} && im < best_im && root.re >= T{0}) {
            best_im = im;
            best_t  = root.re;
        }
    }
    if (best_im == std::numeric_limits<T>::max())
        return std::unexpected(Error{ErrorCode::NoIntersection, "no forward complex pair"});

    return RayProximity<T>{best_t, best_im, ray.origin + ray.direction * best_t};
}

} // namespace spatium::geometry
