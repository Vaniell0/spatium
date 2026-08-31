#pragma once

// Point-mass dynamics on a Riemannian manifold.
//
// `PointOnManifold<M>` is a minimal body type that lives on any
// `Manifold` (Sphere, Hyperbolic, ParametricSurface, …). The free-
// particle path is meant to use `M::exp_map` for position and
// parallel transport for velocity, so the body stays on the
// manifold by construction.
//
// What's currently implemented: analytical `geodesic_step` for
// `Sphere<N>` only — closed-form great-circle motion with the
// velocity rotated in the (position, velocity) plane. Hyperbolic,
// ParametricSurface, and arbitrary `Surface`-only manifolds need
// the generic `parallel_transport`-on-`exp_map` path
// (`mesh/transport.hpp` Schild's ladder is the candidate); that
// generalisation is part of the geometric continuum slice and is
// not in this header yet.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <spatium/spaces/sphere.hpp>
#  include <cmath>
#endif

SPATIUM_EXPORT namespace spatium::physics::mechanics {

template<typename M>
    requires Manifold<M>
struct PointOnManifold {
    using Manifold_  = M;
    using Point      = typename M::PointType;
    using Tangent    = typename M::TangentVector;
    using ScalarType = typename M::ScalarType;

    ScalarType mass{ScalarType{1}};
    Point      position{};
    Tangent    velocity{};

    constexpr PointOnManifold() = default;
    constexpr PointOnManifold(ScalarType m, Point p, Tangent v)
        : mass(m), position(std::move(p)), velocity(std::move(v)) {}
};

// ── Analytical free-particle geodesic on Sphere<N> ────────────
//
// For a particle at `p ∈ S^N` (radius r) with tangent velocity `v ⊥ p`,
// free motion is great-circle:
//   γ(t)  = cos(|v|t/r) · p + sin(|v|t/r) · (v·r / |v|)
//   γ'(t) = rotated v in the (p, v) plane  (|γ'| = |v|, preserved)
//
// After the step, `velocity` is still tangent to the new position and has
// the same magnitude — no drift off the manifold, even over many periods.
template<std::size_t N, Scalar T>
void geodesic_step(PointOnManifold<spatium::Sphere<N, T>>& body,
                   const spatium::Sphere<N, T>& S, T dt) {
    using std::sqrt; using std::sin; using std::cos;
    T speed = sqrt(body.velocity.dot(body.velocity));
    if (speed < spatium::epsilon<T>()) return;

    T r      = S.radius;
    T angle  = speed * dt / r;                       // travelled along great circle
    T c      = cos(angle);
    T s      = sin(angle);

    auto new_pos = Vec<T, N + 1>{
        body.position * c + body.velocity * (r * s / speed)};
    // γ'(t) = -(speed/r)·sin(angle)·p + cos(angle)·v, scaled to keep |γ'| = speed.
    auto new_vel = Vec<T, N + 1>{
        body.position * (-speed * s / r) + body.velocity * c};
    body.position = new_pos;
    body.velocity = new_vel;
}

// Speed (= Riemannian norm of velocity) for a body on a sphere.
template<std::size_t N, Scalar T>
T speed_on_sphere(const PointOnManifold<spatium::Sphere<N, T>>& body) {
    using std::sqrt;
    return sqrt(body.velocity.dot(body.velocity));
}

// Check that position and velocity respect the sphere constraints:
//   |position| == r   and   position · velocity == 0.
// Returns max absolute deviation, useful in tests.
template<std::size_t N, Scalar T>
T manifold_constraint_residual(const PointOnManifold<spatium::Sphere<N, T>>& body,
                               const spatium::Sphere<N, T>& S) {
    using std::abs; using std::sqrt;
    T pos_norm  = sqrt(body.position.dot(body.position));
    T pos_err   = abs(pos_norm - S.radius);
    T tan_err   = abs(body.position.dot(body.velocity)) / (S.radius * S.radius);
    return (pos_err > tan_err) ? pos_err : tan_err;
}

} // namespace spatium::physics::mechanics
