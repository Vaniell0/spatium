#pragma once

// Test-only Euler top stepper. Drives `lie_midpoint_step` /
// `lie_rkmk4_cf_step` on SO(3) with a torque-free body angular
// velocity update so the test suite can compare orientation
// integrators on the canonical asymmetric-rigid-body benchmark.
//
// Lives in tests/helpers/ rather than the public header set
// because its inertia callable + ad-hoc midpoint/RK4 update for ω
// is a benchmark fixture, not a production-quality top integrator.
// LGVI in physics/mechanics/lgvi.hpp is the canonical rigid-body
// path; that one is exact-symmetric on the group, momentum-
// conserving by discrete Noether, and is what the engine ships.

#include <spatium/algebra/groups/so3.hpp>
#include <spatium/physics/mechanics/lie_integrator.hpp>

namespace spatium::tests {

struct EulerTopState {
    spatium::algebra::SO3<double>::ElementType R;     // 3×3 orientation matrix
    spatium::algebra::SO3<double>::AlgebraType  omega; // body angular velocity (Vec3)
};

// One time-step. `inertia_inv(L)` returns ω given body angular momentum L = I·ω.
// `inertia(omega)` returns L. Both supplied as callables so the inertia tensor
// is private to the caller.
template<typename InertiaFn, typename InertiaInvFn>
inline EulerTopState euler_top_step(EulerTopState s, double dt,
                                    InertiaFn&& inertia, InertiaInvFn&& inertia_inv)
{
    using Vec3 = spatium::algebra::SO3<double>::AlgebraType;
    using SO3 = spatium::algebra::SO3<double>;

    SO3 group{};
    auto omega_field = [&](const SO3::ElementType& /*R*/) -> Vec3 {
        return s.omega;
    };
    s.R = spatium::physics::mechanics::lie_midpoint_step(
        group, s.R, dt, omega_field);

    auto domega = [&](const Vec3& w) -> Vec3 {
        Vec3 L = inertia(w);
        Vec3 cross = w.cross(L);
        return Vec3{inertia_inv(cross) * (-1.0)};
    };
    Vec3 k1 = domega(s.omega);
    Vec3 w_mid = Vec3{s.omega + k1 * (dt * 0.5)};
    Vec3 k2 = domega(w_mid);
    s.omega = Vec3{s.omega + k2 * dt};
    return s;
}

// RKMK4 variant: commutator-free order-4 for orientation, classical RK4
// for the (flat) so(3) ω update.
template<typename InertiaFn, typename InertiaInvFn>
inline EulerTopState euler_top_step_rkmk4(EulerTopState s, double dt,
                                          InertiaFn&& inertia, InertiaInvFn&& inertia_inv)
{
    using Vec3 = spatium::algebra::SO3<double>::AlgebraType;
    using SO3 = spatium::algebra::SO3<double>;

    SO3 group{};
    auto omega_field = [&](const SO3::ElementType& /*R*/) -> Vec3 {
        return s.omega;
    };
    s.R = spatium::physics::mechanics::lie_rkmk4_cf_step(
        group, s.R, dt, omega_field);

    auto domega = [&](const Vec3& w) -> Vec3 {
        Vec3 L = inertia(w);
        Vec3 cross = w.cross(L);
        return Vec3{inertia_inv(cross) * (-1.0)};
    };
    Vec3 k1 = domega(s.omega);
    Vec3 k2 = domega(Vec3{s.omega + k1 * (dt * 0.5)});
    Vec3 k3 = domega(Vec3{s.omega + k2 * (dt * 0.5)});
    Vec3 k4 = domega(Vec3{s.omega + k3 *  dt});
    s.omega = Vec3{s.omega + (k1 + (k2 + k3) * 2.0 + k4) * (dt / 6.0)};
    return s;
}

} // namespace spatium::tests
