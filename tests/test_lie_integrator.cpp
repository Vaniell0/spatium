// Lie-group integrator. Validate via torque-free Euler top:
// |L| = |I·ω| should be exactly conserved (rotational kinetic energy too,
// up to integrator order).
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <spatium/algebra/groups/so3.hpp>
#include <spatium/physics/mechanics/lie_integrator.hpp>
#include "helpers/euler_top.hpp"
#include <cmath>

using namespace spatium;
using namespace spatium::algebra;
using namespace spatium::physics::mechanics;
using spatium::tests::EulerTopState;
using spatium::tests::euler_top_step;
using Catch::Approx;

TEST_CASE("Lie-Euler step on SO(3): preserves orthogonality", "[lie][so3]") {
    SO3 group{};
    auto R = group.identity();
    Vec3 omega{0.5, 1.0, -0.3};
    auto field = [&](const SO3::ElementType&) { return omega; };

    for (int i = 0; i < 100; ++i)
        R = lie_euler_step(group, R, 0.01, field);

    // R^T R should still be identity (no drift off the manifold — all updates
    // happen via exp on so(3), which is orthogonality-preserving by construction).
    auto RtR = R.transpose() * R;
    auto I = SO3::ElementType::identity();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            REQUIRE(RtR(i, j) == Approx(I(i, j)).margin(1e-10));
}

TEST_CASE("Lie-midpoint step recovers analytical rotation", "[lie][so3]") {
    SO3 group{};
    auto R = group.identity();
    Vec3 omega_z{0.0, 0.0, 1.0};       // 1 rad/s about z
    auto field = [&](const SO3::ElementType&) { return omega_z; };

    constexpr double dt = 0.01;
    constexpr int steps = 100;          // total t = 1 rad
    for (int i = 0; i < steps; ++i)
        R = lie_midpoint_step(group, R, dt, field);

    // Analytical: R(1) = Rz(1) = [[cos, -sin, 0],[sin, cos, 0],[0,0,1]].
    double c = std::cos(1.0), s = std::sin(1.0);
    REQUIRE(R(0, 0) == Approx(c).margin(1e-9));
    REQUIRE(R(1, 0) == Approx(s).margin(1e-9));
    REQUIRE(R(0, 1) == Approx(-s).margin(1e-9));
    REQUIRE(R(2, 2) == Approx(1.0).margin(1e-9));
}

TEST_CASE("Torque-free Euler top: |L| conserved", "[lie][euler-top]") {
    // Inertia (asymmetric for non-trivial dynamics).
    Vec3 I_diag{1.0, 2.0, 3.0};
    auto inertia     = [&](const Vec3& w) {
        return Vec3{w[0] * I_diag[0], w[1] * I_diag[1], w[2] * I_diag[2]};
    };
    auto inertia_inv = [&](const Vec3& v) {
        return Vec3{v[0] / I_diag[0], v[1] / I_diag[1], v[2] / I_diag[2]};
    };

    EulerTopState s{ SO3::ElementType::identity(), Vec3{1.0, 0.5, 0.2} };

    auto L0_body = inertia(s.omega);
    double L0_mag = std::sqrt(L0_body.dot(L0_body));

    constexpr double dt = 1e-3;
    constexpr int steps = 5000;          // 5 s of evolution
    for (int i = 0; i < steps; ++i)
        s = euler_top_step(s, dt, inertia, inertia_inv);

    auto L_body = inertia(s.omega);
    double L_mag = std::sqrt(L_body.dot(L_body));

    // |L| in body frame is conserved exactly by Euler equations; integrator
    // error gives a small drift bounded by the midpoint method's order.
    REQUIRE(std::abs(L_mag - L0_mag) / L0_mag < 1e-3);

    // R should still be orthogonal after thousands of exp updates.
    auto RtR = s.R.transpose() * s.R;
    auto I_mat = SO3::ElementType::identity();
    double orth_err = 0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            orth_err = std::max(orth_err, std::abs(RtR(i, j) - I_mat(i, j)));
    REQUIRE(orth_err < 1e-9);
}
