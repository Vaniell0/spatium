// Discrete-Lagrangian variational integrator on flat Euclidean space. The
// separable midpoint discretisation reduces to Stormer-Verlet
// leapfrog; we verify the symplectic property (bounded energy oscillation)
// and momentum conservation under translation symmetry.
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <spatium/physics/mechanics/variational.hpp>
#include <cmath>
#include <numbers>

using namespace spatium;
using namespace spatium::physics::mechanics;
using Catch::Approx;

TEST_CASE("Variational integrator: SHO bounded energy", "[variational][sho]") {
    // m·ẍ = -k·x → V(x) = ½kx². Initial: x=1, v=0. ω = √(k/m) = 1.
    constexpr double m = 1.0, k = 1.0;
    constexpr double T_period = 2.0 * std::numbers::pi;
    constexpr double dt = T_period / 200.0;

    PointMass<1, double> body{m, Vec<double, 1>{1.0}, Vec<double, 1>{0.0}};
    auto grad_V = [k](const Vec<double, 1>& q) { return Vec<double, 1>{k * q[0]}; };

    const double E0 = 0.5 * k * 1.0;
    double max_drift = 0;
    int steps = static_cast<int>(50 * T_period / dt);
    for (int i = 0; i < steps; ++i) {
        variational_step_separable(body, grad_V, dt);
        double Ef = body.kinetic_energy()
                  + 0.5 * k * body.state.position[0] * body.state.position[0];
        max_drift = std::max(max_drift, std::abs(Ef - E0) / E0);
    }
    // Energy is bounded — symplectic by construction. Magnitude ~Verlet (2nd order).
    REQUIRE(max_drift < 5e-3);
}

TEST_CASE("Variational integrator: free particle conserves momentum exactly",
          "[variational][noether]") {
    // V ≡ 0 → translation symmetry → momentum p = m·v is exactly conserved
    // by the discrete Noether theorem (Marsden-West).
    PointMass<3, double> body{2.0, Vec<double, 3>{0, 0, 0}, Vec<double, 3>{1.5, -0.5, 0.3}};
    auto grad_V_zero = [](const Vec<double, 3>&) { return Vec<double, 3>{}; };

    const auto p0 = body.momentum();
    for (int i = 0; i < 10000; ++i)
        variational_step_separable(body, grad_V_zero, 0.01);

    auto p1 = body.momentum();
    // Free particle: p stays exact to machine precision.
    REQUIRE(std::abs(p1[0] - p0[0]) < 1e-12);
    REQUIRE(std::abs(p1[1] - p0[1]) < 1e-12);
    REQUIRE(std::abs(p1[2] - p0[2]) < 1e-12);

    // After t = 100 with constant velocity, position must equal v*t.
    REQUIRE(std::abs(body.state.position[0] - 0.0 - 1.5 * 100.0) < 1e-9);
}

TEST_CASE("Variational integrator: 2D circular orbit holds shape",
          "[variational][kepler]") {
    // Inverse-square gravity, circular orbit at r=1, v_circ=1.
    constexpr double GM = 1.0;
    auto grad_V = [GM](const Vec<double, 2>& q) {
        double r2 = q.dot(q);
        double r3 = r2 * std::sqrt(r2);
        return Vec<double, 2>{q * (GM / r3)};
    };
    PointMass<2, double> body{1.0, Vec<double, 2>{1.0, 0.0}, Vec<double, 2>{0.0, 1.0}};

    constexpr double T_period = 2.0 * std::numbers::pi;
    constexpr double dt = T_period / 1000.0;
    const double r0 = 1.0;
    double max_radial_err = 0;
    int steps = static_cast<int>(10 * T_period / dt);
    for (int i = 0; i < steps; ++i) {
        variational_step_separable(body, grad_V, dt);
        double r = std::sqrt(body.state.position.dot(body.state.position));
        max_radial_err = std::max(max_radial_err, std::abs(r - r0));
    }
    // Radius stays near 1 over 10 orbits (verlet-class accuracy).
    REQUIRE(max_radial_err < 1e-2);
}

TEST_CASE("DiscreteLagrangian concept: SeparableMidpointLagrangian satisfies it",
          "[variational][concept]") {
    auto kinetic   = [](const Vec<double, 1>& v) { return 0.5 * v[0] * v[0]; };
    auto potential = [](const Vec<double, 1>& q) { return 0.5 * q[0] * q[0]; };
    // Honest gradient ∇V(q) = q (= d/dq · ½q²) — used to live as a stored
    // constant; now flows through the GradPotential callable.
    auto grad_v    = [](const Vec<double, 1>& q) { return Vec<double, 1>{q[0]}; };

    auto ld = make_separable_midpoint_lagrangian<1, double>(
        kinetic, potential, grad_v, 1.0);

    static_assert(DiscreteLagrangian<decltype(ld), 1, double>);

    Vec<double, 1> q0{1.0}, q1{1.1};
    auto S  = ld.action(q0, q1, 0.1);
    auto p1 = ld.dL_d_dq1(q0, q1, 0.1);
    auto p2 = ld.dL_d_dq2(q0, q1, 0.1);
    REQUIRE(std::isfinite(S));
    REQUIRE(std::isfinite(p1[0]));
    REQUIRE(std::isfinite(p2[0]));
}
