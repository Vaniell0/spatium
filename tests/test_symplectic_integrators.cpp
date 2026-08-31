// RKMK4 (commutator-free) on SO(3) outperforming midpoint, and the
// SymplecticManifold concept + CotangentBundle wrapper passing the
// discrete-symplecticity drift check.
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <spatium/algebra/groups/so3.hpp>
#include <spatium/spaces/euclidean.hpp>
#include <spatium/physics/mechanics/lie_integrator.hpp>
#include <spatium/physics/mechanics/symplectic.hpp>
#include <spatium/physics/mechanics/integrator.hpp>
#include <spatium/physics/mechanics/body.hpp>
#include <spatium/physics/mechanics/force.hpp>
#include "helpers/euler_top.hpp"
#include <cmath>

using namespace spatium;
using namespace spatium::algebra;
using namespace spatium::physics::mechanics;
using spatium::tests::EulerTopState;
using spatium::tests::euler_top_step;
using spatium::tests::euler_top_step_rkmk4;
using Catch::Approx;

// ── RKMK4 (commutator-free) ────────────────────────────────────

TEST_CASE("RKMK4 commutator-free: orthogonality preserved over many steps",
          "[lie][rkmk4]") {
    SO3 group{};
    auto R = group.identity();
    Vec3 omega{0.4, -0.6, 0.8};
    auto field = [&](const SO3::ElementType&) { return omega; };

    for (int i = 0; i < 1000; ++i)
        R = lie_rkmk4_cf_step(group, R, 0.01, field);

    // R^T R should remain identity to ≈ 1e-10 even with 1000 exp's.
    auto RtR = R.transpose() * R;
    auto I = SO3::ElementType::identity();
    double err = 0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            err = std::max(err, std::abs(RtR(i, j) - I(i, j)));
    REQUIRE(err < 1e-10);
}

TEST_CASE("RKMK4 vs midpoint on Euler top: |L| drift drops by ≥ 100×",
          "[lie][rkmk4][euler-top]") {
    Vec3 I_diag{1.0, 2.0, 3.0};
    auto inertia     = [&](const Vec3& w) {
        return Vec3{w[0] * I_diag[0], w[1] * I_diag[1], w[2] * I_diag[2]};
    };
    auto inertia_inv = [&](const Vec3& v) {
        return Vec3{v[0] / I_diag[0], v[1] / I_diag[1], v[2] / I_diag[2]};
    };

    auto run = [&](auto&& stepper) {
        EulerTopState s{ SO3::ElementType::identity(), Vec3{1.0, 0.5, 0.2} };
        Vec3 L0_body = inertia(s.omega);
        double L0_mag = std::sqrt(L0_body.dot(L0_body));

        for (int i = 0; i < 5000; ++i)
            s = stepper(s, 1e-3, inertia, inertia_inv);

        Vec3 L = inertia(s.omega);
        double L_mag = std::sqrt(L.dot(L));
        return std::abs(L_mag - L0_mag) / L0_mag;
    };

    double drift_mid   = run([](auto s, double dt, auto& I, auto& Iinv) {
        return euler_top_step(s, dt, I, Iinv);
    });
    double drift_rkmk4 = run([](auto s, double dt, auto& I, auto& Iinv) {
        return euler_top_step_rkmk4(s, dt, I, Iinv);
    });

    // Midpoint: drift ~1e-4, RKMK4: drift ~1e-8 — order-of-magnitude separation.
    REQUIRE(drift_rkmk4 < drift_mid);
    REQUIRE(drift_rkmk4 * 100.0 < drift_mid);
}

// ── SymplecticManifold + CotangentBundle ──────────────────────

TEST_CASE("CotangentBundle satisfies SymplecticManifold concept",
          "[symplectic][concept]") {
    using Phase = CotangentBundle<Euclidean<3, double>>;
    static_assert(SymplecticManifold<Phase>);
    REQUIRE(Phase::dimension == 6);
}

TEST_CASE("CotangentBundle: canonical 2-form ω = q1·p2 - q2·p1",
          "[symplectic][omega]") {
    using Phase = CotangentBundle<Euclidean<3, double>>;
    Vec3 dq1{1, 0, 0}, dp1{0, 0, 0};
    Vec3 dq2{0, 0, 0}, dp2{1, 0, 0};
    // ω((e1, 0), (0, e1)) = e1 · e1 - 0 = 1.
    REQUIRE(Phase::omega(dq1, dp1, dq2, dp2) == Approx(1.0));
    // Antisymmetry: ω(v, w) = -ω(w, v).
    REQUIRE(Phase::omega(dq2, dp2, dq1, dp1) == Approx(-1.0));
}

TEST_CASE("Verlet step is symplectic: drift of ω after one step is small",
          "[symplectic][verlet]") {
    // Apply Verlet on a 1D harmonic oscillator (q' = p/m, p' = -kq).
    // Build the step map and probe symplecticity numerically.
    using Phase = CotangentBundle<Euclidean<2, double>>;

    auto step_map = [](Phase::State s, double dt) -> Phase::State {
        // Force = -k·q with k = 1 in each axis (uncoupled 2D oscillators).
        constexpr double k = 1.0, m = 1.0;
        Vec<double, 2> a0 = Vec<double, 2>{s.q * (-k / m)};
        Vec<double, 2> v_half = Vec<double, 2>{s.p / m + a0 * (dt * 0.5)};
        Vec<double, 2> q_new  = Vec<double, 2>{s.q + v_half * dt};
        Vec<double, 2> a1 = Vec<double, 2>{q_new * (-k / m)};
        Vec<double, 2> v_new  = Vec<double, 2>{v_half + a1 * (dt * 0.5)};
        return Phase::State{q_new, Vec<double, 2>{v_new * m}};
    };

    Phase::State s0{Vec<double, 2>{1.0, 0.5}, Vec<double, 2>{0.0, 0.3}};
    double drift = verify_symplecticity_drift<Phase>(s0, step_map, 1e-4, 0.01);
    // For a true symplectic map the drift is O(eps) — bound at machine-precision-ish.
    REQUIRE(drift < 1e-3);
}
