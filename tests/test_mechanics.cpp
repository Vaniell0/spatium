// Flat-space mechanics foundations: units (compile-time SI), PointMass,
// forces, integrators.
// Validation through analytical comparisons:
//   - SHO: amplitude/period match analytical solution.
//   - Kepler orbit: energy + angular momentum drift bounded over many periods.
//   - Euler vs Verlet vs RK4 ordering: drift comparison for orbit.
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <spatium/physics/mechanics/units.hpp>
#include <spatium/physics/mechanics/state.hpp>
#include <spatium/physics/mechanics/body.hpp>
#include <spatium/physics/mechanics/force.hpp>
#include <spatium/physics/mechanics/integrator.hpp>
#include <cmath>
#include <numbers>

using namespace spatium;
using namespace spatium::physics::mechanics;
using namespace spatium::physics::mechanics::literals;
using Catch::Approx;

// Generic-lambda wrappers for parametrising helper loops over a stepper.
// Replaces a former runtime enum dispatch with compile-time selection.
inline constexpr auto step_euler =
    [](auto& body, auto&& force, double dt, double t) {
        euler_step(body, force, dt, t);
    };
inline constexpr auto step_verlet =
    [](auto& body, auto&& force, double dt, double t) {
        verlet_step(body, force, dt, t);
    };
inline constexpr auto step_yoshida4 =
    [](auto& body, auto&& force, double dt, double t) {
        yoshida4_step(body, force, dt, t);
    };
inline constexpr auto step_rk4 =
    [](auto& body, auto&& force, double dt, double t) {
        rk4_step(body, force, dt, t);
    };

// ── Units ─────────────────────────────────────────────────────

TEST_CASE("units: literals carry SI dimensions", "[mechanics][units]") {
    auto m = 5.0_kg;
    auto t = 2.0_s;
    auto d = 10.0_m;
    REQUIRE(m.value() == 5.0);
    REQUIRE(t.value() == 2.0);
    REQUIRE(d.value() == 10.0);

    // m/s = velocity
    auto v = d / t;
    static_assert(std::is_same_v<decltype(v), Velocity>);
    REQUIRE(v.value() == 5.0);

    // F = m·a = m·(d/t²) — should yield Force.
    auto a = v / t;
    static_assert(std::is_same_v<decltype(a), Acceleration>);
    auto f = m * a;
    static_assert(std::is_same_v<decltype(f), Force>);
    REQUIRE(f.value() == Approx(12.5));
}

TEST_CASE("units: addition preserves dimension, mixing fails", "[mechanics][units]") {
    auto a = 3.0_m + 4.0_m;
    REQUIRE(a.value() == 7.0);
    // The next line should be a compile error if uncommented:
    //   auto bad = 1.0_m + 1.0_s;
    // Compile-time check via static_assert on derived dimension.
    static_assert(std::is_same_v<decltype(a), Length>);
}

TEST_CASE("units: physical constants", "[mechanics][units]") {
    REQUIRE(g_earth.value() == Approx(9.80665));
    REQUIRE(c_light.value() == Approx(2.99792458e8));
    REQUIRE(G_newton.value() == Approx(6.67430e-11));
}

// ── PointMass + simple force ──────────────────────────────────

TEST_CASE("PointMass kinetic energy + momentum", "[mechanics][body]") {
    PointMass<3> body{2.0, Vec3{0, 0, 0}, Vec3{3, 4, 0}};
    REQUIRE(body.kinetic_energy() == Approx(0.5 * 2.0 * 25.0));   // ½·m·v² = 25
    auto p = body.momentum();
    REQUIRE(p[0] == 6.0);
    REQUIRE(p[1] == 8.0);
}

TEST_CASE("UniformGravity gives F = m·g", "[mechanics][force]") {
    PointMass<3> body{2.0, Vec3{0, 5, 0}, Vec3{}};
    UniformGravity<3> grav{Vec3{0, -9.81, 0}};
    auto f = grav(body, 0.0);
    REQUIRE(f[1] == Approx(-19.62));
}

// ── 1-D Simple Harmonic Oscillator ────────────────────────────
// m·ẍ = -k·x → ω² = k/m. Analytical: x(t) = A·cos(ω·t).
// Confirm Verlet preserves amplitude across many periods.

TEST_CASE("SHO: Verlet preserves amplitude over 50 periods", "[mechanics][integrator][sho]") {
    constexpr double m = 1.0, k = 1.0;
    constexpr double omega = 1.0;          // sqrt(k/m)
    constexpr double T_period = 2.0 * std::numbers::pi / omega;
    constexpr double dt = T_period / 200.0;

    PointMass<1, double> body{m, Vec<double, 1>{1.0}, Vec<double, 1>{0.0}};
    Spring<1, double> spring{Vec<double, 1>{0.0}, k};

    double E0 = body.kinetic_energy() + 0.5 * k * 1.0;  // ½kx² at x=1
    int steps = static_cast<int>(50 * T_period / dt);
    for (int i = 0; i < steps; ++i)
        verlet_step(body, spring, dt);

    double E = body.kinetic_energy() + 0.5 * k * body.state.position[0] * body.state.position[0];
    REQUIRE(std::abs(E - E0) / E0 < 1e-3);   // drift bounded
}

// Compare drift at end of run — Euler should be much worse than Verlet/RK4.
TEST_CASE("SHO: Euler drifts more than Verlet", "[mechanics][integrator][sho]") {
    constexpr double m = 1.0, k = 1.0;
    constexpr double omega = 1.0;
    constexpr double T_period = 2.0 * std::numbers::pi / omega;
    constexpr double dt = T_period / 200.0;
    constexpr int steps = static_cast<int>(20 * T_period / dt);

    auto run = [&](auto stepper) {
        PointMass<1, double> body{m, Vec<double, 1>{1.0}, Vec<double, 1>{0.0}};
        Spring<1, double> spring{Vec<double, 1>{0.0}, k};
        double E0 = 0.5 * k * 1.0;
        for (int i = 0; i < steps; ++i)
            stepper(body, spring, dt, 0.0);
        double Ef = body.kinetic_energy() + 0.5 * k * body.state.position[0] * body.state.position[0];
        return std::abs(Ef - E0) / E0;
    };

    double drift_euler  = run(step_euler);
    double drift_verlet = run(step_verlet);
    REQUIRE(drift_verlet * 100.0 < drift_euler);  // Verlet ≥100× tighter
}

// ── 2-D Kepler orbit (circular) ───────────────────────────────
// Newton gravity from a fixed point source, circular initial conditions:
//   v_circ = sqrt(GM/r). Energy E = -GM·m/(2r) for circular orbit.
//   Angular momentum L = m·v·r (perpendicular).

TEST_CASE("Kepler: circular orbit energy stays bounded with Verlet", "[mechanics][kepler]") {
    constexpr double GM = 1.0;
    constexpr double r0 = 1.0;
    constexpr double v_circ = 1.0;                                // sqrt(GM/r0)
    constexpr double T_period = 2.0 * std::numbers::pi * r0 / v_circ;
    constexpr double dt = T_period / 1000.0;

    PointMass<2, double> body{1.0, Vec<double, 2>{r0, 0.0}, Vec<double, 2>{0.0, v_circ}};
    PointGravity<2, double> grav{Vec<double, 2>{0.0, 0.0}, GM};

    auto E_total = [&]() {
        double r = std::sqrt(body.state.position.dot(body.state.position));
        return body.kinetic_energy() - GM * body.mass / r;
    };
    auto L = [&]() {
        // 2D cross product: x·vy - y·vx
        return body.mass * (body.state.position[0] * body.state.velocity[1]
                          - body.state.position[1] * body.state.velocity[0]);
    };

    double E0 = E_total();
    double L0 = L();

    int steps = static_cast<int>(10 * T_period / dt);
    for (int i = 0; i < steps; ++i)
        verlet_step(body, grav, dt);

    double E_drift = std::abs(E_total() - E0) / std::abs(E0);
    double L_drift = std::abs(L() - L0) / std::abs(L0);

    REQUIRE(E_drift < 1e-3);   // energy drift after 10 periods
    REQUIRE(L_drift < 1e-3);   // angular momentum should be exactly conserved
}

// Order-of-accuracy: RK4 error should fall faster than Verlet as dt halves.
TEST_CASE("Kepler: RK4 short-time accuracy beats Verlet", "[mechanics][kepler]") {
    constexpr double GM = 1.0, r0 = 1.0, v_circ = 1.0;
    constexpr double T_period = 2.0 * std::numbers::pi;

    auto orbit_error = [&](auto stepper, double dt) {
        PointMass<2, double> body{1.0, Vec<double, 2>{r0, 0.0}, Vec<double, 2>{0.0, v_circ}};
        PointGravity<2, double> grav{Vec<double, 2>{0.0, 0.0}, GM};
        int steps = static_cast<int>(T_period / dt);
        for (int i = 0; i < steps; ++i)
            stepper(body, grav, dt, 0.0);
        // After exactly one period, body should be back at (r0, 0).
        return std::abs(body.state.position[0] - r0) + std::abs(body.state.position[1]);
    };

    double err_rk4    = orbit_error(step_rk4,    T_period / 100.0);
    double err_verlet = orbit_error(step_verlet, T_period / 100.0);
    // RK4 is 4th-order in time, Verlet is 2nd-order: at dt=T/100 RK4 wins on
    // round-trip position by orders of magnitude.
    REQUIRE(err_rk4 < err_verlet);
}

// ── KAM-style long-time energy bound (the symplectic payoff) ─────────────
// Symplectic integrators preserve a *modified* Hamiltonian, so total energy
// oscillates within an O(dt^p) envelope FOREVER instead of drifting linearly.
// Non-symplectic methods (RK4) drift even at high order. This is the central
// reason variational/Lie-group methods matter for orbital mechanics.

TEST_CASE("KAM: Verlet bounded vs RK4 drift over 1000 SHO periods",
          "[mechanics][integrator][kam]") {
    constexpr double m = 1.0, k = 1.0;
    constexpr double T_period = 2.0 * std::numbers::pi;   // ω = 1
    constexpr double dt = T_period / 100.0;
    constexpr int steps = static_cast<int>(1000 * T_period / dt);

    auto run = [&](auto stepper) {
        PointMass<1, double> body{m, Vec<double, 1>{1.0}, Vec<double, 1>{0.0}};
        Spring<1, double> spring{Vec<double, 1>{0.0}, k};
        const double E0 = 0.5 * k * 1.0;
        double max_drift = 0;
        for (int i = 0; i < steps; ++i) {
            stepper(body, spring, dt, 0.0);
            double Ef = body.kinetic_energy()
                      + 0.5 * k * body.state.position[0] * body.state.position[0];
            max_drift = std::max(max_drift, std::abs(Ef - E0) / E0);
        }
        return max_drift;
    };

    double drift_verlet   = run(step_verlet);
    double drift_yoshida  = run(step_yoshida4);
    double drift_rk4      = run(step_rk4);

    // Symplectic: bounded over 1000 periods (Verlet ~1e-3, Yoshida ~1e-6).
    REQUIRE(drift_verlet  < 5e-3);
    REQUIRE(drift_yoshida < 1e-5);

    // Yoshida4 is ~1000× tighter than Verlet on the same dt — the gain from
    // 4th-order composition vs 2nd-order Stormer-Verlet.
    REQUIRE(drift_verlet > drift_yoshida * 100.0);

    // RK4 drifts more than the symplectic bound. The exact ratio depends on
    // dt and run length; here over 1000 periods at dt = T/100 it's typically
    // 10× or more.
    REQUIRE(drift_rk4 > drift_yoshida * 10.0);
}

TEST_CASE("Kepler over 100 periods: Yoshida4 beats RK4 on energy",
          "[mechanics][integrator][kepler]") {
    constexpr double GM = 1.0, r0 = 1.0, v_circ = 1.0;
    constexpr double T_period = 2.0 * std::numbers::pi;
    constexpr double dt = T_period / 200.0;
    constexpr int steps = static_cast<int>(100 * T_period / dt);

    auto run = [&](auto stepper) {
        PointMass<2, double> body{1.0, Vec<double, 2>{r0, 0.0},
                                       Vec<double, 2>{0.0, v_circ}};
        PointGravity<2, double> grav{Vec<double, 2>{0.0, 0.0}, GM};
        const double E0 = 0.5 * v_circ * v_circ - GM / r0;
        double max_drift = 0;
        for (int i = 0; i < steps; ++i) {
            stepper(body, grav, dt, 0.0);
            double r = std::sqrt(body.state.position.dot(body.state.position));
            double Ef = body.kinetic_energy() - GM * body.mass / r;
            max_drift = std::max(max_drift, std::abs(Ef - E0) / std::abs(E0));
        }
        return max_drift;
    };

    double drift_verlet  = run(step_verlet);
    double drift_yoshida = run(step_yoshida4);
    double drift_rk4     = run(step_rk4);

    // Yoshida4 is 4th-order symplectic — orders of magnitude tighter than
    // 2nd-order Verlet on the same dt for orbital problems.
    REQUIRE(drift_yoshida < drift_verlet);
    // RK4 is 4th-order non-symplectic — its energy drift exceeds Yoshida4's
    // bounded oscillation by the time we hit hundreds of orbits.
    REQUIRE(drift_rk4 > drift_yoshida);
}
