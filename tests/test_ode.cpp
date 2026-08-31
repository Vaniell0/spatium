#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/algebra/ode.hpp>
#include <cmath>

using namespace spatium;
using Catch::Matchers::WithinAbs;

TEST_CASE("euler_step matches exponential decay at a small step size", "[ode]") {
    double k = 1.5;
    auto f = [k](double, Vec<double, 1> y) { return Vec<double, 1>{-k * y[0]}; };
    Vec<double, 1> y0{1.0};
    double dt = 1e-4;
    int steps = 10000; // t = 1.0
    Vec<double, 1> y = integrate_fixed(f, 0.0, y0, dt, steps,
                                        [](auto&& f_, double t, const Vec<double, 1>& y_, double dt_) {
                                            return euler_step(f_, t, y_, dt_);
                                        });
    double expected = std::exp(-k * 1.0);
    CHECK_THAT(y[0], WithinAbs(expected, 1e-3)); // Euler is only O(dt) accurate
}

TEST_CASE("rk4_step matches exponential decay far more tightly than Euler at the same step size",
          "[ode]") {
    double k = 1.5;
    auto f = [k](double, Vec<double, 1> y) { return Vec<double, 1>{-k * y[0]}; };
    Vec<double, 1> y0{1.0};
    double dt = 1e-2; // coarser step than the Euler test above
    int steps = 100;  // t = 1.0
    Vec<double, 1> y = integrate_fixed(f, 0.0, y0, dt, steps,
                                        [](auto&& f_, double t, const Vec<double, 1>& y_, double dt_) {
                                            return rk4_step(f_, t, y_, dt_);
                                        });
    double expected = std::exp(-k * 1.0);
    CHECK_THAT(y[0], WithinAbs(expected, 1e-8)); // O(dt^4): far tighter at a coarser step
}

TEST_CASE("rk4_step matches the harmonic oscillator's closed form", "[ode]") {
    // State = [position, velocity]; y' = [v, -w^2*x].
    double w = 2.0;
    auto f = [w](double, Vec<double, 2> y) { return Vec<double, 2>{y[1], -w * w * y[0]}; };
    Vec<double, 2> y0{1.0, 0.0};
    double dt = 1e-3;
    int steps = 1000; // t = 1.0
    Vec<double, 2> y = integrate_fixed(f, 0.0, y0, dt, steps,
                                        [](auto&& f_, double t, const Vec<double, 2>& y_, double dt_) {
                                            return rk4_step(f_, t, y_, dt_);
                                        });
    CHECK_THAT(y[0], WithinAbs(std::cos(w * 1.0), 1e-6));
    CHECK_THAT(y[1], WithinAbs(-w * std::sin(w * 1.0), 1e-6));
}

TEST_CASE("Euler's instability past its stability boundary is real, not just slow convergence",
          "[ode]") {
    // For y'=-k*y, explicit Euler is only bounded while k*dt <= 2 -- past
    // that the numerical solution diverges even though the true solution
    // decays to zero. This is the exact signal rsc/include/ode_task.hpp
    // dispatches on, checked here directly against Spatium's own stepper,
    // not assumed from numerical-analysis theory alone.
    // Each step multiplies by (1-k*dt) = -1.5, so magnitude grows as
    // 1.5^n -- geometric, not runaway-fast, so enough steps are needed to
    // see it clearly (1.5^20 ~ 3325, 1.5^40 ~ 1.1e7).
    double k = 1.0;
    double dt = 2.5; // k*dt = 2.5 > 2: past Euler's stability boundary
    auto f = [k](double, Vec<double, 1> y) { return Vec<double, 1>{-k * y[0]}; };
    Vec<double, 1> y = Vec<double, 1>{1.0};
    for (int i = 0; i < 40; ++i) y = euler_step(f, 0.0, y, dt);
    CHECK(std::abs(y[0]) > 1e6); // diverged, while exp(-k*t) -> 0
}
