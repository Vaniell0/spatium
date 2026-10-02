#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/algebra/ode.hpp>
#include <spatium/algebra/dual.hpp>
#include <spatium/core/precision.hpp>
#include <algorithm>
#include <cmath>
#include <utility>

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

// ── Extrapolation: the tolerance follows the scalar ──────────

namespace {
// y'' = -y from (1, 0): y(t) = cos t, y'(t) = -sin t. The error against
// the closed form, and the right-hand side evaluations it cost.
template<class T>
std::pair<double, int> oscillator_error(double t_end, double tol) {
    int calls = 0;
    auto f = [&calls](T, const Vec<T, 2>& y) { ++calls; return Vec<T, 2>{y[1], -y[0]}; };
    const auto r = integrate_extrapolated(f, T(0), Vec<T, 2>{T(1), T(0)}, T(t_end), tol);
    REQUIRE(r.has_value());
    using std::cos; using std::sin;
    const double e0 = primal_double(T((*r)[0] - cos(T(t_end))));
    const double e1 = primal_double(T((*r)[1] + sin(T(t_end))));
    return {std::max(std::abs(e0), std::abs(e1)), calls};
}
}  // namespace

TEST_CASE("integrate_extrapolated meets a tolerance that follows the scalar", "[ode]") {
    // A tolerance a power of the arithmetic's own epsilon, and the answer
    // within a small factor of it: float, double, long double alike.
    for (const double tol : {1e-5, 1e-8, 1e-11}) {
        CHECK(oscillator_error<double>(3.0, tol).first < 20 * tol);
    }
    CHECK(oscillator_error<float>(3.0, std::pow(machine_epsilon<float>(), 2.0 / 3.0)).first < 1e-3);
    CHECK(oscillator_error<long double>(3.0, std::pow(machine_epsilon<long double>(), 2.0 / 3.0)).first < 1e-11);
}

TEST_CASE("integrate_extrapolated beats RK4 at the same cost on a smooth problem", "[ode]") {
    // The comparison that justifies it: to 1e-11, extrapolation spends ~140
    // right-hand-side evaluations where RK4 needs thousands of steps.
    const auto [err, calls] = oscillator_error<double>(3.0, 1e-11);
    CHECK(err < 1e-10);
    CHECK(calls < 400);
}

TEST_CASE("integrate_extrapolated runs backwards and over a zero span", "[ode]") {
    auto f = [](double, const Vec<double, 2>& y) { return Vec<double, 2>{y[1], -y[0]}; };
    const auto back = integrate_extrapolated(f, 3.0, Vec<double, 2>{std::cos(3.0), -std::sin(3.0)}, 0.0, 1e-11);
    REQUIRE(back.has_value());
    CHECK_THAT((*back)[0], WithinAbs(1.0, 1e-9));
    CHECK_THAT((*back)[1], WithinAbs(0.0, 1e-9));

    const Vec<double, 2> y0{0.25, -1.5};
    const auto same = integrate_extrapolated(f, 1.0, y0, 1.0, 1e-11);
    REQUIRE(same.has_value());
    CHECK((*same)[0] == 0.25);
}

TEST_CASE("integrate_extrapolated says so when it cannot reach the tolerance", "[ode]") {
    // y' = y^2, y(0) = 1 blows up at t = 1: integrating to 2 must report it,
    // not return a number.
    auto blow = [](double, const Vec<double, 1>& y) { return Vec<double, 1>{y[0] * y[0]}; };
    const auto r = integrate_extrapolated(blow, 0.0, Vec<double, 1>{1.0}, 2.0, 1e-10);
    CHECK_FALSE(r.has_value());
    if (!r) CHECK(r.error().code == ErrorCode::NotConverged);

    // A right-hand side that is not finite is the same.
    auto nan = [](double, const Vec<double, 1>&) { return Vec<double, 1>{std::nan("")}; };
    CHECK_FALSE(integrate_extrapolated(nan, 0.0, Vec<double, 1>{1.0}, 1.0, 1e-10).has_value());
}

TEST_CASE("integrate_extrapolated carries a Dual's derivative through the integration", "[ode]") {
    // d/dk of y(1) for y' = -k y is -exp(-k): AD through the adaptive
    // extrapolation, against the closed form.
    using D = Dual<double>;
    const D k = D::variable(1.5);
    auto f = [k](D, const Vec<D, 1>& y) { return Vec<D, 1>{D(-1.0) * k * y[0]}; };
    const auto r = integrate_extrapolated(f, D(0.0), Vec<D, 1>{D(1.0)}, D(1.0), 1e-11);
    REQUIRE(r.has_value());
    CHECK_THAT((*r)[0].value, WithinAbs(std::exp(-1.5), 1e-10));
    CHECK_THAT((*r)[0].deriv, WithinAbs(-std::exp(-1.5), 1e-9));
}

#if SPATIUM_HAS_BOOST_MULTIPRECISION
TEST_CASE("integrate_extrapolated reaches fifty-digit accuracy on Real50", "[ode][precision]") {
    // The case a fixed step cannot serve: 1e-30 where double asks 1e-11.
    const auto [err, calls] = oscillator_error<Real50>(3.0, 1e-30);
    CHECK(err < 1e-28);
}
#endif
