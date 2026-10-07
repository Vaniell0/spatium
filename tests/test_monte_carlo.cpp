// Monte Carlo and quasi-Monte Carlo over a box (algebra/monte_carlo.hpp).
//
// Each estimate is held against a value it does not share a path with: a closed form where
// there is one (a product integrand, the volume of a 5-ball by an indicator, a torus's area
// through its metric), and nested quadrature where the dimension still allows it. What is
// asserted of a sampling estimate is statistical: within a few of its own standard errors,
// which is what its error_estimate claims to be.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/algebra/dual.hpp>
#include <spatium/algebra/monte_carlo.hpp>
#include <spatium/algebra/quadrature.hpp>
#include <spatium/spaces/metric_chart.hpp>
#include <cmath>
#include <numbers>

using namespace spatium;
using Catch::Matchers::WithinAbs;

namespace {
constexpr double kPi = std::numbers::pi;

template<std::size_t K>
Vec<double, K> filled(double v) { Vec<double, K> x{}; for (std::size_t i = 0; i < K; ++i) x[i] = v; return x; }

// torus metric g = diag((R + r cos v)^2, r^2)
struct TorusMetric {
    double R, r;
    template<class S>
    Matrix<S, 2, 2> operator()(const Vec<S, 2>& x) const {
        using std::cos;
        Matrix<S, 2, 2> g{};
        const S ring = S(R) + S(r) * cos(x[1]);
        g(0, 0) = ring * ring;
        g(1, 1) = S(r * r);
        return g;
    }
};
}  // namespace

TEST_CASE("an eight-dimensional product, where nested quadrature cannot go", "[monte_carlo]") {
    // the integral of the product of 2 x_i over [0,1]^8 is 1
    const auto f = [](const Vec<double, 8>& x) { double p = 1; for (std::size_t i = 0; i < 8; ++i) p *= 2 * x[i]; return p; };
    const auto lo = filled<8>(0.0), hi = filled<8>(1.0);
    const auto mc = monte_carlo(f, lo, hi);
    const auto qmc = quasi_monte_carlo(f, lo, hi);
    CHECK(mc.trusted());
    CHECK(qmc.trusted());
    CHECK(std::abs(mc.value - 1.0) < 4 * mc.error_estimate);
    CHECK(std::abs(qmc.value - 1.0) < 4 * qmc.error_estimate + 1e-12);
    // the low-discrepancy sequence is better at the same n, if not by much in eight dimensions on
    // 4096 points per shift (measured: a standard error 0.59 of the random one's); in two dimensions
    // below it is far better
    CHECK(qmc.error_estimate < mc.error_estimate);
    CHECK(mc.evaluations == 65536);
}

TEST_CASE("the volume of a 5-ball by an indicator", "[monte_carlo]") {
    // pi^(5/2) / Gamma(7/2) = 8 pi^2 / 15
    const double exact = 8 * kPi * kPi / 15;
    const auto inside = [](const Vec<double, 5>& x) { double r2 = 0; for (std::size_t i = 0; i < 5; ++i) r2 += x[i] * x[i]; return r2 <= 1.0 ? 1.0 : 0.0; };
    const auto lo = filled<5>(-1.0), hi = filled<5>(1.0);
    const auto mc = monte_carlo(inside, lo, hi, MonteCarloOptions{.samples = 1 << 20});
    CHECK(std::abs(mc.value - exact) < 4 * mc.error_estimate);
    CHECK(mc.error_estimate < 0.05);
}

TEST_CASE("sampling against nested quadrature, where both can run", "[monte_carlo][symmetry]") {
    // the integral of sin(x y) over [0,1]^2: two paths that share nothing
    const auto f2 = [](const Vec<double, 2>& x) { return std::sin(x[0] * x[1]); };
    const auto nested = integrate([&](double x) {
        return integrate([&](double y) { return f2(Vec<double, 2>{x, y}); }, Finite<double>{0.0, 1.0}).value;
    }, Finite<double>{0.0, 1.0});
    const auto mc = monte_carlo(f2, filled<2>(0.0), filled<2>(1.0), MonteCarloOptions{.samples = 1 << 18});
    const auto qmc = quasi_monte_carlo(f2, filled<2>(0.0), filled<2>(1.0), MonteCarloOptions{.samples = 1 << 18});
    CHECK(std::abs(mc.value - nested.value) < 4 * mc.error_estimate);
    CHECK(std::abs(qmc.value - nested.value) < 4 * qmc.error_estimate + 1e-12);
    CHECK(qmc.error_estimate < mc.error_estimate / 5);      // in two dimensions the sequence is far better
}

TEST_CASE("the area of a torus by sampling its metric's volume element", "[monte_carlo][symmetry]") {
    const double R = 2.0, r = 0.7;
    const TorusMetric g{R, r};
    const auto area = [&](const Vec<double, 2>& x) { return spaces::volume_element(g, x); };
    const Vec<double, 2> lo{0.0, 0.0}, hi{2 * kPi, 2 * kPi};
    const auto mc = monte_carlo(area, lo, hi, MonteCarloOptions{.samples = 1 << 18});
    const auto qmc = quasi_monte_carlo(area, lo, hi, MonteCarloOptions{.samples = 1 << 18});
    const double exact = 4 * kPi * kPi * R * r;
    CHECK(std::abs(mc.value - exact) < 4 * mc.error_estimate);
    CHECK(std::abs(qmc.value - exact) < 4 * qmc.error_estimate + 1e-9);
}

TEST_CASE("reproducible by seed; the status says what it can", "[monte_carlo]") {
    const auto f = [](const Vec<double, 3>& x) { return x[0] * x[1] + x[2]; };
    const auto lo = filled<3>(0.0), hi = filled<3>(1.0);
    const auto a = monte_carlo(f, lo, hi, MonteCarloOptions{.samples = 1000, .seed = 7});
    const auto b = monte_carlo(f, lo, hi, MonteCarloOptions{.samples = 1000, .seed = 7});
    const auto c = monte_carlo(f, lo, hi, MonteCarloOptions{.samples = 1000, .seed = 8});
    CHECK(a.value == b.value);
    CHECK(a.value != c.value);
    // a non-finite value is a failure, not an average
    CHECK(monte_carlo([](const Vec<double, 1>& x) { return x[0] < 0.5 ? std::nan("") : 1.0; }, filled<1>(0.0), filled<1>(1.0)).status
          == IntegralStatus::Failed);
    // too few points for a standard error to mean anything
    CHECK(monte_carlo(f, lo, hi, MonteCarloOptions{.samples = 10}).status == IntegralStatus::Suspicious);
}

TEST_CASE("a Dual passes through a sampled integral", "[monte_carlo][dual]") {
    using D = Dual<double>;
    // d/dp of the integral of p x over [0,1] is the integral of x, 1/2
    const D p = D::variable(3.0);
    const auto r = monte_carlo([p](const Vec<D, 1>& x) { return p * x[0]; }, Vec<D, 1>{D(0.0)}, Vec<D, 1>{D(1.0)},
                               MonteCarloOptions{.samples = 1 << 16});
    CHECK(std::abs(r.value.deriv - 0.5) < 4 * r.error_estimate.value / 3.0 + 1e-3);
    CHECK_THAT(r.value.value, WithinAbs(1.5, 0.05));
}
