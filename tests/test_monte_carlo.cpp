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
#include <array>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <type_traits>
#include <vector>

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

// ── the Sobol sequence ─────────────────────────────────────────

namespace {
// The raw (unshifted) Sobol points of Gray-code order, n = 0 .. count-1, as fractions.
std::vector<std::vector<double>> sobol_points(std::size_t dims, std::size_t count) {
    std::vector<std::array<std::uint32_t, 33>> v;
    for (std::size_t d = 0; d < dims; ++d) v.push_back(monte_carlo_detail::sobol_directions(d));
    std::vector<std::vector<double>> out;
    std::vector<std::uint32_t> state(dims, 0);
    for (std::size_t n = 0; n < count; ++n) {
        if (n > 0) {
            std::uint64_t value = n - 1;
            int c = 1;
            while (value & 1u) { value >>= 1; ++c; }
            for (std::size_t d = 0; d < dims; ++d) state[d] ^= v[d][c];
        }
        std::vector<double> p(dims);
        for (std::size_t d = 0; d < dims; ++d) p[d] = static_cast<double>(state[d]) / 4294967296.0;
        out.push_back(p);
    }
    return out;
}
}  // namespace

TEST_CASE("the first Sobol points are the known ones", "[monte_carlo][sobol]") {
    const auto p = sobol_points(2, 5);
    const double d1[] = {0.0, 0.5, 0.75, 0.25, 0.375};          // van der Corput in Gray-code order
    const double d2[] = {0.0, 0.5, 0.25, 0.75, 0.375};
    for (std::size_t n = 0; n < 5; ++n) {
        CHECK(p[n][0] == d1[n]);
        CHECK(p[n][1] == d2[n]);
    }
}

TEST_CASE("Sobol's first two coordinates are a (0, m, 2)-net: one point in every elementary box", "[monte_carlo][sobol][symmetry]") {
    // 64 points, every box of width 2^-a by 2^-b with a + b = 6 holds exactly one: a property of
    // the construction that a wrong recurrence or direction number would break
    const auto p = sobol_points(2, 64);
    for (int a = 0; a <= 6; ++a) {
        const int b = 6 - a;
        std::vector<int> count(64, 0);
        for (const auto& q : p) {
            const std::size_t i = static_cast<std::size_t>(q[0] * (1 << a)), j = static_cast<std::size_t>(q[1] * (1 << b));
            ++count[i * (1u << b) + j];
        }
        for (int c : count) CHECK(c == 1);
    }
}

TEST_CASE("every one of the 40 coordinates is stratified: 2^10 points, one per interval", "[monte_carlo][sobol][symmetry]") {
    // the one-dimensional projection of the first 2^m Sobol points is a permutation of the grid
    // i / 2^m -- for each coordinate, which holds the whole direction-number table to account
    const auto p = sobol_points(40, 1024);
    for (std::size_t d = 0; d < 40; ++d) {
        std::vector<int> seen(1024, 0);
        for (const auto& q : p) ++seen[static_cast<std::size_t>(q[d] * 1024)];
        INFO("coordinate " << d);
        for (int c : seen) CHECK(c == 1);
    }
}

TEST_CASE("Sobol against Halton and random points in high dimensions", "[monte_carlo][sobol]") {
    // a smooth integrand of low effective dimension: the product over i of the integral of exp(x_i / K)
    const auto make = [](auto K_tag) {
        constexpr std::size_t K = decltype(K_tag)::value;
        return [](const Vec<double, K>& x) { double s = 0; for (std::size_t i = 0; i < K; ++i) s += x[i]; return std::exp(s / static_cast<double>(K)); };
    };
    const auto check = [&](auto K_tag, bool against_halton) {
        constexpr std::size_t K = decltype(K_tag)::value;
        const auto f = make(K_tag);
        const double exact = std::pow(static_cast<double>(K) * (std::exp(1.0 / static_cast<double>(K)) - 1.0), static_cast<double>(K));
        const auto lo = filled<K>(0.0), hi = filled<K>(1.0);
        MonteCarloOptions o;
        const auto mc = monte_carlo(f, lo, hi, o);
        const auto so = quasi_monte_carlo(f, lo, hi, o);
        INFO("K = " << K);
        CHECK(std::abs(so.value - exact) < 4 * so.error_estimate + 1e-12);
        CHECK(so.error_estimate < mc.error_estimate / 20);          // measured 65x at K = 40
        if (against_halton) {
            o.sequence = Sequence::Halton;
            const auto ha = quasi_monte_carlo(f, lo, hi, o);
            CHECK(so.error_estimate < ha.error_estimate / 4);       // measured 8x at K = 16, 18x at K = 32
        }
    };
    check(std::integral_constant<std::size_t, 8>{}, false);
    check(std::integral_constant<std::size_t, 16>{}, true);
    check(std::integral_constant<std::size_t, 32>{}, true);
    check(std::integral_constant<std::size_t, 40>{}, false);
}

TEST_CASE("more dimensions than a sequence has is a failure, not a quiet wrong answer", "[monte_carlo][sobol]") {
    const auto f = [](const Vec<double, 41>&) { return 1.0; };
    CHECK(quasi_monte_carlo(f, filled<41>(0.0), filled<41>(1.0)).status == IntegralStatus::Failed);
    const auto g = [](const Vec<double, 33>&) { return 1.0; };
    CHECK(quasi_monte_carlo(g, filled<33>(0.0), filled<33>(1.0), MonteCarloOptions{.sequence = Sequence::Halton}).status
          == IntegralStatus::Failed);
    CHECK(quasi_monte_carlo(g, filled<33>(0.0), filled<33>(1.0)).trusted());   // Sobol has 40
}
