// The chance that a random walk on Z^3 ever comes back, for walks whose steps are one orbit
// of the cube's symmetries.
//
// A walk that steps along the orbit of v under the signed permutations (algebra::orbit over
// SignedPermutationAction) comes back with probability p = 1 - 1/W, where W, the mean number
// of visits to the start, is Watson's integral
//
//   W = (1/pi^3) int over [0, pi]^3 of dk / (1 - phi(k)),   phi(k) = mean over the steps s of cos(k.s)
//
// singular like 1/|k|^2 at the corners of the cube where phi = 1. Three orbits are three
// lattices with a closed form in Gamma functions: (1,0,0) the simple cubic, (1,1,1) the
// body-centred, (1,1,0) the face-centred. The integral is taken by three nested
// integrate() calls, whose doubly exponential rule meets the corner singularity at an end
// of each interval.
//
// The denominator is written as the mean of 2 sin^2(k.s / 2): 1 - cos is that exactly, and
// computed as 1 - cos it cancels to 0 for k below about 1e-8, where the integrand then is
// infinite. The quadrature is held to refuse that form rather than return a number.
//
// For the simple cubic walk a second route is the series W = sum_n u_2n of return
// probabilities after 2n steps, u_2n = C(2n,n) S(n) / 36^n with S(n) the sum of squared
// trinomial coefficients, itself checked in two forms. The terms fall like n^(-3/2), so a
// partial sum to N is off by about N^(-1/2). The leading asymptotic 2 (3/(4 pi n))^(3/2)
// summed exactly through zeta(3/2) leaves N^(-3/2), and Richardson in N leaves N^(-5/2):
// the test holds both orders, measured by doubling N. The library has no general summation
// of slowly falling series; this is that border, measured on one series: at N = 1000 the
// series is good to about 3e-9, where the integral reaches 1e-15.

#include <spatium/algebra/groups/action.hpp>
#include <spatium/algebra/groups/hyperoctahedral.hpp>
#include <spatium/algebra/quadrature.hpp>

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <numbers>
#include <vector>

using namespace spatium;

namespace {

constexpr double pi = std::numbers::pi;

struct Lattice {
    const char* name;
    IntVec<3> step;
    std::size_t steps;
    double closed; // W in Gamma functions
};

std::vector<Lattice> lattices() {
    const double g24 = std::tgamma(1.0 / 24) * std::tgamma(5.0 / 24) * std::tgamma(7.0 / 24) *
                       std::tgamma(11.0 / 24);
    return {
        {"simple cubic", {1, 0, 0}, 6, std::sqrt(6.0) / (32 * pi * pi * pi) * g24},
        {"body-centred", {1, 1, 1}, 8, std::pow(std::tgamma(0.25), 4) / (4 * pi * pi * pi)},
        // Watson's 3 Gamma(1/3)^6 / (2^(14/3) pi^4) is the integral of 1/(3 - sum cos cos); with
        // 1 - phi = (3 - sum cos cos) / 3 the visits are three times it
        {"face-centred", {1, 1, 0}, 12,
         9 * std::pow(std::tgamma(1.0 / 3), 6) / (std::pow(2.0, 14.0 / 3) * std::pow(pi, 4))},
    };
}

struct Watson {
    double value;
    double estimate;
    bool trusted;
};

template <class Inverse> Watson watson(Inverse&& inverse, double tolerance) {
    const IntegralOptions o{.tolerance = tolerance};
    const auto r = integrate(
        [&](double a) {
            return integrate(
                       [&](double b) {
                           return integrate([&](double c) { return inverse(a, b, c); },
                                            Finite<double>{0, pi}, o)
                               .value;
                       },
                       Finite<double>{0, pi}, o)
                .value;
        },
        Finite<double>{0, pi}, o);
    const double cube = pi * pi * pi;
    return {r.value / cube, r.error_estimate / cube, r.trusted()};
}

// 1 / (1 - phi(k)) without cancellation: 1 - cos x = 2 sin^2(x / 2).
auto inverse_of(const std::vector<IntVec<3>>& steps) {
    return [&steps](double a, double b, double c) {
        double gap = 0;
        for (const auto& s : steps) {
            const double h = std::sin((a * s[0] + b * s[1] + c * s[2]) / 2);
            gap += 2 * h * h;
        }
        return static_cast<double>(steps.size()) / gap;
    };
}

double log_choose(int n, int k) {
    return std::lgamma(n + 1.0) - std::lgamma(k + 1.0) - std::lgamma(n - k + 1.0);
}

double log_sum_exp(const std::vector<double>& t) {
    const double m = *std::ranges::max_element(t);
    double s = 0;
    for (double x : t)
        s += std::exp(x - m);
    return m + std::log(s);
}

// log of S(n) = sum over j + k + l = n of (n! / (j! k! l!))^2, as the sum itself...
double log_s_by_triples(int n) {
    std::vector<double> t;
    for (int j = 0; j <= n; ++j)
        for (int k = 0; j + k <= n; ++k)
            t.push_back(2 * (std::lgamma(n + 1.0) - std::lgamma(j + 1.0) - std::lgamma(k + 1.0) -
                             std::lgamma(n - j - k + 1.0)));
    return log_sum_exp(t);
}

// ... and by the identity S(n) = sum_k C(n,k)^2 C(2k,k), a single sum.
double log_s(int n) {
    std::vector<double> t;
    for (int k = 0; k <= n; ++k)
        t.push_back(2 * log_choose(n, k) + log_choose(2 * k, k));
    return log_sum_exp(t);
}

struct SeriesSums {
    double partial, with_tail;
};

// sum_{n <= N} u_2n, and the same plus the leading asymptotic's tail past N, summed exactly
// as K (zeta(3/2) - sum_{n <= N} n^(-3/2)).
SeriesSums series_to(int N) {
    constexpr double zeta_3_2 = 2.612375348685488343;
    const double K = 2 * std::pow(3 / (4 * pi), 1.5);
    double sum = 0, harmonic = 0;
    for (int n = 0; n <= N; ++n) {
        sum += std::exp(log_choose(2 * n, n) + log_s(n) - n * std::log(36.0));
        if (n > 0)
            harmonic += std::pow(n, -1.5);
    }
    return {sum, sum + K * (zeta_3_2 - harmonic)};
}

} // namespace

TEST_CASE("Watson's integral for three orbit walks matches its Gamma form", "[algebra][quadrature][return-probability]") {
    for (const Lattice& l : lattices()) {
        INFO(l.name);
        const auto steps = orbit(HyperoctahedralGroup<3>{}, SignedPermutationAction<3>{}, l.step);
        REQUIRE(steps.size() == l.steps);
        const Watson w = watson(inverse_of(steps), 1e-3);
        CHECK(w.trusted);
        // measured: relative errors 4e-7..9e-7 under estimates of 2e-5..5e-5
        CHECK(std::abs(w.value - l.closed) <= w.estimate);
        CHECK(std::abs(w.value - l.closed) <= 2e-6 * l.closed);
    }
    // the return probabilities: 0.3405 (Polya's), 0.2822, 0.2563
    const auto all = lattices();
    CHECK(std::abs((1 - 1 / all[0].closed) - 0.340537329551) < 1e-12);
    CHECK(std::abs((1 - 1 / all[1].closed) - 0.282229988954) < 1e-12);
    CHECK(std::abs((1 - 1 / all[2].closed) - 0.256318236505) < 1e-12);
}

TEST_CASE("Watson's integral to the last digits", "[.deep][algebra][quadrature][return-probability]") {
    // about 55 s in all at -O2; measured relative errors 3e-16, 6e-15, 5e-15
    for (const Lattice& l : lattices()) {
        INFO(l.name);
        const auto steps = orbit(HyperoctahedralGroup<3>{}, SignedPermutationAction<3>{}, l.step);
        const Watson w = watson(inverse_of(steps), 0);
        CHECK(w.trusted);
        CHECK(std::abs(w.value - l.closed) <= 1e-13 * l.closed);
    }
}

TEST_CASE("the integrand written as 1 - cos is refused, not summed", "[algebra][quadrature][return-probability]") {
    const auto naive = [](double a, double b, double c) {
        return 1 / (1 - (std::cos(a) + std::cos(b) + std::cos(c)) / 3);
    };
    const Watson w = watson(naive, 1e-3);
    CHECK_FALSE(w.trusted);
    CHECK(std::isnan(w.value));
}

TEST_CASE("the series of return probabilities, with its tail and Richardson", "[algebra][return-probability]") {
    for (int n = 0; n <= 40; ++n) {
        INFO("n = " << n);
        CHECK(std::abs(log_s(n) - log_s_by_triples(n)) < 1e-12);
    }

    const double W = lattices()[0].closed;
    const SeriesSums at250 = series_to(250), at500 = series_to(500), at1000 = series_to(1000);
    const double r = std::pow(2.0, 1.5);
    const double rich500 = (r * at500.with_tail - at250.with_tail) / (r - 1);
    const double rich1000 = (r * at1000.with_tail - at500.with_tail) / (r - 1);

    // the bare partial sum is far off: about N^(-1/2)
    CHECK(std::abs(at1000.partial - W) > 1e-2);
    // the tail leaves N^(-3/2): doubling N divides the error by about 2^1.5 = 2.83
    const double tail_ratio = std::abs(at500.with_tail - W) / std::abs(at1000.with_tail - W);
    CHECK(tail_ratio > 2.5);
    CHECK(tail_ratio < 3.2);
    // Richardson leaves N^(-5/2): by about 2^2.5 = 5.66
    const double rich_ratio = std::abs(rich500 - W) / std::abs(rich1000 - W);
    CHECK(rich_ratio > 5.0);
    CHECK(rich_ratio < 6.4);
    // measured 2.6e-9 at N = 1000
    CHECK(std::abs(rich1000 - W) < 5e-9);
}
