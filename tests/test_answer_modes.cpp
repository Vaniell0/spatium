// Damped sinusoids in noisy data, counted by description length. Held on a ringdown-shaped signal -- 4096 samples
// a second, 30 ms, modes near 250 and 330 Hz with decay times of a few milliseconds -- with white noise, and
// with the correlated noise of a covariance: the count is the true one, the recovered frequencies and decay
// times lie within a few of their standard errors of the truth, and data of noise alone gives no mode. The cases
// it must fail: a covariance that is not positive definite, a grid with no frequency in it.

#include <catch2/catch_test_macros.hpp>
#include <spatium/answer/modes.hpp>

#include <cmath>
#include <cstdint>
#include <numbers>
#include <random>
#include <vector>

namespace {

namespace sa = spatium::answer;

struct Truth {
    double f, tau, a, b;
};

std::vector<double> signal(const std::vector<double>& t, const std::vector<Truth>& modes) {
    std::vector<double> y(t.size(), 0.0);
    for (const auto& m : modes)
        for (std::size_t i = 0; i < t.size(); ++i)
            y[i] += std::exp(-t[i] / m.tau) * (m.a * std::cos(2 * std::numbers::pi * m.f * t[i]) + m.b * std::sin(2 * std::numbers::pi * m.f * t[i]));
    return y;
}

sa::ModeSpec grid(std::size_t modes) {
    sa::ModeSpec s;
    s.max_modes = modes;
    s.frequency_lo = 120, s.frequency_hi = 450, s.frequency_step = 6;
    s.decay_lo = 0.5e-3, s.decay_hi = 20e-3, s.decay_count = 20;
    return s;
}

const std::size_t n = 123;   // 30 ms at 4096 a second

std::vector<double> times() {
    std::vector<double> t(n);
    for (std::size_t i = 0; i < n; ++i) t[i] = double(i) / 4096.0;
    return t;
}

} // namespace

TEST_CASE("one mode in white noise is counted as one, with errors that cover the truth", "[answer][modes]") {
    const auto t = times();
    const std::vector<Truth> truth{{250.0, 4e-3, 1.0, 0.4}};
    int right_count = 0, covered_f = 0, covered_tau = 0;
    const int trials = 20;
    for (int seed = 0; seed < trials; ++seed) {
        std::mt19937_64 gen(1234 + seed);
        std::normal_distribution<double> unit(0.0, 1.0);
        auto y = signal(t, truth);
        for (auto& v : y) v += 0.05 * unit(gen);
        const sa::Modes m = sa::synthesize_modes(t, y, std::vector<double>(n, 0.05), grid(2));
        REQUIRE(m.lengths.size() == 3);
        if (m.modes.size() != 1) continue;
        ++right_count;
        covered_f += std::abs(m.modes[0].frequency - 250.0) < 3 * m.modes[0].frequency_error;
        covered_tau += std::abs(m.modes[0].decay - 4e-3) < 3 * m.modes[0].decay_error;
    }
    INFO("count right " << right_count << ", f covered " << covered_f << ", tau covered " << covered_tau << " of " << trials);
    CHECK(right_count >= 18);
    CHECK(covered_f >= right_count - 2);
    CHECK(covered_tau >= right_count - 2);
}

TEST_CASE("two modes are counted as two, noise alone as none", "[answer][modes]") {
    const auto t = times();
    const std::vector<Truth> truth{{250.0, 4e-3, 1.0, 0.4}, {340.0, 2e-3, -0.8, 0.5}};
    std::mt19937_64 gen(7);
    std::normal_distribution<double> unit(0.0, 1.0);
    auto y = signal(t, truth);
    for (auto& v : y) v += 0.03 * unit(gen);
    const sa::Modes two = sa::synthesize_modes(t, y, std::vector<double>(n, 0.03), grid(3));
    REQUIRE(two.modes.size() == 2);
    // the length falls from none to one to two modes, and the two are the ones put in, in some order
    CHECK(two.lengths[1] < two.lengths[0]);
    CHECK(two.lengths[2] < two.lengths[1]);
    const auto& m0 = two.modes[0].frequency < two.modes[1].frequency ? two.modes[0] : two.modes[1];
    const auto& m1 = two.modes[0].frequency < two.modes[1].frequency ? two.modes[1] : two.modes[0];
    CHECK(std::abs(m0.frequency - 250.0) < 4 * m0.frequency_error);
    CHECK(std::abs(m1.frequency - 340.0) < 4 * m1.frequency_error);
    CHECK(two.rms < 1.3);

    std::vector<double> noise(n);
    for (auto& v : noise) v = 0.03 * unit(gen);
    const sa::Modes none = sa::synthesize_modes(t, noise, std::vector<double>(n, 0.03), grid(2));
    CHECK(none.modes.empty());
    CHECK(none.lengths.size() == 3);
}

TEST_CASE("correlated noise given as a covariance is whitened, and a bad covariance or a bad grid gives nothing", "[answer][modes]") {
    const auto t = times();
    // AR(1) noise of correlation 0.8 and standard deviation 0.05: the covariance is sigma^2 rho^|i - j|
    const double rho = 0.8, sd = 0.05;
    std::vector<double> cov(n * n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) cov[i * n + j] = sd * sd * std::pow(rho, std::abs(double(i) - double(j)));
    const std::vector<Truth> truth{{250.0, 4e-3, 1.0, 0.4}};
    int one = 0, none_in_noise = 0;
    const int trials = 20;
    for (int seed = 0; seed < trials; ++seed) {
        std::mt19937_64 gen(99 + seed);
        std::normal_distribution<double> unit(0.0, 1.0);
        std::vector<double> e(n);
        e[0] = sd * unit(gen);
        for (std::size_t i = 1; i < n; ++i) e[i] = rho * e[i - 1] + sd * std::sqrt(1 - rho * rho) * unit(gen);
        auto y = signal(t, truth);
        for (std::size_t i = 0; i < n; ++i) y[i] += e[i];
        sa::ModeSpec spec = grid(2);
        spec.covariance = cov;
        one += sa::synthesize_modes(t, y, {}, spec).modes.size() == 1;
        none_in_noise += sa::synthesize_modes(t, e, {}, spec).modes.empty();
    }
    INFO("one mode " << one << ", none in noise alone " << none_in_noise << " of " << trials);
    CHECK(one >= 18);
    CHECK(none_in_noise >= 18);

    sa::ModeSpec bad = grid(2);
    auto notpd = cov;
    notpd[1] = notpd[n] = 2 * sd * sd;   // |rho| > 1
    bad.covariance = notpd;
    CHECK(sa::synthesize_modes(t, signal(t, truth), {}, bad).lengths.empty());
    sa::ModeSpec empty = grid(2);
    empty.frequency_hi = empty.frequency_lo - 1;
    CHECK(sa::synthesize_modes(t, signal(t, truth), std::vector<double>(n, 0.05), empty).lengths.empty());
}
