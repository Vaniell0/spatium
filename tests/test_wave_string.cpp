#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/physics/mechanics/wave_string.hpp>
#include <cmath>
#include <numbers>
#include <vector>

using namespace spatium;
using namespace spatium::physics::mechanics;
using Catch::Matchers::WithinRel;

namespace {

// Small hand-rolled Goertzel-algorithm power estimator -- exactly the
// "small DFT, doesn't need a full spectral library" the task calls for:
// O(N) per candidate frequency, no FFT machinery, evaluated at a scan of
// candidate frequencies to find the dominant one.
double goertzel_power(const std::vector<double>& samples, double sample_rate, double freq) {
    double omega = 2.0 * std::numbers::pi * freq / sample_rate;
    double coeff = 2.0 * std::cos(omega);
    double s0 = 0.0, s1 = 0.0, s2 = 0.0;
    for (double x : samples) {
        s0 = x + coeff * s1 - s2;
        s2 = s1;
        s1 = s0;
    }
    return s1 * s1 + s2 * s2 - coeff * s1 * s2;
}

// Scans [f_min, f_max] in steps of f_step and returns the frequency with
// the largest Goertzel power -- the dominant spectral component.
double find_dominant_frequency(const std::vector<double>& samples, double sample_rate,
                                double f_min, double f_max, double f_step) {
    double best_f = f_min;
    double best_power = -1.0;
    for (double f = f_min; f <= f_max; f += f_step) {
        double p = goertzel_power(samples, sample_rate, f);
        if (p > best_power) {
            best_power = p;
            best_f = f;
        }
    }
    return best_f;
}

} // namespace

TEST_CASE("VibratingString's first step matches the exact zero-velocity Taylor formula", "[wave_string]") {
    // u^1 = u^0 + (r2/2)*laplacian(u^0) is the well-known closed-form first
    // step for zero initial velocity -- checked directly against
    // make_plucked_string's phantom-previous-state construction plus one
    // step(), independent of any later frequency-domain check.
    std::size_t node_count = 11;
    double length = 1.0;
    double c = 2.0;
    double dt = 0.01;
    auto s = make_plucked_string(node_count, length, c, dt, 0.3, 1.0, 0.0);

    std::vector<double> u0 = s.current;
    double h = s.dx();
    double r2 = (c * dt / h) * (c * dt / h);

    s.step(dt);

    for (std::size_t i = 1; i + 1 < node_count; ++i) {
        double laplacian = u0[i + 1] - 2.0 * u0[i] + u0[i - 1];
        double expected = u0[i] + 0.5 * r2 * laplacian;
        CHECK_THAT(s.current[i], WithinRel(expected, 1e-10));
    }
    CHECK(s.current[0] == 0.0);
    CHECK(s.current[node_count - 1] == 0.0);
}

TEST_CASE("make_plucked_string builds the triangular pluck shape and fixes both endpoints", "[wave_string]") {
    std::size_t node_count = 21;
    double length = 2.0;
    double c = 1.0;
    double dt = 0.001;
    double pluck_position = 0.25;
    double pluck_height = 0.5;

    auto s = make_plucked_string(node_count, length, c, dt, pluck_position, pluck_height, 0.0);

    CHECK(s.current[0] == 0.0);
    CHECK(s.current[node_count - 1] == 0.0);

    double h = s.dx();
    std::size_t peak_index = static_cast<std::size_t>(std::round(pluck_position * length / h));
    double peak = s.current[peak_index];
    for (double v : s.current)
        CHECK(v <= peak + 1e-12);
    CHECK_THAT(peak, WithinRel(pluck_height, 1e-2));
}

TEST_CASE("Plucked string's simulated fundamental frequency matches the closed form", "[wave_string]") {
    // Parameters in the ballpark of a guitar's low-E string (real low-E:
    // ~82 Hz, ~0.648 m scale length) -- not tuned to hit that value
    // exactly, just physically plausible order-of-magnitude inputs.
    double tension         = 71.0;    // N
    double linear_density  = 0.0059;  // kg/m
    double length          = 0.648;   // m
    double sample_rate     = 44100.0;
    double dt              = 1.0 / sample_rate;

    double wave_speed = std::sqrt(tension / linear_density);
    double expected_f1 = string_fundamental_frequency(tension, linear_density, length);

    // Choose a spatial grid fine enough to stay safely under the CFL bound
    // at this dt (see wave_string.hpp's courant_number()).
    double courant_margin = 0.99;
    double dx_min = wave_speed * dt / courant_margin;
    auto node_count = static_cast<std::size_t>(std::ceil(length / dx_min)) + 1;
    if (node_count < 5) node_count = 5;

    double dx = length / double(node_count - 1);
    REQUIRE(courant_number(wave_speed, dt, dx) <= 1.0);

    auto str = make_plucked_string(node_count, length, wave_speed, dt,
                                    /*pluck_position=*/0.3, /*pluck_height=*/0.01,
                                    /*damping=*/0.0);

    // Sample displacement at an off-center "pickup" point over a 0.3s
    // window -- long enough for several periods of the ~85 Hz fundamental,
    // short enough that the Goertzel scan below stays fast.
    auto pickup = static_cast<std::size_t>(std::round(0.15 * double(node_count - 1)));
    pickup = std::max<std::size_t>(1, std::min(pickup, node_count - 2));

    auto num_samples = static_cast<std::size_t>(0.3 * sample_rate);
    std::vector<double> samples;
    samples.reserve(num_samples);
    samples.push_back(str.current[pickup]);
    for (std::size_t i = 1; i < num_samples; ++i) {
        str.step(dt);
        samples.push_back(str.current[pickup]);
    }

    double measured_f1 = find_dominant_frequency(samples, sample_rate,
                                                  expected_f1 * 0.5, expected_f1 * 2.0, 0.5);

    // 1% tolerance covers the explicit central-difference scheme's own
    // numerical dispersion (phase velocity is slightly below the ideal
    // continuum c at finite grid resolution) plus the 0.5 Hz scan step,
    // while still being a real, tight, closed-form-vs-simulation check.
    CHECK_THAT(measured_f1, WithinRel(expected_f1, 0.01));
}

TEST_CASE("Doubling tension raises the simulated fundamental by sqrt(2), as the closed form predicts",
          "[wave_string]") {
    double linear_density = 0.005;
    double length = 0.5;
    double sample_rate = 44100.0;
    double dt = 1.0 / sample_rate;

    auto simulate_f1 = [&](double tension) {
        double wave_speed = std::sqrt(tension / linear_density);
        double dx_min = wave_speed * dt / 0.99;
        auto node_count = static_cast<std::size_t>(std::ceil(length / dx_min)) + 1;
        if (node_count < 5) node_count = 5;

        auto str = make_plucked_string(node_count, length, wave_speed, dt, 0.3, 0.01, 0.0);
        auto pickup = static_cast<std::size_t>(std::round(0.2 * double(node_count - 1)));
        pickup = std::max<std::size_t>(1, std::min(pickup, node_count - 2));

        auto num_samples = static_cast<std::size_t>(0.2 * sample_rate);
        std::vector<double> samples;
        samples.reserve(num_samples);
        samples.push_back(str.current[pickup]);
        for (std::size_t i = 1; i < num_samples; ++i) {
            str.step(dt);
            samples.push_back(str.current[pickup]);
        }

        double f1 = string_fundamental_frequency(tension, linear_density, length);
        return find_dominant_frequency(samples, sample_rate, f1 * 0.5, f1 * 2.0, 0.5);
    };

    double f_low  = simulate_f1(80.0);
    double f_high = simulate_f1(160.0);

    CHECK_THAT(f_high / f_low, WithinRel(std::sqrt(2.0), 0.02));
}
