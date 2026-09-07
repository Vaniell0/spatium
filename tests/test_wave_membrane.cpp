#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/physics/mechanics/wave_membrane.hpp>
#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

using namespace spatium;
using namespace spatium::physics::mechanics;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

// Same hand-rolled Goertzel dominant-frequency finder as
// test_wave_string.cpp (kept local to each test file rather than shared,
// since it's a ~15-line test-only helper, not library code).
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

TEST_CASE("VibratingMembrane's first step matches the exact zero-displacement Taylor formula",
          "[wave_membrane]") {
    // For displacement=0, velocity=v0 initial data, u^1 = dt*v0 exactly
    // (no diffusion term survives the very first step, since the
    // Laplacian of an identically-zero field is zero) -- see
    // make_struck_membrane's own file comment for the derivation.
    std::size_t nx = 9, ny = 9;
    double lx = 1.0, ly = 1.0;
    double c = 2.0;
    double dt = 0.001;

    auto m = make_struck_membrane(nx, ny, lx, ly, c, dt,
                                   /*strike_x=*/0.5, /*strike_y=*/0.5,
                                   /*strike_radius=*/0.3, /*strike_height=*/5.0,
                                   /*damping=*/0.0);

    std::vector<double> v0(nx * ny);
    for (std::size_t j = 0; j < ny; ++j)
        for (std::size_t i = 0; i < nx; ++i)
            v0[j * nx + i] = -m.at(m.previous, i, j) / dt; // recover v0 from previous = -dt*v0

    for (auto& u : m.current) CHECK(u == 0.0);

    m.step(dt);

    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            double expected = dt * v0[j * nx + i];
            CHECK_THAT(m.at(m.current, i, j), WithinAbs(expected, 1e-12));
        }
    }
}

TEST_CASE("Struck rectangular membrane's simulated fundamental matches the closed form",
          "[wave_membrane]") {
    // Illustrative small-drum-ish parameters: surface tension and areal
    // density order-of-magnitude typical of a taut membrane, not tuned to
    // a specific real instrument.
    double surface_tension = 3000.0; // N/m
    double areal_density   = 0.2;    // kg/m^2
    double lx = 0.3, ly = 0.3;       // m
    double sample_rate = 44100.0;
    double dt = 1.0 / sample_rate;

    double wave_speed = std::sqrt(surface_tension / areal_density);
    double expected_f11 = membrane_mode_frequency(surface_tension, areal_density, lx, ly, 1, 1);

    double courant_margin = 0.99;
    // 2D CFL: c*dt*sqrt(1/dx^2+1/dy^2) <= margin: with dx=dy=h,
    // c*dt*sqrt(2)/h <= margin  =>  h >= c*dt*sqrt(2)/margin.
    double h_min = wave_speed * dt * std::sqrt(2.0) / courant_margin;
    auto n = static_cast<std::size_t>(std::floor(lx / h_min)) + 1;
    if (n < 9) n = 9;
    if (n % 2 == 0) n += 1; // odd node count -> an exact grid node at the
                            // geometric center, matching strike/pickup position

    double h = lx / double(n - 1);
    REQUIRE(courant_number_2d(wave_speed, dt, h, h) <= 1.0);

    // Strike (and later sample) at the center: the (1,1) mode's antinode,
    // and a node of every (even, n) / (m, even) mode, so the fundamental
    // dominates what's excited and measured.
    auto mem = make_struck_membrane(n, n, lx, ly, wave_speed, dt,
                                     /*strike_x=*/0.5, /*strike_y=*/0.5,
                                     /*strike_radius=*/0.05, /*strike_height=*/50.0,
                                     /*damping=*/0.0);

    std::size_t center = n / 2;
    auto num_samples = static_cast<std::size_t>(0.2 * sample_rate);
    std::vector<double> samples;
    samples.reserve(num_samples);
    samples.push_back(mem.at(mem.current, center, center));
    for (std::size_t i = 1; i < num_samples; ++i) {
        mem.step(dt);
        samples.push_back(mem.at(mem.current, center, center));
    }

    double measured_f11 = find_dominant_frequency(samples, sample_rate,
                                                   expected_f11 * 0.5, expected_f11 * 2.0, 0.5);

    // Slightly looser than the 1D string check (2D numerical dispersion
    // on a coarser relative grid), still a real closed-form comparison.
    CHECK_THAT(measured_f11, WithinRel(expected_f11, 0.02));
}
