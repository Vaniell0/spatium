// Physically-based sound synthesis: a plucked guitar-like string and a
// struck rectangular drumhead, simulated directly from Spatium's own PDE
// time-stepping (physics/mechanics/wave_string.hpp,
// physics/mechanics/wave_membrane.hpp) and exported through the new
// zero-dependency io/wav.hpp writer -- no audio library anywhere in the
// pipeline, from differential equation to playable file. See
// docs/ROADMAP.md's "Sound synthesis" backlog entry and this feature's PR
// description for what's shipped here vs. deliberately deferred
// (geometric room acoustics via BVH ray tracing, mesh modal analysis via
// the cotangent Laplacian's eigenvalues).
//
// Outputs (in the working directory):
//   plucked_string_low.wav   -- guitar-low-E-ish pluck
//   plucked_string_high.wav  -- same code, higher string tension -> higher pitch
//   struck_drum.wav          -- rectangular membrane struck at its center
//
// Console output reports each simulation's exact closed-form predicted
// fundamental next to the physical parameters used, as a human-readable
// sanity check alongside tests/test_wave_string.cpp's and
// tests/test_wave_membrane.cpp's automated dominant-frequency checks.

#include "io_helpers.hpp"

#include <spatium/io/wav.hpp>
#include <spatium/physics/mechanics/wave_membrane.hpp>
#include <spatium/physics/mechanics/wave_string.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <print>
#include <string>
#include <string_view>
#include <vector>

using namespace spatium;
using namespace spatium::physics::mechanics;

namespace {

void print_usage() {
    std::cout <<
        "Usage: sound_synthesis_demo [--help] [--force] [--out-prefix PREFIX]\n"
        "  Simulates a plucked string (1D wave equation) and a struck\n"
        "  rectangular drumhead (2D wave equation) via explicit finite-\n"
        "  difference time-stepping, and writes the results as real WAV\n"
        "  files via the zero-dependency io/wav.hpp writer.\n"
        "  --help, -h         show this message\n"
        "  --force            overwrite existing output files\n"
        "  --out-prefix PFX   prefix prepended to every output filename\n"
        "  Outputs:           <prefix>plucked_string_low.wav\n"
        "                     <prefix>plucked_string_high.wav\n"
        "                     <prefix>struck_drum.wav\n";
}

constexpr double SAMPLE_RATE     = 44100.0;
constexpr double DT              = 1.0 / SAMPLE_RATE;
constexpr double COURANT_MARGIN  = 0.99; // safety margin under the CFL stability bound

// Normalizes a buffer to `peak` absolute amplitude in place. A simulated
// displacement field's raw units (meters) have nothing to do with a WAV
// file's [-1,1] convention -- this is a mixing-stage concern, deliberately
// kept out of wave_string.hpp/wave_membrane.hpp's own physics.
void normalize(std::vector<double>& samples, double peak = 0.8) {
    double max_abs = 0.0;
    for (double s : samples) max_abs = std::max(max_abs, std::abs(s));
    if (max_abs <= 0.0) return;
    double scale = peak / max_abs;
    for (double& s : samples) s *= scale;
}

// Simulates a plucked string and returns its displacement at an
// off-center "pickup" point, sampled at SAMPLE_RATE for `duration` seconds.
std::vector<double> render_plucked_string(double tension, double linear_density,
                                           double length, double duration, double damping) {
    double wave_speed = std::sqrt(tension / linear_density);

    // Fine enough spatial grid to stay safely under the CFL bound at this
    // dt (see wave_string.hpp's courant_number()).
    double dx_min = wave_speed * DT / COURANT_MARGIN;
    auto node_count = static_cast<std::size_t>(std::ceil(length / dx_min)) + 1;
    node_count = std::max<std::size_t>(node_count, 5);

    auto str = make_plucked_string(node_count, length, wave_speed, DT,
                                    /*pluck_position=*/0.3, /*pluck_height=*/0.01, damping);

    auto pickup = static_cast<std::size_t>(std::round(0.15 * double(node_count - 1)));
    pickup = std::clamp<std::size_t>(pickup, 1, node_count - 2);

    auto num_samples = static_cast<std::size_t>(duration * SAMPLE_RATE);
    std::vector<double> samples;
    samples.reserve(num_samples);
    samples.push_back(str.current[pickup]);
    for (std::size_t i = 1; i < num_samples; ++i) {
        str.step(DT);
        samples.push_back(str.current[pickup]);
    }
    normalize(samples);
    return samples;
}

// Simulates a struck rectangular membrane, sampled at its center (the
// (1,1) mode's antinode, and a node of every mode with an even index --
// see membrane_mode_frequency's caller below).
std::vector<double> render_struck_membrane(double surface_tension, double areal_density,
                                            double lx, double ly, double duration, double damping) {
    double wave_speed = std::sqrt(surface_tension / areal_density);

    double h_min = wave_speed * DT * std::sqrt(2.0) / COURANT_MARGIN;
    auto n = static_cast<std::size_t>(std::floor(lx / h_min)) + 1;
    n = std::max<std::size_t>(n, 9);
    if (n % 2 == 0) n += 1; // odd node count -> exact grid node at the center

    auto mem = make_struck_membrane(n, n, lx, ly, wave_speed, DT,
                                     /*strike_x=*/0.5, /*strike_y=*/0.5,
                                     /*strike_radius=*/0.05, /*strike_height=*/50.0, damping);

    std::size_t center = n / 2;
    auto num_samples = static_cast<std::size_t>(duration * SAMPLE_RATE);
    std::vector<double> samples;
    samples.reserve(num_samples);
    samples.push_back(mem.at(mem.current, center, center));
    for (std::size_t i = 1; i < num_samples; ++i) {
        mem.step(DT);
        samples.push_back(mem.at(mem.current, center, center));
    }
    normalize(samples);
    return samples;
}

bool write_note(const std::string& path, const std::vector<double>& samples, bool force) {
    if (!spatium::examples::confirm_overwrite(path, force)) return false;
    auto result = io::save_wav(samples, path);
    if (!result) {
        std::cerr << "error writing " << path << ": " << result.error().to_string() << "\n";
        return false;
    }
    return true;
}

} // namespace

int main(int argc, char** argv) {
    bool force = false;
    std::string out_prefix;
    for (int i = 1; i < argc; ++i) {
        std::string_view a = argv[i];
        if (a == "--help" || a == "-h") { print_usage(); return 0; }
        if (a == "--force") { force = true; continue; }
        if (a == "--out-prefix" && i + 1 < argc) { out_prefix = argv[++i]; continue; }
        std::cerr << "unknown option: " << a << "\n";
        return 1;
    }

    std::println("=== Physically-based sound synthesis: PDEs in, WAV out, zero audio deps ===\n");

    // ── Plucked string (1D wave equation) ─────────────────────────────
    // Guitar-low-E-ish parameters (a real low E string: ~82 Hz, ~0.648 m
    // scale length) -- physically plausible order-of-magnitude inputs,
    // not tuned to hit that exact value.
    {
        double tension = 71.0, linear_density = 0.0059, length = 0.648;
        double f1 = string_fundamental_frequency(tension, linear_density, length);
        std::println("Plucked string (low): tension={:.1f} N, linear_density={:.4f} kg/m, "
                      "length={:.3f} m -> predicted fundamental {:.2f} Hz",
                      tension, linear_density, length, f1);
        auto samples = render_plucked_string(tension, linear_density, length, 3.0, /*damping (1/s)=*/1.2);
        if (write_note(out_prefix + "plucked_string_low.wav", samples, force))
            std::println("  wrote {}plucked_string_low.wav ({} samples, {:.2f}s)",
                          out_prefix, samples.size(), samples.size() / SAMPLE_RATE);
    }

    // Same code path, higher tension -- same physical model, different note.
    {
        double tension = 220.0, linear_density = 0.0059, length = 0.648;
        double f1 = string_fundamental_frequency(tension, linear_density, length);
        std::println("\nPlucked string (high): tension={:.1f} N, linear_density={:.4f} kg/m, "
                      "length={:.3f} m -> predicted fundamental {:.2f} Hz",
                      tension, linear_density, length, f1);
        auto samples = render_plucked_string(tension, linear_density, length, 3.0, /*damping (1/s)=*/1.2);
        if (write_note(out_prefix + "plucked_string_high.wav", samples, force))
            std::println("  wrote {}plucked_string_high.wav ({} samples, {:.2f}s)",
                          out_prefix, samples.size(), samples.size() / SAMPLE_RATE);
    }

    // ── Struck drumhead (2D wave equation) ────────────────────────────
    {
        double surface_tension = 3000.0, areal_density = 0.2, lx = 0.3, ly = 0.3;
        double f11 = membrane_mode_frequency(surface_tension, areal_density, lx, ly, 1, 1);
        std::println("\nStruck drum: surface_tension={:.0f} N/m, areal_density={:.2f} kg/m^2, "
                      "{:.2f}x{:.2f} m -> predicted fundamental (1,1) mode {:.2f} Hz",
                      surface_tension, areal_density, lx, ly, f11);
        auto samples = render_struck_membrane(surface_tension, areal_density, lx, ly, 2.0, /*damping (1/s)=*/2.5);
        if (write_note(out_prefix + "struck_drum.wav", samples, force))
            std::println("  wrote {}struck_drum.wav ({} samples, {:.2f}s)",
                          out_prefix, samples.size(), samples.size() / SAMPLE_RATE);
    }

    std::println("\nDone -- open the .wav files in any player to listen.");
    return 0;
}
