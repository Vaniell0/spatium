#pragma once

// Procedural sky for CPU raytracers: point stars + a handful of spiral-
// galaxy/gas-cloud structures + a background gradient, evaluated live by
// sample_sky()/sample_sky_color(), nothing pre-rasterized.
//
// Promoted out of examples/io_helpers.hpp (2026-08-28) once a 4th caller
// (geodesic_procgen_demo.cpp) needed the color/star-scatter pieces but
// not the whole-sky gradient -- see `wide_sky` below -- the same "found
// duplicated/needed by another caller -> promote to the engine" pattern
// camera.hpp/parallel_for_rows.hpp/supersample.hpp already went through.
//
// No texture buffer anywhere in this file, deliberately -- a first
// version rendered stars/nebulae into a 2048x1024 RGB raster and read it
// back via sample_sky()'s nearest-neighbor lookup. That's a real,
// reported bug, not a style choice: baking a star into a handful of
// texture pixels and then nearest-neighbor-sampling it blows up into a
// large blocky diamond wherever the screen-to-texture mapping is coarse
// (which gravitational lensing guarantees happens somewhere in every
// render) -- no texture resolution fixes that, since nearest-neighbor
// sampling has no notion of a feature smaller than one texel to shrink
// gracefully. `Sky` instead stores the stars/nebulae themselves (unit
// directions + closed-form falloffs), and `sample_sky()` evaluates them
// directly against the EXACT queried direction -- resolution-independent
// by construction, sharp at any zoom or lensing stretch.
//
// Stars are bucketed on a coarse (theta, phi) grid, queried at sample
// time, so a render (which can easily call sample_sky() tens of millions
// of times) only checks nearby stars per call (O(stars-per-bucket), not
// O(all stars)). Nebulae/galaxies are few enough (single digits) to just
// check all of them per call.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/render/color.hpp>
#  include <algorithm>
#  include <cmath>
#  include <cstdint>
#  include <numbers>
#  include <random>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::render {

struct Star {
    Vec<double, 3> dir;
    double ang_size, brightness;
};
struct Spiral {
    Vec<double, 3> center, u, v;
    double ang_radius, tightness;
    int arms;
    Vec<double, 3> tint;
};
struct Cloud {
    Vec<double, 3> center;
    double ang_radius, freq, phase;
    Vec<double, 3> tint;
};

struct Sky {
    Vec<double, 3> tint;
    // true: whole-sky brightness gradient, brighter near mid-latitude and
    // dark at the poles -- correct when most/all of the sky sphere is
    // visible in frame (the wide-FOV GR lensing views this was designed
    // for). false: flat background at `tint` -- correct for a close,
    // narrow-FOV single-object shot, where only a small patch of that
    // gradient would ever be visible and reads as an unexplained dark
    // smear rather than a recognizable sky (found by
    // geodesic_procgen_demo.cpp's planet framing, 2026-08-28).
    bool wide_sky = true;
    std::vector<Star> stars;
    std::vector<Spiral> spirals;
    std::vector<Cloud> clouds;

    static constexpr int kThetaBuckets = 48, kPhiBuckets = 96;
    std::vector<std::vector<std::uint32_t>> star_buckets;  // size kThetaBuckets*kPhiBuckets, indices into stars
};

// Uniform-on-sphere direction, in sample_sky()'s own convention (y =
// "up" = cos(theta), x/z span the phi=0 plane).
inline Vec<double, 3> random_sky_dir(std::mt19937& rng) {
    constexpr double pi = std::numbers::pi;
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    double cz = 2.0 * u01(rng) - 1.0;
    double phi = 2.0 * pi * u01(rng);
    double s = std::sqrt(std::max(0.0, 1.0 - cz * cz));
    return Vec<double, 3>{s * std::cos(phi), cz, s * std::sin(phi)};
}

// Deep gradient background (tinted by `tint`) when `wide_sky` (default),
// with a few procedural spiral-galaxy/gas-cloud structures and analytic
// point stars -- all evaluated live by sample_sky(). Different
// seeds/tints give visually distinct skies.
//
// Spiral arms: a logarithmic-spiral density field, density = 0.5 +
// 0.5*cos(arms*local_angle - tightness*ln(ang_dist)) -- the standard
// closed-form parametrization of a logarithmic spiral, windowed by a
// broad gaussian-ish falloff in angular distance from the center (not a
// bare exponential, which fades from the very center instead of staying
// bright across most of the structure's footprint).
//
// Gas clouds: a handful of summed sinusoids over the direction vector's
// own components -- a cheap value-noise substitute (no permutation
// tables, no external dependency), windowed the same way.
inline Sky make_starfield(int n_stars, std::uint32_t seed = 42,
                           Vec<double, 3> tint = {8.0, 6.0, 20.0}, bool wide_sky = true) {
    constexpr double pi = std::numbers::pi;
    Sky sky;
    sky.tint = tint;
    sky.wide_sky = wide_sky;
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u01(0.0, 1.0);

    auto tangent_frame = [](const Vec<double, 3>& center, Vec<double, 3>& u, Vec<double, 3>& v) {
        Vec<double, 3> arbitrary = std::abs(center[1]) < 0.9 ? Vec<double, 3>{0.0, 1.0, 0.0}
                                                              : Vec<double, 3>{1.0, 0.0, 0.0};
        u = Vec<double, 3>{arbitrary.cross(center).normalized()};
        v = Vec<double, 3>{center.cross(u).normalized()};
    };

    // Real, reported overcorrection from an earlier version, not a
    // subjective nitpick: a single ang_radius ~0.55-0.7pi "hero" galaxy
    // covers so much of the sphere that from any reasonably close
    // vantage point the camera sits INSIDE its footprint, seeing only a
    // near-uniform wash of its own tint -- reads as "space is too
    // bright" and "the nebula has no shape", not as a recognizable
    // galaxy. Also every structure used the same narrow high-R,G,B
    // "pale lavender" tint range regardless of RNG draw, so different
    // structures were indistinguishable in hue -- fixed by drawing an
    // actual random hue per structure (hsv_to_rgb255()) instead.
    // Smaller footprints (0.10-0.30 pi, comparable to how large the Moon
    // or a real nebula photo subtends) plus a lower opacity ceiling
    // (0.5x) keep dark space reading as dark space, with recognizable
    // colored structures in it rather than a colored sky replacing it.
    // arms=2 spirals read as a yin-yang-like swirl -- fine as an
    // occasional accent, cluttered/busy if several overlap in the same
    // view at once (blackhole_demo's wide lensed FOV wraps much more of
    // the sky into one frame than wormhole_demo's narrow one, so it
    // shows this first). Kept modest (2-3, not 3-6) for that reason.
    // tightness also lowered (was 2.5-5.5): high angular frequency in
    // the arm pattern aliases into small fragmented-looking blobs right
    // where gravitational lensing is most extreme (near the photon ring
    // -- adjacent screen pixels there can map to wildly different sky
    // directions), a real reported artifact ("наплывы... поломанней
    // ближе к краю"), not a subjective critique. A gentler winding stays
    // recognizable as a spiral without needing supersampling to survive
    // that magnification.
    for (int i = 0, n = 2 + static_cast<int>(2.0 * u01(rng)); i < n; ++i) {
        auto c = random_sky_dir(rng);
        Vec<double, 3> u, v;
        tangent_frame(c, u, v);
        double hue = u01(rng);
        sky.spirals.push_back({c, u, v, (0.10 + 0.20 * u01(rng)) * pi, 1.0 + 1.0 * u01(rng),
                                u01(rng) < 0.5 ? 2 : 3,
                                hsv_to_rgb255(hue, 0.55 + 0.3 * u01(rng), 0.85 + 0.15 * u01(rng))});
    }
    for (int i = 0, n = 4 + static_cast<int>(3.0 * u01(rng)); i < n; ++i) {
        double hue = u01(rng);
        sky.clouds.push_back({random_sky_dir(rng), (0.08 + 0.14 * u01(rng)) * pi,
                               4.0 + 4.0 * u01(rng), u01(rng) * 2.0 * pi,
                               hsv_to_rgb255(hue, 0.45 + 0.35 * u01(rng), 0.75 + 0.2 * u01(rng))});
    }

    // Stars: small angular size relative to a typical render's own
    // pixel resolution (checked directly against blackhole_demo's
    // 1920x1080 @ FOV_HALF=0.45rad -- about 0.00047 rad/pixel -- so a
    // star's radius stays sub-pixel except where lensing magnifies it,
    // which is the one case it SHOULD grow, correctly).
    sky.stars.resize(static_cast<std::size_t>(n_stars));
    sky.star_buckets.resize(static_cast<std::size_t>(Sky::kThetaBuckets) * Sky::kPhiBuckets);
    for (int i = 0; i < n_stars; ++i) {
        Vec<double, 3> d = random_sky_dir(rng);
        double brightness = 0.35 + 0.65 * u01(rng);
        double ang_size = 0.00025 + 0.00035 * u01(rng) * (0.6 + 0.4 * brightness);
        sky.stars[static_cast<std::size_t>(i)] = {d, ang_size, brightness};

        double theta = std::acos(std::clamp(d[1], -1.0, 1.0));
        double phi = std::atan2(d[2], d[0]);
        int tb = std::clamp(static_cast<int>(theta / pi * Sky::kThetaBuckets), 0,
                             Sky::kThetaBuckets - 1);
        int pb = std::clamp(static_cast<int>((phi + pi) / (2.0 * pi) * Sky::kPhiBuckets), 0,
                             Sky::kPhiBuckets - 1);
        sky.star_buckets[static_cast<std::size_t>(tb) * Sky::kPhiBuckets + pb]
            .push_back(static_cast<std::uint32_t>(i));
    }
    return sky;
}

inline Vec<double, 3> sample_sky_color(const Sky& sky, const Vec<double, 3>& dir_in) {
    double r = dir_in.norm();
    if (r < 1e-12) return {0.0, 0.0, 0.0};
    Vec<double, 3> dir{dir_in / r};
    constexpr double pi = std::numbers::pi;

    double theta = std::acos(std::clamp(dir[1], -1.0, 1.0));
    double phi = std::atan2(dir[2], dir[0]);

    // Background: whole-sky gradient (brighter near mid-latitude, same
    // shape the old per-pixel texture fill used) when wide_sky, else a
    // flat fill at `tint` -- see Sky::wide_sky's doc comment for why a
    // close, narrow-FOV shot needs the flat variant.
    Vec<double, 3> color = sky.tint;
    if (sky.wide_sky) {
        double vgrad = std::abs(theta / pi - 0.5) * 2.0;
        double bg_bright = std::pow(1.0 - vgrad, 1.5);
        color = Vec<double, 3>{sky.tint[0] + 12.0 * bg_bright, sky.tint[1] + 8.0 * bg_bright,
                                sky.tint[2] + 22.0 * bg_bright};
    }

    for (auto& s : sky.spirals) {
        double ang = std::acos(std::clamp(dir.dot(s.center), -1.0, 1.0));
        if (ang > s.ang_radius) continue;
        double local_angle = std::atan2(dir.dot(s.v), dir.dot(s.u));
        double arm_phase = local_angle * s.arms - s.tightness * std::log(ang * 3.0 + 1.0);
        double density = 0.5 + 0.5 * std::cos(arm_phase);
        double u_norm = ang / s.ang_radius;
        double falloff = std::exp(-u_norm * u_norm * 2.0);
        color = Vec<double, 3>{color + s.tint * (density * falloff)};
    }
    for (auto& c : sky.clouds) {
        double ang = std::acos(std::clamp(dir.dot(c.center), -1.0, 1.0));
        if (ang > c.ang_radius) continue;
        double noise =
            std::sin(dir[0] * c.freq + c.phase) * std::sin(dir[2] * c.freq * 1.3 - c.phase) +
            0.5 * std::sin((dir[0] + dir[1]) * c.freq * 0.7 + c.phase * 0.6) +
            0.3 * std::sin((dir[1] - dir[2]) * c.freq * 0.9 - c.phase * 0.4);
        double density = std::clamp(0.5 + 0.35 * noise, 0.0, 1.0);
        double u_norm = ang / c.ang_radius;
        double falloff = std::exp(-u_norm * u_norm * 1.5);
        color = Vec<double, 3>{color + c.tint * (density * falloff * 0.7)};
    }

    if (!sky.star_buckets.empty()) {
        int tb = std::clamp(static_cast<int>(theta / pi * Sky::kThetaBuckets), 0,
                             Sky::kThetaBuckets - 1);
        int pb = static_cast<int>((phi + pi) / (2.0 * pi) * Sky::kPhiBuckets);
        pb = ((pb % Sky::kPhiBuckets) + Sky::kPhiBuckets) % Sky::kPhiBuckets;

        double best = 0.0;
        for (int dtb = -1; dtb <= 1; ++dtb) {
            int tbb = tb + dtb;
            if (tbb < 0 || tbb >= Sky::kThetaBuckets) continue;
            for (int dpb = -1; dpb <= 1; ++dpb) {
                int pbb = ((pb + dpb) % Sky::kPhiBuckets + Sky::kPhiBuckets) % Sky::kPhiBuckets;
                for (std::uint32_t idx : sky.star_buckets[static_cast<std::size_t>(tbb) * Sky::kPhiBuckets +
                                                           pbb]) {
                    const Star& s = sky.stars[idx];
                    double ang = std::acos(std::clamp(dir.dot(s.dir), -1.0, 1.0));
                    if (ang > s.ang_size * 3.0) continue;
                    double falloff = std::exp(-(ang * ang) / (s.ang_size * s.ang_size));
                    best = std::max(best, s.brightness * falloff);
                }
            }
        }
        if (best > 1e-3) {
            double v = best * 255.0;
            color = Vec<double, 3>{std::max(color[0], v), std::max(color[1], v),
                                    std::max(color[2], std::min(255.0, v * 1.05))};
        }
    }

    return color;
}

inline void sample_sky(const Sky& sky, const Vec<double, 3>& dir_in, std::uint8_t out[3]) {
    Vec<double, 3> color = sample_sky_color(sky, dir_in);
    out[0] = static_cast<std::uint8_t>(std::clamp(color[0], 0.0, 255.0));
    out[1] = static_cast<std::uint8_t>(std::clamp(color[1], 0.0, 255.0));
    out[2] = static_cast<std::uint8_t>(std::clamp(color[2], 0.0, 255.0));
}

} // namespace spatium::render
