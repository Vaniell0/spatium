// Phase 4: the actual final-video render, on GPU. Matches
// examples/blackhole_gr_demo.cpp's flyby choreography exactly (same
// camera_position_for_t/camera_for_frame, including the departure
// look-blend fix so the hole actually leaves frame by the end) plus
// the SPATIUM star-text watermark, so this is the CPU pipeline's exact
// visual design, just executed per-ray on the T4 instead of the CPU.
//
// GPU does the expensive part (geodesic integration + volumetric disk
// emission-absorption, one thread per ray, AAxAA supersampled). The
// host does the cheap O(W*H) part exactly like main()'s Kerr path
// does: sky lookup (examples/io_helpers.hpp's real procedural Sky,
// unmodified), apply_bloom, tone_map_to_8bit -- transcribed from
// blackhole_gr_demo.cpp (not shared headers) to match this session's
// version exactly. Usage: render_4k_frame WIDTH HEIGHT [start_frame]
// [end_frame_exclusive], writing frames_gpu/frame_%04d.png per frame.

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../examples/io_helpers.hpp"

#include "disk_physics.hpp"
#include "kerr_geometry.hpp"
#include "render_kernel.h"

#include <spatium/algebra/vector.hpp>
#include <spatium/physics/relativity/geodesic.hpp>
#include <spatium/physics/relativity/kerr.hpp>
#include <spatium/render/camera.hpp>
#include <spatium/render/write_image.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <numbers>
#include <string>
#include <vector>

using namespace spatium;
using namespace spatium::physics::relativity;
using namespace spatium::render;
using namespace spatium::examples;

namespace {

// --- choreography constants, exactly matching blackhole_gr_demo.cpp's
// defaults (M_BH=1, SPIN=0.9, PROGRADE=true, N_FRAMES=150) ---
constexpr double kMass = 1.0;
constexpr double kDefaultSpin = 0.9;
constexpr bool kPrograde = true;
constexpr double kDApproachStart = 45.0, kDPeak = 14.0, kDDepartEnd = 30.0;
constexpr double kElevStartDeg = 12.0, kElevPeakDeg = 30.0, kAzimuthTotalDeg = 80.0;
constexpr double kFocusBlend = 0.22, kFocusFreezeT = 0.12;
constexpr double kFovDeg = 40.0;
constexpr int kNFrames = 150;
constexpr double kDiskTimePerFrame = 2.0;  // see blackhole_gr_demo.cpp's own comment
constexpr double kRDiskOuter = 16.0;
constexpr double kDiskAspectRatio = 0.22, kDiskDensityPower = 2.0;
constexpr double kDiskT0 = 9000.0, kDiskTRef = 6500.0, kDopplerDamping = 0.6;
constexpr double kAbsorption = 0.7, kEmissionScale = 3.0, kTransmittanceCutoff = 1e-4;
constexpr double kBaseDlambda = 0.15, kDlambdaFloorFrac = 0.10;
constexpr int kMaxSteps = 3000;

double smoothstep01(double t) {
    t = std::clamp(t, 0.0, 1.0);
    return t * t * (3.0 - 2.0 * t);
}

// Symmetric smooth "flash" envelope: 0 at |t-center|>=half_width, easing
// via smoothstep up to 1 exactly at center. Used to make the hero
// SPATIUM text (only genuinely legible in ONE target frame -- see
// add_text_stars_screen_space's comment) visibly ignite and fade around
// that frame instead of sitting at constant brightness for the whole
// sequence, where in every other frame those same stars are just
// scattered (not shaped into text) by construction.
double flash_pulse(double t, double center, double half_width) {
    double d = std::abs(t - center);
    double x = 1.0 - std::clamp(d / half_width, 0.0, 1.0);
    return smoothstep01(x);
}

// Truly off outside the flash window -- not a dim baseline (0.3 was
// still visibly "on", same order as ordinary background stars).
constexpr double kHeroMinBrightness = 0.0, kHeroMaxBrightness = 2.5;
// Tight on purpose: the hero text is only geometrically legible in ONE
// target frame (see add_text_stars_screen_space) -- everywhere else
// those stars are just scattered, so the flash should read as a single
// moment igniting at that frame, not a multi-second glow surrounding
// it. ~0.01 normalized time is ~7-8 frames each side at 750 frames
// (under half a second at 30fps).
constexpr double kHeroFlashHalfWidth = 0.01;

// Side/legible copy: fades in from 0 over the first kSideFadeInEndT of
// the sequence (approach only, not a flash -- it stays lit afterward),
// so the frame-0-through-early-approach shots (where this sky direction
// is heavily lensed -- see the fade-in comment at its call site) don't
// show it at full strength while still badly distorted.
// Not capped to ordinary-star range: invisible at the start (brightness
// 0, before the fade-in) is what matters, not dimming it once lit --
// at full brightness it should read clearly as text, same as it always
// did.
constexpr double kSideMaxBrightness = 1.5;
constexpr double kSideFadeInEndT = 0.08;

Vec<double, 3> camera_position_for_t(double t) {
    double ease = smoothstep01(t);
    double D = t < 0.4 ? kDApproachStart + (kDPeak - kDApproachStart) * smoothstep01(t / 0.4)
                       : kDPeak + (kDDepartEnd - kDPeak) * smoothstep01((t - 0.4) / 0.6);
    double elev = (kElevStartDeg + (kElevPeakDeg - kElevStartDeg) * ease) * std::numbers::pi / 180.0;
    double az = kAzimuthTotalDeg * ease * std::numbers::pi / 180.0;
    return Vec<double, 3>{D * std::cos(elev) * std::cos(az), D * std::cos(elev) * std::sin(az),
                           D * std::sin(elev)};
}

Vec<double, 3> motion_direction_for_t(double t) {
    constexpr double eps = 1e-3;
    Vec<double, 3> p0 = camera_position_for_t(std::clamp(t - eps, 0.0, 1.0));
    Vec<double, 3> p1 = camera_position_for_t(std::clamp(t + eps, 0.0, 1.0));
    return Vec<double, 3>{(p1 - p0).normalized()};
}

constexpr double kDepartureLookStartT = 0.6;

Camera<double> camera_for_frame(int frame, int n_frames) {
    double t = double(frame) / double(n_frames - 1);
    Vec<double, 3> pos = camera_position_for_t(t);
    Vec<double, 3> frozen_pos = camera_position_for_t(kFocusFreezeT);
    Vec<double, 3> frozen_fwd{(-frozen_pos).normalized()};
    Vec<double, 3> live_fwd{(-pos).normalized()};
    Vec<double, 3> focus_fwd{
        Vec<double, 3>{live_fwd * kFocusBlend + frozen_fwd * (1.0 - kFocusBlend)}.normalized()};

    double departure = smoothstep01((t - kDepartureLookStartT) / (1.0 - kDepartureLookStartT));
    Vec<double, 3> motion_fwd = motion_direction_for_t(t);
    Vec<double, 3> fwd{
        Vec<double, 3>{focus_fwd * (1.0 - departure) + motion_fwd * departure}.normalized()};
    Vec<double, 3> target{pos + fwd * pos.norm()};

    // Tried tilting `up` off the exact spin axis (1deg->30deg over the
    // sequence) on the theory that up.e_phi==0's exact-alignment coincidence
    // (see kerr_render_kernel.cu's degeneracy comment) causes the visible
    // vertical seam -- tested directly at 1deg (frame 0) and ~11deg (frame
    // 300): the seam is UNCHANGED at both. Disproven, not just unhelpful --
    // reverted rather than left in, since it also changes the whole
    // sequence's framing (a camera roll) for zero benefit. The seam's real
    // mechanism is still not understood (this is the fifth ruled-out
    // hypothesis: aliasing, disk texture, max-steps fallback, n_ph epsilon,
    // now this).
    return Camera<double>{.position = pos, .target = target, .up = {0.0, 0.0, 1.0},
                           .fov_deg = kFovDeg};
}

struct SpheroidalPoint {
    double r, theta, phi;
};

SpheroidalPoint to_spheroidal(const Vec<double, 3>& pos, double a) {
    double X = pos[0], Y = pos[1], Z = pos[2];
    double r2sum = X * X + Y * Y + Z * Z;
    double a2 = a * a;
    double r2 = 0.5 * (r2sum - a2) + std::sqrt(0.25 * (r2sum - a2) * (r2sum - a2) + a2 * Z * Z);
    double r = std::sqrt(std::max(r2, 1e-12));
    double theta = std::acos(std::clamp(Z / r, -1.0, 1.0));
    double phi = std::atan2(Y, X);
    return {r, theta, phi};
}

void apply_bloom(std::vector<double>& img, int w, int h) {
    constexpr double THRESHOLD = 190.0;  // see blackhole_gr_demo.cpp's own comment
    constexpr int RADIUS = 6;
    constexpr double STRENGTH = 0.35;

    std::vector<double> bright(img.size(), 0.0);
    for (std::size_t i = 0; i < img.size() / 3; ++i) {
        double r = img[3 * i], g = img[3 * i + 1], b = img[3 * i + 2];
        double lum = 0.299 * r + 0.587 * g + 0.114 * b;
        double scale = std::max(0.0, lum - THRESHOLD) / (lum + 1e-6);
        bright[3 * i] = r * scale; bright[3 * i + 1] = g * scale; bright[3 * i + 2] = b * scale;
    }

    std::vector<double> tmp(bright.size());
    auto box_blur_pass = [&] {
        for (int y = 0; y < h; ++y)
            for (int x = 0; x < w; ++x) {
                double sum[3] = {0.0, 0.0, 0.0};
                int count = 0;
                for (int dx = -RADIUS; dx <= RADIUS; ++dx) {
                    int xx = x + dx;
                    if (xx < 0 || xx >= w) continue;
                    std::size_t idx = 3 * (static_cast<std::size_t>(y) * w + xx);
                    sum[0] += bright[idx]; sum[1] += bright[idx + 1]; sum[2] += bright[idx + 2];
                    ++count;
                }
                std::size_t idx = 3 * (static_cast<std::size_t>(y) * w + x);
                tmp[idx] = sum[0] / count; tmp[idx + 1] = sum[1] / count; tmp[idx + 2] = sum[2] / count;
            }
        std::swap(bright, tmp);
        for (int x = 0; x < w; ++x)
            for (int y = 0; y < h; ++y) {
                double sum[3] = {0.0, 0.0, 0.0};
                int count = 0;
                for (int dy = -RADIUS; dy <= RADIUS; ++dy) {
                    int yy = y + dy;
                    if (yy < 0 || yy >= h) continue;
                    std::size_t idx = 3 * (static_cast<std::size_t>(yy) * w + x);
                    sum[0] += bright[idx]; sum[1] += bright[idx + 1]; sum[2] += bright[idx + 2];
                    ++count;
                }
                std::size_t idx = 3 * (static_cast<std::size_t>(y) * w + x);
                tmp[idx] = sum[0] / count; tmp[idx + 1] = sum[1] / count; tmp[idx + 2] = sum[2] / count;
            }
        std::swap(bright, tmp);
    };
    box_blur_pass();

    for (std::size_t i = 0; i < img.size(); ++i) img[i] += STRENGTH * bright[i];
}

std::vector<std::uint8_t> tone_map_to_8bit(const std::vector<double>& hdr) {
    constexpr double KNEE = 200.0;
    constexpr double HEADROOM = 255.0 - KNEE;
    constexpr double SCALE = 400.0;
    std::vector<std::uint8_t> out(hdr.size());
    for (std::size_t i = 0; i < hdr.size(); ++i) {
        double c = hdr[i];
        double mapped = c <= KNEE ? c : KNEE + HEADROOM * (1.0 - std::exp(-(c - KNEE) / SCALE));
        out[i] = static_cast<std::uint8_t>(std::clamp(mapped, 0.0, 255.0));
    }
    return out;
}

// "SPATIUM" watermark -- transcribed from blackhole_gr_demo.cpp (not a
// shared header there either). See that file for the full rationale;
// unchanged here.
struct Glyph5x7 { std::uint8_t rows[7]; };

const Glyph5x7* find_glyph(char c) {
    static const Glyph5x7 kS{{0b01111, 0b10000, 0b10000, 0b01110, 0b00001, 0b00001, 0b11110}};
    static const Glyph5x7 kP{{0b11110, 0b10001, 0b10001, 0b11110, 0b10000, 0b10000, 0b10000}};
    static const Glyph5x7 kA{{0b01110, 0b10001, 0b10001, 0b11111, 0b10001, 0b10001, 0b10001}};
    static const Glyph5x7 kT{{0b11111, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100}};
    static const Glyph5x7 kI{{0b11111, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100, 0b11111}};
    static const Glyph5x7 kU{{0b10001, 0b10001, 0b10001, 0b10001, 0b10001, 0b10001, 0b01110}};
    static const Glyph5x7 kM{{0b10001, 0b11011, 0b10101, 0b10001, 0b10001, 0b10001, 0b10001}};
    switch (c) {
        case 'S': return &kS;
        case 'P': return &kP;
        case 'A': return &kA;
        case 'T': return &kT;
        case 'I': return &kI;
        case 'U': return &kU;
        case 'M': return &kM;
        default: return nullptr;
    }
}

void add_star_to_sky(Sky& sky, const Vec<double, 3>& dir, double ang_size, double brightness) {
    sky.stars.push_back({dir, ang_size, brightness});
    int idx = static_cast<int>(sky.stars.size()) - 1;
    constexpr double pi = std::numbers::pi;
    double theta = std::acos(std::clamp(dir[1], -1.0, 1.0));
    double phi = std::atan2(dir[2], dir[0]);
    int tb = std::clamp(static_cast<int>(theta / pi * Sky::kThetaBuckets), 0, Sky::kThetaBuckets - 1);
    int pb = std::clamp(static_cast<int>((phi + pi) / (2.0 * pi) * Sky::kPhiBuckets), 0,
                         Sky::kPhiBuckets - 1);
    sky.star_buckets[static_cast<std::size_t>(tb) * Sky::kPhiBuckets + static_cast<std::size_t>(pb)]
        .push_back(idx);
}

void add_text_stars(Sky& sky, const std::string& text, const Vec<double, 3>& center_dir,
                     double cell_ang, double brightness) {
    Vec<double, 3> c{center_dir.normalized()};
    Vec<double, 3> arbitrary =
        std::abs(c[1]) < 0.9 ? Vec<double, 3>{0.0, 1.0, 0.0} : Vec<double, 3>{1.0, 0.0, 0.0};
    Vec<double, 3> u{arbitrary.cross(c).normalized()};
    Vec<double, 3> v{c.cross(u).normalized()};

    double total_w = static_cast<double>(text.size()) * 6.0 - 1.0;
    double x0 = -total_w * 0.5;
    for (std::size_t li = 0; li < text.size(); ++li) {
        const Glyph5x7* g = find_glyph(text[li]);
        if (!g) continue;
        for (int row = 0; row < 7; ++row) {
            for (int col = 0; col < 5; ++col) {
                if (!((g->rows[row] >> (4 - col)) & 1)) continue;
                double dx = (x0 + static_cast<double>(li) * 6.0 + col) * cell_ang;
                double dy = (3.0 - row) * cell_ang;
                Vec<double, 3> raw{c + u * dx + v * dy};
                add_star_to_sky(sky, raw.normalized(), cell_ang * 0.55, brightness);
            }
        }
    }
}

// Solves, for ONE specific screen-space ray in ONE specific target
// frame, the actual sky direction that ray's geodesic reaches on
// escape -- using the real Dual<T>-exact CPU engine (geodesic_step()),
// not the closed-form GPU kernel's own math, though Phase 2/3 already
// verified those two agree to ~1e-13-1e-16 relative, so a star placed
// exactly at this solved direction reproduces this same screen pixel
// when the GPU kernel later retraces the identical ray for real. This
// is the inverse of what add_text_stars() does: that function commits
// to a shape in the SOURCE sky and accepts whatever the real lensing
// map does to it (arcs, wrap-around, whatever falls out); this
// function commits to a shape on the TARGET SCREEN for one frame and
// solves backward for where each dot must actually live in the sky --
// legible in that one frame by construction, not by accident.
// Returns false (skip that dot) if the ray is captured or doesn't
// resolve within the step budget.
bool solve_kerr_exit_direction(const Camera<double>& cam, const CameraBasis<double>& basis,
                                double sx, double sy, double mass, double spin,
                                Vec<double, 3>& out_dir) {
    Vec<double, 3> dir_world = camera_ray_dir(basis, sx, sy);

    SpheroidalPoint cam_sp = to_spheroidal(cam.position, spin);
    auto tet = spatium::gpu::spheroidal_tetrad(cam_sp.r, cam_sp.theta, cam_sp.phi, spin);
    Vec<double, 3> tet_e_r{tet.e_r.x, tet.e_r.y, tet.e_r.z};
    Vec<double, 3> tet_e_theta{tet.e_theta.x, tet.e_theta.y, tet.e_theta.z};
    Vec<double, 3> tet_e_phi{tet.e_phi.x, tet.e_phi.y, tet.e_phi.z};
    double g_rr, g_thth, g_phph;
    spatium::gpu::kerr_metric_diag(mass, spin, cam_sp.r, cam_sp.theta, &g_rr, &g_thth, &g_phph);
    double sigma =
        cam_sp.r * cam_sp.r + spin * spin * std::cos(cam_sp.theta) * std::cos(cam_sp.theta);
    double g_tt = -(1.0 - 2.0 * mass * cam_sp.r / sigma);
    double ut_cam = 1.0 / std::sqrt(-g_tt);

    double n_r = dir_world.dot(tet_e_r), n_th = dir_world.dot(tet_e_theta),
           n_ph = dir_world.dot(tet_e_phi);
    // See kerr_render_kernel.cu's matching comment: up.e_phi==0
    // identically, so there's always one exact-n_ph==0 screen column
    // (a real degenerate-orbit separatrix, not a bug) that renders
    // badly under finite step budgets/precision. Nudged on n_ph (an
    // O(1) unit-vector dot product, distance-independent), not on
    // uphi itself (uphi's own scale shrinks with camera distance, so a
    // fixed epsilon there clamped far more than just the degenerate
    // column at wide shots -- caught before the long render committed
    // to it).
    constexpr double kDegeneracyEpsilon = 1e-4;
    if (std::abs(n_ph) < kDegeneracyEpsilon) n_ph = kDegeneracyEpsilon;
    double ur = n_r / std::sqrt(g_rr), utheta = n_th / std::sqrt(g_thth),
           uphi = n_ph / std::sqrt(g_phph);

    KerrMetric<double> metric{mass, spin};
    Vec<double, 8> state{0.0, cam_sp.r, cam_sp.theta, cam_sp.phi, ut_cam, ur, utheta, uphi};
    double D_cam = cam.position.norm();
    double r_horizon = kerr_outer_horizon_radius(mass, spin);
    double r_capture = r_horizon * 1.15, r_escape = 40.0;

    for (int step = 0; step < 3000; ++step) {
        double r = state[1];
        double dl = 0.15 * std::max(0.10, r / D_cam);
        state = geodesic_step(metric, state, dl);
        double r_new = state[1];
        if (!std::isfinite(r_new) || r_new <= r_capture) return false;
        if (r_new >= r_escape && state[5] > 0.0) {
            double r_e = state[1], theta_e = state[2], phi_e = state[3];
            double g_rr_e, g_thth_e, g_phph_e;
            spatium::gpu::kerr_metric_diag(mass, spin, r_e, theta_e, &g_rr_e, &g_thth_e, &g_phph_e);
            double nr_e = state[5] * std::sqrt(g_rr_e), nth_e = state[6] * std::sqrt(g_thth_e),
                   nph_e = state[7] * std::sqrt(g_phph_e);
            double mag = std::sqrt(nr_e * nr_e + nth_e * nth_e + nph_e * nph_e);
            nr_e /= mag; nth_e /= mag; nph_e /= mag;
            auto tet_e = spatium::gpu::spheroidal_tetrad(r_e, theta_e, phi_e, spin);
            out_dir = Vec<double, 3>{tet_e.e_r.x * nr_e + tet_e.e_theta.x * nth_e + tet_e.e_phi.x * nph_e,
                                      tet_e.e_r.y * nr_e + tet_e.e_theta.y * nth_e + tet_e.e_phi.y * nph_e,
                                      tet_e.e_r.z * nr_e + tet_e.e_theta.z * nth_e + tet_e.e_phi.z * nph_e};
            return true;
        }
    }
    return false;
}

// Places `text` as a dot-matrix directly in SCREEN space for
// `target_frame`, solving each on dot backward to its real sky
// position via solve_kerr_exit_direction() -- see that function's own
// comment for why this is the inverse of add_text_stars(). (sx0, sy0)
// is the block's screen-space center in camera_ray_dir()'s own NDC
// convention; cell_ndc is one font cell's width in that same space.
int add_text_stars_screen_space(Sky& sky, const std::string& text, int target_frame, int n_frames,
                                 double mass, double spin, double sx0, double sy0, double cell_ndc,
                                 double brightness) {
    Camera<double> cam = camera_for_frame(target_frame, n_frames);
    CameraBasis<double> basis = make_camera_basis(cam);

    double total_w = static_cast<double>(text.size()) * 6.0 - 1.0;
    double x0 = -total_w * 0.5;
    int placed = 0;
    for (std::size_t li = 0; li < text.size(); ++li) {
        const Glyph5x7* g = find_glyph(text[li]);
        if (!g) continue;
        for (int row = 0; row < 7; ++row) {
            for (int col = 0; col < 5; ++col) {
                if (!((g->rows[row] >> (4 - col)) & 1)) continue;
                double sx = sx0 + (x0 + static_cast<double>(li) * 6.0 + col) * cell_ndc;
                double sy = sy0 + (3.0 - row) * cell_ndc;
                Vec<double, 3> star_dir;
                if (!solve_kerr_exit_direction(cam, basis, sx, sy, mass, spin, star_dir)) continue;
                // cell_ndc is a tan()-scaled screen offset (same
                // convention as camera_ray_dir's sx/sy args), so for
                // the small values used here it's already a good
                // small-angle approximation of the actual angular
                // size in radians -- no extra tan_half factor needed.
                add_star_to_sky(sky, star_dir.normalized(), cell_ndc * 0.6, brightness);
                ++placed;
            }
        }
    }
    return placed;
}

std::vector<spatium::gpu::DiscParticleParams> make_disc_particles(double r_disk_inner) {
    constexpr std::size_t N = 16;
    std::vector<spatium::gpu::DiscParticleParams> particles(N);
    for (std::size_t i = 0; i < N; ++i) {
        double t = double(i) / double(N);
        double r_center = r_disk_inner + (kRDiskOuter - r_disk_inner) * (0.1 + 0.8 * t);
        double ecc = 0.10 + 0.20 * std::abs(std::sin(t * 11.0));
        double r_min = r_center * (1.0 - ecc);
        double r_max = r_center * (1.0 + ecc);
        particles[i] = {1.0 / r_max, 1.0 / r_min, t * 2.0 * std::numbers::pi * 5.0,
                         1.5 + 3.0 * std::fmod(t * 7.0, 1.0)};
    }
    return particles;
}

}  // namespace

int main(int argc, char** argv) {
    int W = argc > 1 ? std::atoi(argv[1]) : 1920;
    int H = argc > 2 ? std::atoi(argv[2]) : 1080;
    int start_frame = argc > 3 ? std::atoi(argv[3]) : 0;
    int n_frames = argc > 4 ? std::atoi(argv[4]) : kNFrames;  // total sequence length
    int end_frame = argc > 5 ? std::atoi(argv[5]) : n_frames;  // exclusive
    double spin = argc > 6 ? std::atof(argv[6]) : kDefaultSpin;

    double r_horizon0 = kerr_outer_horizon_radius(kMass, spin);
    double r_disk_inner0 = kerr_isco_radius(kMass, spin, kPrograde);
    std::vector<spatium::gpu::DiscParticleParams> particles = make_disc_particles(r_disk_inner0);

    // 40000 (up from 10000) -- denser background shrinks the gaps
    // between adjacent stars' lensed images near the photon sphere
    // (found, this session, to be a real discreteness artifact, not
    // aliasing or a disk effect -- confirmed by a control render with
    // the disk's emission/absorption both zeroed, where the same thin
    // dark line persisted unchanged).
    //
    // Sky (incl. SPATIUM watermark) is built once and reused across
    // every frame, exactly matching main()'s own construction in
    // blackhole_gr_demo.cpp.
    Sky sky = make_starfield(40000, 7, {8.0, 6.0, 18.0});
    // Set below; consumed per-frame in the main loop to fade this copy in
    // from 0 over the first several frames of the approach, rather than
    // sitting at full brightness from frame 0 -- at D_cam=45 the lensing
    // near this sky direction stretches the text across a huge chunk of
    // the frame (found by inspecting frame 0 directly), so showing it at
    // full strength before the camera has closed in reads as a giant
    // scattered blob, not a legible word.
    std::size_t side_star_begin = 0, side_star_end = 0;
    {
        Camera<double> cam0 = camera_for_frame(0, n_frames);
        Vec<double, 3> fwd0{Vec<double, 3>{cam0.target - cam0.position}.normalized()};
        Vec<double, 3> arbitrary0 =
            std::abs(fwd0[1]) < 0.9 ? Vec<double, 3>{0.0, 1.0, 0.0} : Vec<double, 3>{1.0, 0.0, 0.0};
        Vec<double, 3> u0{arbitrary0.cross(fwd0).normalized()};
        Vec<double, 3> v0{fwd0.cross(u0).normalized()};
        Vec<double, 3> side_dir{Vec<double, 3>{fwd0 + u0 * 0.22 + v0 * 0.05}.normalized()};
        side_star_begin = sky.stars.size();
        add_text_stars(sky, "SPATIUM", side_dir, 0.0008, kSideMaxBrightness);
        side_star_end = sky.stars.size();
    }
    // Set by the hero-text block below; consumed per-frame in the main
    // loop to flash these specific stars' brightness up around the one
    // frame they're actually shaped into legible text for.
    std::size_t hero_star_begin = 0, hero_star_end = 0;
    double hero_target_t = 0.0;
    {
        // Hero copy: placed in SCREEN space for t~=0.75 (well into
        // departure -- the hole is still in frame but off-center,
        // departure look-blend ~=0.32 -- rather than the peak t=0.4)
        // via solve_kerr_exit_direction()'s inverse lensing -- see that
        // function's own comment. Legible in that exact frame by
        // construction; unlike the forward-placed legible copy above,
        // this one is solved, not guessed. Frame index computed from a
        // fixed t fraction (not a hardcoded frame number) so it lands
        // on the same point of the choreography regardless of n_frames.
        int target_frame = static_cast<int>(std::round(0.6 * (n_frames - 1)));
        hero_target_t = n_frames > 1 ? double(target_frame) / double(n_frames - 1) : 0.0;
        hero_star_begin = sky.stars.size();
        int n_placed = add_text_stars_screen_space(sky, "SPATIUM", target_frame, n_frames, kMass,
                                                     spin, 0.15, 0.15, 0.012, kHeroMaxBrightness);
        hero_star_end = sky.stars.size();
        std::printf("hero SPATIUM: frame %d/%d, %d letter-dots resolved to real sky positions\n",
                    target_frame, n_frames, n_placed);
    }

    std::filesystem::create_directories("frames_gpu");

    for (int frame = start_frame; frame < end_frame; ++frame) {
    Camera<double> cam = camera_for_frame(frame, n_frames);
    CameraBasis<double> basis = make_camera_basis(cam);
    double D_cam = cam.position.norm();
    // Scaled by normalized progress (not frame*kDiskTimePerFrame)
    // so the disk's own orbital advance over the FULL sequence stays
    // the originally-tuned ~90 time units (see kDiskTimePerFrame's own
    // comment: "half an orbit over the full sequence") regardless of
    // how many frames that sequence is stretched across.
    double t_norm = n_frames > 1 ? double(frame) / double(n_frames - 1) : 0.0;
    double disk_time = t_norm * kDiskTimePerFrame * (kNFrames - 1);

    // Flash the hero SPATIUM stars up from a faint baseline to full
    // brightness right as this frame crosses their one legible target
    // frame, then back down -- see flash_pulse's comment. sky is shared
    // across the whole render, so this mutates brightness in place
    // before sample_sky_color() reads it later this iteration.
    if (hero_star_end > hero_star_begin) {
        double flash = flash_pulse(t_norm, hero_target_t, kHeroFlashHalfWidth);
        double hero_brightness = kHeroMinBrightness + (kHeroMaxBrightness - kHeroMinBrightness) * flash;
        for (std::size_t si = hero_star_begin; si < hero_star_end; ++si)
            sky.stars[si].brightness = hero_brightness;
    }
    if (side_star_end > side_star_begin) {
        double side_brightness = kSideMaxBrightness * smoothstep01(t_norm / kSideFadeInEndT);
        for (std::size_t si = side_star_begin; si < side_star_end; ++si)
            sky.stars[si].brightness = side_brightness;
    }

    SpheroidalPoint cam_sp = to_spheroidal(cam.position, spin);
    auto tet = spatium::gpu::spheroidal_tetrad(cam_sp.r, cam_sp.theta, cam_sp.phi, spin);
    double g_rr, g_thth, g_phph;
    spatium::gpu::kerr_metric_diag(kMass, spin, cam_sp.r, cam_sp.theta, &g_rr, &g_thth, &g_phph);
    double sigma = cam_sp.r * cam_sp.r + spin * spin * std::cos(cam_sp.theta) * std::cos(cam_sp.theta);
    double g_tt = -(1.0 - 2.0 * kMass * cam_sp.r / sigma);
    double g_tphi =
        -2.0 * kMass * cam_sp.r * spin * std::sin(cam_sp.theta) * std::sin(cam_sp.theta) / sigma;
    double ut_cam = 1.0 / std::sqrt(-g_tt);

    double r_capture = r_horizon0 * 1.15, r_escape = 40.0;

    std::printf("frame=%d/%d t=%.4f D_cam=%.3f cam_sp=(r=%.3f,theta=%.3f,phi=%.3f)\n", frame, n_frames,
                t_norm, D_cam, cam_sp.r, cam_sp.theta, cam_sp.phi);

    // AAxAA jittered-grid supersampling per pixel, matching
    // render::supersample_pixel_hdr()'s own convention -- the earlier
    // 1-ray-per-pixel test showed real aliasing (a jagged dotted line)
    // right where the image changes fastest per pixel, near the
    // critical curve. Each pixel's AA*AA subrays are laid out
    // contiguously so the compositing loop below can average them
    // after sky lookup, exactly like the CPU path composites full
    // colors per subray before averaging, not raw geodesic state.
    constexpr int AA = 2;
    long n_pixels = static_cast<long>(W) * H;
    long n_rays = n_pixels * AA * AA;
    std::vector<double> dirs(static_cast<std::size_t>(n_rays) * 3);
    for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x) {
            long pixel_idx = static_cast<long>(y) * W + x;
            for (int j = 0; j < AA; ++j)
                for (int i = 0; i < AA; ++i) {
                    double jx = (i + 0.5) / AA, jy = (j + 0.5) / AA;
                    Vec<double, 3> d = camera_pixel_dir(cam, basis, x, y, W, H, jx, jy);
                    long idx = pixel_idx * AA * AA + j * AA + i;
                    dirs[idx * 3 + 0] = d[0]; dirs[idx * 3 + 1] = d[1]; dirs[idx * 3 + 2] = d[2];
                }
        }

    KerrRenderParams p{};
    p.mass = kMass; p.spin = spin; p.D_cam = D_cam;
    p.cam_r = cam_sp.r; p.cam_theta = cam_sp.theta; p.cam_phi = cam_sp.phi;
    p.e_r_x = tet.e_r.x; p.e_r_y = tet.e_r.y; p.e_r_z = tet.e_r.z;
    p.e_theta_x = tet.e_theta.x; p.e_theta_y = tet.e_theta.y; p.e_theta_z = tet.e_theta.z;
    p.e_phi_x = tet.e_phi.x; p.e_phi_y = tet.e_phi.y; p.e_phi_z = tet.e_phi.z;
    p.ut_cam = ut_cam; p.sqrt_grr_cam = std::sqrt(g_rr); p.sqrt_gthth_cam = std::sqrt(g_thth);
    p.sqrt_gphph_cam = std::sqrt(g_phph);
    p.g_cam_00 = g_tt; p.g_cam_03 = g_tphi; p.g_cam_33 = g_phph;
    p.prograde = kPrograde ? 1 : 0;
    p.disk_time = disk_time;
    p.r_horizon = r_horizon0; p.r_capture = r_capture; p.r_escape = r_escape;
    p.r_disk_inner = r_disk_inner0; p.r_disk_outer = kRDiskOuter;
    p.disk_aspect_ratio = kDiskAspectRatio; p.disk_density_power = kDiskDensityPower;
    p.disk_t0 = kDiskT0; p.disk_t_ref = kDiskTRef; p.doppler_damping = kDopplerDamping;
    p.absorption = kAbsorption; p.emission_scale = kEmissionScale;
    p.transmittance_cutoff = kTransmittanceCutoff;
    p.base_dlambda = kBaseDlambda; p.dlambda_floor_frac = kDlambdaFloorFrac; p.max_steps = kMaxSteps;

    std::vector<double> gpu_out(static_cast<std::size_t>(n_rays) * 8);

    auto t0 = std::chrono::steady_clock::now();
    spatium_gpu_render_kerr(p, particles.data(), static_cast<int>(particles.size()), dirs.data(),
                             gpu_out.data(), static_cast<int>(n_rays));
    auto t1 = std::chrono::steady_clock::now();
    double gpu_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::printf("GPU kernel: %.1f ms for %ld rays (%.3f us/ray)\n", gpu_ms, n_rays,
                1000.0 * gpu_ms / double(n_rays));

    std::vector<double> img(static_cast<std::size_t>(n_pixels) * 3);
    for (long px = 0; px < n_pixels; ++px) {
        double sum[3] = {0.0, 0.0, 0.0};
        for (int s = 0; s < AA * AA; ++s) {
            long i = px * AA * AA + s;
            double ax = gpu_out[i * 8 + 0], ay = gpu_out[i * 8 + 1], az = gpu_out[i * 8 + 2];
            double transmittance = gpu_out[i * 8 + 3];
            bool captured = gpu_out[i * 8 + 7] > 0.5;
            if (!captured) {
                Vec<double, 3> exit_dir{gpu_out[i * 8 + 4], gpu_out[i * 8 + 5], gpu_out[i * 8 + 6]};
                Vec<double, 3> sky_color = sample_sky_color(sky, exit_dir);
                ax += sky_color[0] * transmittance;
                ay += sky_color[1] * transmittance;
                az += sky_color[2] * transmittance;
            }
            sum[0] += ax; sum[1] += ay; sum[2] += az;
        }
        double inv_n = 1.0 / (AA * AA);
        img[px * 3 + 0] = sum[0] * inv_n; img[px * 3 + 1] = sum[1] * inv_n;
        img[px * 3 + 2] = sum[2] * inv_n;
    }

    // The vertical degenerate-orbit seam (see kerr_render_kernel.cu's
    // comment) is left as-is here on purpose -- kept, not patched.
    apply_bloom(img, W, H);
    std::vector<std::uint8_t> img8 = tone_map_to_8bit(img);
    char path[64];
    std::snprintf(path, sizeof(path), "frames_gpu/frame_%04d.png", frame);
    write_png_rgb(path, W, H, img8);
    std::printf("wrote %s (%dx%d)\n", path, W, H);
    }  // for frame
    return 0;
}
