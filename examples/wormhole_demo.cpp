// Ellis/Morris-Thorne wormhole ray-tracer — camera flythrough.
//
// Null-geodesic backward integration against the "zero-tidal-force"
// traversable-wormhole metric — spherically symmetric, the same
// r^2*dOmega^2 angular structure as Schwarzschild, just with
// r(l) = sqrt(b0^2 + l^2) as the "radius" instead of r itself, so the
// same equatorial-plane restriction blackhole_demo.cpp relies on
// applies unchanged (angular momentum is still a conserved 3-vector,
// confining a photon started in-plane to that plane forever):
//
//     ds^2 = -dt^2 + dl^2 + r(l)^2 (dtheta^2 + sin^2(theta) dphi^2)
//
// Integrated by the AFFINE parameter lambda (proper path length, E=1),
// not by phi the way blackhole_demo integrates its u(phi) -- a real bug
// found and fixed while building this, not assumed away: the natural
// phi-parametrized second derivative,
//     d^2 l / dphi^2 = l * (2 r(l)^2 - b^2) / b^2,
// is regular at a turning point but NOT well-behaved far from the
// throat -- dl/dphi ~ l^2/b as l -> infinity (dphi/dlambda = b/r^2 -> 0
// there, so phi barely advances while l races away), which made a fixed
// dphi step wildly inaccurate exactly where the render needs to be
// clean (the "far, nearly straight line" case) and showed up as
// garbled, streaked starfields even at weak lensing. Fixed the same way
// blackhole_demo avoids the analogous problem with r itself (it
// integrates the *bounded* u=1/r, never r directly): use lambda, under
// which both l and phi stay bounded per step everywhere:
//     dl/dlambda = +-sqrt(1 - b^2/r(l)^2)      (bounded in [-1, 1])
//     dphi/dlambda = b / r(l)^2                (bounded, r(l) >= b0)
// Differentiating the first removes its own sign ambiguity at a turning
// point the same way the phi-parametrized derivation did, for ANY
// spherically-symmetric r(l), not just the sqrt form (r_of_l() below is
// a shape-adjustable generalization, dr/dl taken numerically so the
// formula never needs re-deriving by hand when the shape changes):
//     d^2 l / dlambda^2 = (b^2 / r(l)^3) * dr/dl                     (*)
// — regular everywhere, and well-behaved at large |l| too this time
// (r(l)~|l| there, dr/dl~1, so (*) ~ b^2/l^3 -> 0: free-streaming at
// infinity, exactly the expected physics). b < b0 (b0 = throat radius): the
// sqrt's argument stays positive for every l, so the ray just carries
// straight through to the far side. b >= b0: there IS a turning point
// (where r(l)=b) and the integrator decelerates l to zero and reverses
// on its own via (*) -- no manual sign-flip logic needed. phi is
// integrated alongside l via simple trapezoidal accumulation of
// dphi/dlambda, same velocity-Verlet step.
//
// Camera-ray geometry is exact here, not a "camera far away" flat-space
// approximation the way blackhole_demo's b=D*sin(alpha) is: the metric
// has a perfectly regular local orthonormal frame (e_l = d/dl,
// e_phi = (1/r(l)) d/dphi) at ANY l — including l=0, the throat itself,
// where the camera needs to be valid mid-flythrough, unlike Schwarzschild
// where r=0 isn't even reachable. b = r(l_cam)*sin(alpha),
// dl/dlambda(0) = cos(alpha), alpha = angle between the ray direction
// and the l-axis, hold exactly everywhere since there's no coordinate
// singularity to approximate around (r(l) >= b0 > 0 always).
//
// The escaping ray's direction is reconstructed in full 3D by rotating
// the pixel's own transverse basis vector by the total angle swept
// during integration (Rodrigues' rotation around the l-axis) — skipping
// this would sample the wrong point on each side's starfield for any
// ray that circles the throat before escaping, which is exactly the
// near-critical-impact-parameter "Einstein ring" case that makes this
// worth rendering correctly in the first place.
//
// Output: wormhole_frames/frame_%04d.png — assemble into video
// separately (deliberately no ffmpeg subprocess call from here):
//   ffmpeg -framerate 60 -i wormhole_frames/frame_%04d.png -pix_fmt yuv420p wormhole.mp4

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <spatium/vendor/stb_image_write.h>

#include "io_helpers.hpp"

#include <spatium/algebra/vector.hpp>
#include <spatium/render/sky.hpp>
#include <spatium/render/supersample.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <numbers>
#include <print>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

using spatium::Vec;
using spatium::render::make_starfield;
using spatium::render::sample_sky_color;
using spatium::render::Sky;

namespace {

constexpr int    W          = 960;
constexpr int    H          = 540;
constexpr double B0         = 1.0;   // throat radius
constexpr double FOV_HALF   = 0.55;  // half-angle in radians
constexpr int    MAX_STEPS  = 6000;
constexpr double DLAMBDA    = 0.05;  // affine-parameter step
constexpr double L_ESCAPE   = 60.0;  // far-field cutoff, either side
constexpr int    N_FRAMES   = 240;  // finer steps between frames for frame-by-frame viewing
constexpr double L_START    = 12.0;  // camera path: +L_START -> -L_START

// r(l)^4 = b0^4 + l^4 -- a sharper-shouldered generalization of the
// plain Ellis form (r(l)^2 = b0^2+l^2, n=2): flatter near the throat,
// steeper on the sides, a visually sharper "neck" at the same b0.
// Hardcoded to n=4 and worked in r^4 terms deliberately, not via
// std::pow(x, n) with a runtime n -- a first version used
// std::pow()+a numeric central-difference dr/dl (correct, and general
// for any n), but that runs ~24 std::pow() calls per Verlet step; measured
// directly, it turned a ~6s five-frame render into one that didn't
// finish in 10x that. r^4 = b0^4+l^4 needs only multiplies, and
// dr/dl = l^3/r^3 (the n=4 case of the general dr/dl = l^(n-1)/r^(n-1)
// result) needs at most one sqrt -- see file header for the general
// derivation this specializes.
double r4_of_l(double l) {
    double l2 = l * l;
    constexpr double b0_4 = B0 * B0 * B0 * B0;
    return b0_4 + l2 * l2;
}
double r_of_l(double l) { double r4 = r4_of_l(l); return std::sqrt(std::sqrt(r4)); }

// d^2 l/dlambda^2 = (b^2/r^3) * dr/dl = b^2*l^3/r^6 = b^2*l^3/(r4*sqrt(r4)).
double lpp_of(double l, double b) {
    double r4 = r4_of_l(l);
    double l3 = l * l * l;
    return (b * b * l3) / (r4 * std::sqrt(r4));
}

// dphi/dlambda = b/r(l)^2 = b/sqrt(r4) -- integrated alongside l by
// simple trapezoidal accumulation each Verlet step.
double phip_of(double l, double b) {
    double r4 = r4_of_l(l);
    return b / std::sqrt(r4);
}

constexpr double DISK_HALF_WIDTH = 0.12;  // thin band of l around the throat

enum class TraceResult { Escaped, HitDisk, Unresolved };

// Velocity-Verlet in the affine parameter (see file header for why
// lambda, not phi). Three outcomes: Escaped (|l| reaches L_ESCAPE on
// either side -- a traversable wormhole has no horizon, nothing is ever
// "captured" the way a black hole geodesic can be); HitDisk (crossed
// into the thin band of l around the throat -- a real object along the
// ray's own path, not a texture, so it gets properly gravitationally
// lensed for free); Unresolved if MAX_STEPS runs out first (a numerical
// safety net -- expected only within the near-critical impact-parameter
// ring, b close to b0, where the turning point sits near an unstable
// equilibrium and genuinely takes long to leave).
//
// Honest about what this actually is, not oversold as a literal
// accretion disk: every pixel's own ray stays confined to ITS OWN
// orbital plane (the (l_hat, that pixel's t_hat) plane -- what makes
// the whole per-pixel-plane technique work at all for a spherically
// symmetric metric), and "|l| < DISK_HALF_WIDTH" is checked in that same
// per-pixel plane. A real accretion disk sits in ONE fixed plane shared
// by every ray, which would need actually tracking theta(lambda) instead
// of assuming it stays at pi/2 for every pixel -- not implemented here.
// What this renders is a spherically-symmetric glowing shell hugging
// r~b0, seen from any angle as a genuine, correctly-lensed bright ring
// around the throat (visually close to the classic shot this was asked
// for), not a tilted elliptical disk.
//
// Disk detection only arms once the ray starts outside the band --
// otherwise a camera positioned inside it (mid-flythrough, near the
// throat) would trigger a false hit on its very first step.
TraceResult trace_geodesic(double l0, double lp0, double b, double& l_out, double& phi_out) {
    double l = l0, lp = lp0, phi = 0.0;
    double lpp = lpp_of(l, b);
    bool disk_armed = std::abs(l0) > DISK_HALF_WIDTH * 1.5;
    for (int i = 0; i < MAX_STEPS; ++i) {
        double l_new   = l + DLAMBDA * lp + 0.5 * DLAMBDA * DLAMBDA * lpp;
        double lpp_new = lpp_of(l_new, b);
        double lp_new  = lp + 0.5 * DLAMBDA * (lpp + lpp_new);
        double phi_new = phi + 0.5 * DLAMBDA * (phip_of(l, b) + phip_of(l_new, b));
        l = l_new; lp = lp_new; lpp = lpp_new; phi = phi_new;
        if (disk_armed && std::abs(l) <= DISK_HALF_WIDTH) {
            l_out = l; phi_out = phi;
            return TraceResult::HitDisk;
        }
        if (std::abs(l) >= L_ESCAPE) { l_out = l; phi_out = phi; return TraceResult::Escaped; }
    }
    l_out = l;
    phi_out = phi;
    return TraceResult::Unresolved;
}

void print_usage() {
    std::print("Usage: wormhole_demo [--mode flythrough|ring] [--frames N] [--force] [--help]\n"
               "  Ellis/Morris-Thorne wormhole null-geodesic renderer.\n"
               "  --mode flythrough (default)  camera flies through the throat\n"
               "  --mode ring                  single external Einstein-ring shot\n"
               "  --frames N   flythrough frame count (default 120)\n"
               "  --force      overwrite existing output files\n"
               "  --help       show this message\n"
               "  Output:      wormhole_frames/frame_%04d.png (960x540 RGB), or\n"
               "               einstein_ring.png (1600x900 RGB) in --mode ring\n");
}

// Renders one frame at an explicit camera state (l_cam, tilt away from
// -l_hat, azimuthal precession of the tilt direction) -- shared by the
// flythrough loop (which animates all three) and the standalone
// Einstein-ring shot (which just picks one good fixed vantage point).
// Parallelized the same way blackhole_demo.cpp splits rows across
// std::jthread workers.
void render_wormhole_frame(double l_cam, double tilt, double precess, const Vec<double, 3>& l_hat,
                            const Sky& sky_a, const Sky& sky_b, int width, int height,
                            double fov_half, std::vector<std::uint8_t>& img) {
    double r_cam = r_of_l(l_cam);
    double aspect = double(width) / height;
    double tan_fov = std::tan(fov_half);

    Vec<double, 3> tilt_axis{0.0, std::cos(precess), std::sin(precess)};  // always perp to l_hat
    Vec<double, 3> base_fwd{-1.0, 0.0, 0.0};
    Vec<double, 3> cam_fwd{(base_fwd * std::cos(tilt) +
                             Vec<double, 3>{tilt_axis.cross(base_fwd)} * std::sin(tilt))
                                .normalized()};
    Vec<double, 3> cam_right{cam_fwd.cross(Vec<double, 3>{0.0, 1.0, 0.0}).normalized()};
    Vec<double, 3> cam_up{cam_right.cross(cam_fwd).normalized()};

    img.assign(3 * static_cast<std::size_t>(width) * height, 0);

    auto render_rows = [&](int y0, int y1) {
        for (int y = y0; y < y1; ++y) {
            for (int x = 0; x < width; ++x) {
                std::uint8_t* px = &img[3 * (y * width + x)];

                // Per-ray color, called aa*aa times per pixel by
                // supersample_pixel (spatium/render/supersample.hpp) at
                // sub-pixel (sx, sy) -- smooths the disk edge and the
                // b=b0 transmit/reflect boundary automatically instead
                // of leaving them hard-aliased.
                auto ray_color = [&](double sx, double sy) -> Vec<double, 3> {
                    Vec<double, 3> dir = (cam_fwd + cam_right * sx + cam_up * sy).normalized();

                    double dir_dot_l = dir.dot(l_hat);
                    Vec<double, 3> t_hat = dir - l_hat * dir_dot_l;
                    double t_norm = t_hat.norm();

                    if (t_norm < 1e-9) {
                        // Radial ray: travels straight along the axis
                        // forever (b=0, no deflection) -- pick the side it
                        // points toward directly.
                        return sample_sky_color(dir_dot_l > 0.0 ? sky_a : sky_b, dir);
                    }
                    t_hat = t_hat / t_norm;

                    double b = r_cam * t_norm;      // r(l_cam) * sin(alpha)
                    double lp0 = dir_dot_l;         // dl/dlambda(0) = cos(alpha)

                    double l_final, phi_out;
                    auto result = trace_geodesic(l_cam, lp0, b, l_final, phi_out);
                    if (result == TraceResult::Unresolved) {
                        return {255.0, 0.0, 255.0};  // flags an unresolved ray
                    }
                    if (result == TraceResult::HitDisk) {
                        // Hot plasma gradient, gently textured by the
                        // crossing angle -- a first version used 8 sharp
                        // alternating light/dark bands (sin(phi_out*8)),
                        // which combined with the naturally concentric
                        // rings different impact parameters produce here
                        // into a hard yellow/black bullseye, not a disk (a
                        // real, reported "looks like a bee" problem, not a
                        // subjective nitpick). One slow, gentle wave and a
                        // brightness floor well above zero keeps this a
                        // warm gradient with no near-black rings anywhere.
                        double t = 0.5 + 0.5 * std::sin(phi_out * 1.5 + l_final * 4.0);
                        double heat = 0.7 + 0.3 * t;  // stays in [0.7, 1.0], never dips dark
                        return Vec<double, 3>{std::clamp(255.0 * heat, 0.0, 255.0),
                                               std::clamp(150.0 + 95.0 * t, 0.0, 255.0),
                                               std::clamp(55.0 + 55.0 * t, 0.0, 255.0)};
                    }

                    double r_final = r_of_l(l_final);
                    double v_phi = b / r_final;
                    double v_l = (l_final > 0.0 ? 1.0 : -1.0) *
                                 std::sqrt(std::max(0.0, 1.0 - v_phi * v_phi));

                    // Rotate the transverse basis by the total angle swept
                    // (Rodrigues' formula around l_hat) -- t_hat2 completes
                    // the local orthonormal triad.
                    Vec<double, 3> t_hat2{l_hat.cross(t_hat)};
                    Vec<double, 3> t_rot{t_hat * std::cos(phi_out) + t_hat2 * std::sin(phi_out)};
                    Vec<double, 3> outdir{l_hat * v_l + t_rot * v_phi};

                    return sample_sky_color(l_final > 0.0 ? sky_a : sky_b, outdir);
                };

                spatium::render::supersample_pixel(x, y, width, height, tan_fov, aspect, ray_color,
                                                    px);
            }
        }
    };

    unsigned nthreads = std::max(1u, std::thread::hardware_concurrency());
    nthreads = std::min(nthreads, static_cast<unsigned>(height));
    int rows_per_thread = (height + static_cast<int>(nthreads) - 1) / static_cast<int>(nthreads);
    {
        std::vector<std::jthread> workers;
        workers.reserve(nthreads);
        for (unsigned t = 0; t < nthreads; ++t) {
            int y0 = static_cast<int>(t) * rows_per_thread;
            int y1 = std::min(height, y0 + rows_per_thread);
            if (y0 >= y1) continue;
            workers.emplace_back([&render_rows, y0, y1] { render_rows(y0, y1); });
        }
    }  // jthreads join here
}

}  // namespace

int main(int argc, char** argv) {
    bool force = false;
    int n_frames = N_FRAMES;
    std::string mode = "flythrough";
    for (int i = 1; i < argc; ++i) {
        std::string_view a = argv[i];
        if (a == "--help" || a == "-h") { print_usage(); return 0; }
        if (a == "--force") { force = true; continue; }
        if (a == "--frames" && i + 1 < argc) { n_frames = std::atoi(argv[++i]); continue; }
        if (a == "--mode" && i + 1 < argc) { mode = argv[++i]; continue; }
        std::print(stderr, "unknown option: {}\n", a);
        return 1;
    }

    // Two visually distinct starfields -- crossing the throat should be
    // unmistakable. Side A = l > 0 (deep blue, matches blackhole_demo's
    // own palette), side B = l < 0 (warm amber).
    Sky sky_a = make_starfield(12000, /*seed=*/42, {8.0, 6.0, 20.0});
    Sky sky_b = make_starfield(9000, /*seed=*/7, {22.0, 12.0, 6.0});
    const Vec<double, 3> l_hat{1.0, 0.0, 0.0};  // wormhole's fixed symmetry axis

    if (mode == "ring") {
        // The classic external shot: camera far from the throat (deep
        // in the asymptotically-flat region, where the on-axis camera
        // machinery above is exact, not an approximation -- no separate
        // off-axis derivation needed, see the note on DISK_HALF_WIDTH's
        // TraceResult for why a genuinely general off-axis camera isn't
        // a small extension here: this wormhole's spatial topology is
        // R x S^2 -- an l-axis crossed with a full 2-sphere at every l,
        // not an ordinary spherically-symmetric-around-a-point space --
        // so "move the camera to the side in flat 3D" doesn't have a
        // single unambiguous meaning the way it does for a black hole).
        // A large fixed tilt breaks the perfect symmetry a dead-on shot
        // would have, framing the throat as a distinct ring off-center
        // with the lensed starfield wrapped around it, rather than
        // through it.
        constexpr double kRingLCam = 20.0;
        constexpr double kRingTilt = 0.42;
        constexpr int kRingW = 1600, kRingH = 900;
        std::vector<std::uint8_t> img;
        auto t0 = std::chrono::steady_clock::now();
        render_wormhole_frame(kRingLCam, kRingTilt, /*precess=*/0.0, l_hat, sky_a, sky_b, kRingW,
                               kRingH, FOV_HALF, img);
        double ms = std::chrono::duration<double, std::milli>(
                        std::chrono::steady_clock::now() - t0)
                        .count();
        std::print("wormhole_demo: einstein_ring.png {}x{}, {:.0f} ms\n", kRingW, kRingH, ms);
        if (spatium::examples::confirm_overwrite("einstein_ring.png", force))
            stbi_write_png("einstein_ring.png", kRingW, kRingH, 3, img.data(), kRingW * 3);
        return 0;
    }

    std::error_code ec;
    std::filesystem::create_directories("wormhole_frames", ec);

    auto t0 = std::chrono::steady_clock::now();

    for (int frame = 0; frame < n_frames; ++frame) {
        double t01 = n_frames > 1 ? double(frame) / (n_frames - 1) : 0.0;
        double l_cam = L_START * (1.0 - 2.0 * t01);  // +L_START -> -L_START through the throat

        // Tilted, animated camera: looking straight down -l_hat the
        // whole way keeps every frame perfectly axially symmetric, which
        // is exactly why the mid-flythrough frames looked flat/
        // featureless -- tilting away from the axis by an angle that
        // ramps up approaching the throat and precesses in azimuth over
        // the flythrough makes the lensing visibly asymmetric (some
        // pixels closer to the critical impact parameter than others)
        // and makes background stars visibly slide across the frame as
        // the tilt direction rotates, per the user's own request.
        constexpr double kTiltMax = 0.35;      // radians, peak deviation from the axis
        constexpr double kPrecessCycles = 1.5; // full azimuthal turns over the flythrough
        double tilt = kTiltMax * std::sin(std::numbers::pi * t01);
        double precess = kPrecessCycles * 2.0 * std::numbers::pi * t01;

        std::vector<std::uint8_t> img;
        render_wormhole_frame(l_cam, tilt, precess, l_hat, sky_a, sky_b, W, H, FOV_HALF, img);

        char path[64];
        std::snprintf(path, sizeof(path), "wormhole_frames/frame_%04d.png", frame);
        if (spatium::examples::confirm_overwrite(path, force))
            stbi_write_png(path, W, H, 3, img.data(), W * 3);

        std::print("\r  frame {}/{}  l_cam={:+.2f}", frame + 1, n_frames, l_cam);
    }

    auto t1 = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::print("\rwormhole_demo: {} frames at {}x{}, {:.0f} ms ({:.0f} ms/frame)\n", n_frames, W, H,
               ms, ms / n_frames);
    std::print("Assemble with:\n"
               "  ffmpeg -framerate 30 -i wormhole_frames/frame_%04d.png "
               "-pix_fmt yuv420p wormhole.mp4\n");
    return 0;
}
