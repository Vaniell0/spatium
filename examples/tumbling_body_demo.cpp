// Dzhanibekov effect (intermediate-axis theorem) — torque-free rigid
// body tumble, static camera.
//
// A free rigid body spun almost exactly about its INTERMEDIATE moment-
// of-inertia axis is a textbook instability (Euler's equations have a
// saddle point there): a vanishingly small perturbation grows until the
// body flips its spin axis, settles briefly, then flips back —
// indefinitely, not just once. This is the "twisting wrench in orbit"
// video (Dzhanibekov 1985) and the classroom "tennis racket theorem."
// It's exactly solvable (Jacobi elliptic functions) for the ideal free
// body, so "unpredictable-looking" here means genuinely striking motion
// from deterministic physics, not chaos in the sensitive-to-initial-
// conditions sense.
//
// Physics: `physics/mechanics/lgvi.hpp`'s `lgvi_rigid_body_step()` —
// already built and tested this project, unmodified here. It's a
// Lie-group variational integrator on SO(3): orientation stays exactly
// on the group (no renormalization drift, ever) and spatial angular
// momentum is conserved to round-off by discrete Noether — verified
// directly with a standalone probe before committing to this demo's
// parameters: over a t=0..200 run at h=1e-3, |Π| stayed at 1.00010 the
// entire time and 4 clean flips appeared roughly every 50-52 time units,
// with the state at t=100 landing almost exactly back on the t=0 start —
// one full flip-and-recover cycle, the natural loop this demo renders
// (see the frame-count constants below for the exact numbers, fixed
// once after the first cut turned out to spin far too fast to actually
// watch). A naive Euler/RK4 integrator run this long would visibly drift
// and eventually make the flip look wrong; this is the actual point of
// using LGVI here, not an arbitrary implementation choice.
//
// Geometry: the tumbling body is a solid ellipsoid with three distinct
// semi-axes (a < b < c) — for a uniform ellipsoid the principal moments
// are I_x=(b²+c²)/5, I_y=(a²+c²)/5, I_z=(a²+b²)/5 (mass=1), so I_x is
// largest (excludes the shortest semi-axis) and I_z is smallest
// (excludes the longest) — the INTERMEDIATE moment I_y always belongs to
// the intermediate-LENGTH semi-axis b, for any ellipsoid. J_diag below
// is computed from the same (a,b,c) used to render the body, not chosen
// independently — the physics and the geometry are the same object.
//
// Rendering: the ellipsoid stays fixed and axis-aligned in its own body
// frame; each frame's camera ray is transformed into body space by R^T
// (R = current orientation, an SO(3) matrix — R^T = R^-1 for a rotation)
// before calling `ray_quadric()` — a closed-form, always-fast
// intersection (`Quadric::ellipsoid()`, ~25ns/ray, see README's raycast
// benchmark table) rather than Newton-UV ray marching, so a rotating
// object costs nothing extra per frame versus a static one. The hit
// normal is rotated back to world space (R · n) for shading. Grid-line
// UV coordinates are derived from the body-space hit point after the
// fact (u = atan2(y/b, x/a), v = acos(z/c)) purely for surface markings
// — decoupled from the intersection method, unlike the parametric-
// surface demos where UV comes from the Newton solve itself.
//
// Built on the render/ engine primitives added the same session this
// demo showcases: `render::Camera`/`make_camera_basis`/`camera_pixel_dir`
// (camera.hpp), `render::parallel_for_rows` (work-stealing row
// parallelism), `render::supersample_pixel` (2x2 antialiasing), and
// `render::write_png_rgb` — no per-file reimplementation of any of them,
// which is the entire point of this demo: it's meant to be the pattern
// future demos (and the existing ones, revisited later) follow.
//
// Output: tumbling_frames/frame_%04d.png — assemble separately:
//   ffmpeg -framerate 60 -i tumbling_frames/frame_%04d.png -pix_fmt yuv420p tumbling_body.mp4

// Only `#define` here -- `render/write_image.hpp` below is the one
// place that includes <spatium/vendor/stb_image_write.h> in this file.
// Including it a second time directly (the old per-demo pattern) while
// the macro is still defined recompiles the implementation twice and
// fails with redefinition errors.
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "io_helpers.hpp"

#include <spatium/algebra/groups/so3.hpp>
#include <spatium/algebra/vector.hpp>
#include <spatium/geometry/ray_surface.hpp>
#include <spatium/physics/mechanics/lgvi.hpp>
#include <spatium/render/camera.hpp>
#include <spatium/render/color.hpp>
#include <spatium/render/parallel_for_rows.hpp>
#include <spatium/render/sky.hpp>
#include <spatium/render/supersample.hpp>
#include <spatium/render/write_image.hpp>

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
#include <vector>

using spatium::Vec;
using spatium::geometry::Quadric;
using spatium::geometry::ray_quadric;
using spatium::geometry::Ray;
using spatium::physics::mechanics::LGVIRigidBodyState;
using spatium::physics::mechanics::lgvi_rigid_body_step;
using SO3 = spatium::algebra::SO3;
using spatium::render::hsv_to_rgb255;
using spatium::render::make_starfield;
using spatium::render::sample_sky_color;
using spatium::render::Sky;
using spatium::render::Camera;
using spatium::render::make_camera_basis;
using spatium::render::camera_ray_dir;
using spatium::render::parallel_for_rows;
using spatium::render::supersample_pixel;
using spatium::render::write_png_rgb;

namespace {

constexpr int    W = 960;
constexpr int    H = 540;

// Ellipsoid semi-axes, a < b < c -- see file header for why the
// intermediate-inertia axis this depends on is always the intermediate-
// length one.
constexpr double A_AXIS = 0.5;
constexpr double B_AXIS = 0.85;
constexpr double C_AXIS = 1.6;

// Verified via a standalone probe (2026-08-26, not guessed): h=1e-3 with
// this perturbation gives 4 clean flips over ~200 time units, momentum
// magnitude constant at 1.00010 throughout, and the state at t=100
// (Pi=(+0.0066,+1.0000,+0.0088)) lands almost exactly back on the start
// (Pi=(0.01,1.0,0.01)) -- one full flip-and-recover cycle, a natural
// loop point.
//
// Frame pacing, fixed after the first render looked far too fast to
// actually see the tumble (a real, reported problem, not a style
// preference): base spin rate omega = Pi_y/J_y ~ 1/0.562 ~ 1.78 rad per
// time unit. The first cut used 333 substeps/frame = 0.333 time units/
// frame ~ 34 degrees/frame at 60fps -- far too much angular travel per
// frame to read as smooth rotation. 40 substeps/frame = 0.04 time units/
// frame ~ 4 degrees/frame is comfortably watchable, and 2500 frames
// covers exactly the t=0..100 loop above.
constexpr double DT_PHYSICS       = 1e-3;
constexpr int    SUBSTEPS_PER_FRAME = 40;    // ~4 deg/frame at 60fps
constexpr int    N_FRAMES         = 2500;    // covers one full t=0..100 cycle
constexpr double PI_MAGNITUDE     = 1.0;
constexpr double PI_PERTURBATION  = 0.01;    // tiny push off the unstable axis

constexpr double FOV_DEG   = 32.0;
constexpr double GRID_SPACING_U = std::numbers::pi / 6.0;   // 12 lines around
constexpr double GRID_SPACING_V = std::numbers::pi / 6.0;
constexpr double GRID_HALF_WIDTH = 0.02;
constexpr double GRID_DARKEN     = 0.15;

Vec<double, 3> ellipsoid_j_diag() {
    double a2 = A_AXIS * A_AXIS, b2 = B_AXIS * B_AXIS, c2 = C_AXIS * C_AXIS;
    return {(b2 + c2) / 5.0, (a2 + c2) / 5.0, (a2 + b2) / 5.0};
}

bool near_periodic_line(double coord, double spacing, double half_width) {
    double m = std::fmod(coord, spacing);
    if (m < 0.0) m += spacing;
    return m < half_width || m > spacing - half_width;
}

void print_usage() {
    std::print("Usage: tumbling_body_demo [--frames N] [--force] [--help]\n"
               "  Dzhanibekov-effect rigid-body tumble via LGVI (physics/mechanics/lgvi.hpp).\n"
               "  --frames N   frame count (default 600)\n"
               "  --force      overwrite existing output files\n"
               "  --help       show this message\n"
               "  Output:      tumbling_frames/frame_%04d.png (960x540 RGB)\n");
}

}  // namespace

int main(int argc, char** argv) {
    bool force = false;
    int n_frames = N_FRAMES;
    for (int i = 1; i < argc; ++i) {
        std::string_view a = argv[i];
        if (a == "--help" || a == "-h") { print_usage(); return 0; }
        if (a == "--force") { force = true; continue; }
        if (a == "--frames" && i + 1 < argc) { n_frames = std::atoi(argv[++i]); continue; }
        std::print(stderr, "unknown option: {}\n", a);
        return 1;
    }

    const Vec<double, 3> j_diag = ellipsoid_j_diag();
    const Quadric<double> ellipsoid = Quadric<double>::ellipsoid(A_AXIS, B_AXIS, C_AXIS);
    Sky sky = make_starfield(9000, /*seed=*/3, /*tint=*/{6.0, 6.0, 10.0});

    LGVIRigidBodyState body;
    body.R = SO3::ElementType::identity();
    body.Pi = SO3::AlgebraType{PI_PERTURBATION, PI_MAGNITUDE, PI_PERTURBATION};

    const Camera<double> cam{
        .position = {3.6, 2.4, 1.8}, .target = {0.0, 0.0, 0.0}, .up = {0.0, 0.0, 1.0},
        .fov_deg = FOV_DEG};
    const auto basis = make_camera_basis(cam);
    const Vec<double, 3> light = Vec<double, 3>{Vec<double, 3>{1.0, -0.4, 1.2}.normalized()};

    std::error_code ec;
    std::filesystem::create_directories("tumbling_frames", ec);

    auto t0 = std::chrono::steady_clock::now();

    for (int frame = 0; frame < n_frames; ++frame) {
        for (int s = 0; s < SUBSTEPS_PER_FRAME; ++s)
            body = lgvi_rigid_body_step(body, j_diag, DT_PHYSICS);

        // R^T rotates world-space vectors into the ellipsoid's own fixed
        // body frame -- the object itself never moves, the camera ray
        // does, which is equivalent and avoids rebuilding anything
        // per-frame beyond this one 3x3 transpose.
        SO3::ElementType Rt = body.R.transpose();

        std::vector<std::uint8_t> img(3 * static_cast<std::size_t>(W) * H, 0);
        Vec<double, 3> origin_body = Rt * cam.position;

        parallel_for_rows(H, [&](int y) {
            for (int x = 0; x < W; ++x) {
                std::uint8_t* px = &img[3 * (static_cast<std::size_t>(y) * W + x)];

                auto ray_color = [&](double sx, double sy) -> Vec<double, 3> {
                    // sx, sy arrive already in the final tan(fov/2)-scaled
                    // NDC convention (supersample_pixel computed them the
                    // same way camera_pixel_dir would for this pixel), so
                    // camera_ray_dir consumes them directly.
                    Vec<double, 3> dir_world = camera_ray_dir(basis, sx, sy);
                    Vec<double, 3> dir_body = Rt * dir_world;

                    auto hits = ray_quadric(Ray<3, double>{origin_body, dir_body}, ellipsoid);
                    if (hits.empty()) return sample_sky_color(sky, dir_world);

                    const auto& hit = hits.front();
                    Vec<double, 3> n_body = hit.normal;
                    if (n_body.dot(dir_body) > 0.0) n_body = Vec<double, 3>{-n_body};
                    Vec<double, 3> n_world = Vec<double, 3>{body.R * n_body};

                    double diff = std::max(0.0, n_world.dot(light));
                    double shaded = 0.25 + 0.75 * diff;

                    double u = std::atan2(hit.point[1] / B_AXIS, hit.point[0] / A_AXIS);
                    double v = std::acos(std::clamp(hit.point[2] / C_AXIS, -1.0, 1.0));
                    double hue = (u + std::numbers::pi) / (2.0 * std::numbers::pi);
                    Vec<double, 3> base = hsv_to_rgb255(hue, 0.55, 0.95);
                    if (near_periodic_line(u + std::numbers::pi, GRID_SPACING_U, GRID_HALF_WIDTH) ||
                        near_periodic_line(v, GRID_SPACING_V, GRID_HALF_WIDTH))
                        shaded *= GRID_DARKEN;

                    return Vec<double, 3>{base * shaded};
                };

                supersample_pixel(x, y, W, H, basis.tan_half,
                                   static_cast<double>(W) / H, ray_color, px);
            }
        });

        char path[64];
        std::snprintf(path, sizeof(path), "tumbling_frames/frame_%04d.png", frame);
        if (spatium::examples::confirm_overwrite(path, force)) write_png_rgb(path, W, H, img);

        std::print("\r  frame {}/{}", frame + 1, n_frames);
    }

    auto t1 = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::print("\rtumbling_body_demo: {} frames at {}x{}, {:.0f} ms ({:.0f} ms/frame)\n", n_frames,
               W, H, ms, ms / n_frames);
    std::print("Assemble with:\n"
               "  ffmpeg -framerate 60 -i tumbling_frames/frame_%04d.png "
               "-pix_fmt yuv420p tumbling_body.mp4\n");
    return 0;
}
