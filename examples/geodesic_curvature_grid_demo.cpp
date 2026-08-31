// A 3x3 grid of parallel light rays crossing near a Schwarzschild black
// hole, traced as real geodesics (the same geodesic.hpp engine
// blackhole_gr_demo.cpp uses -- no approximation), rendered as bent
// curves in 3D rather than composited into a shaded scene. This is the
// classic "light bending near a massive object" diagram, done with an
// actual numerical geodesic integrator instead of a hand-drawn
// illustration: initially parallel rays visibly converge, some spiral
// in and are captured, others deflect and escape -- spacetime
// curvature made visible through the ray geometry itself.
//
// Each ray starts far from the hole moving in +x, at a grid of (y,z)
// offsets (the impact-parameter plane). Schwarzschild's full spherical
// symmetry means each ray, taken alone, stays confined to its own
// plane through the origin -- but with rays at many different (y,z)
// offsets this file doesn't exploit that per-ray-plane trick (unlike
// blackhole_gr_demo.cpp's camera path): it integrates the full
// (t,r,theta,phi) state directly from each ray's own actual spherical
// coordinates, since the point here is recording and drawing the real
// 3D path, not just an exit direction.

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "io_helpers.hpp"

#include <spatium/algebra/vector.hpp>
#include <spatium/physics/relativity/geodesic.hpp>
#include <spatium/physics/relativity/schwarzschild.hpp>
#include <spatium/render/camera.hpp>
#include <spatium/render/write_image.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <numbers>
#include <string_view>
#include <vector>

using namespace spatium;
using namespace spatium::physics::relativity;
using namespace spatium::render;
using namespace spatium::examples;

namespace {

// Standard spherical-polar orthonormal tetrad in the ambient Cartesian
// embedding -- Schwarzschild's a=0 special case of the oblate-
// spheroidal tetrad blackhole_gr_demo.cpp's Kerr path uses (rho=r when
// a=0), rewritten here directly since a=0 lets it simplify to the
// textbook spherical tetrad.
struct Tetrad { Vec<double, 3> e_r, e_theta, e_phi; };

Tetrad spherical_tetrad(double theta, double phi) {
    double s = std::sin(theta), c = std::cos(theta);
    double cp = std::cos(phi), sp = std::sin(phi);
    Vec<double, 3> e_r{s * cp, s * sp, c};
    Vec<double, 3> e_th{c * cp, c * sp, -s};
    Vec<double, 3> e_ph{-sp, cp, 0.0};
    return {e_r, e_th, e_ph};
}

struct RayPath {
    std::vector<Vec<double, 3>> points;
    bool captured;
};

// Builds the initial 4-velocity for a photon at Cartesian `start` with
// LOCAL flat-space direction `dir` (unit vector), via the same static-
// observer tetrad decomposition blackhole_gr_demo.cpp's Schwarzschild
// fast path uses -- exact at any r outside the horizon, not just a
// large-r approximation.
RayPath trace_ray(const SchwarzschildMetric<double>& metric, double mass,
                   const Vec<double, 3>& start, const Vec<double, 3>& dir, double r_capture,
                   double r_escape, long max_steps) {
    double r0 = start.norm();
    double theta0 = std::acos(std::clamp(start[2] / r0, -1.0, 1.0));
    double phi0 = std::atan2(start[1], start[0]);
    Tetrad tet = spherical_tetrad(theta0, phi0);

    double n_r = dir.dot(tet.e_r), n_th = dir.dot(tet.e_theta), n_ph = dir.dot(tet.e_phi);
    double f0 = 1.0 - 2.0 * mass / r0;
    double ut = 1.0 / std::sqrt(f0);
    double ur = n_r * std::sqrt(f0);
    double utheta = n_th / r0;
    double uphi = n_ph / (r0 * std::sin(theta0));

    Vec<double, 8> state{0.0, r0, theta0, phi0, ut, ur, utheta, uphi};
    RayPath path;
    path.points.push_back(start);
    path.captured = false;

    for (long step = 0; step < max_steps; ++step) {
        double r = state[1];
        double dl = 0.1 * std::max(0.08, r / 20.0);
        state = geodesic_step(metric, state, dl);
        double r_new = state[1];
        if (!std::isfinite(r_new) || r_new <= r_capture) {
            path.captured = true;
            break;
        }
        double theta = state[2], phi = state[3];
        double s = std::sin(theta);
        path.points.push_back(
            Vec<double, 3>{r_new * s * std::cos(phi), r_new * s * std::sin(phi), r_new * std::cos(theta)});
        if (r_new >= r_escape && state[5] > 0.0) break;
    }
    return path;
}

// Draws a soft anti-aliased line segment into an HDR RGB buffer via a
// simple 1-pixel-radius splat along the segment -- more than adequate
// for thin bright curves on a dark background, no need for a full
// scanline rasterizer here.
void draw_line(std::vector<double>& img, int w, int h, double x0, double y0, double x1, double y1,
               Vec<double, 3> color) {
    double dx = x1 - x0, dy = y1 - y0;
    double len = std::sqrt(dx * dx + dy * dy);
    int n = std::max(1, static_cast<int>(len * 2.0));
    for (int i = 0; i <= n; ++i) {
        double t = static_cast<double>(i) / n;
        double x = x0 + dx * t, y = y0 + dy * t;
        int xi = static_cast<int>(std::floor(x)), yi = static_cast<int>(std::floor(y));
        for (int oy = 0; oy <= 1; ++oy)
            for (int ox = 0; ox <= 1; ++ox) {
                int px = xi + ox, py = yi + oy;
                if (px < 0 || px >= w || py < 0 || py >= h) continue;
                double wx = 1.0 - std::abs(x - (xi + ox));
                double wy = 1.0 - std::abs(y - (yi + oy));
                double weight = wx * wy;
                std::size_t idx = 3 * (static_cast<std::size_t>(py) * w + px);
                img[idx] += color[0] * weight;
                img[idx + 1] += color[1] * weight;
                img[idx + 2] += color[2] * weight;
            }
    }
}

// Projects a 3D world point to 2D screen pixel coordinates under
// `cam`/`basis`; returns false if the point is behind the camera.
bool project(const Camera<double>& cam, const CameraBasis<double>& basis, int w, int h,
             const Vec<double, 3>& point, double& px, double& py) {
    Vec<double, 3> local{point - cam.position};
    double depth = local.dot(basis.fwd);
    if (depth <= 1e-6) return false;
    double aspect = static_cast<double>(w) / h;
    double sx = local.dot(basis.right) / (depth * basis.tan_half * aspect);
    double sy = local.dot(basis.up) / (depth * basis.tan_half);
    px = (sx + 1.0) * 0.5 * w;
    py = (1.0 - sy) * 0.5 * h;
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    int W = argc > 1 ? std::atoi(argv[1]) : 1280;
    int H = argc > 2 ? std::atoi(argv[2]) : 720;
    bool force = argc > 3 && std::string_view(argv[3]) == "--force";

    if (!confirm_overwrite("geodesic_curvature_grid.png", force)) return 0;

    double M = 1.0;
    SchwarzschildMetric<double> metric{M};
    double r_horizon = schwarzschild_horizon_radius(M);
    double r_capture = r_horizon * 1.02;
    double r_escape = 60.0;
    double b_crit = schwarzschild_critical_impact_parameter(M);

    // Rays start at x=-30M moving in +x, on a grid of (y,z) impact-
    // parameter offsets spanning well inside the critical impact
    // parameter (captured) to well outside it (barely deflected) --
    // the full spectrum the classic diagram needs to make its point.
    constexpr int N = 15;
    constexpr double kSpan = 2.2;  // grid half-extent, in units of b_crit
    std::vector<RayPath> paths;
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            double y = (-kSpan + 2.0 * kSpan * i / (N - 1)) * b_crit;
            double z = (-kSpan + 2.0 * kSpan * j / (N - 1)) * b_crit;
            Vec<double, 3> start{-30.0, y, z};
            Vec<double, 3> dir{1.0, 0.0, 0.0};
            paths.push_back(trace_ray(metric, M, start, dir, r_capture, r_escape, 4000));
        }
    }

    long n_captured = std::count_if(paths.begin(), paths.end(), [](const RayPath& p) { return p.captured; });
    std::print("{} rays: {} captured, {} escaped\n", paths.size(), n_captured,
               static_cast<long>(paths.size()) - n_captured);

    // Oblique viewing camera -- off to the side, slightly above the
    // ray bundle's own plane of travel, so the 3D bending (rays that
    // started on a flat grid ending up on a curved, non-planar surface
    // after passing the hole) is visible rather than looking like a
    // flat 2D diagram.
    Camera<double> cam{.position = {8.0, -38.0, 16.0}, .target = {0.0, 0.0, 0.0},
                        .up = {0.0, 0.0, 1.0}, .fov_deg = 42.0};
    CameraBasis<double> basis = make_camera_basis(cam);

    std::vector<double> img(3 * static_cast<std::size_t>(W) * H, 0.0);

    // Horizon silhouette: approximate as a flat dark disk facing the
    // camera at the origin, radius = horizon radius -- a diagram aid,
    // not a lensed render (the whole point here is the ray paths).
    {
        Vec<double, 3> to_cam{(cam.position).normalized()};
        Vec<double, 3> up_ish = std::abs(to_cam[2]) < 0.9 ? Vec<double, 3>{0, 0, 1} : Vec<double, 3>{1, 0, 0};
        Vec<double, 3> u{up_ish.cross(to_cam).normalized()};
        Vec<double, 3> v{to_cam.cross(u).normalized()};
        constexpr int kDiskSegs = 96;
        double px_prev = 0, py_prev = 0;
        bool have_prev = false;
        for (int i = 0; i <= kDiskSegs; ++i) {
            double ang = 2.0 * std::numbers::pi * i / kDiskSegs;
            Vec<double, 3> pt{Vec<double, 3>{u * (r_horizon * std::cos(ang)) +
                                              v * (r_horizon * std::sin(ang))}};
            double px, py;
            if (project(cam, basis, W, H, pt, px, py)) {
                if (have_prev) draw_line(img, W, H, px_prev, py_prev, px, py, {40.0, 30.0, 20.0});
                px_prev = px; py_prev = py; have_prev = true;
            } else {
                have_prev = false;
            }
        }
    }

    for (const auto& path : paths) {
        Vec<double, 3> color = path.captured ? Vec<double, 3>{180.0, 60.0, 40.0}
                                              : Vec<double, 3>{120.0, 200.0, 255.0};
        double px_prev = 0, py_prev = 0;
        bool have_prev = false;
        for (const auto& pt : path.points) {
            double px, py;
            if (project(cam, basis, W, H, pt, px, py)) {
                if (have_prev) draw_line(img, W, H, px_prev, py_prev, px, py, color);
                px_prev = px; py_prev = py; have_prev = true;
            } else {
                have_prev = false;
            }
        }
    }

    std::vector<std::uint8_t> img8(img.size());
    for (std::size_t i = 0; i < img.size(); ++i)
        img8[i] = static_cast<std::uint8_t>(std::clamp(img[i], 0.0, 255.0));
    write_png_rgb("geodesic_curvature_grid.png", W, H, img8);
    std::print("wrote geodesic_curvature_grid.png ({}x{})\n", W, H);
    return 0;
}
