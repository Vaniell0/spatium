// Schwarzschild black hole ray-tracer.
//
// Backward null-geodesic tracing in the equatorial plane of each pixel ray.
// Schwarzschild photon equation in u = 1/r form:
//
//     d²u/dφ² + u = 3 M u²              (1)
//     (du/dφ)²   = 1/b² - u² + 2 M u³   (first integral, sets initial u')
//
// where b is the impact parameter and M is the geometric-units mass.
// Each pixel ray defines a 2-plane spanned by (camera position, ray direction);
// φ is measured inside that plane from the camera's outward radial.
// We integrate (1) with velocity-Verlet until either u reaches the horizon
// (ray captured → black) or u falls below a far-field cutoff (ray escapes →
// sample starfield at the outgoing direction).
//
// Output: blackhole.png (1920x1080 RGB).

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
#include <filesystem>
#include <numbers>
#include <print>
#include <thread>
#include <random>
#include <string_view>
#include <vector>

using spatium::Vec;

namespace {

constexpr int    W = 1920;
constexpr int    H = 1080;
constexpr double M_BH       = 1.0;        // Schwarzschild mass (geometric units)
constexpr double R_HORIZON  = 2.0 * M_BH; // event horizon
constexpr double D_CAM      = 30.0 * M_BH;// camera distance from BH
constexpr double FOV_HALF   = 0.45;       // half-angle in radians (~26°)
constexpr int    MAX_STEPS  = 8000;
constexpr double DPHI       = 0.005;      // φ-step
constexpr double U_CAPTURE  = 1.0 / (R_HORIZON * 1.001); // just outside horizon
constexpr double U_ESCAPE   = 1.0 / (300.0 * M_BH);      // far-field cutoff

using spatium::render::make_starfield;
using spatium::render::sample_sky;
using spatium::render::Sky;

// Returns true if the ray escapes; sets phi_out to the total angular sweep.
// Returns false if u ≥ horizon (captured).
bool trace_geodesic(double u0, double up0, double& phi_out) {
    double u = u0;
    double up = up0;
    double phi = 0.0;
    auto upp_of = [](double uu) { return -uu + 3.0 * M_BH * uu * uu; };
    double upp = upp_of(u);
    for (int i = 0; i < MAX_STEPS; ++i) {
        // Velocity-Verlet on u'' = -u + 3 M u²
        double u_new   = u + DPHI * up + 0.5 * DPHI * DPHI * upp;
        double upp_new = upp_of(u_new);
        double up_new  = up + 0.5 * DPHI * (upp + upp_new);
        u = u_new; up = up_new; upp = upp_new;
        phi += DPHI;
        if (u >= U_CAPTURE) return false;
        if (u <= U_ESCAPE && up < 0.0) { phi_out = phi; return true; }
    }
    phi_out = phi;
    return true;
}

void print_usage() {
    std::print("Usage: blackhole_demo [--force] [--help]\n"
               "  Schwarzschild null-geodesic ray-tracer.\n"
               "  --force   overwrite existing output file\n"
               "  --help    show this message\n"
               "  Output:   blackhole.png (1920x1080 RGB)\n");
}

}  // namespace

int main(int argc, char** argv) {
    bool force = false;
    for (int i = 1; i < argc; ++i) {
        std::string_view a = argv[i];
        if (a == "--help" || a == "-h") { print_usage(); return 0; }
        if (a == "--force") { force = true; continue; }
        std::print(stderr, "unknown option: {}\n", a);
        return 1;
    }

    auto t0 = std::chrono::steady_clock::now();
    Sky sky = make_starfield(12000);

    std::vector<std::uint8_t> img(3 * W * H, 0);

    // Camera at (D, 0, 0), looking toward origin (BH).
    const Vec<double, 3> cam_pos  {D_CAM, 0.0, 0.0};
    const Vec<double, 3> cam_fwd  {-1.0,  0.0, 0.0};
    const Vec<double, 3> cam_right{ 0.0,  0.0, 1.0};  // image x → world +z
    const Vec<double, 3> cam_up   { 0.0,  1.0, 0.0};  // image y → world +y
    const double aspect = double(W) / H;
    const double tan_fov = std::tan(FOV_HALF);

    auto render_rows = [&](int y0, int y1) -> int {
        int local_captured = 0;
        for (int y = y0; y < y1; ++y) {
            for (int x = 0; x < W; ++x) {
                std::uint8_t* px = &img[3 * (y * W + x)];
                bool any_captured = false;

                // Per-ray color, called aa*aa times per pixel by
                // supersample_pixel (spatium/render/supersample.hpp) at
                // sub-pixel (sx, sy) -- averaging smooths the horizon
                // silhouette and photon-ring edge automatically instead
                // of leaving them hard-aliased.
                auto ray_color = [&](double sx, double sy) -> Vec<double, 3> {
                    Vec<double, 3> dir = (cam_fwd + cam_right * sx + cam_up * sy).normalized();

                    // Geodesic plane: spanned by r_hat (outward radial) and tangential
                    // component of dir. φ measured from r_hat toward t_hat.
                    Vec<double, 3> r_hat = (cam_pos / cam_pos.norm());
                    double dir_dot_r = dir.dot(r_hat);
                    Vec<double, 3> t_hat = dir - r_hat * dir_dot_r;
                    double t_norm = t_hat.norm();

                    if (t_norm < 1e-9) {
                        // Pure-radial ray: hits BH if pointing inward.
                        if (dir_dot_r < 0) { any_captured = true; return {0.0, 0.0, 0.0}; }
                        return spatium::render::sample_sky_color(sky, dir);
                    }
                    t_hat = t_hat / t_norm;

                    // Initial conditions for u(φ):
                    //   u₀  = 1/D
                    //   sin α = |t-component of dir| = t_norm
                    //   b    = D * sin α  (asymptotic; exact at large D/M)
                    //   (u')² = 1/b² - u² + 2Mu³, sign chosen by radial direction.
                    const double D = cam_pos.norm();
                    const double u0 = 1.0 / D;
                    const double sin_a = t_norm;
                    const double b = D * sin_a;
                    double rhs = 1.0 / (b * b) - u0 * u0 + 2.0 * M_BH * u0 * u0 * u0;
                    if (rhs < 0) rhs = 0;
                    double up0 = std::sqrt(rhs);
                    // dir_dot_r > 0  ⇒ pointing outward ⇒ r grows ⇒ u shrinks ⇒ u' < 0.
                    // dir_dot_r < 0  ⇒ pointing inward  ⇒ u grows ⇒ u' > 0.
                    if (dir_dot_r > 0) up0 = -up0;

                    double phi_out = 0.0;
                    if (!trace_geodesic(u0, up0, phi_out)) {
                        any_captured = true;
                        return {0.0, 0.0, 0.0};
                    }
                    // Escape direction: the ray's outgoing tangent in world coords.
                    // The outward radial after sweep φ_out is r̂(φ) = cos(φ)·r̂ + sin(φ)·t̂.
                    // For a near-asymptotic geodesic the outgoing direction aligns with
                    // that radial — accurate enough for starfield sampling.
                    double cph = std::cos(phi_out), sph = std::sin(phi_out);
                    Vec<double, 3> outdir = r_hat * cph + t_hat * sph;
                    Vec<double, 3> color = spatium::render::sample_sky_color(sky, outdir);

                    // Photon-ring glow: a ray that swept phi_out well past
                    // half a turn grazed close to the photon sphere
                    // (b ~ 3*sqrt(3)*M) before escaping -- this code already
                    // computes that, but until now only used it for the
                    // escape *direction*, discarding how close the ray came
                    // to capture. The GR-native analog of
                    // ray_quadric_proximity()'s complex-root glow
                    // (parametric_analytical_demo.cpp's glow_sphere.png):
                    // same idea (render a "how close was this near-miss"
                    // quantity instead of throwing it away), different
                    // physics providing the number.
                    if (phi_out > std::numbers::pi) {
                        double glow = std::min(1.0, (phi_out - std::numbers::pi) / (2.0 * std::numbers::pi));
                        color = Vec<double, 3>{std::min(255.0, color[0] + glow * 220.0),
                                                std::min(255.0, color[1] + glow * 160.0),
                                                std::min(255.0, color[2] + glow * 60.0)};
                    }
                    return color;
                };

                spatium::render::supersample_pixel(x, y, W, H, tan_fov, aspect, ray_color, px);
                if (any_captured) ++local_captured;
            }
        }
        return local_captured;
    };

    // Rows split across std::jthread workers -- each pixel is fully
    // independent (a fresh geodesic integration into a disjoint slice of
    // img), so no synchronization is needed beyond the implicit join
    // when the workers go out of scope.
    unsigned nthreads = std::max(1u, std::thread::hardware_concurrency());
    nthreads = std::min(nthreads, static_cast<unsigned>(H));
    std::vector<int> captured_per_thread(nthreads, 0);
    int rows_per_thread = (H + static_cast<int>(nthreads) - 1) / static_cast<int>(nthreads);
    {
        std::vector<std::jthread> workers;
        workers.reserve(nthreads);
        for (unsigned t = 0; t < nthreads; ++t) {
            int y0 = static_cast<int>(t) * rows_per_thread;
            int y1 = std::min(H, y0 + rows_per_thread);
            if (y0 >= y1) continue;
            workers.emplace_back([&render_rows, &captured_per_thread, t, y0, y1] {
                captured_per_thread[t] = render_rows(y0, y1);
            });
        }
    }  // jthreads join here
    int captured = 0;
    for (int c : captured_per_thread) captured += c;

    auto t1 = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::print("\rblackhole_demo: {}x{}, captured {} pixels ({:.1f}%), {:.0f} ms\n",
               W, H, captured, 100.0 * captured / (W * H), ms);

    if (!spatium::examples::confirm_overwrite("blackhole.png", force)) return 0;
    if (!stbi_write_png("blackhole.png", W, H, 3, img.data(), W * 3)) {
        std::print(stderr, "failed to write blackhole.png\n");
        return 1;
    }
    std::print("wrote blackhole.png\n");
    return 0;
}
