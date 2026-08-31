// Wave equation on a sphere via finite-difference cellular automaton.
//
// State lives on a regular (θ, φ) grid; the wave PDE u_tt = c² Δu is stepped
// with explicit central differences:
//
//   u^{n+1}_{ij} = 2 u^n_{ij} − u^{n−1}_{ij}
//                + (cΔt/h)² ( Σ_neighbours u^n − k u^n_{ij} )
//
// The discrete laplacian uses 4 neighbours with periodic wrap on φ and pole
// folding on θ. The grid is sampled via bilinear interpolation by a
// ParametricSurface whose radial displacement equals 1 + ε·u(θ, φ), so the
// scene is rendered exactly through the existing Newton-UV ray-parametric
// solver — no rasterizer required.
//
// Output: wave_NNNN.png frames (256x192 RGB).

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <spatium/vendor/stb_image_write.h>

#include "io_helpers.hpp"

#include <spatium/algebra/vector.hpp>
#include <spatium/geometry/line.hpp>
#include <spatium/geometry/make.hpp>
#include <spatium/geometry/ray_parametric.hpp>
#include <spatium/spaces/parametric.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <numbers>
#include <print>
#include <string>
#include <string_view>
#include <vector>

using spatium::Vec;
using spatium::geometry::ray;
using spatium::geometry::ray_parametric_first;
using spatium::geometry::RayParametricConfig;
using spatium::ParametricSurface;
using spatium::parametric;
using std::numbers::pi;

namespace {

constexpr int    NTHETA = 96;        // grid rows (θ in (0, π))
constexpr int    NPHI   = 192;       // grid cols (φ in [0, 2π), periodic)
constexpr int    FRAMES = 96;
constexpr int    SUBSTEPS = 3;       // wave steps per render frame
constexpr double C_WAVE = 1.0;
constexpr double DT = 0.018;
constexpr double DISPLACE = 0.18;    // radial amplitude
constexpr int    W_IMG = 320;
constexpr int    H_IMG = 240;

struct Field {
    int rows, cols;
    std::vector<double> data;
    Field(int r, int c) : rows(r), cols(c), data(r * c, 0.0) {}
    double& at(int i, int j) { return data[i * cols + j]; }
    double  at(int i, int j) const { return data[i * cols + j]; }
};

double wrap_phi(double phi) {
    constexpr double TWOPI = 2.0 * pi;
    phi = std::fmod(phi, TWOPI);
    if (phi < 0) phi += TWOPI;
    return phi;
}

// Bilinear sampling, periodic in φ, clamp at θ poles.
double sample(const Field& f, double theta, double phi) {
    phi = wrap_phi(phi);
    double tt = std::clamp(theta / pi, 0.0, 1.0) * (f.rows - 1);
    double pp = phi / (2.0 * pi) * f.cols;
    int i0 = int(std::floor(tt));
    int i1 = std::min(i0 + 1, f.rows - 1);
    double a = tt - i0;
    int j0 = int(std::floor(pp)) % f.cols;
    if (j0 < 0) j0 += f.cols;
    int j1 = (j0 + 1) % f.cols;
    double b = pp - std::floor(pp);
    double v00 = f.at(i0, j0), v01 = f.at(i0, j1);
    double v10 = f.at(i1, j0), v11 = f.at(i1, j1);
    double v0 = v00 * (1 - b) + v01 * b;
    double v1 = v10 * (1 - b) + v11 * b;
    return v0 * (1 - a) + v1 * a;
}

// One wave step on the (θ, φ) grid.
void step_wave(const Field& u_prev, const Field& u_now, Field& u_next,
               double c, double dt, double dtheta, double dphi) {
    const double a_theta = (c * dt / dtheta) * (c * dt / dtheta);
    const double a_phi   = (c * dt / dphi  ) * (c * dt / dphi);
    for (int i = 0; i < u_now.rows; ++i) {
        // sin(θ) factor for spherical laplacian on regular grid.
        // Δ u = (1/sin θ) ∂_θ(sin θ ∂_θ u) + (1/sin² θ) ∂²_φ u
        // We use the conservative-symmetric form:
        //   Δu ≈ (u_north sin_n + u_south sin_s − u(sin_n+sin_s)) / (sin_i Δθ²)
        //      + (u_east + u_west − 2u) / (sin² θ Δφ²)
        double theta = (i + 0.5) / u_now.rows * pi;
        double sin_i = std::max(std::sin(theta), 1e-3);
        double theta_n = std::max(theta - dtheta, 1e-3);
        double theta_s = std::min(theta + dtheta, pi - 1e-3);
        double sin_n = std::sin(theta_n);
        double sin_s = std::sin(theta_s);
        for (int j = 0; j < u_now.cols; ++j) {
            int jl = (j - 1 + u_now.cols) % u_now.cols;
            int jr = (j + 1) % u_now.cols;
            // θ-neighbours fold across pole: when stepping past pole, j → (j + N/2) % N.
            int i_n = i - 1, j_n = j;
            if (i_n < 0)               { i_n = 0;            j_n = (j + u_now.cols / 2) % u_now.cols; }
            int i_s = i + 1, j_s = j;
            if (i_s >= u_now.rows)     { i_s = u_now.rows-1; j_s = (j + u_now.cols / 2) % u_now.cols; }

            double u_c = u_now.at(i, j);
            double u_n = u_now.at(i_n, j_n);
            double u_s = u_now.at(i_s, j_s);
            double u_e = u_now.at(i, jr);
            double u_w = u_now.at(i, jl);

            double lap_theta = (u_n * sin_n + u_s * sin_s - u_c * (sin_n + sin_s)) / sin_i;
            double lap_phi   = (u_e + u_w - 2.0 * u_c) / (sin_i * sin_i);
            double delta = a_theta * lap_theta + a_phi * lap_phi;
            u_next.at(i, j) = 2.0 * u_c - u_prev.at(i, j) + delta;
        }
    }
    // Mild damping at the boundary in θ to keep the simulation stable
    // long-term — purely cosmetic.
    for (int j = 0; j < u_now.cols; ++j) {
        u_next.at(0,             j) *= 0.999;
        u_next.at(u_now.rows - 1, j) *= 0.999;
    }
}

// Initial condition: gaussian bump in front of the camera.
void seed_pulse(Field& u) {
    const double theta0 = pi * 0.5;
    const double phi0   = 0.0;
    const double sigma  = 0.18;
    for (int i = 0; i < u.rows; ++i) {
        double theta = (i + 0.5) / u.rows * pi;
        for (int j = 0; j < u.cols; ++j) {
            double phi = (j + 0.5) / u.cols * 2.0 * pi;
            // Geodesic distance on unit sphere via dot product.
            double cos_d = std::sin(theta) * std::sin(theta0) * std::cos(phi - phi0)
                         + std::cos(theta) * std::cos(theta0);
            cos_d = std::clamp(cos_d, -1.0, 1.0);
            double d = std::acos(cos_d);
            u.at(i, j) = std::exp(-(d * d) / (sigma * sigma));
        }
    }
}

void print_usage() {
    std::print("Usage: wave_ca_demo [--force] [--help]\n"
               "  Wave equation on a sphere via finite-difference cellular automaton.\n"
               "  --force   overwrite existing wave_NNNN.png frames\n"
               "  --help    show this message\n"
               "  Outputs:  wave_0000.png ... wave_{:04d}.png "
               "({}x{} RGB, {} frames)\n",
               FRAMES - 1, W_IMG, H_IMG, FRAMES);
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

    Field u_prev(NTHETA, NPHI);
    Field u_now (NTHETA, NPHI);
    Field u_next(NTHETA, NPHI);
    seed_pulse(u_now);

    const double dtheta = pi / NTHETA;
    const double dphi   = 2.0 * pi / NPHI;
    const double cfl = std::min(dtheta, dphi) / C_WAVE;
    if (DT > 0.5 * cfl) {
        std::print("warning: dt={} above CFL ~{}\n", DT, 0.5 * cfl);
    }

    // Surface: r = 1 + DISPLACE · u(θ, φ).
    const Field* sample_field = &u_now;  // updated each frame
    auto wave_surf = parametric<double>(
        [&sample_field](double theta, double phi) -> Vec<double, 3> {
            double u = sample(*sample_field, theta, phi);
            double r = 1.0 + DISPLACE * u;
            return {r * std::sin(theta) * std::cos(phi),
                    r * std::sin(theta) * std::sin(phi),
                    r * std::cos(theta)};
        },
        typename ParametricSurface<double>::Domain{0.0, pi, 0.0, 2.0 * pi},
        false, true);

    // Camera + lighting.
    Vec<double, 3> cam_pos{2.7, 2.0, 1.6};
    Vec<double, 3> cam_target{0, 0, 0};
    Vec<double, 3> cam_up{0, 0, 1};
    Vec<double, 3> light_dir = Vec<double, 3>{1.0, -0.6, 1.4}.normalized();

    auto fwd = (cam_target - cam_pos).normalized();
    auto right = fwd.cross(cam_up).normalized();
    auto up    = right.cross(fwd).normalized();
    const double tan_half = std::tan(35.0 * pi / 360.0);
    const double aspect = double(W_IMG) / H_IMG;

    RayParametricConfig<double> cfg;  // defaults

    auto t0 = std::chrono::steady_clock::now();
    std::vector<std::uint8_t> img(3 * W_IMG * H_IMG, 0);

    for (int frame = 0; frame < FRAMES; ++frame) {
        // Advance physics.
        for (int s = 0; s < SUBSTEPS; ++s) {
            step_wave(u_prev, u_now, u_next, C_WAVE, DT, dtheta, dphi);
            std::swap(u_prev.data, u_now.data);
            std::swap(u_now.data,  u_next.data);
        }
        sample_field = &u_now;

        // Render via Newton-UV ray-parametric.
        std::fill(img.begin(), img.end(), std::uint8_t(0));
        for (int py = 0; py < H_IMG; ++py) {
            for (int px = 0; px < W_IMG; ++px) {
                double nx = (2.0 * (px + 0.5) / W_IMG  - 1.0) * aspect * tan_half;
                double ny = (1.0 - 2.0 * (py + 0.5) / H_IMG)         * tan_half;
                Vec<double, 3> dir = (fwd + right * nx + up * ny).normalized();
                auto r = spatium::unwrap(ray(cam_pos, dir));
                auto hit = ray_parametric_first(r, wave_surf, cfg);

                std::uint8_t* p = &img[3 * (py * W_IMG + px)];
                if (hit) {
                    auto n = hit->normal;
                    if (n.dot(dir) > 0) n = Vec<double, 3>{-n};
                    double diff = std::max(0.0, n.dot(light_dir));
                    // Color tinted by displacement value at hit (red = crest, blue = trough).
                    double u_at = sample(u_now, hit->u, hit->v);
                    double shade = 0.20 + 0.80 * diff;
                    double mix = std::clamp(u_at, 0.0, 1.0);
                    Vec<double, 3> warm{0.95, 0.45, 0.25};
                    Vec<double, 3> cool{0.30, 0.50, 0.85};
                    Vec<double, 3> base = u_at > 0
                        ? Vec<double, 3>{warm * mix + cool * (1.0 - mix)}
                        : cool;
                    p[0] = std::uint8_t(std::clamp(base[0] * shade, 0.0, 1.0) * 255);
                    p[1] = std::uint8_t(std::clamp(base[1] * shade, 0.0, 1.0) * 255);
                    p[2] = std::uint8_t(std::clamp(base[2] * shade, 0.0, 1.0) * 255);
                } else {
                    double v01 = 0.5 * (dir[2] + 1.0);
                    p[0] = std::uint8_t( 14 + (int)(40 * v01));
                    p[1] = std::uint8_t( 16 + (int)(45 * v01));
                    p[2] = std::uint8_t( 28 + (int)(60 * v01));
                }
            }
        }

        char path[64];
        std::snprintf(path, sizeof(path), "wave_%04d.png", frame);
        if (spatium::examples::confirm_overwrite(path, force))
            stbi_write_png(path, W_IMG, H_IMG, 3, img.data(), W_IMG * 3);
        if ((frame % 8) == 0)
            std::print("\r  frame {}/{}", frame, FRAMES);
    }
    auto t1 = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::print("\rwave_ca_demo: {} frames, {}x{}, {:.1f} s ({:.0f} ms/frame)\n",
               FRAMES, W_IMG, H_IMG, ms / 1000.0, ms / FRAMES);
    return 0;
}
