// Collatz orbits as a 4D parametric surface, projected into 3D.
//
// We treat (n, k) as parameters: n is the starting integer (1 ≤ n ≤ N),
// k is the iteration step (0 ≤ k ≤ K). Apply T(x) = x/2 if even, 3x+1 else.
// Each (n, k) gives a value v(n, k); the orbit lives in 4D as
//
//     ψ(n, k) = ( n, k, log₂(v + 1), parity(v) ) ∈ R⁴
//
// We render two views:
//   collatz_height.png  — z = log₂(v+1), shaded by parity (even=cool, odd=warm).
//   collatz_4d.png      — full ℝ⁴ → ℝ³ via stereographic projection that
//                         pulls the 4-th axis (parity, normalised to [-1,1])
//                         into the 3D scene as bend along z.
//
// Both views go through the existing CPU mesh renderer (custom raster, tiny
// painter's-algorithm pipeline — kept self-contained so the demo stays
// dependency-free beyond Spatium's algebra header and stb_image_write).

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <spatium/vendor/stb_image_write.h>

#include "io_helpers.hpp"

#include <spatium/algebra/vector.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <numbers>
#include <print>
#include <string>
#include <string_view>
#include <vector>

using spatium::Vec;
using std::numbers::pi;

namespace {

constexpr int N_MAX = 192;   // starting integers 1..N
constexpr int K_MAX = 96;    // steps
constexpr int W_IMG = 1024;
constexpr int H_IMG = 768;

// Returns next Collatz step; saturates at 1 (we just stay there).
std::uint64_t collatz_step(std::uint64_t v) {
    if (v <= 1) return 1;
    return (v % 2 == 0) ? (v / 2) : (3 * v + 1);
}

struct OrbitGrid {
    int N, K;
    std::vector<double> log2_value;  // size N*K, log₂(v+1)
    std::vector<int>    parity;      // size N*K, 0 even / 1 odd

    OrbitGrid(int n, int k)
        : N(n), K(k), log2_value(n * k, 0), parity(n * k, 0) {}

    void at_set(int i, int j, double lv, int p) {
        log2_value[i * K + j] = lv;
        parity    [i * K + j] = p;
    }
    double lv(int i, int j) const { return log2_value[i * K + j]; }
    int    pa(int i, int j) const { return parity[i * K + j]; }
};

OrbitGrid build_orbits(int N, int K) {
    OrbitGrid g(N, K);
    for (int i = 0; i < N; ++i) {
        std::uint64_t v = std::uint64_t(i + 1);
        for (int j = 0; j < K; ++j) {
            g.at_set(i, j, std::log2(double(v) + 1.0), int(v & 1));
            v = collatz_step(v);
        }
    }
    return g;
}

// ── Tiny CPU rasteriser ──────────────────────────────────────────

struct Frame {
    int w, h;
    std::vector<std::uint8_t> rgb;   // 3*w*h
    std::vector<double>       depth; // w*h, +∞ = empty
    Frame(int W, int H)
        : w(W), h(H),
          rgb(3 * W * H, 0),
          depth(W * H, std::numeric_limits<double>::infinity()) {
        // Sky gradient as background.
        for (int y = 0; y < H; ++y) {
            double t = double(y) / (H - 1);
            std::uint8_t r = std::uint8_t(15 + 14 * (1 - t));
            std::uint8_t g = std::uint8_t(18 + 18 * (1 - t));
            std::uint8_t b = std::uint8_t(28 + 28 * (1 - t));
            for (int x = 0; x < W; ++x) {
                rgb[3 * (y * W + x) + 0] = r;
                rgb[3 * (y * W + x) + 1] = g;
                rgb[3 * (y * W + x) + 2] = b;
            }
        }
    }
};

void put_px(Frame& f, int x, int y, double z, const Vec<double, 3>& c) {
    if (x < 0 || y < 0 || x >= f.w || y >= f.h) return;
    int idx = y * f.w + x;
    if (z >= f.depth[idx]) return;
    f.depth[idx] = z;
    f.rgb[3 * idx + 0] = std::uint8_t(std::clamp(c[0], 0.0, 1.0) * 255);
    f.rgb[3 * idx + 1] = std::uint8_t(std::clamp(c[1], 0.0, 1.0) * 255);
    f.rgb[3 * idx + 2] = std::uint8_t(std::clamp(c[2], 0.0, 1.0) * 255);
}

struct Vert3 {
    Vec<double, 3> world;
    Vec<double, 3> color;
};

// Project to screen with simple perspective.
struct Camera {
    Vec<double, 3> eye, target, up;
    double fov_deg, near_z, far_z;
};

void build_basis(const Camera& c, Vec<double, 3>& fwd,
                 Vec<double, 3>& right, Vec<double, 3>& up_o) {
    fwd   = (c.target - c.eye).normalized();
    right = fwd.cross(c.up).normalized();
    up_o  = right.cross(fwd).normalized();
}

bool project(const Camera& c, const Vec<double, 3>& fwd,
             const Vec<double, 3>& right, const Vec<double, 3>& up_o,
             const Vec<double, 3>& p, int W, int H,
             double& sx, double& sy, double& z) {
    Vec<double, 3> v = p - c.eye;
    z = v.dot(fwd);
    if (z < c.near_z || z > c.far_z) return false;
    double x = v.dot(right);
    double y = v.dot(up_o);
    double tan_half = std::tan(c.fov_deg * pi / 360.0);
    double aspect = double(W) / H;
    sx = (x / z) / (aspect * tan_half);
    sy = (y / z) / tan_half;
    if (std::abs(sx) > 1.5 || std::abs(sy) > 1.5) return false;
    sx = (sx + 1.0) * 0.5 * W;
    sy = (1.0 - sy) * 0.5 * H;
    return true;
}

// Filled triangle with linear color & depth interpolation (barycentric).
void raster_tri(Frame& f, const Camera& cam, const Vec<double,3>& fwd,
                const Vec<double,3>& right, const Vec<double,3>& upo,
                const Vert3& A, const Vert3& B, const Vert3& C) {
    double ax, ay, az, bx, by, bz, cx, cy, cz;
    if (!project(cam, fwd, right, upo, A.world, f.w, f.h, ax, ay, az)) return;
    if (!project(cam, fwd, right, upo, B.world, f.w, f.h, bx, by, bz)) return;
    if (!project(cam, fwd, right, upo, C.world, f.w, f.h, cx, cy, cz)) return;
    int xmin = std::max(0,           int(std::floor(std::min({ax, bx, cx}))));
    int xmax = std::min(f.w - 1,     int(std::ceil (std::max({ax, bx, cx}))));
    int ymin = std::max(0,           int(std::floor(std::min({ay, by, cy}))));
    int ymax = std::min(f.h - 1,     int(std::ceil (std::max({ay, by, cy}))));
    if (xmax < xmin || ymax < ymin) return;
    double denom = (by - cy)*(ax - cx) + (cx - bx)*(ay - cy);
    if (std::abs(denom) < 1e-9) return;
    for (int y = ymin; y <= ymax; ++y) {
        for (int x = xmin; x <= xmax; ++x) {
            double w0 = ((by - cy)*(x - cx) + (cx - bx)*(y - cy)) / denom;
            double w1 = ((cy - ay)*(x - cx) + (ax - cx)*(y - cy)) / denom;
            double w2 = 1.0 - w0 - w1;
            if (w0 < 0 || w1 < 0 || w2 < 0) continue;
            double z = w0 * az + w1 * bz + w2 * cz;
            Vec<double, 3> col = A.color * w0 + B.color * w1 + C.color * w2;
            put_px(f, x, y, z, col);
        }
    }
}

// ── Surface meshes ───────────────────────────────────────────────

Vec<double, 3> parity_color(int parity, double height01) {
    // even = cool blue/teal, odd = warm orange/red.
    Vec<double, 3> cool {0.25, 0.55, 0.95};
    Vec<double, 3> warm {0.95, 0.40, 0.20};
    Vec<double, 3> base = parity ? warm : cool;
    double shade = 0.45 + 0.55 * height01;
    return base * shade;
}

// View 1: (n, k, height). Direct (n, k) → 3D.
std::vector<Vert3> build_height_view(const OrbitGrid& g, double xy_scale,
                                     double z_scale) {
    std::vector<Vert3> vs;
    vs.reserve(g.N * g.K);
    double maxlv = 1e-9;
    for (auto v : g.log2_value) maxlv = std::max(maxlv, v);
    for (int i = 0; i < g.N; ++i) {
        for (int j = 0; j < g.K; ++j) {
            double x = (i / double(g.N - 1) - 0.5) * xy_scale;
            double y = (j / double(g.K - 1) - 0.5) * xy_scale * 0.7;
            double h01 = g.lv(i, j) / maxlv;
            double z = h01 * z_scale;
            Vert3 v{Vec<double,3>{x, y, z},
                    parity_color(g.pa(i, j), h01)};
            vs.push_back(v);
        }
    }
    return vs;
}

// View 2: ℝ⁴ → ℝ³ stereographic. We treat parity (∈{0,1}) → w ∈ {-1, +1}
// and embed (x, y, z, w) in S³, then stereo-project from the south pole.
std::vector<Vert3> build_4d_view(const OrbitGrid& g, double xy_scale,
                                 double z_scale, double w_amp) {
    std::vector<Vert3> vs;
    vs.reserve(g.N * g.K);
    double maxlv = 1e-9;
    for (auto v : g.log2_value) maxlv = std::max(maxlv, v);
    for (int i = 0; i < g.N; ++i) {
        for (int j = 0; j < g.K; ++j) {
            double x = (i / double(g.N - 1) - 0.5) * xy_scale;
            double y = (j / double(g.K - 1) - 0.5) * xy_scale * 0.7;
            double h01 = g.lv(i, j) / maxlv;
            double z = h01 * z_scale;
            double w = (g.pa(i, j) ? +1.0 : -1.0) * w_amp;
            // Stereographic: (X, Y, Z, W) → (X, Y, Z) / (1 − W/R) with R=2.
            double scale = 1.0 / std::max(1e-3, 1.0 - w / 2.0);
            Vert3 v{Vec<double,3>{x * scale, y * scale, z * scale},
                    parity_color(g.pa(i, j), h01)};
            vs.push_back(v);
        }
    }
    return vs;
}

void render_grid(const std::vector<Vert3>& vs, int N, int K,
                 const Camera& cam, Frame& f) {
    Vec<double, 3> fwd, right, upo;
    build_basis(cam, fwd, right, upo);
    auto ix = [K](int i, int j) { return i * K + j; };
    for (int i = 0; i < N - 1; ++i) {
        for (int j = 0; j < K - 1; ++j) {
            const auto& A = vs[ix(i,   j)];
            const auto& B = vs[ix(i+1, j)];
            const auto& C = vs[ix(i+1, j+1)];
            const auto& D = vs[ix(i,   j+1)];
            raster_tri(f, cam, fwd, right, upo, A, B, C);
            raster_tri(f, cam, fwd, right, upo, A, C, D);
        }
    }
}

}  // namespace

namespace {

void print_usage() {
    std::print("Usage: collatz_demo [--force] [--help]\n"
               "  Renders Collatz orbits as a 4D parametric surface projected into 3D.\n"
               "  --force   overwrite existing output files\n"
               "  --help    show this message\n"
               "  Outputs:  collatz_height.png, collatz_4d.png\n");
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
    OrbitGrid g = build_orbits(N_MAX, K_MAX);

    {
        Frame f(W_IMG, H_IMG);
        Camera cam{
            .eye    = { 4.0, -8.5,  6.0},
            .target = { 0.0,  0.0,  2.0},
            .up     = { 0.0,  0.0,  1.0},
            .fov_deg= 32.0, .near_z = 0.1, .far_z = 200.0
        };
        auto vs = build_height_view(g, /*xy_scale*/ 8.0, /*z_scale*/ 5.0);
        render_grid(vs, g.N, g.K, cam, f);
        if (spatium::examples::confirm_overwrite("collatz_height.png", force)) {
            stbi_write_png("collatz_height.png", f.w, f.h, 3, f.rgb.data(), f.w * 3);
            std::print("wrote collatz_height.png\n");
        }
    }

    {
        Frame f(W_IMG, H_IMG);
        Camera cam{
            .eye    = { 5.5, -10.5,  6.5},
            .target = { 0.0,   0.0,  2.0},
            .up     = { 0.0,   0.0,  1.0},
            .fov_deg= 32.0, .near_z = 0.1, .far_z = 200.0
        };
        auto vs = build_4d_view(g, /*xy_scale*/ 8.0, /*z_scale*/ 5.0,
                                   /*w_amp*/ 0.8);
        render_grid(vs, g.N, g.K, cam, f);
        if (spatium::examples::confirm_overwrite("collatz_4d.png", force)) {
            stbi_write_png("collatz_4d.png", f.w, f.h, 3, f.rgb.data(), f.w * 3);
            std::print("wrote collatz_4d.png\n");
        }
    }

    auto t1 = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::print("collatz_demo: N={} starts × K={} steps = {} cells, {:.0f} ms\n",
               N_MAX, K_MAX, N_MAX * K_MAX, ms);
    return 0;
}
