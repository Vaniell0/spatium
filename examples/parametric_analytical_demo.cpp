// Analytical rendering of ParametricSurface via Newton UV.
//
// No tesselation, no BVH — each pixel's ray is intersected with the
// surface by solving f(u,v) = o + t·d for (u, v, t). Normals come from
// the Jacobian's cross product at the converged (u, v), so shading is
// exact at sub-pixel resolution.
//
// Exports:
//   klein_analytical.png     — classic 3D immersion of the Klein bottle
//   mobius_analytical.png    — non-orientable Möbius strip
//   bumpy_analytical.png     — displaced sphere (arbitrary f(u, v))
//   parametric_gallery.png   — 2x2 grid: analytic torus + three above
//   schwarzschild_grid.png   — Flamm paraboloid, spacetime curvature grid
//   wormhole_grid.png        — Ellis wormhole embedding, curvature grid
//   curvature_grids.png      — the two grids above, side by side
//   glow_sphere.png          — sphere with a ray_quadric_proximity() glow

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <spatium/vendor/stb_image_write.h>

#include "io_helpers.hpp"

#include <spatium/algebra/vector.hpp>
#include <spatium/geometry/line.hpp>
#include <spatium/geometry/make.hpp>
#include <spatium/geometry/ray_parametric.hpp>
#include <spatium/geometry/ray_surface.hpp>
#include <spatium/render/supersample.hpp>
#include <spatium/spaces/parametric.hpp>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <functional>
#include <numbers>
#include <optional>
#include <print>
#include <string>
#include <thread>
#include <utility>
#include <vector>

using namespace spatium;
using namespace spatium::geometry;
using Clock = std::chrono::steady_clock;

namespace {

struct Camera {
    Vec<double, 3> position;
    Vec<double, 3> target;
    Vec<double, 3> up{0, 0, 1};
    double fov_deg = 38.0;
};

struct ShadingConfig {
    Vec<double, 3> base_color{0.85, 0.55, 0.30};
    Vec<double, 3> light_dir{1.0, -0.6, 1.4};  // normalized on use
    Vec<double, 3> sky_top{0.18, 0.22, 0.30};
    Vec<double, 3> sky_bottom{0.08, 0.10, 0.14};
    double ambient = 0.20;
    bool two_sided = true;
    bool tint_by_uv = false;

    // "Spacetime curvature grid": darken a band around evenly-spaced UV
    // lines so the coordinate mesh is visibly drawn on the curved
    // surface, same idea as the classic embedding-diagram illustration
    // but computed from the real metric instead of an ad hoc paraboloid.
    bool grid_lines = false;
    double grid_spacing_u = std::numbers::pi / 6.0;  // 12 lines around phi
    double grid_spacing_v = 1.0;
    double grid_half_width = 0.025;  // in UV units
    double grid_darken = 0.12;       // multiplier applied on a gridline
};

bool near_periodic_line(double coord, double spacing, double half_width) {
    double m = std::fmod(coord, spacing);
    if (m < 0.0) m += spacing;
    return m < half_width || m > spacing - half_width;
}

// HSV → RGB for UV tinting (shows parameterization orientation).
Vec<double, 3> hsv(double h, double s, double v) {
    double c = v * s;
    double hp = h * 6.0;
    double x = c * (1.0 - std::abs(std::fmod(hp, 2.0) - 1.0));
    Vec<double, 3> rgb{0, 0, 0};
    if      (hp < 1.0) rgb = {c, x, 0};
    else if (hp < 2.0) rgb = {x, c, 0};
    else if (hp < 3.0) rgb = {0, c, x};
    else if (hp < 4.0) rgb = {0, x, c};
    else if (hp < 5.0) rgb = {x, 0, c};
    else               rgb = {c, 0, x};
    double m = v - c;
    return Vec<double, 3>{rgb + Vec<double, 3>{m, m, m}};
}

struct RenderStats {
    std::size_t hits = 0;
    std::size_t pixels = 0;
    double elapsed_ms = 0.0;
};

// CPU raytrace a ParametricSurface into an RGBA image tile.
// Returns the tile plus stats.
//
// Two real, measured fixes applied here (2026-08-25) after a full
// 384x384 gallery render took 2m22s and per-surface stats showed why:
// Klein bottle/Möbius/bumpy-sphere ran at 242,000-367,000 ns/ray vs. the
// closed-form quartic torus's 213 ns/ray:
//   1. cell_radius no longer recomputed per ray (estimate_cell_radius()
//      is a 25-surface-eval, ray-independent cost that ray_parametric()
//      used to pay on every single pixel) -- computed once here instead.
//   2. Each scanline carries the previous pixel's converged (u, v) as a
//      uv_hint: adjacent pixels almost always hit near the same surface
//      point, so most interior pixels now cost one Newton solve instead
//      of up to grid_u*grid_v=144. Reset at the start of each row (a
//      hint only helps against its own immediate neighbor).
//   3. Rows are split across std::jthread workers -- no shared mutable
//      state between them beyond each thread's own disjoint row range
//      of the output buffer, so no synchronization is needed beyond the
//      implicit join when the workers go out of scope.
RenderStats render_parametric(const ParametricSurface<double>& surf,
                              const Camera& cam,
                              const ShadingConfig& shade,
                              int width, int height,
                              std::vector<std::uint8_t>& rgba,
                              int offset_x = 0, int offset_y = 0,
                              int stride = 0) {
    if (stride == 0) stride = width;
    if (static_cast<int>(rgba.size()) < stride * (offset_y + height) * 4)
        rgba.resize(stride * (offset_y + height) * 4);

    auto fwd = Vec<double, 3>{(cam.target - cam.position).normalized()};
    auto right = Vec<double, 3>{fwd.cross(cam.up).normalized()};
    auto up    = Vec<double, 3>{right.cross(fwd).normalized()};
    double tan_half = std::tan(cam.fov_deg * std::numbers::pi / 360.0);
    double aspect = static_cast<double>(width) / height;
    auto light = Vec<double, 3>{shade.light_dir.normalized()};
    auto dom = surf.domain();
    double u_range = dom.u_max - dom.u_min;
    double v_range = dom.v_max - dom.v_min;

    RayParametricConfig<double> cfg;  // defaults: 12x12 seed grid, 24 iterations
    std::optional<double> cell_radius = estimate_cell_radius(surf, cfg);

    RenderStats stats{};
    stats.pixels = static_cast<std::size_t>(width) * height;

    auto render_rows = [&](int y0, int y1) -> std::size_t {
        std::size_t local_hits = 0;
        for (int py = y0; py < y1; ++py) {
            std::optional<std::pair<double, double>> hint;
            for (int px = 0; px < width; ++px) {
                std::uint8_t* pixel = &rgba[((offset_y + py) * stride + (offset_x + px)) * 4];
                bool any_hit = false;

                // Per-ray color, called aa*aa times per pixel by
                // supersample_pixel (spatium/render/supersample.hpp) at
                // sub-pixel (nx, ny). `hint` is shared across every
                // subsample call in scanline order (not reset per
                // pixel), so the Newton-UV coherence fast path still
                // applies within a pixel's own subsamples, not just
                // between pixels.
                auto ray_color = [&](double nx, double ny) -> Vec<double, 3> {
                    auto dir = Vec<double, 3>{(fwd + right * nx + up * ny).normalized()};
                    auto r = unwrap(ray(cam.position, dir));

                    auto hit = ray_parametric_first(r, surf, cfg, hint, cell_radius);
                    if (hit) {
                        any_hit = true;
                        hint = std::make_pair(hit->u, hit->v);
                        auto n = hit->normal;
                        // Two-sided shading: orient normal toward the camera
                        if (shade.two_sided && n.dot(dir) > 0.0) n = Vec<double, 3>{-n};
                        double diff = std::max(0.0, n.dot(light));
                        double shaded = shade.ambient + (1.0 - shade.ambient) * diff;

                        Vec<double, 3> base = shade.base_color;
                        if (shade.tint_by_uv) {
                            double u01 = (hit->u - dom.u_min) / u_range;
                            double v01 = (hit->v - dom.v_min) / v_range;
                            base = hsv(std::fmod(u01 + 0.5 * v01, 1.0), 0.55, 0.85);
                        }
                        if (shade.grid_lines &&
                            (near_periodic_line(hit->u - dom.u_min, shade.grid_spacing_u,
                                                 shade.grid_half_width) ||
                             near_periodic_line(hit->v - dom.v_min, shade.grid_spacing_v,
                                                 shade.grid_half_width)))
                            shaded *= shade.grid_darken;
                        return Vec<double, 3>{base * shaded * 255.0};
                    }
                    hint.reset();
                    double v01 = 0.5 * (dir[2] + 1.0);
                    Vec<double, 3> sky{shade.sky_bottom + (shade.sky_top - shade.sky_bottom) * v01};
                    return Vec<double, 3>{sky * 255.0};
                };

                std::uint8_t rgb[3];
                spatium::render::supersample_pixel(px, py, width, height, tan_half, aspect,
                                                      ray_color, rgb);
                pixel[0] = rgb[0];
                pixel[1] = rgb[1];
                pixel[2] = rgb[2];
                pixel[3] = 255;
                if (any_hit) ++local_hits;
            }
        }
        return local_hits;
    };

    unsigned nthreads = std::max(1u, std::thread::hardware_concurrency());
    nthreads = std::min(nthreads, static_cast<unsigned>(std::max(1, height)));
    std::vector<std::size_t> hits_per_thread(nthreads, 0);
    int rows_per_thread = (height + static_cast<int>(nthreads) - 1) / static_cast<int>(nthreads);

    auto start = Clock::now();
    {
        std::vector<std::jthread> workers;
        workers.reserve(nthreads);
        for (unsigned t = 0; t < nthreads; ++t) {
            int y0 = static_cast<int>(t) * rows_per_thread;
            int y1 = std::min(height, y0 + rows_per_thread);
            if (y0 >= y1) continue;
            workers.emplace_back([&render_rows, &hits_per_thread, t, y0, y1] {
                hits_per_thread[t] = render_rows(y0, y1);
            });
        }
    }  // jthreads join here
    stats.elapsed_ms = std::chrono::duration<double, std::milli>(Clock::now() - start).count();
    for (auto h : hits_per_thread) stats.hits += h;
    return stats;
}

// Analytic torus tile (quartic closed-form — complements the parametric path).
RenderStats render_torus_quartic(const Torus<double>& torus,
                                  const Camera& cam,
                                  const ShadingConfig& shade,
                                  int width, int height,
                                  std::vector<std::uint8_t>& rgba,
                                  int offset_x, int offset_y, int stride) {
    auto fwd = Vec<double, 3>{(cam.target - cam.position).normalized()};
    auto right = Vec<double, 3>{fwd.cross(cam.up).normalized()};
    auto up    = Vec<double, 3>{right.cross(fwd).normalized()};
    double tan_half = std::tan(cam.fov_deg * std::numbers::pi / 360.0);
    double aspect = static_cast<double>(width) / height;
    auto light = Vec<double, 3>{shade.light_dir.normalized()};

    RenderStats stats{};
    stats.pixels = static_cast<std::size_t>(width) * height;

    auto render_rows = [&](int y0, int y1) -> std::size_t {
        std::size_t local_hits = 0;
        for (int py = y0; py < y1; ++py) {
            for (int px = 0; px < width; ++px) {
                std::uint8_t* pixel = &rgba[((offset_y + py) * stride + (offset_x + px)) * 4];
                bool any_hit = false;

                auto ray_color = [&](double nx, double ny) -> Vec<double, 3> {
                    auto dir = Vec<double, 3>{(fwd + right * nx + up * ny).normalized()};
                    auto r = unwrap(ray(cam.position, dir));
                    auto hits = ray_torus(r, torus);

                    if (!hits.empty()) {
                        any_hit = true;
                        auto& h = hits.front();
                        auto n = h.normal;
                        if (shade.two_sided && n.dot(dir) > 0.0) n = Vec<double, 3>{-n};
                        double diff = std::max(0.0, n.dot(light));
                        double shaded = shade.ambient + (1.0 - shade.ambient) * diff;
                        return Vec<double, 3>{shade.base_color * shaded * 255.0};
                    }
                    double v01 = 0.5 * (dir[2] + 1.0);
                    Vec<double, 3> sky{shade.sky_bottom + (shade.sky_top - shade.sky_bottom) * v01};
                    return Vec<double, 3>{sky * 255.0};
                };

                std::uint8_t rgb[3];
                spatium::render::supersample_pixel(px, py, width, height, tan_half, aspect,
                                                      ray_color, rgb);
                pixel[0] = rgb[0];
                pixel[1] = rgb[1];
                pixel[2] = rgb[2];
                pixel[3] = 255;
                if (any_hit) ++local_hits;
            }
        }
        return local_hits;
    };

    unsigned nthreads = std::max(1u, std::thread::hardware_concurrency());
    nthreads = std::min(nthreads, static_cast<unsigned>(std::max(1, height)));
    std::vector<std::size_t> hits_per_thread(nthreads, 0);
    int rows_per_thread = (height + static_cast<int>(nthreads) - 1) / static_cast<int>(nthreads);

    auto start = Clock::now();
    {
        std::vector<std::jthread> workers;
        workers.reserve(nthreads);
        for (unsigned t = 0; t < nthreads; ++t) {
            int y0 = static_cast<int>(t) * rows_per_thread;
            int y1 = std::min(height, y0 + rows_per_thread);
            if (y0 >= y1) continue;
            workers.emplace_back([&render_rows, &hits_per_thread, t, y0, y1] {
                hits_per_thread[t] = render_rows(y0, y1);
            });
        }
    }  // jthreads join here
    stats.elapsed_ms = std::chrono::duration<double, std::milli>(Clock::now() - start).count();
    for (auto h : hits_per_thread) stats.hits += h;
    return stats;
}

// ── Surfaces ──────────────────────────────────────────────────

// Classic 3D immersion of the Klein bottle (bottle form).
ParametricSurface<double> make_klein_bottle() {
    return parametric<double>(
        [](double u, double v) -> Vec<double, 3> {
            const double pi = std::numbers::pi;
            double cu = std::cos(u), su = std::sin(u);
            double cv = std::cos(v), sv = std::sin(v);
            double r = 4.0 * (1.0 - cu * 0.5);
            double x, y;
            if (u < pi) {
                x = 6.0 * cu * (1.0 + su) + r * cu * cv;
                y = 16.0 * su + r * su * cv;
            } else {
                x = 6.0 * cu * (1.0 + su) + r * std::cos(v + pi);
                y = 16.0 * su;
            }
            double z = r * sv;
            return {x * 0.08, y * 0.08, z * 0.08};
        },
        typename ParametricSurface<double>::Domain{0.0, 2.0 * std::numbers::pi,
                                                    0.0, 2.0 * std::numbers::pi},
        true, true);
}

// Displaced sphere with bumps — showcases arbitrary f(u, v).
ParametricSurface<double> make_bumpy_sphere() {
    return parametric<double>(
        [](double u, double v) -> Vec<double, 3> {
            double r = 1.0 + 0.12 * std::sin(5.0 * u) * std::sin(4.0 * v);
            return {r * std::sin(v) * std::cos(u),
                    r * std::sin(v) * std::sin(u),
                    r * std::cos(v)};
        },
        typename ParametricSurface<double>::Domain{0.0, 2.0 * std::numbers::pi,
                                                    0.02, std::numbers::pi - 0.02},
        true, false);
}

// ── Spacetime curvature grids (embedding diagrams) ──────────────
//
// The classic "rubber sheet" illustration, done from the real metric
// instead of an ad hoc paraboloid: isometrically embed the equatorial
// spatial slice as a surface of revolution (r(v) cos u, r(v) sin u,
// z(v)) in flat R^3, so ambient distances on the surface match the
// metric's own proper distances. z(v) is the profile curve that makes
// this true: (dz/dv)^2 + (dr/dv)^2 = 1 (unit-speed radial direction).

// Schwarzschild: dz/dr = sqrt(2M/(r-2M)) (standard derivation from
// ds^2 = dr^2/(1-2M/r) + r^2 dphi^2), integrates to the closed-form
// Flamm paraboloid z(r) = 2*sqrt(2M(r-2M)).
ParametricSurface<double> make_schwarzschild_embedding(double M = 1.0, double r_max = 8.0) {
    return parametric<double>(
        [M](double u, double v) -> Vec<double, 3> {
            double r = v;
            double z = 2.0 * std::sqrt(2.0 * M * std::max(0.0, r - 2.0 * M));
            return {r * std::cos(u), r * std::sin(u), z};
        },
        typename ParametricSurface<double>::Domain{0.0, 2.0 * std::numbers::pi, 2.0 * M, r_max},
        true, false);
}

// Ellis/Morris-Thorne wormhole (same metric wormhole_demo.cpp integrates
// null geodesics against): ds^2 = dl^2 + r(l)^2 dphi^2 with
// r(l)=sqrt(b0^2+l^2) -- dl's own coefficient is already 1, so
// dz/dl = sqrt(1-(dr/dl)^2) = b0/r(l), integrating to the closed form
// z(l) = b0*asinh(l/b0). Symmetric double funnel (an hourglass), unlike
// Schwarzschild's single one -- there's no horizon, both sides connect
// through the throat at l=0.
ParametricSurface<double> make_wormhole_embedding(double b0 = 1.0, double l_max = 8.0) {
    return parametric<double>(
        [b0](double u, double v) -> Vec<double, 3> {
            double l = v;
            double r = std::sqrt(b0 * b0 + l * l);
            double z = b0 * std::asinh(l / b0);
            return {r * std::cos(u), r * std::sin(u), z};
        },
        typename ParametricSurface<double>::Domain{0.0, 2.0 * std::numbers::pi, -l_max, l_max},
        true, false);
}

void save_png(const char* path, int w, int h, const std::vector<std::uint8_t>& rgba,
              bool force) {
    if (!spatium::examples::confirm_overwrite(path, force)) return;
    stbi_write_png(path, w, h, 4, rgba.data(), w * 4);
}

void print_stats(const char* label, const RenderStats& s, int w, int h) {
    double ns_per_ray = s.pixels ? (s.elapsed_ms * 1e6 / static_cast<double>(s.pixels)) : 0.0;
    std::println("  {:22} {:>4}x{:<4}  {:>7} hits / {:>7}   {:>7.1f} ms   {:>7.0f} ns/ray",
                 label, w, h, s.hits, s.pixels, s.elapsed_ms, ns_per_ray);
}

int render_single(int width, int height, const std::string& out_prefix, bool force) {
    std::vector<std::uint8_t> buf;

    // Klein bottle
    {
        auto surf = make_klein_bottle();
        Camera cam{.position = {3.5, 3.0, 2.2}, .target = {0, 0, 0}};
        ShadingConfig shade;
        shade.base_color = {0.40, 0.75, 0.85};
        shade.tint_by_uv = true;
        buf.assign(static_cast<std::size_t>(width) * height * 4, 0);
        auto s = render_parametric(surf, cam, shade, width, height, buf);
        print_stats("klein bottle", s, width, height);
        save_png((out_prefix + "klein_analytical.png").c_str(), width, height, buf, force);
    }

    // Möbius strip
    {
        auto surf = make_mobius<double>(2.0, 0.6);
        Camera cam{.position = {3.2, 3.2, 2.6}, .target = {0, 0, 0}};
        ShadingConfig shade;
        shade.base_color = {0.90, 0.55, 0.30};
        shade.tint_by_uv = true;
        buf.assign(static_cast<std::size_t>(width) * height * 4, 0);
        auto s = render_parametric(surf, cam, shade, width, height, buf);
        print_stats("möbius strip", s, width, height);
        save_png((out_prefix + "mobius_analytical.png").c_str(), width, height, buf, force);
    }

    // Bumpy sphere
    {
        auto surf = make_bumpy_sphere();
        Camera cam{.position = {2.8, 2.2, 1.8}, .target = {0, 0, 0}};
        ShadingConfig shade;
        shade.base_color = {0.70, 0.65, 0.95};
        buf.assign(static_cast<std::size_t>(width) * height * 4, 0);
        auto s = render_parametric(surf, cam, shade, width, height, buf);
        print_stats("bumpy sphere", s, width, height);
        save_png((out_prefix + "bumpy_analytical.png").c_str(), width, height, buf, force);
    }

    return 0;
}

int render_gallery(int tile_w, int tile_h, const std::string& out_prefix, bool force) {
    int w = tile_w * 2, h = tile_h * 2;
    std::vector<std::uint8_t> buf(static_cast<std::size_t>(w) * h * 4, 0);

    // Top-left: analytic torus (quartic)
    {
        Torus<double> t{.major_radius = 1.0, .minor_radius = 0.3};
        Camera cam{.position = {3.5, 3.5, 2.8}, .target = {0, 0, 0}};
        ShadingConfig shade;
        shade.base_color = {0.95, 0.75, 0.35};
        auto s = render_torus_quartic(t, cam, shade, tile_w, tile_h, buf, 0, 0, w);
        print_stats("torus (quartic)", s, tile_w, tile_h);
    }
    // Top-right: Klein bottle
    {
        auto surf = make_klein_bottle();
        Camera cam{.position = {3.5, 3.0, 2.2}, .target = {0, 0, 0}};
        ShadingConfig shade;
        shade.base_color = {0.40, 0.75, 0.85};
        shade.tint_by_uv = true;
        auto s = render_parametric(surf, cam, shade, tile_w, tile_h, buf, tile_w, 0, w);
        print_stats("klein bottle", s, tile_w, tile_h);
    }
    // Bottom-left: Möbius strip
    {
        auto surf = make_mobius<double>(2.0, 0.6);
        Camera cam{.position = {3.2, 3.2, 2.6}, .target = {0, 0, 0}};
        ShadingConfig shade;
        shade.base_color = {0.90, 0.55, 0.30};
        shade.tint_by_uv = true;
        auto s = render_parametric(surf, cam, shade, tile_w, tile_h, buf, 0, tile_h, w);
        print_stats("möbius strip", s, tile_w, tile_h);
    }
    // Bottom-right: Bumpy sphere
    {
        auto surf = make_bumpy_sphere();
        Camera cam{.position = {2.8, 2.2, 1.8}, .target = {0, 0, 0}};
        ShadingConfig shade;
        shade.base_color = {0.70, 0.65, 0.95};
        auto s = render_parametric(surf, cam, shade, tile_w, tile_h, buf, tile_w, tile_h, w);
        print_stats("bumpy sphere", s, tile_w, tile_h);
    }

    save_png((out_prefix + "parametric_gallery.png").c_str(), w, h, buf, force);
    std::println("  saved gallery {}x{}", w, h);
    return 0;
}

// Spacetime curvature grids: Schwarzschild's Flamm paraboloid and the
// wormhole's symmetric double funnel, each with its coordinate mesh
// drawn on via ShadingConfig::grid_lines -- "the grid of curvature" as
// an actual rendered image, not a schematic illustration.
int render_grids(int tile_w, int tile_h, const std::string& out_prefix, bool force) {
    // Individual PNGs first.
    {
        auto surf = make_schwarzschild_embedding();
        Camera cam{.position = {13.0, 13.0, 9.0}, .target = {0.0, 0.0, 2.0}};
        ShadingConfig shade;
        shade.base_color = {0.35, 0.55, 0.95};
        shade.grid_lines = true;
        shade.grid_spacing_v = 1.0;  // one ring per unit of r
        std::vector<std::uint8_t> buf(static_cast<std::size_t>(tile_w) * tile_h * 4, 0);
        auto s = render_parametric(surf, cam, shade, tile_w, tile_h, buf);
        print_stats("schwarzschild grid", s, tile_w, tile_h);
        save_png((out_prefix + "schwarzschild_grid.png").c_str(), tile_w, tile_h, buf, force);
    }
    {
        auto surf = make_wormhole_embedding();
        Camera cam{.position = {11.0, 11.0, 6.0}, .target = {0.0, 0.0, 0.0}};
        ShadingConfig shade;
        shade.base_color = {0.85, 0.45, 0.65};
        shade.grid_lines = true;
        shade.grid_spacing_v = 1.0;  // one ring per unit of l
        std::vector<std::uint8_t> buf(static_cast<std::size_t>(tile_w) * tile_h * 4, 0);
        auto s = render_parametric(surf, cam, shade, tile_w, tile_h, buf);
        print_stats("wormhole grid", s, tile_w, tile_h);
        save_png((out_prefix + "wormhole_grid.png").c_str(), tile_w, tile_h, buf, force);
    }

    // Side-by-side composite.
    int w = tile_w * 2, h = tile_h;
    std::vector<std::uint8_t> buf(static_cast<std::size_t>(w) * h * 4, 0);
    {
        auto surf = make_schwarzschild_embedding();
        Camera cam{.position = {13.0, 13.0, 9.0}, .target = {0.0, 0.0, 2.0}};
        ShadingConfig shade;
        shade.base_color = {0.35, 0.55, 0.95};
        shade.grid_lines = true;
        shade.grid_spacing_v = 1.0;
        render_parametric(surf, cam, shade, tile_w, tile_h, buf, 0, 0, w);
    }
    {
        auto surf = make_wormhole_embedding();
        Camera cam{.position = {11.0, 11.0, 6.0}, .target = {0.0, 0.0, 0.0}};
        ShadingConfig shade;
        shade.base_color = {0.85, 0.45, 0.65};
        shade.grid_lines = true;
        shade.grid_spacing_v = 1.0;
        render_parametric(surf, cam, shade, tile_w, tile_h, buf, tile_w, 0, w);
    }
    save_png((out_prefix + "curvature_grids.png").c_str(), w, h, buf, force);
    std::println("  saved curvature grids {}x{}", w, h);
    return 0;
}

// Complex-root-into-glow: ray_quadric_proximity() (spatium/geometry/
// ray_surface.hpp) already computes, for every ray that MISSES a
// Quadric, how close it came via the complex roots' imaginary part --
// an ordinary renderer just discards that (no hit = background) the
// same way ray_quadric() itself does. Render it instead: a glow around
// the sphere's silhouette, brighter the smaller miss_distance is --
// literally turning the discarded complex part of a "no intersection"
// root into light.
RenderStats render_glow_sphere(const Quadric<double>& q, const Camera& cam,
                                const ShadingConfig& shade, double glow_scale, int width,
                                int height, std::vector<std::uint8_t>& rgba) {
    rgba.assign(static_cast<std::size_t>(width) * height * 4, 0);

    auto fwd = Vec<double, 3>{(cam.target - cam.position).normalized()};
    auto right = Vec<double, 3>{fwd.cross(cam.up).normalized()};
    auto up    = Vec<double, 3>{right.cross(fwd).normalized()};
    double tan_half = std::tan(cam.fov_deg * std::numbers::pi / 360.0);
    double aspect = static_cast<double>(width) / height;
    auto light = Vec<double, 3>{shade.light_dir.normalized()};

    RenderStats stats{};
    stats.pixels = static_cast<std::size_t>(width) * height;

    auto render_rows = [&](int y0, int y1) -> std::size_t {
        std::size_t local_hits = 0;
        for (int y = y0; y < y1; ++y) {
            for (int x = 0; x < width; ++x) {
                std::uint8_t* pixel = &rgba[(y * width + x) * 4];
                bool any_hit = false;

                auto ray_color = [&](double nx, double ny) -> Vec<double, 3> {
                    auto dir = Vec<double, 3>{(fwd + right * nx + up * ny).normalized()};
                    auto r = unwrap(ray(cam.position, dir));
                    auto hits = ray_quadric(r, q);

                    if (!hits.empty()) {
                        any_hit = true;
                        auto& h = hits.front();
                        auto n = h.normal;
                        if (shade.two_sided && n.dot(dir) > 0.0) n = Vec<double, 3>{-n};
                        double diff = std::max(0.0, n.dot(light));
                        double shaded = shade.ambient + (1.0 - shade.ambient) * diff;
                        return Vec<double, 3>{shade.base_color * shaded * 255.0};
                    }
                    double v01 = 0.5 * (dir[2] + 1.0);
                    Vec<double, 3> sky{shade.sky_bottom + (shade.sky_top - shade.sky_bottom) * v01};
                    if (auto prox = ray_quadric_proximity(r, q)) {
                        double glow = std::min(3.0, glow_scale / (prox->miss_distance + 0.05));
                        sky = Vec<double, 3>{sky + shade.base_color * glow};
                    }
                    return Vec<double, 3>{sky * 255.0};
                };

                std::uint8_t rgb[3];
                spatium::render::supersample_pixel(x, y, width, height, tan_half, aspect,
                                                      ray_color, rgb);
                pixel[0] = rgb[0];
                pixel[1] = rgb[1];
                pixel[2] = rgb[2];
                pixel[3] = 255;
                if (any_hit) ++local_hits;
            }
        }
        return local_hits;
    };

    unsigned nthreads = std::max(1u, std::thread::hardware_concurrency());
    nthreads = std::min(nthreads, static_cast<unsigned>(std::max(1, height)));
    std::vector<std::size_t> hits_per_thread(nthreads, 0);
    int rows_per_thread = (height + static_cast<int>(nthreads) - 1) / static_cast<int>(nthreads);

    auto start = Clock::now();
    {
        std::vector<std::jthread> workers;
        workers.reserve(nthreads);
        for (unsigned t = 0; t < nthreads; ++t) {
            int y0 = static_cast<int>(t) * rows_per_thread;
            int y1 = std::min(height, y0 + rows_per_thread);
            if (y0 >= y1) continue;
            workers.emplace_back([&render_rows, &hits_per_thread, t, y0, y1] {
                hits_per_thread[t] = render_rows(y0, y1);
            });
        }
    }  // jthreads join here
    stats.elapsed_ms = std::chrono::duration<double, std::milli>(Clock::now() - start).count();
    for (auto h : hits_per_thread) stats.hits += h;
    return stats;
}

int render_glow(int width, int height, const std::string& out_prefix, bool force) {
    auto q = Quadric<double>::sphere(1.0);
    Camera cam{.position = {3.2, 2.4, 1.8}, .target = {0.0, 0.0, 0.0}};
    ShadingConfig shade;
    shade.base_color = {0.55, 0.80, 1.0};  // also the glow's tint
    std::vector<std::uint8_t> buf;
    auto s = render_glow_sphere(q, cam, shade, /*glow_scale=*/0.18, width, height, buf);
    print_stats("glow sphere", s, width, height);
    save_png((out_prefix + "glow_sphere.png").c_str(), width, height, buf, force);
    return 0;
}

void print_usage() {
    std::println("Usage: parametric_analytical_demo [--mode single|gallery|grids] [--width N]"
                 " [--height N] [--out-prefix PATH] [--force]");
    std::println("  --mode single (default)  emit individual PNGs");
    std::println("  --mode gallery           emit one 2x2 composite PNG");
    std::println("  --mode grids             emit Schwarzschild/wormhole curvature-grid PNGs");
    std::println("  --mode glow              emit a sphere with complex-root-proximity glow");
    std::println("  --width, --height        pixels per image / tile (default 960)");
    std::println("  --out-prefix             output filename prefix (default \"\")");
    std::println("  --force                  overwrite existing output files");
}

} // namespace

int main(int argc, char* argv[]) {
    std::string mode = "single";
    std::string out_prefix;
    int width = 960, height = 960;
    bool force = false;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if      (a == "--help" || a == "-h") { print_usage(); return 0; }
        else if (a == "--mode" && i + 1 < argc) mode = argv[++i];
        else if (a == "--width" && i + 1 < argc) width = std::atoi(argv[++i]);
        else if (a == "--height" && i + 1 < argc) height = std::atoi(argv[++i]);
        else if (a == "--out-prefix" && i + 1 < argc) out_prefix = argv[++i];
        else if (a == "--force") force = true;
        else { std::println("Unknown argument: {}", a); print_usage(); return 1; }
    }

    std::println("Analytical parametric render  mode={}  tile={}x{}", mode, width, height);

    if (mode == "single")  return render_single(width, height, out_prefix, force);
    if (mode == "gallery") return render_gallery(width, height, out_prefix, force);
    if (mode == "grids")   return render_grids(width, height, out_prefix, force);
    if (mode == "glow")    return render_glow(width, height, out_prefix, force);

    std::println("Unknown mode: {}", mode);
    print_usage();
    return 1;
}
