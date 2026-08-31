// Unified Vulkan viewer demo — multi-scene dispatcher.
//
// Scenes:
//   primitives (default) — unified scene of the library's core primitives
//                          (sphere, box, torus, cylinder, cone, ellipsoid,
//                          triangle) with BVH raycast visualization and a
//                          console benchmark vs brute-force / analytical
//                          ray_quadric.
//   torus      — Clifford torus S¹×S¹ ⊂ S³ stereographically projected
//                into R³, colored by x₄.
//   klein      — classic Klein bottle R⁴→R³ immersion, animated 4D
//                rotation separates / merges the self-intersection.

#include "io_helpers.hpp"

#include <spatium/spatium.hpp>
#include <spatium/spatial/bvh.hpp>
#include <spatium/mesh/primitives.hpp>
#include <spatium/mesh/operations.hpp>
#include <spatium/viewer/app.hpp>
#include <spatium/io/table.hpp>
#include <spatium/geometry/ray_surface.hpp>
#include <spatium/vendor/stb_image_write.h>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <print>
#include <random>
#include <string>
#include <vector>

#if SPATIUM_HAS_IMGUI
#include <imgui.h>
#endif

using namespace spatium;
using namespace spatium::literals;
using namespace spatium::geometry;
using namespace spatium::mesh;
using namespace spatium::spatial;
using namespace spatium::viewer;

namespace {

using Clock = std::chrono::steady_clock;

// ── Shared helpers ────────────────────────────────────────────

Mesh<Euclidean<3>> shifted(Mesh<Euclidean<3>> m, Vec3 offset) {
    for (auto& v : m.vertices) v = Vec3{v + offset};
    return m;
}

template<Surface S>
Mesh<Euclidean<3>> to_euclidean_mesh(const Mesh<S>& src) {
    Mesh<Euclidean<3>> out;
    out.vertices.reserve(src.vertices.size());
    for (const auto& v : src.vertices) out.vertices.push_back(Vec3{v});
    out.faces = src.faces;
    return out;
}

std::vector<Triangle3> triangles_of(const Mesh<Euclidean<3>>& m) {
    std::vector<Triangle3> out;
    out.reserve(m.face_count());
    for (auto [a, b, c] : m.triangles())
        out.push_back(Triangle3(a, b, c));
    return out;
}

// ── Scene: primitives + raycast ───────────────────────────────

struct Probe {
    Vec3 origin;
    std::vector<Ray<3>> rays;
};

Probe make_probe(std::size_t n, Vec3 origin, Vec3 aim, double spread) {
    Probe p{origin, {}};
    p.rays.reserve(n);

    auto forward = (aim - origin).normalized();
    Vec3 up{0, 0, 1};
    if (std::abs(forward.dot(up)) > 0.95) up = Vec3{0, 1, 0};
    auto right = forward.cross(up).normalized();
    up = right.cross(forward).normalized();

    std::mt19937 rng(1337);
    std::uniform_real_distribution<double> jitter(-spread, spread);

    for (std::size_t i = 0; i < n; ++i) {
        auto d = Vec3{forward + right * jitter(rng) + up * jitter(rng)}.normalized();
        if (auto r = ray(origin, d)) p.rays.push_back(*r);
    }
    return p;
}

int run_primitives(std::size_t ray_count, bool do_viewer, bool force) {
    auto sphere_mesh = shifted(uv_sphere_mesh(24, 16, 0.7),   Vec3{-2.5,  2.0, 0});
    auto box         = shifted(box_mesh(Vec3{0.5, 0.5, 0.5}), Vec3{ 0.0,  2.0, 0});
    auto cube        = shifted(box_mesh(Vec3{0.8, 0.2, 0.3}), Vec3{ 2.5,  2.0, 0});

    auto torus    = shifted(to_euclidean_mesh(parametric_mesh(make_torus(0.8, 0.3), 48, 24)),
                            Vec3{-2.5, -1.0, 0});
    auto cylinder = shifted(to_euclidean_mesh(parametric_mesh(make_cylinder(0.5, 1.5), 32, 8)),
                            Vec3{ 0.0, -1.0, -0.75});
    auto cone     = shifted(to_euclidean_mesh(parametric_mesh(make_cone(0.7, 1.5), 32, 8)),
                            Vec3{ 2.5, -1.0, -0.75});

    auto lod2 = subdivide(icosahedron(), Sphere<2>{}, 2);
    Mesh<Euclidean<3>> ellipsoid;
    ellipsoid.vertices.reserve(lod2.vertex_count());
    for (const auto& v : lod2.vertices)
        ellipsoid.vertices.push_back(Vec3{v[0] * 1.0, v[1] * 0.5, v[2] * 0.3});
    ellipsoid.faces = lod2.faces;
    ellipsoid = shifted(ellipsoid, Vec3{0.0, 0.0, 1.5});

    Mesh<Euclidean<3>> triangle_mesh;
    triangle_mesh.vertices = {Vec3{-1.0, 0.0, -1.8}, Vec3{1.0, 0.0, -1.8}, Vec3{0.0, 0.8, -1.2}};
    triangle_mesh.faces    = {{0, 1, 2}};

    struct Piece { std::string name; Mesh<Euclidean<3>> mesh; Vec4f color; };
    std::vector<Piece> pieces = {
        {"sphere",    sphere_mesh,   {0.45f, 0.75f, 1.00f, 0.90f}},
        {"box",       box,           {1.00f, 0.65f, 0.30f, 0.90f}},
        {"oriented",  cube,          {0.95f, 0.45f, 0.55f, 0.90f}},
        {"torus",     torus,         {0.55f, 0.90f, 0.60f, 0.90f}},
        {"cylinder",  cylinder,      {0.80f, 0.70f, 0.30f, 0.90f}},
        {"cone",      cone,          {0.90f, 0.55f, 0.80f, 0.90f}},
        {"ellipsoid", ellipsoid,     {0.55f, 0.65f, 0.95f, 0.90f}},
        {"triangle",  triangle_mesh, {0.85f, 0.85f, 0.40f, 0.95f}},
    };

    std::vector<Triangle3> tris;
    for (auto& p : pieces) {
        auto t = triangles_of(p.mesh);
        tris.insert(tris.end(), t.begin(), t.end());
    }

    auto build_start = Clock::now();
    auto bvh = BVH<Triangle3>::build(tris);
    auto build_us = std::chrono::duration<double, std::micro>(Clock::now() - build_start).count();

    auto probe = make_probe(ray_count, Vec3{6, 6, 5}, Vec3{0, 0, 0}, 0.35);

    auto bvh_start = Clock::now();
    std::size_t bvh_hits = 0;
    std::vector<Vec3> hit_points;
    hit_points.reserve(probe.rays.size());
    for (const auto& r : probe.rays) {
        if (auto h = bvh.ray_cast(r)) { ++bvh_hits; hit_points.push_back(h->point); }
    }
    auto bvh_us = std::chrono::duration<double, std::micro>(Clock::now() - bvh_start).count();

    auto brute_start = Clock::now();
    std::size_t brute_hits = 0;
    for (const auto& r : probe.rays) {
        double best_t = std::numeric_limits<double>::max();
        bool hit_any = false;
        for (auto& tri : tris) {
            if (auto it = intersect(r, tri)) {
                double t = (*it - r.origin).dot(r.direction);
                if (t >= 0 && t < best_t) { best_t = t; hit_any = true; }
            }
        }
        if (hit_any) ++brute_hits;
    }
    auto brute_us = std::chrono::duration<double, std::micro>(Clock::now() - brute_start).count();

    auto Q_sphere    = Quadric<double>::sphere(1.0);
    auto Q_ellipsoid = Quadric<double>::ellipsoid(1.0, 0.5, 0.3);
    auto Q_cylinder  = Quadric<double>::cylinder_z(0.5);

    auto quad_start = Clock::now();
    std::size_t quad_hits = 0;
    for (const auto& r : probe.rays) {
        if (!ray_quadric(r, Q_sphere).empty())    ++quad_hits;
        if (!ray_quadric(r, Q_ellipsoid).empty()) ++quad_hits;
        if (!ray_quadric(r, Q_cylinder).empty())  ++quad_hits;
    }
    auto quad_us = std::chrono::duration<double, std::micro>(Clock::now() - quad_start).count();

    // Analytical torus (quartic) — placed where the mesh torus sits in the scene
    Torus<double> scene_torus{
        .center = Vec3{-2.5, -1.0, 0.0},
        .axis = Vec3{0, 0, 1},
        .major_radius = 0.8,
        .minor_radius = 0.3
    };
    auto torus_start = Clock::now();
    std::size_t torus_hits = 0;
    for (const auto& r : probe.rays) {
        if (!ray_torus(r, scene_torus).empty()) ++torus_hits;
    }
    auto torus_us = std::chrono::duration<double, std::micro>(Clock::now() - torus_start).count();

    section("Scene");
    {
        Table tbl("Piece", "Verts  Faces");
        for (auto& p : pieces)
            tbl.row(p.name, std::format("{:>5} / {:>5}", p.mesh.vertex_count(), p.mesh.face_count()));
        tbl.row("TOTAL tris", std::format("{}", tris.size())).print();
    }

    section("Raycast");
    Table("Method", "Hits / Rays    Time")
        .row("BVH triangle",     std::format("{:>5} / {:>5}    {:>8.1f} us ({:.2f} ns/ray)",
                                             bvh_hits, ray_count, bvh_us,
                                             bvh_us * 1000.0 / static_cast<double>(ray_count)))
        .row("Brute triangle",   std::format("{:>5} / {:>5}    {:>8.1f} us ({:.2f} ns/ray)",
                                             brute_hits, ray_count, brute_us,
                                             brute_us * 1000.0 / static_cast<double>(ray_count)))
        .row("Analytic quadric", std::format("{:>5} / {:>5}    {:>8.1f} us ({:.2f} ns/ray, 3 quadrics)",
                                             quad_hits, ray_count * 3, quad_us,
                                             quad_us * 1000.0 / static_cast<double>(ray_count * 3)))
        .row("Analytic torus",   std::format("{:>5} / {:>5}    {:>8.1f} us ({:.2f} ns/ray, quartic)",
                                             torus_hits, ray_count, torus_us,
                                             torus_us * 1000.0 / static_cast<double>(ray_count)))
        .print();

    std::println("  BVH build: {:.1f} us  ({} nodes, {} triangles)", build_us,
                 bvh.node_count(), tris.size());
    std::println("  BVH speedup over brute: {:.1f}x", brute_us / bvh_us);

    // ── Analytical PNG preview ────────────────────────────────
    // CPU raytrace a clean torus with Lambert shading + ray_torus analytic normals.
    {
        constexpr int W = 512, H = 512;
        Torus<double> t{.major_radius = 1.0, .minor_radius = 0.3};

        Vec3 cam_pos{4.0, 4.0, 3.2};
        Vec3 target{0, 0, 0};
        Vec3 up{0, 0, 1};
        auto fwd = (target - cam_pos).normalized();
        auto right = fwd.cross(up).normalized();
        up = right.cross(fwd).normalized();

        double fov = 32.0 * std::numbers::pi / 180.0;
        double tan_half = std::tan(fov * 0.5);
        double aspect = static_cast<double>(W) / H;

        auto light = Vec3{1.0, -0.6, 1.4}.normalized();

        std::vector<std::uint8_t> rgba(W * H * 4);
        auto png_start = Clock::now();
        std::size_t pixel_hits = 0;
        for (int y = 0; y < H; ++y) {
            for (int x = 0; x < W; ++x) {
                double nx = (2.0 * (x + 0.5) / W - 1.0) * aspect * tan_half;
                double ny = (1.0 - 2.0 * (y + 0.5) / H) * tan_half;
                auto dir = (fwd + right * nx + up * ny).normalized();
                Ray<3> r{cam_pos, dir};

                auto hits = ray_torus(r, t);
                std::uint8_t* px = &rgba[(y * W + x) * 4];
                if (!hits.empty()) {
                    ++pixel_hits;
                    auto& h = hits.front();
                    double diff = std::max(0.0, h.normal.dot(light));
                    // Simple ambient + diffuse shading with a warm base
                    double base_r = 0.90, base_g = 0.55, base_b = 0.30;
                    double amb = 0.15;
                    double shade = amb + (1.0 - amb) * diff;
                    px[0] = static_cast<std::uint8_t>(std::clamp(base_r * shade, 0.0, 1.0) * 255);
                    px[1] = static_cast<std::uint8_t>(std::clamp(base_g * shade, 0.0, 1.0) * 255);
                    px[2] = static_cast<std::uint8_t>(std::clamp(base_b * shade, 0.0, 1.0) * 255);
                    px[3] = 255;
                } else {
                    // Gradient sky
                    double v = 0.5 * (dir[2] + 1.0);
                    px[0] = static_cast<std::uint8_t>((0.15 + 0.10 * v) * 255);
                    px[1] = static_cast<std::uint8_t>((0.18 + 0.15 * v) * 255);
                    px[2] = static_cast<std::uint8_t>((0.25 + 0.25 * v) * 255);
                    px[3] = 255;
                }
            }
        }
        auto png_us = std::chrono::duration<double, std::micro>(Clock::now() - png_start).count();

        const char* png_path = "torus_analytical.png";
        if (spatium::examples::confirm_overwrite(png_path, force)) {
            stbi_write_png(png_path, W, H, 4, rgba.data(), W * 4);
            std::println("\n  Analytical torus render: {} ({}x{}, {} hits / {} rays, {:.1f} ms, {:.1f} ns/ray)",
                         png_path, W, H, pixel_hits, W * H, png_us / 1000.0,
                         png_us * 1000.0 / static_cast<double>(W * H));
        }
    }

    if (!do_viewer) return 0;

    App app("Spatium \u2014 Primitives + Raycast", 1280, 720);
    for (auto& p : pieces) app.add_mesh(p.mesh, p.color);

    app.add_mesh(shifted(box_mesh(Vec3{0.1, 0.1, 0.1}), probe.origin),
                 Vec4f{1.0f, 0.95f, 0.25f, 1.0f});

    PointCloudData hits;
    hits.positions.reserve(hit_points.size() * 3);
    for (const auto& h : hit_points) {
        hits.positions.push_back(static_cast<float>(h[0]));
        hits.positions.push_back(static_cast<float>(h[1]));
        hits.positions.push_back(static_cast<float>(h[2]));
    }
    hits.point_count = static_cast<uint32_t>(hit_points.size());
    app.add_point_cloud(std::move(hits), Vec4f{1.0f, 0.25f, 0.25f, 1.0f});

    app.camera.distance = 14.0f;
    app.point_size_ = 3.5f;

    std::println("\nControls: drag=orbit, right-drag=pan, scroll=zoom, W=wire, P/O=point size, S=screenshot, Q=quit");
    app.run();
    return 0;
}

// ── Scene: Clifford torus (S³ → R³) with 4D rotation ───────────

// Clifford torus S¹×S¹ ⊂ S³ at radius 1/√2. Parameterized by (θ, φ):
//   x1 = r cos θ, x2 = r sin θ, x3 = r cos φ, x4 = r sin φ
// A 4D rotation in the chosen plane is applied before stereographic
// projection from the north pole (x4 = 1) into R³.
namespace clifford_torus {

constexpr double r4 = 0.70710678118654752;  // 1/√2

enum class Plane { X1X3, X2X4, X3X4 };

Vec3 eval(double theta, double phi, double alpha, Plane plane) {
    double x1 = r4 * std::cos(theta);
    double x2 = r4 * std::sin(theta);
    double x3 = r4 * std::cos(phi);
    double x4 = r4 * std::sin(phi);

    double c = std::cos(alpha), s = std::sin(alpha);
    switch (plane) {
        case Plane::X1X3: { double x1r = x1 * c - x3 * s; double x3r = x1 * s + x3 * c; x1 = x1r; x3 = x3r; break; }
        case Plane::X2X4: { double x2r = x2 * c - x4 * s; double x4r = x2 * s + x4 * c; x2 = x2r; x4 = x4r; break; }
        case Plane::X3X4: { double x3r = x3 * c - x4 * s; double x4r = x3 * s + x4 * c; x3 = x3r; x4 = x4r; break; }
    }

    // Stereographic projection from x4 = 1
    double d = 1.0 / std::max(1.0 - x4, 1e-4);
    return {x1 * d, x2 * d, x3 * d};
}

} // namespace clifford_torus

int run_torus(std::size_t resolution, bool do_viewer) {
    using clifford_torus::Plane;

    int res            = static_cast<int>(resolution);
    double alpha       = 0.0;
    double speed       = 0.35;
    bool animating     = true;
    int plane_choice   = static_cast<int>(Plane::X3X4);

    auto build_surface = [&](double a) {
        auto fn = [=](double theta, double phi) {
            return clifford_torus::eval(theta, phi, a, static_cast<Plane>(plane_choice));
        };
        return parametric(fn, periodic(0.0, 2_pi, 0.0, 2_pi));
    };

    std::println("Clifford torus: S^1 x S^1 in S^3 -> R^3 (stereographic)");
    std::println("Resolution: {} x {}", res, res);

    auto mesh = parametric_mesh(build_surface(alpha), res, res);
    std::println("Mesh: {} vertices, {} faces", mesh.vertex_count(), mesh.face_count());

    double area = 0;
    for (auto [a, b, c] : mesh.triangles())
        area += cross(Vec3{b - a}, Vec3{c - a}).norm() * 0.5;
    std::println("Surface area at α=0: {:.4f} (intrinsic: 2*pi^2 = {:.4f})", area, 2_pi * std::numbers::pi);

    if (!do_viewer) return 0;

    App app("Spatium \u2014 Clifford Torus (S\u00b3 \u2192 R\u00b3)", 1280, 720);
    app.add_mesh(mesh, Vec4f{0.45f, 0.75f, 0.95f, 0.9f});
    app.fit_to(mesh);

    auto last_time = Clock::now();
    int current_res = res;

    app.set_key_callback([&](int key, int action, int /*mods*/) {
        if (action == 0) return;  // GLFW_RELEASE
        if (key == ' ') animating = !animating;
        if (key == '-' || key == '_') speed = std::max(speed - 0.1, 0.0);
        if (key == '=' || key == '+') speed += 0.1;
        if (key == 'R' || key == 'r') alpha = 0.0;
    });

    app.set_frame_callback([&]() {
        auto now = Clock::now();
        double dt = std::chrono::duration<double>(now - last_time).count();
        last_time = now;

        bool dirty = false;
        if (animating) { alpha += speed * dt; dirty = true; }
        if (current_res != res) { current_res = res; dirty = true; }
        if (dirty) {
            auto m = parametric_mesh(build_surface(alpha), current_res, current_res);
            app.update_mesh_vertices(0, mesh_to_render_data(m));
        }
    });

#if SPATIUM_HAS_IMGUI
    app.enable_imgui();
    app.set_gui_callback([&]() {
        ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(320, 0), ImGuiCond_FirstUseEver);
        ImGui::Begin("Clifford Torus");

        ImGui::Text("S\u00b9\u00d7S\u00b9 \u2282 S\u00b3 \u2192 R\u00b3");
        ImGui::Separator();

        ImGui::Checkbox("Animate", &animating);
        ImGui::SameLine();
        if (ImGui::Button("Reset \u03b1")) alpha = 0.0;

        ImGui::SliderFloat("Speed", (float*)&speed, 0.0f, 2.0f, "%.2f rad/s");

        float af = static_cast<float>(alpha);
        if (ImGui::SliderFloat("\u03b1 (manual)", &af, 0.0f, 2.0f * 3.14159265f, "%.3f rad")) {
            alpha = af;
            animating = false;
        }

        ImGui::Separator();
        ImGui::Text("Rotation plane in R\u2074");
        ImGui::RadioButton("x\u2081-x\u2083", &plane_choice, 0); ImGui::SameLine();
        ImGui::RadioButton("x\u2082-x\u2084", &plane_choice, 1); ImGui::SameLine();
        ImGui::RadioButton("x\u2083-x\u2084", &plane_choice, 2);

        ImGui::Separator();
        if (ImGui::SliderInt("Resolution", &res, 16, 192)) {
            if (res < 8) res = 8;
        }

        ImGui::Separator();
        ImGui::Text("Keys: SPACE pause | +/- speed | R reset | W wire");
        ImGui::End();
    });
#endif

    std::println("\nControls: SPACE pause, +/- speed, R reset \u03b1, W wire, drag=orbit, scroll=zoom");
#if SPATIUM_HAS_IMGUI
    std::println("ImGui panel: rotation plane, manual \u03b1 slider, live resolution.");
#endif

    app.run();
    return 0;
}

// ── Scene: Klein bottle (R⁴ → R³), animated ───────────────────

struct KleinBottle {
    // The two branches must agree not just in value at u=pi (they did) but
    // in how the v-parametrized cross-section is embedded, or the mesh
    // twists right at the seam where the neck re-enters the body: the
    // u>=pi branch's y dropped the r*sin(u)*cos(v) term the u<pi branch
    // has, so the cross-section's v-shape flipped from "varies with v in
    // both x and y" to "only in x" discontinuously at the seam, even
    // though position itself matched there. Restoring the term makes y a
    // single smooth expression across both branches (sin(u) -> 0 at u=pi
    // from both sides regardless, so this doesn't change position at the
    // seam -- only the shape approaching it).
    Vec3 eval_r3(double u, double v) const {
        double r = 4.0 * (1.0 - std::cos(u) / 2.0);
        double x, y;
        if (u < std::numbers::pi) {
            x = 6.0 * std::cos(u) * (1.0 + std::sin(u)) + r * std::cos(u) * std::cos(v);
        } else {
            x = 6.0 * std::cos(u) * (1.0 + std::sin(u)) + r * std::cos(v + std::numbers::pi);
        }
        y = 16.0 * std::sin(u) + r * std::sin(u) * std::cos(v);
        return {x, y, r * std::sin(v)};
    }

    double eval_w(double u, double v) const {
        double r = 4.0 * (1.0 - std::cos(u) / 2.0);
        return r * std::sin(v) * std::sin(u / 2.0);
    }

    Vec3 project(double u, double v, double alpha) const {
        auto p = eval_r3(u, v);
        double w = eval_w(u, v);

        double z = p[2] * std::cos(alpha) - w * std::sin(alpha);
        double w_rot = p[2] * std::sin(alpha) + w * std::cos(alpha);
        double s = 0.06 / (1.0 - w_rot * 0.015);

        return {p[0] * s, p[1] * s, z * s};
    }
};

int run_klein(std::size_t resolution, bool do_viewer) {
    KleinBottle klein;
    double alpha = 0.0;

    std::println("Klein bottle: classic immersion R^4 -> R^3");
    std::println("Resolution: {} x {}", resolution, resolution);

    auto surface = parametric([&](double u, double v) { return klein.project(u, v, alpha); },
                              periodic(0.0, 2_pi, 0.0, 2_pi));
    auto mesh = parametric_mesh(surface, resolution, resolution);
    std::println("Mesh: {} vertices, {} faces", mesh.vertex_count(), mesh.face_count());

    if (!do_viewer) return 0;

    App app("Spatium \u2014 Klein Bottle (4D \u2192 3D)", 1280, 720);
    app.add_mesh(mesh, Vec4f{0.3f, 0.75f, 0.85f, 0.9f});
    app.fit_to(mesh);

    double rotation_speed = 0.3;
    bool animating = true;
    auto last_time = Clock::now();

    app.set_key_callback([&](int key, int /*action*/, int /*mods*/) {
        if (key == ' ') animating = !animating;
        if (key == '-' || key == '_') rotation_speed = std::max(rotation_speed - 0.1, 0.0);
        if (key == '=' || key == '+') rotation_speed += 0.1;
    });

    app.set_frame_callback([&]() {
        auto now = Clock::now();
        double dt = std::chrono::duration<double>(now - last_time).count();
        last_time = now;

        if (!animating) return;
        alpha += rotation_speed * dt;

        auto surf = parametric([&](double u, double v) { return klein.project(u, v, alpha); },
                               periodic(0.0, 2_pi, 0.0, 2_pi));
        auto m = parametric_mesh(surf, resolution, resolution);
        app.update_mesh_vertices(0, mesh_to_render_data(m));
    });

    std::println("\nAnimated 4D rotation — self-intersection separates/merges");
    std::println("Controls: drag=orbit, right-drag=pan, scroll=zoom");
    std::println("  SPACE=pause, -/+=speed, W=wireframe, S=screenshot, Q=quit");
    app.run();
    return 0;
}

// ── Dispatch ─────────────────────────────────────────────────

void print_usage() {
    std::println("Usage: primitives_demo [--scene <name>] [options]");
    std::println("Scenes:");
    std::println("  primitives (default) — unified primitives + BVH raycast visualization");
    std::println("  torus                — Clifford torus S\u00b3 \u2192 R\u00b3");
    std::println("  klein                — Klein bottle R\u2074 \u2192 R\u00b3, animated");
    std::println("Options:");
    std::println("  --rays N         primitives scene: ray count (default 1024)");
    std::println("  --resolution N   torus/klein scene: UV grid resolution (default 64/80)");
    std::println("  --no-viewer      skip Vulkan window, console output only");
    std::println("  --force          overwrite existing output files");
    std::println("  --help, -h       this message");
}

} // namespace

int main(int argc, char* argv[]) {
    std::string scene = "primitives";
    std::size_t ray_count = 1024;
    std::size_t resolution = 0;  // 0 = scene default
    bool do_viewer = true;
    bool force = false;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--help" || a == "-h") { print_usage(); return 0; }
        else if (a == "--no-viewer") do_viewer = false;
        else if (a == "--force") force = true;
        else if (a == "--scene" && i + 1 < argc) scene = argv[++i];
        else if (a == "--rays" && i + 1 < argc) ray_count = std::stoul(argv[++i]);
        else if (a == "--resolution" && i + 1 < argc) resolution = std::stoul(argv[++i]);
        else { std::println("Unknown argument: {}", a); print_usage(); return 1; }
    }

    if (scene == "primitives")
        return run_primitives(ray_count, do_viewer, force);
    if (scene == "torus")
        return run_torus(resolution == 0 ? 64 : resolution, do_viewer);
    if (scene == "klein")
        return run_klein(resolution == 0 ? 80 : resolution, do_viewer);

    std::println("Unknown scene: {}", scene);
    print_usage();
    return 1;
}
