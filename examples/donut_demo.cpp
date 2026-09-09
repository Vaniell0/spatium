// Getting started with Spatium, in the tradition of every 3D beginner's
// first project: a donut. Blender's own famously-loved intro tutorial
// opens with the default cube -- delete it, then build the donut from
// nothing. This is the same idea for spatium::io::build's declarative
// scene DSL -- each step below is one line, and there is no for-loop
// anywhere in this file. See docs/getting-started-dsl.md for the same
// steps as prose.
//
// Step 0: the default cube, then delete it -- here, "delete" means
//         explode: every vertex's displacement is a closed-form function
//         of (its own position, time), built from spatium::algebra::
//         PerlinNoise -- no per-frame simulation state, no velocity
//         accumulated across steps. At t=0 it's a stationary cube; ask
//         for a later t and the same formula places the fragments
//         further out. That's the point of keeping motion declarative:
//         "explode outward with turbulence" is one expression, evaluated
//         wherever you need it, not a loop you have to run.
// Step 1: the dough is a torus -- a real spatium::ParametricSurface,
//         not a mesh.
// Step 2: the icing is an *offset* of the dough's own surface --
//         icing(u,v) = dough(u,v) + thickness * dough.normal_at(u,v) --
//         a new analytic surface, composed from a function, not a
//         separate shape draped/projected onto the first (that was
//         tried and was redundant: a surface built to be projected onto
//         another almost-identical surface teaches nothing offset()
//         doesn't already do directly).
// Step 3: the sprinkles are small cylinders, scattered across the icing
//         by spatium::sample_surface_uniform -- rejection sampling
//         weighted by the icing's own first-fundamental-form area
//         element, so coverage is even *by surface area*, not by
//         parameter or texture noise. No mesh exists anywhere in this
//         pipeline until the very last step, and only because a
//         renderer needs triangles -- the scene itself (`scene`, a
//         build::Trace) stays a flat, inspectable record of these three
//         operations the whole time; the loop over `scene.node(i)`
//         below is walking that record, not building geometry.
//
// Two output modes:
//   (default)     print the trace and materialized stats to the console
//   --photo PATH  also CPU-raytrace it to a PNG (BVH<Triangle3>, same
//                 engine as every other offline demo in this tree --
//                 tessellation happens here, only for display, not as
//                 part of the scene's own representation)

#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "io_helpers.hpp"

#include <spatium/algebra/noise.hpp>
#include <spatium/geometry/triangle.hpp>
#include <spatium/io/build.hpp>
#include <spatium/render/camera.hpp>
#include <spatium/render/parallel_for_rows.hpp>
#include <spatium/render/supersample.hpp>
#include <spatium/render/write_image.hpp>
#include <spatium/spatial/bvh.hpp>

#include <spatium/mesh/primitives.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <fstream>
#include <functional>
#include <numbers>
#include <random>
#include <utility>
#include <print>
#include <string>
#include <string_view>
#include <vector>

using namespace spatium;
namespace bd = spatium::io::build;
using spatium::io::Material;
using spatium::geometry::Ray;
using spatium::geometry::Triangle3;
using spatium::spatial::BVH;
using spatium::render::Camera;
using spatium::render::make_camera_basis;
using spatium::render::camera_ray_dir;
using spatium::render::parallel_for_rows;
using spatium::render::supersample_pixel;
using spatium::render::write_png_rgb;

namespace {

// Build-up timeline, in the same "t" materialize() already takes --
// shared by the dust particles' motion and the donut's grow-in, so
// there's one place that defines "when does what happen."
constexpr double T_EXPLODE = 0.9;   // cube -> dust, flying outward
constexpr double T_CONVERGE = 2.0;  // dust flies from its burst positions to the BOOM letterforms
constexpr double T_HOLD = 2.5;      // BOOM holds
constexpr double T_DISSOLVE = 3.2;  // dust shrinks to nothing
constexpr double T_DONUT_START = 2.9, T_DONUT_END = 3.9; // donut grows in, overlapping the dissolve
constexpr double T_BUILD_END = T_DONUT_END;

double smoothstep01(double t) {
    double e = std::clamp(t, 0.0, 1.0);
    return e * e * (3.0 - 2.0 * e);
}

// One dust particle's *core position* (before adding its own tiny local
// offset -- see particle_motion() below) at time t: annihilation, not a
// flight to a fixed point. Every particle swirls the whole time (three
// independent PerlinNoise channels sampled along its own path, not a
// literal orbit formula, but reads as turbulent wandering) growing out
// of the burst, and every particle is pulled toward its own nearby
// point on the BOOM letterforms as t moves through the converge window
// -- but by how much varies per particle (`pull_strength`, in [0,1]):
// near 1 and it lands and holds there; near 0 and the pull is just
// enough to bend its swirl in close past the letterform on the way by,
// not enough to stop it -- flying past *because* of the shape, not
// independently of it, which is also why they end up close enough to
// pick up color (dust_color below) before drifting off again. Only
// strongly-pulled particles get the extra "hits the surface and
// spreads" splat jitter near arrival -- a near-miss brushing past
// doesn't spread on impact, mid particles get a partial jitter.
Vec<double, 3> dust_core(double t, const Vec<double, 3>& burst_dir, double burst_dist,
                          const Vec<double, 3>& target, double pull_strength,
                          double swirl_seed, const algebra::PerlinNoise& swirl) {
    double explode_e = smoothstep01(t / T_EXPLODE);
    // A floor on the radius, not 0 at explode_e=0: 220 specks literally
    // coincident at the origin for the first several frames reproducibly
    // stalled BVH::build() for whole minutes (independent of the dust
    // shape -- tried both a UV-sphere and an icosahedron, same stall) --
    // this keeps them spread apart from frame 0 on, avoiding whatever
    // degenerate case that many near-identical bounding boxes hits.
    double radius = burst_dist * (0.08 + 0.92 * explode_e);
    double freq = 0.6 + 0.5 * std::fmod(swirl_seed, 1.0);
    double phase = swirl_seed * 11.0;
    Vec<double, 3> turb{
        swirl(burst_dir[0] * 3.0 + phase, burst_dir[1] * 3.0, t * freq),
        swirl(burst_dir[0] * 3.0 + phase + 7.0, burst_dir[1] * 3.0, t * freq),
        swirl(burst_dir[0] * 3.0 + phase, burst_dir[1] * 3.0 + 7.0, t * freq)};
    Vec<double, 3> swirl_pos = Vec<double, 3>{burst_dir * radius + turb * (0.5 + 0.6 * radius)};

    if (t < T_EXPLODE) return swirl_pos;

    double pull = smoothstep01((t - T_EXPLODE) / (T_CONVERGE - T_EXPLODE)) * pull_strength;
    Vec<double, 3> approach = Vec<double, 3>{swirl_pos * (1.0 - pull) + target * pull};

    // Splat: a jitter burst that peaks as the particle arrives (pull
    // near 1) and settles back down through the hold phase, not before.
    // Scales with pull_strength too -- a brush-past shouldn't spread.
    double arriving = pull * pull;
    double settle = t < T_HOLD ? 1.0 : 1.0 - smoothstep01((t - T_HOLD) / (T_DISSOLVE - T_HOLD + 0.001));
    Vec<double, 3> splat{
        swirl(target[0] * 5.0 + phase, target[1] * 5.0, t * 4.0),
        swirl(target[0] * 5.0 + phase + 3.0, target[1] * 5.0, t * 4.0),
        swirl(target[0] * 5.0 + phase, target[1] * 5.0 + 3.0, t * 4.0)};
    return Vec<double, 3>{approach + splat * (0.22 * arriving * settle * pull_strength)};
}

Vec<double, 3> particle_motion(const Vec<double, 3>& local_p, double t,
                                const Vec<double, 3>& burst_dir, double burst_dist,
                                const Vec<double, 3>& target, double pull_strength, double swirl_seed,
                                const algebra::PerlinNoise& swirl) {
    Vec<double, 3> pos = dust_core(t, burst_dir, burst_dist, target, pull_strength, swirl_seed, swirl);
    double scale = 1.0;
    if (t >= T_HOLD) scale = 1.0 - smoothstep01((t - T_HOLD) / (T_DISSOLVE - T_HOLD));
    return Vec<double, 3>{pos + local_p * scale};
}

std::vector<std::pair<double, double>> load_boom_points(const std::string& path) {
    std::vector<std::pair<double, double>> pts;
    std::ifstream in(path);
    double u, v;
    while (in >> u >> v) pts.emplace_back(u, v);
    return pts;
}

const char* kind_name(bd::Kind k) {
    switch (k) {
        case bd::Kind::Space:   return "Space";
        case bd::Kind::Offset:  return "Offset";
        case bd::Kind::Scatter: return "Scatter";
        case bd::Kind::Compose: return "Compose";
        case bd::Kind::Literal: return "Literal";
    }
    return "?";
}

// Icing guide: the same torus formula make_torus() uses, restricted to a
// v-band around the tube's "top" (v=pi/2 sits opposite the hole -- see
// spaces/parametric.hpp's make_torus) -- so offset() only builds icing
// over the dough's upper cap, not the whole ring. Not promoted to the
// library: the exact band is donut-specific set dressing, not a general
// primitive (unlike offset_surface() itself, which this still goes
// through directly -- this is the base being offset, not a second
// surface projected onto the first).
ParametricSurface<double> torus_cap(double major_r, double minor_r) {
    constexpr double pi = std::numbers::pi;
    return ParametricSurface<double>(
        [=](double u, double v) -> Vec<double, 3> {
            return {
                (major_r + minor_r * std::cos(v)) * std::cos(u),
                (major_r + minor_r * std::cos(v)) * std::sin(u),
                minor_r * std::sin(v)};
        },
        {0.0, 2.0 * pi, pi * 0.02, pi * 0.98}, // wide band; the noisy thickness
                                               // falloff below does the actual
                                               // edge shaping, not this domain
        /*periodic_u=*/true, /*periodic_v=*/false);
}

// A tiny round dust speck. mesh::uv_sphere_mesh() was the first attempt
// here and reproducibly stalled BVH::build() for whole minutes on
// specific frames of the dust cloud -- its poles are `slices` triangles
// all sharing one exact vertex (near-zero-area at low slice counts),
// and enough of those degenerate triangles clustered near the origin
// early in the explosion looks like what a naive SAH split can pathologically
// stall on. mesh::icosahedron() has no poles at all (12 uniform
// vertices, 20 uniform-ish triangles, a regular solid) -- same round
// silhouette at this size, none of the degenerate geometry.
mesh::Mesh<Euclidean<3, double>> dust_speck(double radius) {
    auto ico = mesh::icosahedron(Sphere<2, double>{radius});
    mesh::Mesh<Euclidean<3, double>> m;
    m.vertices.assign(ico.vertices.begin(), ico.vertices.end());
    m.faces = ico.faces;
    return m;
}

// A solid, capped cylinder -- spaces/parametric.hpp's make_cylinder() is
// the *lateral surface only* (an open tube, no end disks -- fine for a
// ray-traced pipe, wrong for a sprinkle that should read as solid).
// Built directly as a mesh, not a ParametricSurface: a cylinder's two
// flat end caps don't fit one (u,v) formula seamlessly. Centered at the
// origin along z -- scatter() maps local z to the target's normal, so
// half the cylinder sits below the placement point (embedded in the
// surface) and half above, by construction.
mesh::Mesh<Euclidean<3, double>> solid_cylinder(double radius, double height, int sides = 8) {
    mesh::Mesh<Euclidean<3, double>> m;
    constexpr double pi = std::numbers::pi;
    double half = height * 0.5;
    std::uint32_t top_center = 0, bottom_center = 1;
    m.vertices.push_back({0.0, 0.0, half});
    m.vertices.push_back({0.0, 0.0, -half});
    std::vector<std::uint32_t> top_ring, bottom_ring;
    for (int i = 0; i < sides; ++i) {
        double a = 2.0 * pi * static_cast<double>(i) / static_cast<double>(sides);
        double x = radius * std::cos(a), y = radius * std::sin(a);
        top_ring.push_back(static_cast<std::uint32_t>(m.vertices.size()));
        m.vertices.push_back({x, y, half});
        bottom_ring.push_back(static_cast<std::uint32_t>(m.vertices.size()));
        m.vertices.push_back({x, y, -half});
    }
    for (int i = 0; i < sides; ++i) {
        int j = (i + 1) % sides;
        m.faces.push_back({top_center, top_ring[static_cast<std::size_t>(j)], top_ring[static_cast<std::size_t>(i)]});
        m.faces.push_back({bottom_center, bottom_ring[static_cast<std::size_t>(i)], bottom_ring[static_cast<std::size_t>(j)]});
        m.faces.push_back({top_ring[static_cast<std::size_t>(i)], top_ring[static_cast<std::size_t>(j)], bottom_ring[static_cast<std::size_t>(i)]});
        m.faces.push_back({top_ring[static_cast<std::size_t>(j)], bottom_ring[static_cast<std::size_t>(j)], bottom_ring[static_cast<std::size_t>(i)]});
    }
    return m;
}

// Blender-viewport-style coordinate axes (X red, Y green, Z blue, the
// same convention Blender's own gizmo uses) -- three thin solid
// cylinders through the origin. Shown during the build-up (cube ->
// donut), hidden once the donut is finished and the camera starts to
// orbit, so the final hero shot stays clean.
std::vector<bd::Placed<double>> axes_placed(double length = 3.5, double radius = 0.012) {
    std::vector<bd::Placed<double>> out;
    auto arm = [&](const Vec<double, 3>& color) {
        auto m = solid_cylinder(radius, length, 10);
        return bd::Placed<double>{m, Material<double>{.base_color = color, .roughness = 0.6}};
    };
    auto x_arm = arm({0.85, 0.15, 0.15});
    for (auto& v : x_arm.mesh.vertices) v = Vec<double, 3>{v[2], v[0], v[1]}; // z-axis cylinder -> x
    auto y_arm = arm({0.15, 0.75, 0.15});
    for (auto& v : y_arm.mesh.vertices) v = Vec<double, 3>{v[0], v[2], v[1]}; // z-axis cylinder -> y
    auto z_arm = arm({0.15, 0.25, 0.85});
    out.push_back(x_arm);
    out.push_back(y_arm);
    out.push_back(z_arm);
    return out;
}

// Smooth per-vertex normals (area-weighted face-normal average) -- flat
// per-triangle normals on a coarse tessellation read as faceted plastic;
// this is the standard fix, needs nothing beyond the mesh's own topology.
std::vector<Vec<double, 3>> smooth_normals(const mesh::Mesh<Euclidean<3, double>>& m) {
    std::vector<Vec<double, 3>> n(m.vertex_count(), Vec<double, 3>{0, 0, 0});
    for (const auto& f : m.faces) {
        auto& a = m.vertices[f[0]];
        auto& b = m.vertices[f[1]];
        auto& c = m.vertices[f[2]];
        Vec<double, 3> face_n{Vec<double, 3>{b - a}.cross(Vec<double, 3>{c - a})}; // area-weighted (unnormalized)
        n[f[0]] = Vec<double, 3>{n[f[0]] + face_n};
        n[f[1]] = Vec<double, 3>{n[f[1]] + face_n};
        n[f[2]] = Vec<double, 3>{n[f[2]] + face_n};
    }
    for (auto& v : n) {
        double len = v.norm();
        v = len > 1e-12 ? Vec<double, 3>{v / len} : Vec<double, 3>{0, 0, 1};
    }
    return n;
}

struct Prim {
    Triangle3 tri;
    std::array<Vec<double, 3>, 3> vertex_normals;
    Vec<double, 3> color;
    double roughness;
};

std::vector<std::uint8_t> render_frame(const std::vector<bd::Placed<double>>& scene, const Camera<double>& cam,
                                        int W, int H) {
    std::vector<Triangle3> tris;
    std::vector<Prim> prims;
    for (const auto& obj : scene) {
        auto vn = smooth_normals(obj.mesh);
        for (const auto& f : obj.mesh.faces) {
            Triangle3 t(obj.mesh.vertices[f[0]], obj.mesh.vertices[f[1]], obj.mesh.vertices[f[2]]);
            tris.push_back(t);
            prims.push_back(Prim{t, {vn[f[0]], vn[f[1]], vn[f[2]]}, obj.material.base_color, obj.material.roughness});
        }
    }
    auto bvh = BVH<Triangle3>::build(tris);

    const Vec<double, 3> background{0.55, 0.75, 0.92}; // plain light blue, no starfield
    const auto basis = make_camera_basis(cam);
    const Vec<double, 3> light = Vec<double, 3>{Vec<double, 3>{0.55, -0.4, 1.0}.normalized()};
    constexpr int MAX_DEPTH = 3;

    // trace_ray returns linear [0,1] color throughout -- the single
    // conversion back to [0,255] happens exactly once, at the very end
    // of ray_color below, so reflection blending stays on one scale.
    std::function<Vec<double, 3>(const Ray<3, double>&, int)> trace_ray =
        [&](const Ray<3, double>& ray, int depth) -> Vec<double, 3> {
        auto hit = bvh.ray_cast(ray);
        if (!hit) return background;

        const auto& prim = prims[hit->index];
        // Interpolate the smooth vertex normals via the hit's own
        // barycentric weights, not the triangle's single flat normal.
        double w0 = double{1} - hit->u - hit->v, w1 = hit->u, w2 = hit->v;
        Vec<double, 3> n{Vec<double, 3>{prim.vertex_normals[0] * w0 + prim.vertex_normals[1] * w1 +
                                         prim.vertex_normals[2] * w2}
                              .normalized()};
        if (n.dot(ray.direction) > 0.0) n = Vec<double, 3>{-n};

        double diff = std::max(0.0, n.dot(light));
        double roughness = std::clamp(prim.roughness, 0.0, 1.0);

        Vec<double, 3> view = Vec<double, 3>{-ray.direction};
        Vec<double, 3> half = Vec<double, 3>{Vec<double, 3>{view + light}.normalized()};
        double shininess = 4.0 + 90.0 * (1.0 - roughness); // rough: broad/dull, glossy: tight/bright
        double spec = std::pow(std::max(0.0, n.dot(half)), shininess) * (1.0 - roughness) * 0.6;

        Vec<double, 3> local = Vec<double, 3>{prim.color * (0.22 + 0.78 * diff)};
        Vec<double, 3> color = Vec<double, 3>{local + Vec<double, 3>{1.0, 1.0, 1.0} * spec};

        // A real reflected ray for glossy materials, not just a
        // highlight -- what actually earns the word "raytracing" here.
        if (roughness < 0.6 && depth < MAX_DEPTH) {
            Vec<double, 3> refl_dir = Vec<double, 3>{ray.direction - n * (2.0 * ray.direction.dot(n))};
            Vec<double, 3> refl_origin = Vec<double, 3>{hit->point + n * 1e-4};
            Vec<double, 3> refl_color = trace_ray(Ray<3, double>{refl_origin, refl_dir}, depth + 1);
            double reflectivity = (1.0 - roughness) * 0.28;
            color = Vec<double, 3>{color * (1.0 - reflectivity) + refl_color * reflectivity};
        }
        return color;
    };

    std::vector<std::uint8_t> img(3 * static_cast<std::size_t>(W) * H, 0);
    parallel_for_rows(H, [&](int y) {
        for (int x = 0; x < W; ++x) {
            std::uint8_t* px = &img[3 * (static_cast<std::size_t>(y) * W + x)];
            auto ray_color = [&](double sx, double sy) -> Vec<double, 3> {
                Vec<double, 3> dir = camera_ray_dir(basis, sx, sy);
                Vec<double, 3> c = trace_ray(Ray<3, double>{cam.position, dir}, 0);
                return Vec<double, 3>{
                    std::clamp(c[0], 0.0, 1.0) * 255.0,
                    std::clamp(c[1], 0.0, 1.0) * 255.0,
                    std::clamp(c[2], 0.0, 1.0) * 255.0};
            };
            supersample_pixel(x, y, W, H, basis.tan_half,
                              static_cast<double>(W) / H, ray_color, px);
        }
    });
    return img;
}

Camera<double> hero_camera() {
    return {.position = {5.0, -4.2, 3.6}, .target = {0.0, 0.0, 0.0}, .up = {0.0, 0.0, 1.0}, .fov_deg = 38.0};
}

void render_photo(const std::vector<bd::Placed<double>>& scene, const std::string& out_path, bool force) {
    constexpr int W = 960, H = 720;
    auto img = render_frame(scene, hero_camera(), W, H);
    if (spatium::examples::confirm_overwrite(out_path, force))
        write_png_rgb(out_path, W, H, img);
    std::println("  -> {}", out_path);
}

// The lesson, as a video: build-up (cube deletes, donut grows in its
// place, axes visible like a viewport) then a clean orbit once it's
// finished (axes hidden). Two phases, same materialize() call at
// different t -- no separate animation system, the trace already
// describes motion as a function of time.
void render_video(const bd::Trace<double>& scene, std::size_t lesson_idx, const std::string& dir,
                   bool force, int build_frames, int orbit_frames) {
    namespace fs = std::filesystem;
    fs::create_directories(dir);
    constexpr int W = 960, H = 720;
    Camera<double> cam = hero_camera();
    double orbit_radius = std::sqrt(cam.position[0] * cam.position[0] + cam.position[1] * cam.position[1]);
    double orbit_start_angle = std::atan2(cam.position[1], cam.position[0]);
    int frame = 0;

    auto write_frame = [&](const std::vector<bd::Placed<double>>& objs, const Camera<double>& c) {
        std::string path = std::format("{}/frame_{:04d}.png", dir, frame);
        if (spatium::examples::confirm_overwrite(path, force))
            write_png_rgb(path, W, H, render_frame(objs, c, W, H));
        ++frame;
    };

    for (int i = 0; i < build_frames; ++i) {
        double t = T_BUILD_END * static_cast<double>(i) / static_cast<double>(build_frames - 1);
        auto objs = bd::materialize(scene, lesson_idx, t);
        auto ax = axes_placed();
        objs.insert(objs.end(), ax.begin(), ax.end());
        write_frame(objs, cam);
    }
    auto full = bd::materialize(scene, lesson_idx, T_BUILD_END); // finished donut, axes hidden from here on
    for (int i = 0; i < orbit_frames; ++i) {
        double frac = static_cast<double>(i) / static_cast<double>(orbit_frames);
        double angle = orbit_start_angle + frac * 2.0 * std::numbers::pi;
        Camera<double> orbit_cam = cam;
        orbit_cam.position = {orbit_radius * std::cos(angle), orbit_radius * std::sin(angle), cam.position[2]};
        write_frame(full, orbit_cam);
    }
    std::println("donut_demo: wrote {} frames -> {}/", frame, dir);
    std::println("  ffmpeg -framerate 30 -i {}/frame_%04d.png -pix_fmt yuv420p donut.mp4", dir);
}

} // namespace

int main(int argc, char* argv[]) {
    bool force = false;
    bool photo = false;
    bool video = false;
    std::string out_path = "donut.png";
    std::string video_dir = "donut_frames";
    std::size_t sprinkle_count = 600; // 3x denser
    double t = T_BUILD_END; // how far into the build-up to render (T_BUILD_END = fully formed)
    int build_frames = 75, orbit_frames = 45; // half the sampling density -- ~9x more
                                              // dust triangles from copies_per_point makes
                                              // full 150/90 too costly for today; render at
                                              // 15fps (still smooth) instead of 30fps

    for (int i = 1; i < argc; ++i) {
        std::string_view a = argv[i];
        if (a == "--force") { force = true; continue; }
        if (a == "--photo") { photo = true; if (i + 1 < argc && argv[i + 1][0] != '-') out_path = argv[++i]; continue; }
        if (a == "--video") { video = true; if (i + 1 < argc && argv[i + 1][0] != '-') video_dir = argv[++i]; continue; }
        if (a == "--sprinkles" && i + 1 < argc) { sprinkle_count = static_cast<std::size_t>(std::atoi(argv[++i])); continue; }
        if (a == "--t" && i + 1 < argc) { t = std::atof(argv[++i]); continue; }
        if (a == "--help") {
            std::print("donut_demo [--photo [path]] [--video [dir]] [--sprinkles N] [--t seconds] [--force]\n"
                       "  Spatium's getting-started demo: delete the default cube (explode it,\n"
                       "  --t controls how far), then torus dough, offset-surface icing,\n"
                       "  area-sampled sprinkles -- declarative steps, see the file's own header\n"
                       "  comment. Prints the trace to the console by default; --photo also\n"
                       "  CPU-raytraces it to a PNG; --video renders the whole build-up + orbit\n"
                       "  as a PNG frame sequence (assemble with ffmpeg, printed at the end).\n");
            return 0;
        }
        std::print(stderr, "unknown option: {}\n", a);
        return 1;
    }

    auto t0 = std::chrono::steady_clock::now();

    bd::Trace<double> scene;

    // Shared by dough/icing/sprinkles below: grows from nothing to full
    // size during [T_DONUT_START, T_DONUT_END] (overlapping the dust's
    // dissolve), so the donut visibly *replaces* the dust rather than
    // simply being present the whole time -- still one closed-form
    // function of (point, t), no simulation loop.
    auto grow_scale = [](const Vec<double, 3>& p, double time) -> Vec<double, 3> {
        double e = smoothstep01((time - T_DONUT_START) / (T_DONUT_END - T_DONUT_START));
        return Vec<double, 3>{p * e};
    };

    // Step 0 -- the default cube. Visible briefly, static, then it's
    // gone -- what actually reads as "the cube" from here on is the
    // dust field below (`dust`), not this node's own motion.
    auto cube = scene.cube({0.9, 0.9, 0.9})
                    .colored(Material<double>{.base_color = {0.55, 0.55, 0.58}})
                    .moving([](const Vec<double, 3>& p, double time) -> Vec<double, 3> {
                        return Vec<double, 3>{p * (time < 0.12 ? 1.0 : 0.0)};
                    });

    // Step 0.5 -- delete the cube by *exploding* it: not a shrink this
    // time, real dust -- ~220 tiny cubes flying from the cube's own
    // volume, converging into the shape of the word "BOOM" (points
    // rasterized from a real font via ImageMagick, not hand-placed --
    // see examples/data/boom_points.txt and the shell one-liner that
    // made it), holding briefly, then dissolving before the donut grows
    // in. Each particle is its own tiny Literal node with its own
    // `.moving()` closure -- the DSL's per-node motion slot, just used
    // 220 times instead of once; still no simulation state anywhere.
    auto boom_uv = load_boom_points("examples/data/boom_points.txt");
    if (boom_uv.empty()) boom_uv = load_boom_points("data/boom_points.txt");
    Camera<double> text_cam = hero_camera();
    auto text_basis = make_camera_basis(text_cam);
    constexpr double text_scale = 2.6;
    // Positioned relative to the camera's own basis (a point along its
    // forward axis, offset along its screen-up), not a guessed world
    // coordinate -- guarantees it lands on-screen regardless of exactly
    // where the camera sits, rather than trial-and-error against one
    // specific camera pose.
    const Vec<double, 3> text_center{
        Vec<double, 3>{text_cam.position + text_basis.fwd * 7.0 + text_basis.up * 0.85}};

    std::mt19937 burst_rng(11);
    std::uniform_real_distribution<double> unit(-1.0, 1.0);
    std::uniform_real_distribution<double> dist_amt(1.6, 2.6);
    std::uniform_real_distribution<double> unit01(0.0, 1.0);
    algebra::PerlinNoise swirl_noise(4);

    // Grey far from any target, shifting toward red-yellow the closer a
    // particle's *current* position (recomputed live at whatever t this
    // gets asked about, via color_fn -- not baked in once) is to the
    // letterform point it's headed for -- unconditionally: even a
    // near-miss that only brushes past a letter on a weak pull picks up
    // color as it passes close, then fades back toward grey again as it
    // moves off. There's no separate "these ones never get colored"
    // case anymore; distance alone decides it.
    auto dust_color = [](const Vec<double, 3>& p, const Vec<double, 3>& target) {
        Vec<double, 3> grey{0.45, 0.44, 0.42};
        double d = std::clamp((p - target).norm() / 1.4, 0.0, 1.0);
        double hot = 1.0 - d;
        return Vec<double, 3>{grey * (1.0 - hot) + Vec<double, 3>{0.95, 0.35 + 0.35 * hot, 0.06} * hot};
    };

    // Each BOOM-letterform point gets several independent dust specks
    // converging on it (own burst direction/seed/pull each) rather than
    // regenerating a denser point cloud from the font raster -- denser
    // text *and* a denser swirling cloud together, same 220-point
    // letterform data. pull_strength is continuous, not a landed/missed
    // coin flip: most particles are weakly pulled (curve in close past a
    // letter, swept along by the shape rather than drifting
    // independently, then carried on past by the swirl) and a smaller
    // fraction pull strongly enough to actually land and hold.
    constexpr std::size_t copies_per_point = 90;
    std::vector<bd::Handle<double>> dust;
    dust.reserve(boom_uv.size() * copies_per_point);
    for (std::size_t i = 0; i < boom_uv.size() * copies_per_point; ++i) {
        Vec<double, 3> burst_dir = Vec<double, 3>{Vec<double, 3>{unit(burst_rng), unit(burst_rng), unit(burst_rng)}.normalized()};
        double burst_dist = dist_amt(burst_rng);
        double swirl_seed = unit01(burst_rng);
        double roll = unit01(burst_rng);
        // A real, if small, positional pull on the rest blurred the
        // letterforms into a haze (every gap between letters had
        // weakly-pulled specks drifting into it) -- pure swirl for
        // those keeps the shape crisp; they still pick up color via
        // dust_color's live distance check on whatever near passes
        // their own turbulent wandering happens to bring them, which is
        // the "just gets colored flying past" part without needing an
        // explicit pull to manufacture it.
        double pull_strength = roll < 0.3 ? (0.85 + 0.15 * unit01(burst_rng)) : 0.0; // ~30% land firmly
        auto [u, v] = boom_uv[i % boom_uv.size()];
        Vec<double, 3> target = Vec<double, 3>{
            text_center + text_basis.right * (u * text_scale) + text_basis.up * (v * text_scale)};

        dust.push_back(
            scene.literal(dust_speck(0.014))
                .colored(std::function<Vec<double, 3>(const Vec<double, 3>&, double)>{
                    [target, dust_color](const Vec<double, 3>& p, double) { return dust_color(p, target); }})
                .moving([burst_dir, burst_dist, target, pull_strength, swirl_seed, swirl_noise](
                            const Vec<double, 3>& p, double time) {
                    return particle_motion(p, time, burst_dir, burst_dist, target, pull_strength, swirl_seed, swirl_noise);
                }));
    }

    // Step 1 -- the dough is a torus, offset by a fine noise bump so it
    // actually has bready surface texture (not just a rough *shading*
    // response on a mathematically perfect torus -- real geometry, so
    // the bumps show up in the silhouette and catch light like bumps
    // should, not a flat-shaded illusion of them).
    // Calibrated against a real reference photo (Evan-Amos, "a pink,
    // frosted doughnut... Dunkin' Donuts", public domain, Wikimedia
    // Commons File:Pink-Frosted-Donut.jpg) rather than guessed -- two
    // things it corrected: real dough is close to smooth (a matte
    // *shading* response, not fine bumpy geometry -- the fine-noise
    // texture from an earlier pass was reading as noise, not "bread");
    // real icing covers almost the entire top uniformly, edge mostly
    // clean along the torus's natural equator, with a *few* localized
    // drips past it, not an all-over ragged boundary.
    algebra::PerlinNoise surface_noise(2);
    auto dough_bump = std::function<double(double, double)>{[surface_noise](double u, double v) {
        return 0.020 * surface_noise(u * 3.0, v * 3.0, 0.0); // visible but still broad, not fine-grain
    }};
    auto dough_base = scene.torus(2.0, 1.0, 160, 80);
    auto dough = scene.offset(dough_base, dough_bump)
                     .colored(Material<double>{.base_color = {0.87, 0.58, 0.27}, .roughness = 0.92})
                     .moving(grow_scale);

    // Step 2 -- the icing is the dough's own surface, offset outward --
    // over almost the whole top (torus_cap), clean-edged except for a
    // handful of localized drips, gently undulating (glaze pools and
    // settles in broad waves, not fine ripples) rather than a mirror.
    algebra::PerlinNoise icing_noise(3);
    auto icing_base = scene.offset(scene.space(torus_cap(2.0, 1.0), 160, 64), dough_bump);
    auto icing = scene.offset(icing_base, std::function<double(double, double)>{[icing_noise](double u, double v) {
                        constexpr double pi = std::numbers::pi;
                        double edge_dist = std::min(v - pi * 0.02, pi * 0.98 - v); // distance to the band's edge
                        double wobble = icing_noise(std::cos(u) * 2.0, std::sin(u) * 2.0, 0.0) * (pi * 0.03);
                        double drip = std::max(0.0, icing_noise(std::cos(u) * 1.3, std::sin(u) * 1.3, 8.0) - 0.35) * (pi * 0.35);
                        double falloff = std::clamp((edge_dist + wobble + drip) / (pi * 0.09), 0.0, 1.0);
                        falloff = falloff * falloff * (3.0 - 2.0 * falloff); // soft, not torn
                        double pooling = icing_noise(u * 3.0, v * 3.0, 4.0) * 0.045; // broad, gentle waves
                        return (0.10 + pooling) * falloff;
                    }})
                     .colored(Material<double>{.base_color = {0.98, 0.55, 0.68}, .roughness = 0.40})
                     .moving(grow_scale);

    // Step 3 -- sprinkles scatter across the icing, area-weighted, in
    // several colors -- one scatter() call per color (each its own
    // area-weighted placement, not a single-color scatter recolored
    // after the fact) so the mix reads as actual rainbow sprinkles.
    // Flat slivers, not standing cylinders (scatter() maps an item's
    // local z to the target's normal, so a shape *long* in z stands up
    // like a quill; long in x/y and thin in z lies flat instead) --
    // scatter() also gives each instance its own random spin about the
    // normal now, so they no longer all face the same way.
    // solid_cylinder()'s length runs along local z -- exactly the axis
    // scatter() maps onto the target's *normal*, so used as-is it stands
    // the sprinkle straight up again (the quill problem, back). Permute
    // (x,y,z) -> (z,x,y) so the length axis becomes x (a tangent
    // direction instead), then nudge down slightly so it sits sunk into
    // the icing rather than merely resting exactly half-in.
    auto sprinkle_mesh = solid_cylinder(0.024, 0.13, 8); // bigger
    for (auto& v : sprinkle_mesh.vertices) v = Vec<double, 3>{v[2], v[0], v[1] - 0.024};
    auto sprinkle = scene.literal(sprinkle_mesh);
    std::vector<Vec<double, 3>> sprinkle_colors{
        {0.95, 0.20, 0.25}, {0.98, 0.75, 0.15}, {0.25, 0.65, 0.35},
        {0.30, 0.45, 0.90}, {0.85, 0.30, 0.75}, {0.98, 0.98, 0.95}};
    std::vector<bd::Handle<double>> sprinkle_groups;
    std::size_t per_group = sprinkle_count / sprinkle_colors.size();
    for (std::size_t g = 0; g < sprinkle_colors.size(); ++g)
        sprinkle_groups.push_back(
            scene.scatter(sprinkle, icing, per_group, static_cast<std::uint32_t>(g * 97 + 11))
                .colored(Material<double>{.base_color = sprinkle_colors[g], .roughness = 0.35})
                .moving(grow_scale)); // scatter() places against icing's *analytic*, always-full-size
                                      // surface (resolve_surface() doesn't see .moving()) -- this is
                                      // what actually keeps sprinkles in sync with the growing donut

    std::vector<bd::Handle<double>> scene_children{cube, dough, icing};
    scene_children.insert(scene_children.end(), sprinkle_groups.begin(), sprinkle_groups.end());
    scene_children.insert(scene_children.end(), dust.begin(), dust.end());
    auto lesson = scene.compose(scene_children);

    std::println("donut_demo: a Trace is real data -- here it is, {} nodes:", scene.size());
    for (std::size_t i = 0; i < scene.size(); ++i)
        std::println("  [{}] {}", i, kind_name(scene.node(i).kind));

    auto placed = bd::materialize(scene, lesson.index, t);
    std::size_t verts = 0, faces = 0;
    for (auto& obj : placed) { verts += obj.mesh.vertex_count(); faces += obj.mesh.face_count(); }

    auto t1 = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::println("materialized at t={}: exploded cube + dough + icing + {} sprinkles -> {} vertices, {} triangles, {:.1f} ms",
                 t, sprinkle_count, verts, faces, ms);

    if (photo) render_photo(placed, out_path, force);
    if (video) render_video(scene, lesson.index, video_dir, force, build_frames, orbit_frames);
    return 0;
}
