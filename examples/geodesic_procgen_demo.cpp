// Geodesic procedural generation — a "planet" whose surface regions come
// from the mesh's own intrinsic (geodesic) metric, not from texture noise
// laid over flat geometry the way most procgen does it.
//
// Pipeline, entirely composed from existing pieces, no new geometry
// algorithm written for this demo:
//   1. `make_bumpy_sphere()` (spaces/implicit.hpp) — an ImplicitSurface
//      "world", F(x,y,z) = 0.
//   2. `marching_cubes()` (spaces/implicit.hpp) — tessellates it into a
//      Mesh<ImplicitSurface<T>>.
//   3. `MeshTopology::build()` + repeated `geodesic_voronoi()` calls
//      (mesh/topology.hpp, mesh/voronoi.hpp) — farthest-point sampling:
//      each new site is the vertex farthest (by mesh-graph geodesic
//      distance, multi-source Dijkstra) from every site already chosen,
//      which is exactly what geodesic_voronoi()'s own per-vertex distance
//      field already computes, called once per site instead of writing a
//      separate single-source Dijkstra.
//   4. `face_labels()` — per-triangle region id (or "boundary" for a
//      mixed-label face straddling two regions).
//   5. `BVH<Triangle3>` + `render::Camera`/`parallel_for_rows`/
//      `supersample_pixel`/`write_png_rgb` — the same CPU-raytracer
//      engine tumbling_body_demo.cpp established as the pattern every
//      offline demo should follow, pointed at a real mesh via
//      `BVH::ray_cast()` (returns `Hit::index`, the winning triangle's
//      position in the array passed to `BVH::build()` — built here in
//      exactly `mesh.triangles()`'s order, so `hit->index` is directly
//      `face_labels[hit->index]`).
//
// Output: geodesic_world.png

#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "io_helpers.hpp"

#include <spatium/geometry/triangle.hpp>
#include <spatium/mesh/topology.hpp>
#include <spatium/mesh/voronoi.hpp>
#include <spatium/render/camera.hpp>
#include <spatium/render/color.hpp>
#include <spatium/render/parallel_for_rows.hpp>
#include <spatium/render/sky.hpp>
#include <spatium/render/supersample.hpp>
#include <spatium/render/write_image.hpp>
#include <spatium/spaces/implicit.hpp>
#include <spatium/spatial/bvh.hpp>

#include <chrono>
#include <cstdlib>
#include <print>
#include <string>
#include <string_view>
#include <vector>

using spatium::ImplicitSurface;
using spatium::Vec;
using spatium::make_bumpy_sphere;
using spatium::marching_cubes;
using spatium::geometry::Ray;
using spatium::geometry::Triangle3;
using spatium::mesh::face_labels;
using spatium::mesh::geodesic_voronoi;
using spatium::mesh::MeshTopology;
using spatium::mesh::no_vertex;
using spatium::spatial::BVH;
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

using World = ImplicitSurface<double>;

// Farthest-point sampling via repeated geodesic_voronoi(): each new site
// is the vertex with the largest distance to its nearest already-chosen
// site -- exactly the per-vertex field geodesic_voronoi() already
// computes, so no separate single-source Dijkstra is written here.
std::vector<uint32_t> farthest_point_sample(const MeshTopology<World>& topo, const World& space,
                                             uint32_t first, std::size_t count) {
    std::vector<uint32_t> sites{first};
    for (std::size_t i = 1; i < count; ++i) {
        auto vd = geodesic_voronoi(topo, space, sites);
        uint32_t farthest = 0;
        double best = -1.0;
        for (uint32_t v = 0; v < vd.distances.size(); ++v) {
            if (vd.distances[v] > best) { best = vd.distances[v]; farthest = v; }
        }
        sites.push_back(farthest);
    }
    return sites;
}

} // namespace

int main(int argc, char* argv[]) {
    bool force = false;
    int resolution = 64;
    std::size_t regions = 12;
    std::string out_path = "geodesic_world.png";

    for (int i = 1; i < argc; ++i) {
        std::string_view a = argv[i];
        if (a == "--force") { force = true; continue; }
        if (a == "--regions" && i + 1 < argc) { regions = std::atoi(argv[++i]); continue; }
        if (a == "--resolution" && i + 1 < argc) { resolution = std::atoi(argv[++i]); continue; }
        if (a == "--out" && i + 1 < argc) { out_path = argv[++i]; continue; }
        if (a == "--help") {
            std::print("geodesic_procgen_demo [--regions N] [--resolution N] [--out path] [--force]\n"
                       "  Generates a bumpy-sphere world, partitions its surface into N regions\n"
                       "  by geodesic distance (not flat Euclidean texture noise), renders it.\n");
            return 0;
        }
        std::print(stderr, "unknown option: {}\n", a);
        return 1;
    }

    auto t0 = std::chrono::steady_clock::now();

    World world = make_bumpy_sphere<double>(1.0, 0.15);
    auto mesh = marching_cubes(world, static_cast<std::size_t>(resolution));

    auto topo = MeshTopology<World>::build(mesh);
    auto sites = farthest_point_sample(topo, world, 0, regions);
    auto vd = geodesic_voronoi(topo, world, sites);
    auto flabels = face_labels(vd, mesh);

    std::vector<Triangle3> tris;
    tris.reserve(mesh.face_count());
    for (auto [a, b, c] : mesh.triangles())
        tris.push_back(Triangle3(a, b, c));
    auto bvh = BVH<Triangle3>::build(tris);

    std::vector<Vec<double, 3>> palette(regions);
    for (std::size_t i = 0; i < regions; ++i)
        palette[i] = hsv_to_rgb255(static_cast<double>(i) / static_cast<double>(regions), 0.55, 0.92);
    const Vec<double, 3> boundary_color{18.0, 18.0, 22.0};

    constexpr int W = 960, H = 720;
    // wide_sky=false: this is a close, narrow-FOV single-object shot, not
    // a wide-FOV whole-sky view -- see Sky::wide_sky's doc comment.
    Sky sky = make_starfield(2000, /*seed=*/7, /*tint=*/{6.0, 6.0, 14.0}, /*wide_sky=*/false);

    const Camera<double> cam{
        .position = {4.2, 3.4, 2.6}, .target = {0.0, 0.0, 0.0}, .up = {0.0, 0.0, 1.0},
        .fov_deg = 40.0};
    const auto basis = make_camera_basis(cam);
    const Vec<double, 3> light = Vec<double, 3>{Vec<double, 3>{0.6, -0.5, 1.0}.normalized()};

    std::vector<std::uint8_t> img(3 * static_cast<std::size_t>(W) * H, 0);

    parallel_for_rows(H, [&](int y) {
        for (int x = 0; x < W; ++x) {
            std::uint8_t* px = &img[3 * (static_cast<std::size_t>(y) * W + x)];

            auto ray_color = [&](double sx, double sy) -> Vec<double, 3> {
                Vec<double, 3> dir = camera_ray_dir(basis, sx, sy);
                auto hit = bvh.ray_cast(Ray<3, double>{cam.position, dir});
                if (!hit) return sample_sky_color(sky, dir);

                Vec<double, 3> n = hit->normal;
                if (n.dot(dir) > 0.0) n = Vec<double, 3>{-n};
                double diff = std::max(0.0, n.dot(light));
                double shaded = 0.30 + 0.70 * diff;

                auto label = flabels[hit->index];
                Vec<double, 3> base = (label == no_vertex) ? boundary_color : palette[label];
                return Vec<double, 3>{base * shaded};
            };

            supersample_pixel(x, y, W, H, basis.tan_half,
                              static_cast<double>(W) / H, ray_color, px);
        }
    });

    if (spatium::examples::confirm_overwrite(out_path, force))
        write_png_rgb(out_path, W, H, img);

    auto t1 = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::print("geodesic_procgen_demo: {} vertices, {} faces, {} regions, {:.0f} ms -> {}\n",
               mesh.vertex_count(), mesh.face_count(), regions, ms, out_path);
    return 0;
}
