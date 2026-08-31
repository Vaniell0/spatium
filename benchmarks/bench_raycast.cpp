// Raycast throughput benchmarks across BVH, brute force and analytical
// ray-quadric. Uses a single tessellated unit sphere to give BVH vs
// brute-force a fair comparison; quadric path uses the analytical sphere.

#include <benchmark/benchmark.h>
#include <spatium/spatial/bvh.hpp>
#include <spatium/geometry/triangle.hpp>
#include <spatium/geometry/line.hpp>
#include <spatium/geometry/ray_surface.hpp>
#include <spatium/geometry/ray_parametric.hpp>
#include <spatium/geometry/make.hpp>
#include <spatium/mesh/primitives.hpp>
#include <spatium/mesh/subdivision.hpp>
#include <spatium/spaces/sphere.hpp>
#include <spatium/spaces/parametric.hpp>
#include <random>
#include <vector>

using namespace spatium;
using namespace spatium::geometry;
using namespace spatium::mesh;
using namespace spatium::spatial;

namespace {

std::vector<Triangle3> unit_sphere_triangles(std::size_t subdivisions) {
    Sphere<2> s;
    auto mesh = subdivide(icosahedron(s), s, subdivisions);
    std::vector<Triangle3> tris;
    tris.reserve(mesh.face_count());
    for (auto [a, b, c] : mesh.triangles())
        tris.push_back(Triangle3(a, b, c));
    return tris;
}

std::vector<Ray<3>> make_rays(std::size_t count, double spread = 0.6) {
    std::mt19937 rng(2026);
    std::uniform_real_distribution<double> j(-spread, spread);
    std::vector<Ray<3>> out;
    out.reserve(count);
    Vec3 origin{3, 3, 3};
    Vec3 aim{0, 0, 0};
    auto forward = (aim - origin).normalized();
    Vec3 up{0, 0, 1};
    auto right = forward.cross(up).normalized();
    up = right.cross(forward).normalized();
    while (out.size() < count) {
        auto d = Vec3{forward + right * j(rng) + up * j(rng)}.normalized();
        if (auto r = ray(origin, d)) out.push_back(*r);
    }
    return out;
}

} // namespace

// ── BVH ray_cast on tessellated sphere ────────────────────────

static void BM_RayCast_BVH_Sphere(benchmark::State& state) {
    auto tris = unit_sphere_triangles(static_cast<std::size_t>(state.range(0)));
    auto bvh = BVH<Triangle3>::build(tris);
    auto rays = make_rays(256);
    std::size_t i = 0;
    for (auto _ : state) {
        auto hit = bvh.ray_cast(rays[i++ & 255]);
        benchmark::DoNotOptimize(hit);
    }
    state.counters["tris"] = static_cast<double>(tris.size());
}
BENCHMARK(BM_RayCast_BVH_Sphere)->Arg(2)->Arg(4)->Arg(5);

static void BM_RayCast_Brute_Sphere(benchmark::State& state) {
    auto tris = unit_sphere_triangles(static_cast<std::size_t>(state.range(0)));
    auto rays = make_rays(256);
    std::size_t i = 0;
    for (auto _ : state) {
        auto& r = rays[i++ & 255];
        double best_t = 1e30;
        bool any = false;
        for (auto& tri : tris) {
            if (auto it = intersect(r, tri)) {
                double t = (*it - r.origin).dot(r.direction);
                if (t >= 0 && t < best_t) { best_t = t; any = true; }
            }
        }
        benchmark::DoNotOptimize(any);
    }
    state.counters["tris"] = static_cast<double>(tris.size());
}
BENCHMARK(BM_RayCast_Brute_Sphere)->Arg(2)->Arg(4)->Arg(5);

// ── Analytical ray-quadric ────────────────────────────────────

static void BM_RayQuadric_Sphere(benchmark::State& state) {
    auto q = Quadric<double>::sphere(1.0);
    auto rays = make_rays(256);
    std::size_t i = 0;
    for (auto _ : state) {
        auto hits = ray_quadric(rays[i++ & 255], q);
        benchmark::DoNotOptimize(hits);
    }
}
BENCHMARK(BM_RayQuadric_Sphere);

static void BM_RayQuadric_Ellipsoid(benchmark::State& state) {
    auto q = Quadric<double>::ellipsoid(1.0, 0.5, 0.3);
    auto rays = make_rays(256);
    std::size_t i = 0;
    for (auto _ : state) {
        auto hits = ray_quadric(rays[i++ & 255], q);
        benchmark::DoNotOptimize(hits);
    }
}
BENCHMARK(BM_RayQuadric_Ellipsoid);

static void BM_RayQuadric_Cylinder(benchmark::State& state) {
    auto q = Quadric<double>::cylinder_z(0.5);
    auto rays = make_rays(256);
    std::size_t i = 0;
    for (auto _ : state) {
        auto hits = ray_quadric(rays[i++ & 255], q);
        benchmark::DoNotOptimize(hits);
    }
}
BENCHMARK(BM_RayQuadric_Cylinder);

static void BM_RayQuadric_Proximity(benchmark::State& state) {
    auto q = Quadric<double>::sphere(1.0);
    // Rays aimed away from sphere to hit the miss-proximity path
    auto rays = make_rays(256, 1.5);
    std::size_t i = 0;
    for (auto _ : state) {
        auto prox = ray_quadric_proximity(rays[i++ & 255], q);
        benchmark::DoNotOptimize(prox);
    }
}
BENCHMARK(BM_RayQuadric_Proximity);

// ── Analytical ray-torus ─────────────────────────────────────

static void BM_RayTorus(benchmark::State& state) {
    Torus<double> torus{.major_radius = 1.0, .minor_radius = 0.3};
    auto rays = make_rays(256);
    std::size_t i = 0;
    for (auto _ : state) {
        auto hits = ray_torus(rays[i++ & 255], torus);
        benchmark::DoNotOptimize(hits);
    }
}
BENCHMARK(BM_RayTorus);

static void BM_RayTorus_Proximity(benchmark::State& state) {
    Torus<double> torus{.major_radius = 1.0, .minor_radius = 0.3};
    auto rays = make_rays(256, 2.5);
    std::size_t i = 0;
    for (auto _ : state) {
        auto prox = ray_torus_proximity(rays[i++ & 255], torus);
        benchmark::DoNotOptimize(prox);
    }
}
BENCHMARK(BM_RayTorus_Proximity);

// ── Ray-ParametricSurface via Newton UV ──────────────────────
// Compares against BM_RayTorus (closed-form quartic) on the same torus.

static void BM_RayParametric_Torus(benchmark::State& state) {
    auto surf = make_torus<double>(1.0, 0.3);
    auto rays = make_rays(256);
    std::size_t i = 0;
    for (auto _ : state) {
        auto hits = ray_parametric(rays[i++ & 255], surf);
        benchmark::DoNotOptimize(hits);
    }
}
BENCHMARK(BM_RayParametric_Torus);

static void BM_RayParametric_Mobius(benchmark::State& state) {
    auto surf = make_mobius<double>(1.0, 0.3);
    auto rays = make_rays(256);
    std::size_t i = 0;
    for (auto _ : state) {
        auto hits = ray_parametric(rays[i++ & 255], surf);
        benchmark::DoNotOptimize(hits);
    }
}
BENCHMARK(BM_RayParametric_Mobius);

static void BM_RayParametric_TorusFirst(benchmark::State& state) {
    auto surf = make_torus<double>(1.0, 0.3);
    auto rays = make_rays(256);
    std::size_t i = 0;
    for (auto _ : state) {
        auto hit = ray_parametric_first(rays[i++ & 255], surf);
        benchmark::DoNotOptimize(hit);
    }
}
BENCHMARK(BM_RayParametric_TorusFirst);
