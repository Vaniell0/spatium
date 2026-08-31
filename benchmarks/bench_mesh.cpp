#include <benchmark/benchmark.h>
#include <spatium/mesh/mesh.hpp>
#include <spatium/mesh/subdivision.hpp>
#include <spatium/mesh/primitives.hpp>
#include <spatium/spaces/sphere.hpp>

using namespace spatium;
using namespace spatium::mesh;

using Sphere2 = Sphere<2>;

// ── Subdivision ───────────────────────────────────────────────

static void BM_SubdivideL1(benchmark::State& state) {
    Sphere2 sphere;
    auto mesh = icosahedron(sphere);
    for (auto _ : state) {
        auto result = subdivide_once(mesh, sphere);
        benchmark::DoNotOptimize(result);
    }
}
BENCHMARK(BM_SubdivideL1);

static void BM_SubdivideL2(benchmark::State& state) {
    Sphere2 sphere;
    auto mesh = icosahedron(sphere);
    for (auto _ : state) {
        auto result = subdivide(mesh, sphere, 2);
        benchmark::DoNotOptimize(result);
    }
}
BENCHMARK(BM_SubdivideL2);

static void BM_SubdivideL3(benchmark::State& state) {
    Sphere2 sphere;
    auto mesh = icosahedron(sphere);
    for (auto _ : state) {
        auto result = subdivide(mesh, sphere, 3);
        benchmark::DoNotOptimize(result);
    }
}
BENCHMARK(BM_SubdivideL3);

// ── Area computation ──────────────────────────────────────────

static void BM_MeshArea_L1(benchmark::State& state) {
    Sphere2 sphere;
    auto mesh = subdivide(icosahedron(sphere), sphere, 1);
    for (auto _ : state) {
        benchmark::DoNotOptimize(mesh.area(sphere));
    }
}
BENCHMARK(BM_MeshArea_L1);

static void BM_MeshArea_L3(benchmark::State& state) {
    Sphere2 sphere;
    auto mesh = subdivide(icosahedron(sphere), sphere, 3);
    for (auto _ : state) {
        benchmark::DoNotOptimize(mesh.area(sphere));
    }
}
BENCHMARK(BM_MeshArea_L3);
