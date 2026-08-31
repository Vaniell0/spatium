#include <benchmark/benchmark.h>
#include <spatium/geometry/intersection.hpp>
#include <spatium/geometry/triangle.hpp>
#include <spatium/geometry/line.hpp>
#include <spatium/geometry/box.hpp>

using namespace spatium;
using namespace spatium::geometry;

// ── Ray-Triangle ──────────────────────────────────────────────

static void BM_RayTriangleHit(benchmark::State& state) {
    Triangle3 tri({0, 0, 0}, {2, 0, 0}, {0, 2, 0});
    auto ray = *Ray3::from(Vec3{0.5, 0.5, 5.0}, Vec3{0.0, 0.0, -1.0});
    for (auto _ : state) {
        auto r = intersect(ray, tri);
        benchmark::DoNotOptimize(r);
    }
}
BENCHMARK(BM_RayTriangleHit);

static void BM_RayTriangleMiss(benchmark::State& state) {
    Triangle3 tri({0, 0, 0}, {1, 0, 0}, {0, 1, 0});
    auto ray = *Ray3::from(Vec3{5.0, 5.0, 5.0}, Vec3{0.0, 0.0, -1.0});
    for (auto _ : state) {
        auto r = intersect(ray, tri);
        benchmark::DoNotOptimize(r);
    }
}
BENCHMARK(BM_RayTriangleMiss);

// ── Ray-Box ───────────────────────────────────────────────────

static void BM_RayBoxHit(benchmark::State& state) {
    Box3 box{Vec3{-1, -1, -1}, Vec3{1, 1, 1}};
    auto ray = *Ray3::from(Vec3{0, 0, 5}, Vec3{0, 0, -1});
    for (auto _ : state) {
        auto r = intersect_parameters(ray, box);
        benchmark::DoNotOptimize(r);
    }
}
BENCHMARK(BM_RayBoxHit);

static void BM_RayBoxMiss(benchmark::State& state) {
    Box3 box{Vec3{-1, -1, -1}, Vec3{1, 1, 1}};
    auto ray = *Ray3::from(Vec3{5, 5, 5}, Vec3{1, 0, 0});
    for (auto _ : state) {
        auto r = intersect_parameters(ray, box);
        benchmark::DoNotOptimize(r);
    }
}
BENCHMARK(BM_RayBoxMiss);

