#include <benchmark/benchmark.h>
#include <spatium/spatial/bvh.hpp>
#include <spatium/geometry/triangle.hpp>
#include <spatium/geometry/line.hpp>
#include <random>

using namespace spatium;
using namespace spatium::geometry;
using namespace spatium::spatial;

static std::vector<Triangle3> make_grid(int nx, int ny) {
    std::vector<Triangle3> tris;
    for (int y = 0; y < ny; ++y)
        for (int x = 0; x < nx; ++x) {
            double fx = x, fy = y;
            tris.push_back(Triangle3({fx, fy, 0}, {fx + 1, fy, 0}, {fx, fy + 1, 0}));
            tris.push_back(Triangle3({fx + 1, fy, 0}, {fx + 1, fy + 1, 0}, {fx, fy + 1, 0}));
        }
    return tris;
}

// ── Build ─────────────────────────────────────────────────────

static void BM_BVH_Build_1K(benchmark::State& state) {
    auto tris = make_grid(22, 23); // ~1012 triangles
    for (auto _ : state) {
        auto bvh = BVH<Triangle3>::build(tris);
        benchmark::DoNotOptimize(bvh);
    }
}
BENCHMARK(BM_BVH_Build_1K);

static void BM_BVH_Build_10K(benchmark::State& state) {
    auto tris = make_grid(70, 72); // ~10080 triangles
    for (auto _ : state) {
        auto bvh = BVH<Triangle3>::build(tris);
        benchmark::DoNotOptimize(bvh);
    }
}
BENCHMARK(BM_BVH_Build_10K);

// ── Ray cast ──────────────────────────────────────────────────

static void BM_BVH_RayCast_1K(benchmark::State& state) {
    auto tris = make_grid(22, 23);
    auto bvh = BVH<Triangle3>::build(tris);
    auto ray = *Ray3::from(Vec3{11.0, 11.5, 10.0}, Vec3{0.0, 0.0, -1.0});
    for (auto _ : state) {
        auto hit = bvh.ray_cast(ray);
        benchmark::DoNotOptimize(hit);
    }
}
BENCHMARK(BM_BVH_RayCast_1K);

static void BM_BVH_RayCast_10K(benchmark::State& state) {
    auto tris = make_grid(70, 72);
    auto bvh = BVH<Triangle3>::build(tris);
    auto ray = *Ray3::from(Vec3{35.0, 36.0, 10.0}, Vec3{0.0, 0.0, -1.0});
    for (auto _ : state) {
        auto hit = bvh.ray_cast(ray);
        benchmark::DoNotOptimize(hit);
    }
}
BENCHMARK(BM_BVH_RayCast_10K);

// ── Nearest ───────────────────────────────────────────────────

static void BM_BVH_Nearest_1K(benchmark::State& state) {
    auto tris = make_grid(22, 23);
    auto bvh = BVH<Triangle3>::build(tris);
    for (auto _ : state) {
        auto r = bvh.nearest(Vec3{11.0, 11.5, 3.0});
        benchmark::DoNotOptimize(r);
    }
}
BENCHMARK(BM_BVH_Nearest_1K);

static void BM_BVH_Nearest_10K(benchmark::State& state) {
    auto tris = make_grid(70, 72);
    auto bvh = BVH<Triangle3>::build(tris);
    for (auto _ : state) {
        auto r = bvh.nearest(Vec3{35.0, 36.0, 3.0});
        benchmark::DoNotOptimize(r);
    }
}
BENCHMARK(BM_BVH_Nearest_10K);

// ── Brute-force comparison ────────────────────────────────────

static void BM_BruteForce_RayCast_1K(benchmark::State& state) {
    auto tris = make_grid(22, 23);
    auto ray = *Ray3::from(Vec3{11.0, 11.5, 10.0}, Vec3{0.0, 0.0, -1.0});
    for (auto _ : state) {
        double best_t = 1e30;
        Vec3 best_point;
        for (auto& tri : tris) {
            auto r = intersect(ray, tri);
            if (r) {
                auto t = (r.value() - ray.origin).dot(ray.direction);
                if (t >= 0 && t < best_t) {
                    best_t = t;
                    best_point = r.value();
                }
            }
        }
        benchmark::DoNotOptimize(best_point);
    }
}
BENCHMARK(BM_BruteForce_RayCast_1K);

static void BM_BruteForce_RayCast_10K(benchmark::State& state) {
    auto tris = make_grid(70, 72);
    auto ray = *Ray3::from(Vec3{35.0, 36.0, 10.0}, Vec3{0.0, 0.0, -1.0});
    for (auto _ : state) {
        double best_t = 1e30;
        Vec3 best_point;
        for (auto& tri : tris) {
            auto r = intersect(ray, tri);
            if (r) {
                auto t = (r.value() - ray.origin).dot(ray.direction);
                if (t >= 0 && t < best_t) {
                    best_t = t;
                    best_point = r.value();
                }
            }
        }
        benchmark::DoNotOptimize(best_point);
    }
}
BENCHMARK(BM_BruteForce_RayCast_10K);
