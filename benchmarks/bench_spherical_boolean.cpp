// Spherical polygon intersection: closed-form clip-composition vs.
// brute-force Monte Carlo area estimation at comparable accuracy.
// Mirrors bench_raycast.cpp's own closed-form-vs-brute-force framing
// (BM_RayQuadric_* vs BM_RayCast_Brute_Sphere there).

#include <benchmark/benchmark.h>
#include <spatium/geometry/spherical_polygon.hpp>
#include <numbers>
#include <random>

using namespace spatium;
using namespace spatium::geometry;

namespace {

// Two octant-sized spherical triangles overlapping in a known
// pi/4-area region (see test_spherical_polygon.cpp for the derivation
// of that value) -- large enough to require several clip iterations
// (one per edge of the second polygon), not a degenerate corner case.
SphericalPolygon<double> octant_a() {
    return {{Vec3{1, 0, 0}, Vec3{0, 1, 0}, Vec3{0, 0, 1}}, 1.0};
}

SphericalPolygon<double> octant_b_rotated45() {
    constexpr double s = std::numbers::sqrt2 / 2.0;
    return {{Vec3{s, s, 0}, Vec3{-s, s, 0}, Vec3{0, 0, 1}}, 1.0};
}

} // namespace

// ── Closed-form: multi-edge great-circle clip ─────────────────────

static void BM_SphericalIntersection_ClosedForm(benchmark::State& state) {
    auto a = octant_a();
    auto b = octant_b_rotated45();
    for (auto _ : state) {
        auto r = spherical_intersection(a, b);
        benchmark::DoNotOptimize(r);
    }
}
BENCHMARK(BM_SphericalIntersection_ClosedForm);

// ── Brute force: Monte Carlo area estimate at matching accuracy ──
// Samples uniform points on the sphere (normalized Gaussian vectors)
// and counts the fraction landing inside both polygons' contains().
// The `samples` argument controls estimator accuracy -- state.range(0)
// swept from a coarse to a fine estimate so the crossover against the
// closed-form cost above is visible directly in the benchmark table
// rather than asserted separately.

static void BM_SphericalIntersection_MonteCarlo(benchmark::State& state) {
    auto a = octant_a();
    auto b = octant_b_rotated45();
    auto samples = static_cast<std::size_t>(state.range(0));
    std::mt19937 rng(2026);
    std::normal_distribution<double> gauss(0.0, 1.0);
    for (auto _ : state) {
        std::size_t inside = 0;
        for (std::size_t i = 0; i < samples; ++i) {
            Vec3 v{gauss(rng), gauss(rng), gauss(rng)};
            auto p = v.normalized();
            if (a.contains(p) && b.contains(p)) ++inside;
        }
        double area = 4.0 * std::numbers::pi * double(inside) / double(samples);
        benchmark::DoNotOptimize(area);
    }
    state.counters["samples"] = static_cast<double>(samples);
}
// 1e5 samples ~ 1% relative error on this pi/4-out-of-4*pi fraction;
// 1e6 ~ 0.3%; swept up to 1e7 to show the cost of tightening further.
BENCHMARK(BM_SphericalIntersection_MonteCarlo)->Arg(100000)->Arg(1000000)->Arg(10000000);
