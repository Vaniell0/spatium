// Throughput of the Block-D point-to-surface narrow-phase queries
// + the IPC barrier evaluated on top of them. A synthetic cloud of
// query points covers both the in-band region (where the barrier is
// active and force is non-zero) and the far field (early-out path).
//
// Numbers to compare against:
//   - sphere/torus closed form: should land in tens of ns / op,
//   - parametric (Newton on UV): hundreds of µs / op for 8×8 grid +
//     ~5 Newton iters; the price of "any f:(u,v)→R³ becomes contact
//     target" without bespoke per-surface math.

#include <benchmark/benchmark.h>
#include <spatium/physics/mechanics/narrow_phase.hpp>
#include <spatium/spaces/parametric.hpp>
#include <random>
#include <vector>
#include <cmath>
#include <numbers>

using namespace spatium;
using namespace spatium::physics::mechanics;

namespace {

constexpr double kSphereR     = 1.0;
constexpr double kTorusR      = 1.0;
constexpr double kTorusr      = 0.25;
constexpr double kBandHat     = 0.05;

std::vector<Vec<double, 3>> query_cloud(std::size_t n,
                                        double inner, double outer)
{
    std::mt19937 rng(2026);
    std::uniform_real_distribution<double> radial(inner, outer);
    std::uniform_real_distribution<double> ang(0.0, 2.0 * std::numbers::pi);
    std::uniform_real_distribution<double> z(-0.5, 0.5);

    std::vector<Vec<double, 3>> pts;
    pts.reserve(n);
    while (pts.size() < n) {
        double r = radial(rng);
        double a = ang(rng);
        pts.push_back({r * std::cos(a), r * std::sin(a), z(rng)});
    }
    return pts;
}

geometry::Torus<double> make_torus() {
    geometry::Torus<double> t;
    t.major_radius = kTorusR;
    t.minor_radius = kTorusr;
    return t;
}

ParametricSurface<double> make_torus_surface() {
    using std::cos; using std::sin;
    constexpr double pi = std::numbers::pi_v<double>;
    return ParametricSurface<double>(
        [](double u, double v) -> Vec<double, 3> {
            return {(kTorusR + kTorusr * cos(v)) * cos(u),
                    (kTorusR + kTorusr * cos(v)) * sin(u),
                    kTorusr * sin(v)};
        },
        ParametricSurface<double>::Domain{0.0, 2.0 * pi, 0.0, 2.0 * pi},
        true, true);
}

} // namespace

// ── point_to_sphere ────────────────────────────────────────────

static void BM_PointToSphere_OutsideBand(benchmark::State& state) {
    auto pts = query_cloud(2048, 1.5, 3.0);   // far from the surface
    Vec<double, 3> c{};
    std::size_t i = 0;
    for (auto _ : state) {
        auto q = point_to_sphere(pts[i++ & 2047], c, kSphereR);
        benchmark::DoNotOptimize(q);
    }
}
BENCHMARK(BM_PointToSphere_OutsideBand);

static void BM_PointToSphere_InBand(benchmark::State& state) {
    // Points are on a thin shell just outside the unit sphere, all
    // inside the IPC active band.
    auto pts = query_cloud(2048, 1.0 + 1e-3, 1.0 + kBandHat * 0.8);
    Vec<double, 3> c{};
    std::size_t i = 0;
    for (auto _ : state) {
        auto q = point_to_sphere(pts[i++ & 2047], c, kSphereR);
        benchmark::DoNotOptimize(q);
    }
}
BENCHMARK(BM_PointToSphere_InBand);

// ── point_to_torus ─────────────────────────────────────────────

static void BM_PointToTorus_OutsideBand(benchmark::State& state) {
    auto pts = query_cloud(2048, 0.4, 2.0);
    auto t = make_torus();
    std::size_t i = 0;
    for (auto _ : state) {
        auto q = point_to_torus(pts[i++ & 2047], t);
        benchmark::DoNotOptimize(q);
    }
}
BENCHMARK(BM_PointToTorus_OutsideBand);

// ── point_to (Newton via Surface::project) ──────────

static void BM_PointToParametric_TorusSurface(benchmark::State& state) {
    auto pts = query_cloud(256, 1.1, 1.6);    // small cloud — Newton is slow
    auto surf = make_torus_surface();
    std::size_t i = 0;
    for (auto _ : state) {
        auto q = point_to(pts[i++ & 255], surf);
        benchmark::DoNotOptimize(q);
    }
}
BENCHMARK(BM_PointToParametric_TorusSurface);

// ── full IPC pipeline (query + energy + force) ─────────────────

static void BM_IpcContactPipeline_Sphere(benchmark::State& state) {
    auto pts = query_cloud(2048, 1.0 + 1e-3, 1.0 + kBandHat * 0.9);
    Vec<double, 3> c{};
    std::size_t i = 0;
    for (auto _ : state) {
        auto q = point_to_sphere(pts[i++ & 2047], c, kSphereR);
        auto E = ipc_contact_energy(q, kBandHat);
        auto F = ipc_contact_force(q, kBandHat);
        benchmark::DoNotOptimize(E);
        benchmark::DoNotOptimize(F);
    }
}
BENCHMARK(BM_IpcContactPipeline_Sphere);

static void BM_IpcContactPipeline_Torus(benchmark::State& state) {
    auto pts = query_cloud(2048, kTorusR + kTorusr + 1e-3,
                                  kTorusR + kTorusr + kBandHat * 0.9);
    auto t = make_torus();
    std::size_t i = 0;
    for (auto _ : state) {
        auto q = point_to_torus(pts[i++ & 2047], t);
        auto E = ipc_contact_energy(q, kBandHat);
        auto F = ipc_contact_force(q, kBandHat);
        benchmark::DoNotOptimize(E);
        benchmark::DoNotOptimize(F);
    }
}
BENCHMARK(BM_IpcContactPipeline_Torus);
