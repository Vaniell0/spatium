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
#include <spatium/physics/mechanics/rigid_contact.hpp>
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

// ── rigid_contact.hpp: pairwise dynamic-body queries ────────────
// Same query-cost comparison as above, for the closed-form pairwise
// contact queries (discrete sphere-sphere, discrete sphere-AABB, and
// continuous/swept sphere-sphere) that back the XPBD rigid-collision
// path — sit these next to the static-surface numbers above
// (`BM_PointToSphere*`/`BM_PointToTorus*` at 3.7-17.7 ns/op, the full
// point_to+energy+force IPC pipeline at 12.6-25.5 ns/op via
// `BM_IpcContactPipeline_*`) for an honest side-by-side against the
// closed-form work this header adds. The swept query is expected to
// cost more than the discrete one — it runs a full quadratic solve
// (`ray_quadric`) instead of a handful of dot products — so the two
// numbers are reported separately rather than implying "continuous
// detection is free."

static void BM_SphereSphereContact_Separated(benchmark::State& state) {
    auto pts = query_cloud(2048, 3.0, 6.0);   // well outside contact range
    Vec<double, 3> a{};
    std::size_t i = 0;
    for (auto _ : state) {
        auto q = sphere_sphere_contact(a, 1.0, pts[i++ & 2047], 1.0);
        benchmark::DoNotOptimize(q);
    }
}
BENCHMARK(BM_SphereSphereContact_Separated);

static void BM_SphereSphereContact_Overlapping(benchmark::State& state) {
    // Centers within 2·radius of the origin -- guaranteed overlap with
    // a unit sphere fixed at the origin.
    auto pts = query_cloud(2048, 0.1, 1.8);
    Vec<double, 3> a{};
    std::size_t i = 0;
    for (auto _ : state) {
        auto q = sphere_sphere_contact(a, 1.0, pts[i++ & 2047], 1.0);
        benchmark::DoNotOptimize(q);
    }
}
BENCHMARK(BM_SphereSphereContact_Overlapping);

static void BM_SweepSphereSphere_Crossing(benchmark::State& state) {
    // A crosses a stationary unit sphere at the origin every call --
    // the general (non-degenerate, real-root) path through
    // `ray_quadric`'s quadratic solve, the case a moving-body pair
    // actually needs continuous detection for.
    Vec<double, 3> b{};
    Vec<double, 3> disp_b{};
    std::vector<Vec<double, 3>> starts = query_cloud(2048, 3.0, 6.0);
    std::size_t i = 0;
    for (auto _ : state) {
        Vec<double, 3> a0 = starts[i++ & 2047];
        Vec<double, 3> disp_a = Vec<double, 3>{b - a0};   // straight at the target
        auto s = sweep_sphere_sphere(a0, 1.0, disp_a, b, 1.0, disp_b);
        benchmark::DoNotOptimize(s);
    }
}
BENCHMARK(BM_SweepSphereSphere_Crossing);

static void BM_SweepSphereSphere_NoCrossing(benchmark::State& state) {
    // Miss case: motion that passes well clear of the target sphere --
    // still runs the full quadratic solve (a miss is a real-vs-complex
    // root distinction, not an early out), so worth its own number.
    Vec<double, 3> b{};
    Vec<double, 3> disp_b{};
    std::mt19937 rng(99);
    std::uniform_real_distribution<double> off(3.0, 6.0);
    std::vector<Vec<double, 3>> pairs_a, pairs_disp;
    for (int k = 0; k < 2048; ++k) {
        Vec<double, 3> a0{-off(rng), 5.0, 0.0};
        pairs_a.push_back(a0);
        pairs_disp.push_back(Vec<double, 3>{off(rng) * 2.0, 0.0, 0.0});   // sweeps past, not through
    }
    std::size_t i = 0;
    for (auto _ : state) {
        std::size_t k = i++ & 2047;
        auto s = sweep_sphere_sphere(pairs_a[k], 1.0, pairs_disp[k], b, 1.0, disp_b);
        benchmark::DoNotOptimize(s);
    }
}
BENCHMARK(BM_SweepSphereSphere_NoCrossing);

static void BM_SphereAabbContact_Outside(benchmark::State& state) {
    auto pts = query_cloud(2048, 3.0, 6.0);
    Vec<double, 3> box_min{-1, -1, -1}, box_max{1, 1, 1};
    std::size_t i = 0;
    for (auto _ : state) {
        auto q = sphere_aabb_contact(pts[i++ & 2047], 0.5, box_min, box_max);
        benchmark::DoNotOptimize(q);
    }
}
BENCHMARK(BM_SphereAabbContact_Outside);

static void BM_SphereAabbContact_DeepPenetration(benchmark::State& state) {
    // Centers well inside the box on all axes -- exercises the
    // minimum-translation-axis fallback path, not just the cheap clamp.
    auto pts = query_cloud(2048, 0.0, 3.0);
    Vec<double, 3> box_min{-10, -10, -10}, box_max{10, 10, 10};
    std::size_t i = 0;
    for (auto _ : state) {
        auto q = sphere_aabb_contact(pts[i++ & 2047], 0.5, box_min, box_max);
        benchmark::DoNotOptimize(q);
    }
}
BENCHMARK(BM_SphereAabbContact_DeepPenetration);

// ── broad phase: O(n^2) AABB sweep at demo scale ────────────────
// This header targets "tens, not thousands" of bodies (see
// rigid_contact.hpp's file header for why that rules out spatial/
// bvh.hpp here) -- benchmark exactly that regime instead of an
// unrealistically large n that would make the O(n^2) choice look
// worse than it is at the scale it's actually meant for.
static void BM_BroadPhaseAabbPairs(benchmark::State& state) {
    auto pts = query_cloud(static_cast<std::size_t>(state.range(0)), 0.0, 10.0);
    std::vector<geometry::Box<3, double>> boxes;
    boxes.reserve(pts.size());
    for (auto& p : pts) boxes.push_back(sphere_aabb(p, 0.3));

    for (auto _ : state) {
        auto pairs = broad_phase_aabb_pairs(boxes);
        benchmark::DoNotOptimize(pairs);
    }
}
BENCHMARK(BM_BroadPhaseAabbPairs)->Arg(16)->Arg(32)->Arg(64);
