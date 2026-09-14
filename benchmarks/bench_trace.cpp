// Cost of the DSL's per-node type erasure, measured rather than assumed.
//
// spatium::io::build::TraceNode holds its motion and color hooks as
// std::function (build.hpp's `transform`, `color_fn`, `thickness`).
// materialize_mesh() then applies `transform` once per *vertex*:
//
//     if (n.transform)
//         for (auto& v : out.vertices) v = n.transform(v, t);
//
// examples/donut_demo.cpp builds 220 letterform points x 90 copies =
// 19 800 Literal nodes, each a 12-vertex icosahedron, each carrying its
// own .moving() and .colored() closure -- so one frame is ~19 800
// indirect color calls and ~237 600 indirect vertex calls. That is the
// shape benchmarked here.
//
// The interesting number is not the total, it is how much of the total
// the indirection actually accounts for, so the same work is measured
// three ways:
//
//   - through std::function, as the DSL does today,
//   - through the identical lambda called directly, so the compiler can
//     inline it (what a lowered/monomorphised trace could reach),
//   - with no motion hook at all, to show the floor that neither form
//     can go below (mesh copy + assembly).
//
// Both a trivial callee and a realistic one (Perlin swirl, matching the
// donut's own particle motion) are measured: a cheap callee makes the
// indirect call look disproportionately expensive, an expensive one
// amortises it, and the honest answer lies between the two.

#include <benchmark/benchmark.h>
#include <spatium/io/build.hpp>
#include <spatium/algebra/noise.hpp>
#include <spatium/mesh/primitives.hpp>
#include <spatium/spaces/sphere.hpp>
#include <functional>
#include <random>
#include <vector>

using namespace spatium;
using namespace spatium::io::build;

namespace {

using V3 = Vec<double, 3>;

// The donut demo's own dust speck: an icosahedron, 12 vertices.
mesh::Mesh<Euclidean<3, double>> dust_speck(double radius) {
    auto ico = mesh::icosahedron(Sphere<2, double>{radius});
    mesh::Mesh<Euclidean<3, double>> m;
    m.vertices.assign(ico.vertices.begin(), ico.vertices.end());
    m.faces = ico.faces;
    return m;
}

// Trivial callee: a time-scaled translation. Nothing to amortise an
// indirect call against.
struct CheapMotion {
    V3 dir;
    V3 operator()(const V3& p, double t) const { return V3{p + dir * t}; }
};

// Realistic callee, shaped like the donut demo's particle_motion():
// three noise lookups plus the arithmetic around them.
struct SwirlMotion {
    const algebra::PerlinNoise* noise;
    V3 dir;
    double dist;
    V3 operator()(const V3& p, double t) const {
        V3 base{p + dir * (dist * t)};
        double s = 0.7;
        V3 swirl{
            (*noise)(base[0] * s, base[1] * s, t),
            (*noise)(base[1] * s, base[2] * s, t),
            (*noise)(base[2] * s, base[0] * s, t)};
        return V3{base + swirl * 0.35};
    }
};

// Build a dust-shaped scene: `count` Literal nodes composed into one
// root, each carrying its own motion hook when `with_motion`.
template<typename Motion>
std::pair<Trace<double>, std::size_t>
build_dust(std::size_t count, Motion motion, bool with_motion) {
    Trace<double> trace;
    auto speck = dust_speck(0.014);
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> unit(-1.0, 1.0);

    std::vector<Handle<double>> nodes;
    nodes.reserve(count);
    for (std::size_t i = 0; i < count; ++i) {
        auto h = trace.literal(speck);
        if (with_motion) {
            V3 dir{unit(rng), unit(rng), unit(rng)};
            auto m = motion;
            m.dir = dir;
            h = h.moving([m](const V3& p, double t) { return m(p, t); });
        }
        nodes.push_back(h);
    }
    auto root = trace.compose(nodes);
    return {std::move(trace), root.index};
}

constexpr std::size_t kDustCount = 19800; // donut_demo: 220 * 90

} // namespace

// ── Whole-scene materialize, the real operation ───────────────────

static void BM_Materialize_NoMotion(benchmark::State& state) {
    auto [trace, root] = build_dust(kDustCount, CheapMotion{}, false);
    for (auto _ : state)
        benchmark::DoNotOptimize(materialize(trace, root, 0.5));
    state.SetItemsProcessed(state.iterations() * kDustCount);
}
BENCHMARK(BM_Materialize_NoMotion)->Unit(benchmark::kMillisecond);

static void BM_Materialize_CheapMotion(benchmark::State& state) {
    auto [trace, root] = build_dust(kDustCount, CheapMotion{}, true);
    for (auto _ : state)
        benchmark::DoNotOptimize(materialize(trace, root, 0.5));
    state.SetItemsProcessed(state.iterations() * kDustCount);
}
BENCHMARK(BM_Materialize_CheapMotion)->Unit(benchmark::kMillisecond);

static void BM_Materialize_SwirlMotion(benchmark::State& state) {
    static algebra::PerlinNoise noise(3);
    auto [trace, root] = build_dust(kDustCount, SwirlMotion{&noise, {}, 1.0}, true);
    for (auto _ : state)
        benchmark::DoNotOptimize(materialize(trace, root, 0.5));
    state.SetItemsProcessed(state.iterations() * kDustCount);
}
BENCHMARK(BM_Materialize_SwirlMotion)->Unit(benchmark::kMillisecond);

// ── The indirection in isolation ──────────────────────────────────
//
// Same vertices, same callee, same arithmetic; the only difference is
// whether the call goes through std::function or is visible to the
// optimiser. This is the part a lowered trace could actually recover.

std::vector<V3> make_vertices(std::size_t n) {
    std::vector<V3> v;
    v.reserve(n);
    std::mt19937 rng(11);
    std::uniform_real_distribution<double> unit(-1.0, 1.0);
    for (std::size_t i = 0; i < n; ++i) v.push_back(V3{unit(rng), unit(rng), unit(rng)});
    return v;
}

template<typename Motion>
static void run_direct(benchmark::State& state, Motion motion) {
    auto verts = make_vertices(static_cast<std::size_t>(state.range(0)));
    for (auto _ : state) {
        for (auto& v : verts) v = motion(v, 0.5);
        benchmark::DoNotOptimize(verts.data());
    }
    state.SetItemsProcessed(state.iterations() * state.range(0));
}

template<typename Motion>
static void run_erased(benchmark::State& state, Motion motion) {
    auto verts = make_vertices(static_cast<std::size_t>(state.range(0)));
    std::function<V3(const V3&, double)> f =
        [motion](const V3& p, double t) { return motion(p, t); };
    for (auto _ : state) {
        for (auto& v : verts) v = f(v, 0.5);
        benchmark::DoNotOptimize(verts.data());
    }
    state.SetItemsProcessed(state.iterations() * state.range(0));
}

static void BM_Vertices_Cheap_Direct(benchmark::State& s) {
    run_direct(s, CheapMotion{V3{0.1, 0.2, 0.3}});
}
static void BM_Vertices_Cheap_StdFunction(benchmark::State& s) {
    run_erased(s, CheapMotion{V3{0.1, 0.2, 0.3}});
}
static void BM_Vertices_Swirl_Direct(benchmark::State& s) {
    static algebra::PerlinNoise noise(3);
    run_direct(s, SwirlMotion{&noise, V3{0.1, 0.2, 0.3}, 1.0});
}
static void BM_Vertices_Swirl_StdFunction(benchmark::State& s) {
    static algebra::PerlinNoise noise(3);
    run_erased(s, SwirlMotion{&noise, V3{0.1, 0.2, 0.3}, 1.0});
}

// 237 600 = the donut demo's real per-frame vertex count (19 800 x 12).
BENCHMARK(BM_Vertices_Cheap_Direct)->Arg(237600);
BENCHMARK(BM_Vertices_Cheap_StdFunction)->Arg(237600);
BENCHMARK(BM_Vertices_Swirl_Direct)->Arg(237600);
BENCHMARK(BM_Vertices_Swirl_StdFunction)->Arg(237600);
