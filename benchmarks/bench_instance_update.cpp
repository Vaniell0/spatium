// Following a changing instance set with BVH::insert / remove, against
// building the tree again.
//
// The instances are the real thing: cylinders scattered over a torus by
// the scene DSL, cooked and laid out. Three cases per (N, k):
//
//   Rebuild  -- BVH::build over the N instances, what a regenerated scene
//               pays today for its instance tree.
//   Edit     -- k removals and k insertions on a tree of N, the tree
//               operations alone (no diff): a random k leave, k new ones
//               come in from a second scatter.
//   Diff     -- render::detail::update_instance_tree, which also finds
//               which instances changed (hash, sort, compare) in two
//               N-sized lists that differ in k places, flipping between
//               them.
//
// Then the whole regeneration of a scene, one step of a count slider: a
// scatter of N copies becomes one of N + k (the DSL scene of
// examples/scatter_regen.cpp).
//
//   RegenCookCold   -- cook() of the new scene from scratch.
//   RegenCookResume -- cook() given the previous cook: the copies of an
//                      unchanged scatter are taken, only the new ones placed.
//   RegenUpdate     -- lay_out(previous scene, new cook): the instance tree
//                      updated, with the instances found by hashing every one
//                      (the first measurement of this file) or by comparing
//                      the scatter's generator (now).
//   RegenFresh      -- lay_out(new cook) from scratch.
//
// Run on the efficiency cores with medians of five:
//   taskset -c 4-11 spatium_bench --benchmark_filter='Instance|Regen' --benchmark_repetitions=5
//   (add --benchmark_report_aggregates_only=true for the medians alone)

#include <benchmark/benchmark.h>
#include <spatium/io/build.hpp>
#include <spatium/render/cooked_scene.hpp>
#include <map>
#include <optional>
#include <random>

using namespace spatium;
namespace bd = spatium::io::build;

namespace {

using Inst = geometry::Instanced<geometry::BoundedQuadric<double>>;
using Tree = spatial::BVH<Inst>;

struct Set {
    render::CookedScene<double> scene;   // owns the quadric the instances point into
    std::vector<Inst> a, pool;           // N instances; N more from another seed
    std::vector<render::InstanceLook<double>> looks_a, looks_pool;
};

render::CookedScene<double> scatter(std::size_t n, std::uint32_t seed) {
    bd::Trace<double> tr;
    auto torus = tr.torus(2.0, 1.0, 48, 24);
    auto item = tr.cylinder(0.011, 0.062, 8, 2);
    auto sc = tr.scatter(item, torus, n, seed, 0.45, bd::SeatAxis::X);
    auto root = tr.compose({torus, sc});
    return render::lay_out(bd::cook(tr, root.index, 0.0));
}

// Both lists must point into one quadric table: the second scene's
// instances are re-pointed at the first's.
const Set& set_of(std::size_t n) {
    static std::map<std::size_t, Set> sets;
    auto it = sets.find(n);
    if (it != sets.end()) return it->second;
    Set s;
    s.scene = scatter(n, 42);
    auto other = scatter(n, 43);
    s.a = s.scene.instances.shapes();
    s.pool = other.instances.shapes();
    for (auto& p : s.pool) p.shape = s.scene.quadrics.data();
    s.looks_a = s.scene.looks;
    s.looks_pool = other.looks;
    return sets.emplace(n, std::move(s)).first->second;
}

}  // namespace

static void BM_InstanceRebuild(benchmark::State& state) {
    const auto& s = set_of(static_cast<std::size_t>(state.range(0)));
    for (auto _ : state) {
        auto tree = Tree::build(s.a);
        benchmark::DoNotOptimize(tree);
    }
}

static void BM_InstanceEdit(benchmark::State& state) {
    const auto& s = set_of(static_cast<std::size_t>(state.range(0)));
    const auto k = static_cast<std::size_t>(state.range(1));
    auto tree = Tree::build(s.a);
    std::mt19937_64 g(1);
    std::size_t next = 0;
    for (auto _ : state) {
        for (std::size_t i = 0; i < k; ++i) tree.remove(g() % tree.size());
        for (std::size_t i = 0; i < k; ++i) tree.insert(s.pool[next++ % s.pool.size()]);
    }
    state.counters["rebuilds"] = static_cast<double>(tree.rebuilds());
}

static void BM_InstanceDiff(benchmark::State& state) {
    const auto& s = set_of(static_cast<std::size_t>(state.range(0)));
    const auto k = static_cast<std::size_t>(state.range(1));
    // b is a with its first k replaced by pool entries.
    auto b = s.a;
    auto looks_b = s.looks_a;
    for (std::size_t i = 0; i < k; ++i) {
        b[i] = s.pool[i];
        looks_b[i] = s.looks_pool[i];
    }
    auto tree = Tree::build(s.a);
    auto looks = s.looks_a;
    bool to_b = true;
    for (auto _ : state) {
        // The N-sized copies the call consumes are not what is measured.
        state.PauseTiming();
        auto next = to_b ? b : s.a;
        auto next_looks = to_b ? looks_b : s.looks_a;
        state.ResumeTiming();
        render::detail::update_instance_tree<double>(tree, looks, std::move(next), std::move(next_looks), 1.0);
        to_b = !to_b;
    }
}

#define SPATIUM_INSTANCE_ARGS \
    ->Args({1000, 1})->Args({1000, 10})->Args({1000, 100}) \
    ->Args({10000, 1})->Args({10000, 10})->Args({10000, 100}) \
    ->Args({100000, 1})->Args({100000, 10})->Args({100000, 100}) \
    ->Unit(benchmark::kMicrosecond)

BENCHMARK(BM_InstanceRebuild)->Arg(1000)->Arg(10000)->Arg(100000)->Unit(benchmark::kMicrosecond);
BENCHMARK(BM_InstanceEdit) SPATIUM_INSTANCE_ARGS;
BENCHMARK(BM_InstanceDiff) SPATIUM_INSTANCE_ARGS;

// ── Regenerating a whole scene ───────────────────────────────────────────

namespace {

struct Regen {
    bd::Trace<double> trace;
    std::size_t root = 0;
    Regen(std::size_t count, std::uint32_t seed = 42) {
        auto torus = trace.torus(2.0, 1.0, 240, 120);
        auto item = trace.cylinder(0.011, 0.062, 8, 2);
        auto sc = trace.scatter(item, torus, count, seed, 0.45, bd::SeatAxis::X);
        root = trace.compose({torus, sc}).index;
    }
};

// The step: N copies to N + k.
struct Step {
    bd::Cooked<double> a, b;
    explicit Step(std::size_t n, std::size_t k) {
        Regen ra(n), rb(n + k);
        a = bd::cook(ra.trace, ra.root, 0.0);
        b = bd::cook(rb.trace, rb.root, 0.0);
    }
};

}  // namespace

static void BM_RegenCookCold(benchmark::State& state) {
    const auto n = static_cast<std::size_t>(state.range(0)), k = static_cast<std::size_t>(state.range(1));
    Regen r(n + k);
    for (auto _ : state) {
        std::optional<bd::Cooked<double>> c;
        c.emplace(bd::cook(r.trace, r.root, 0.0));
        benchmark::DoNotOptimize(c);
        state.PauseTiming();
        c.reset();
        state.ResumeTiming();
    }
}

#ifndef SPATIUM_BENCH_WITHOUT_COOK_RESUME
static void BM_RegenCookResume(benchmark::State& state) {
    const auto n = static_cast<std::size_t>(state.range(0)), k = static_cast<std::size_t>(state.range(1));
    Regen ra(n), rb(n + k);
    const auto before = bd::cook(ra.trace, ra.root, 0.0);
    for (auto _ : state) {
        std::optional<bd::Cooked<double>> c;
        c.emplace(bd::cook(rb.trace, rb.root, 0.0, &before));
        benchmark::DoNotOptimize(c);
        state.PauseTiming();
        c.reset();
        state.ResumeTiming();
    }
}
#endif

static void BM_RegenUpdate(benchmark::State& state) {
    const Step s(static_cast<std::size_t>(state.range(0)), static_cast<std::size_t>(state.range(1)));
    render::StaticCache<double> cache;
    for (auto _ : state) {
        state.PauseTiming();
        std::optional<render::CookedScene<double>> prev, next;
        prev.emplace(render::lay_out(s.a, {}, nullptr, &cache));
        state.ResumeTiming();
        next.emplace(render::lay_out(std::move(*prev), s.b, {}, nullptr, &cache));
        benchmark::DoNotOptimize(next);
        state.PauseTiming();
        next.reset();
        prev.reset();
        state.ResumeTiming();
    }
}

static void BM_RegenFresh(benchmark::State& state) {
    const Step s(static_cast<std::size_t>(state.range(0)), static_cast<std::size_t>(state.range(1)));
    render::StaticCache<double> cache;
    for (auto _ : state) {
        std::optional<render::CookedScene<double>> next;
        next.emplace(render::lay_out(s.b, {}, nullptr, &cache));
        benchmark::DoNotOptimize(next);
        state.PauseTiming();
        next.reset();
        state.ResumeTiming();
    }
}

#define SPATIUM_REGEN_ARGS SPATIUM_INSTANCE_ARGS->Iterations(5)

BENCHMARK(BM_RegenCookCold) SPATIUM_REGEN_ARGS;
#ifndef SPATIUM_BENCH_WITHOUT_COOK_RESUME
BENCHMARK(BM_RegenCookResume) SPATIUM_REGEN_ARGS;
#endif
BENCHMARK(BM_RegenUpdate) SPATIUM_REGEN_ARGS;
BENCHMARK(BM_RegenFresh) SPATIUM_REGEN_ARGS;
