#include <benchmark/benchmark.h>
#include <spatium/physics/atomic/orbital.hpp>
#include <spatium/physics/atomic/atom_model.hpp>
#include <spatium/spaces/implicit.hpp>

using namespace spatium;
using namespace spatium::physics::atomic;

// ── Rejection sampling ───────────────────────────────────────

static void BM_SampleOrbital_1s_1k(benchmark::State& state) {
    for (auto _ : state)
        benchmark::DoNotOptimize(sample_orbital_points(1, 0, 0, 1000));
}
BENCHMARK(BM_SampleOrbital_1s_1k);

static void BM_SampleOrbital_1s_10k(benchmark::State& state) {
    for (auto _ : state)
        benchmark::DoNotOptimize(sample_orbital_points(1, 0, 0, 10000));
}
BENCHMARK(BM_SampleOrbital_1s_10k);

static void BM_SampleOrbital_2p_1k(benchmark::State& state) {
    for (auto _ : state)
        benchmark::DoNotOptimize(sample_orbital_points(2, 1, 0, 1000));
}
BENCHMARK(BM_SampleOrbital_2p_1k);

static void BM_SampleOrbital_3d_1k(benchmark::State& state) {
    for (auto _ : state)
        benchmark::DoNotOptimize(sample_orbital_points(3, 2, 0, 1000));
}
BENCHMARK(BM_SampleOrbital_3d_1k);

// ── Marching cubes ───────────────────────────────────────────

static void BM_MarchingCubes_1s(benchmark::State& state) {
    auto surf = make_orbital(1, 0, 0, 0.005);
    for (auto _ : state)
        benchmark::DoNotOptimize(marching_cubes(surf, 32));
}
BENCHMARK(BM_MarchingCubes_1s);

static void BM_MarchingCubes_2p(benchmark::State& state) {
    auto surf = make_orbital(2, 1, 0, 0.00125);
    for (auto _ : state)
        benchmark::DoNotOptimize(marching_cubes(surf, 32));
}
BENCHMARK(BM_MarchingCubes_2p);

// ── AtomModel build ──────────────────────────────────────────

static void BM_AtomModel_H(benchmark::State& state) {
    for (auto _ : state)
        benchmark::DoNotOptimize(AtomModel<>::build(1, 10000));
}
BENCHMARK(BM_AtomModel_H);

static void BM_AtomModel_C(benchmark::State& state) {
    for (auto _ : state)
        benchmark::DoNotOptimize(AtomModel<>::build(6, 5000));
}
BENCHMARK(BM_AtomModel_C);

static void BM_AtomModel_Fe(benchmark::State& state) {
    for (auto _ : state)
        benchmark::DoNotOptimize(AtomModel<>::build(26, 1000));
}
BENCHMARK(BM_AtomModel_Fe)->Iterations(1);
