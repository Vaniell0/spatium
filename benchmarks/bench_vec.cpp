#include <benchmark/benchmark.h>
#include <spatium/algebra/vector.hpp>
#include <random>

using namespace spatium;

template<typename T, std::size_t N>
static Vec<T, N> random_vec(std::mt19937& rng) {
    std::uniform_real_distribution<T> dist(T(-1), T(1));
    Vec<T, N> v;
    for (std::size_t i = 0; i < N; ++i) v[i] = dist(rng);
    return v;
}

// ── Dot product ───────────────────────────────────────────────

template<typename T, std::size_t N>
static void BM_VecDot(benchmark::State& state) {
    std::mt19937 rng(42);
    auto a = random_vec<T, N>(rng), b = random_vec<T, N>(rng);
    for (auto _ : state) {
        benchmark::DoNotOptimize(a.dot(b));
        a[0] += T(1e-15);
    }
}

BENCHMARK(BM_VecDot<float, 3>);
BENCHMARK(BM_VecDot<float, 4>);
BENCHMARK(BM_VecDot<double, 3>);
BENCHMARK(BM_VecDot<double, 4>);

// ── Norm ──────────────────────────────────────────────────────

template<typename T, std::size_t N>
static void BM_VecNorm(benchmark::State& state) {
    std::mt19937 rng(42);
    auto a = random_vec<T, N>(rng);
    for (auto _ : state) {
        benchmark::DoNotOptimize(a.norm());
        a[0] += T(1e-15);
    }
}

BENCHMARK(BM_VecNorm<float, 3>);
BENCHMARK(BM_VecNorm<float, 4>);
BENCHMARK(BM_VecNorm<double, 3>);
BENCHMARK(BM_VecNorm<double, 4>);

// ── Cross product ─────────────────────────────────────────────

template<typename T>
static void BM_VecCross(benchmark::State& state) {
    std::mt19937 rng(42);
    auto a = random_vec<T, 3>(rng), b = random_vec<T, 3>(rng);
    for (auto _ : state) {
        Vec<T, 3> c = a.cross(b);
        benchmark::DoNotOptimize(c);
        a[0] += T(1e-15);
    }
}

BENCHMARK(BM_VecCross<float>);
BENCHMARK(BM_VecCross<double>);

// ── Normalize ─────────────────────────────────────────────────

template<typename T, std::size_t N>
static void BM_VecNormalize(benchmark::State& state) {
    std::mt19937 rng(42);
    auto a = random_vec<T, N>(rng);
    for (auto _ : state) {
        Vec<T, N> n = a.normalized();
        benchmark::DoNotOptimize(n);
        a[0] += T(1e-15);
    }
}

BENCHMARK(BM_VecNormalize<float, 3>);
BENCHMARK(BM_VecNormalize<double, 3>);

// ── Expression Template chain: a + b + c + d ──────────────────

template<typename T, std::size_t N>
static void BM_VecAddChain4(benchmark::State& state) {
    std::mt19937 rng(42);
    auto a = random_vec<T, N>(rng), b = random_vec<T, N>(rng);
    auto c = random_vec<T, N>(rng), d = random_vec<T, N>(rng);
    for (auto _ : state) {
        Vec<T, N> r = a + b + c + d;
        benchmark::DoNotOptimize(r);
        a[0] += T(1e-15);
    }
}

BENCHMARK(BM_VecAddChain4<float, 3>);
BENCHMARK(BM_VecAddChain4<float, 4>);
BENCHMARK(BM_VecAddChain4<double, 3>);
BENCHMARK(BM_VecAddChain4<double, 4>);

// ── Mixed expression: (a - b).dot(c) ─────────────────────────

template<typename T, std::size_t N>
static void BM_VecExprDot(benchmark::State& state) {
    std::mt19937 rng(42);
    auto a = random_vec<T, N>(rng), b = random_vec<T, N>(rng), c = random_vec<T, N>(rng);
    for (auto _ : state) {
        benchmark::DoNotOptimize((a - b).dot(c));
        a[0] += T(1e-15);
    }
}

BENCHMARK(BM_VecExprDot<float, 3>);
BENCHMARK(BM_VecExprDot<double, 3>);
