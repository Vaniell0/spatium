// Benchmarks for the native cyclic-Jacobi eigendecomposition
// (algebra/eigen_decomp.hpp) and the A^T A-based SVD (algebra/svd.hpp),
// compared against Eigen's SelfAdjointEigenSolver / JacobiSVD at matching
// sizes when SPATIUM_EIGEN is available. This project measures and
// discloses real comparative numbers rather than asserting "native is
// fine" without one -- see bench_raycast.cpp's BVH-vs-brute-force-vs-
// analytical style for the house convention this follows.

#include <benchmark/benchmark.h>
#include <spatium/algebra/eigen_decomp.hpp>
#include <spatium/algebra/svd.hpp>
#include <random>

#if SPATIUM_HAS_EIGEN
#  include <spatium/algebra/eigen_interop.hpp>
#  include <Eigen/Eigenvalues>
#  include <Eigen/SVD>
#endif

using namespace spatium;

namespace {

template<std::size_t N>
Matrix<double, N, N> random_symmetric(unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> d(-5.0, 5.0);
    Matrix<double, N, N> S;
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = i; j < N; ++j) {
            double v = d(rng);
            S(i, j) = v;
            S(j, i) = v;
        }
    return S;
}

template<std::size_t M, std::size_t N>
Matrix<double, M, N> random_matrix(unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> d(-5.0, 5.0);
    Matrix<double, M, N> A;
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t j = 0; j < N; ++j)
            A(i, j) = d(rng);
    return A;
}

} // namespace

// ── Native eigen_decompose() ──────────────────────────────────

static void BM_EigenDecompose_Native_4(benchmark::State& state) {
    auto S = random_symmetric<4>(2026);
    for (auto _ : state) {
        auto eig = eigen_decompose(S);
        benchmark::DoNotOptimize(eig);
    }
}
BENCHMARK(BM_EigenDecompose_Native_4);

static void BM_EigenDecompose_Native_8(benchmark::State& state) {
    auto S = random_symmetric<8>(2026);
    for (auto _ : state) {
        auto eig = eigen_decompose(S);
        benchmark::DoNotOptimize(eig);
    }
}
BENCHMARK(BM_EigenDecompose_Native_8);

static void BM_EigenDecompose_Native_16(benchmark::State& state) {
    auto S = random_symmetric<16>(2026);
    for (auto _ : state) {
        auto eig = eigen_decompose(S);
        benchmark::DoNotOptimize(eig);
    }
}
BENCHMARK(BM_EigenDecompose_Native_16);

// ── Native svd() ───────────────────────────────────────────────

static void BM_SVD_Native_8x4(benchmark::State& state) {
    auto A = random_matrix<8, 4>(2026);
    for (auto _ : state) {
        auto r = svd(A);
        benchmark::DoNotOptimize(r);
    }
}
BENCHMARK(BM_SVD_Native_8x4);

static void BM_SVD_Native_16x8(benchmark::State& state) {
    auto A = random_matrix<16, 8>(2026);
    for (auto _ : state) {
        auto r = svd(A);
        benchmark::DoNotOptimize(r);
    }
}
BENCHMARK(BM_SVD_Native_16x8);

#if SPATIUM_HAS_EIGEN

// ── Eigen::SelfAdjointEigenSolver, same sizes ─────────────────

static void BM_EigenDecompose_Eigen_4(benchmark::State& state) {
    auto em = to_eigen(random_symmetric<4>(2026));
    for (auto _ : state) {
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 4, 4>> solver(em);
        auto evals = solver.eigenvalues();
        benchmark::DoNotOptimize(evals);
    }
}
BENCHMARK(BM_EigenDecompose_Eigen_4);

static void BM_EigenDecompose_Eigen_8(benchmark::State& state) {
    auto em = to_eigen(random_symmetric<8>(2026));
    for (auto _ : state) {
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 8, 8>> solver(em);
        auto evals = solver.eigenvalues();
        benchmark::DoNotOptimize(evals);
    }
}
BENCHMARK(BM_EigenDecompose_Eigen_8);

static void BM_EigenDecompose_Eigen_16(benchmark::State& state) {
    auto em = to_eigen(random_symmetric<16>(2026));
    for (auto _ : state) {
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 16, 16>> solver(em);
        auto evals = solver.eigenvalues();
        benchmark::DoNotOptimize(evals);
    }
}
BENCHMARK(BM_EigenDecompose_Eigen_16);

// ── Eigen::JacobiSVD, same sizes ───────────────────────────────

static void BM_SVD_Eigen_8x4(benchmark::State& state) {
    auto em = to_eigen(random_matrix<8, 4>(2026));
    for (auto _ : state) {
        Eigen::JacobiSVD<Eigen::Matrix<double, 8, 4>> solver(
            em, Eigen::ComputeThinU | Eigen::ComputeThinV);
        auto svals = solver.singularValues();
        benchmark::DoNotOptimize(svals);
    }
}
BENCHMARK(BM_SVD_Eigen_8x4);

static void BM_SVD_Eigen_16x8(benchmark::State& state) {
    auto em = to_eigen(random_matrix<16, 8>(2026));
    for (auto _ : state) {
        Eigen::JacobiSVD<Eigen::Matrix<double, 16, 8>> solver(
            em, Eigen::ComputeThinU | Eigen::ComputeThinV);
        auto svals = solver.singularValues();
        benchmark::DoNotOptimize(svals);
    }
}
BENCHMARK(BM_SVD_Eigen_16x8);

#endif // SPATIUM_HAS_EIGEN
