#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/algebra/eigen_decomp.hpp>
#include <spatium/core/precision.hpp>
#include <array>
#include <cmath>

using namespace spatium;
using Catch::Matchers::WithinAbs;

namespace {

template<std::size_t N>
void check_matrix_close(const Matrix<double, N, N>& a, const Matrix<double, N, N>& b,
                         double tol = 1e-9) {
    for (std::size_t i = 0; i < a.data.size(); ++i)
        CHECK_THAT(a.data[i], WithinAbs(b.data[i], tol));
}

template<std::size_t N>
void check_orthonormal(const Matrix<double, N, N>& V, double tol = 1e-9) {
    check_matrix_close(V.transpose() * V, Matrix<double, N, N>::identity(), tol);
}

template<std::size_t N>
void check_ascending(const Vec<double, N>& values) {
    for (std::size_t i = 1; i < N; ++i) CHECK(values[i - 1] <= values[i] + 1e-9);
}

} // namespace

TEST_CASE("eigen_decompose: diagonal matrix returns the diagonal, sorted ascending", "[eigen_decomp]") {
    Matrix<double, 4, 4> S;
    S(0, 0) = 4; S(1, 1) = 1; S(2, 2) = 3; S(3, 3) = 2;

    auto eig = eigen_decompose(S);
    check_ascending(eig.values);
    CHECK_THAT(eig.values[0], WithinAbs(1.0, 1e-9));
    CHECK_THAT(eig.values[1], WithinAbs(2.0, 1e-9));
    CHECK_THAT(eig.values[2], WithinAbs(3.0, 1e-9));
    CHECK_THAT(eig.values[3], WithinAbs(4.0, 1e-9));
    check_orthonormal(eig.vectors);
    check_matrix_close(reconstruct(eig), S);
}

TEST_CASE("eigen_decompose: 3x3 rotated block, known eigenvalues {4,9,25}", "[eigen_decomp]") {
    // S = R diag(4,9,25) R^T for R a 45-degree rotation on axes (0,1),
    // axis 2 untouched -- worked out by hand, same construction as
    // spaces/spd.hpp's own eigen_sym test for the identical matrix:
    // S(0,0)=S(1,1)=(4+9)/2=6.5, S(0,1)=S(1,0)=(4-9)/2=-2.5.
    Matrix<double, 3, 3> S;
    S(0, 0) = 6.5; S(1, 1) = 6.5; S(2, 2) = 25.0;
    S(0, 1) = -2.5; S(1, 0) = -2.5;

    auto eig = eigen_decompose(S);
    check_ascending(eig.values);
    CHECK_THAT(eig.values[0], WithinAbs(4.0, 1e-9));
    CHECK_THAT(eig.values[1], WithinAbs(9.0, 1e-9));
    CHECK_THAT(eig.values[2], WithinAbs(25.0, 1e-9));
    check_orthonormal(eig.vectors);
    check_matrix_close(reconstruct(eig), S);
}

TEST_CASE("eigen_decompose: 4x4 two independent 45-degree blocks, known eigenvalues {2,8,10,20}", "[eigen_decomp]") {
    // Block on axes (0,1): eigenvalues {2,8} -> diag=(2+8)/2=5, off-diag=(2-8)/2=-3.
    // Block on axes (2,3): eigenvalues {10,20} -> diag=(10+20)/2=15, off-diag=(10-20)/2=-5.
    // Blocks are independent (zero coupling), so the full spectrum is the
    // union of both: {2,8,10,20}.
    Matrix<double, 4, 4> S;
    S(0, 0) = 5.0; S(1, 1) = 5.0; S(0, 1) = -3.0; S(1, 0) = -3.0;
    S(2, 2) = 15.0; S(3, 3) = 15.0; S(2, 3) = -5.0; S(3, 2) = -5.0;

    auto eig = eigen_decompose(S);
    check_ascending(eig.values);
    CHECK_THAT(eig.values[0], WithinAbs(2.0, 1e-9));
    CHECK_THAT(eig.values[1], WithinAbs(8.0, 1e-9));
    CHECK_THAT(eig.values[2], WithinAbs(10.0, 1e-9));
    CHECK_THAT(eig.values[3], WithinAbs(20.0, 1e-9));
    check_orthonormal(eig.vectors);
    check_matrix_close(reconstruct(eig), S);
}

TEST_CASE("eigen_decompose: repeated eigenvalue diag(3,3,10)", "[eigen_decomp]") {
    // Same shape as the axisymmetric-inertia-tensor case that broke
    // spaces/spd.hpp's detail::eigen_sym (see docs/ROADMAP.md's "SPD(n):
    // affine-invariant metric" entry) -- Jacobi doesn't share that failure
    // mode (it never solves a null space, only accumulates rotations), but
    // this is still the canonical repeated-eigenvalue regression case to
    // check orthonormality/reconstruction stay correct under it.
    Matrix<double, 3, 3> S;
    S(0, 0) = 3.0; S(1, 1) = 3.0; S(2, 2) = 10.0;

    auto eig = eigen_decompose(S);
    CHECK_THAT(eig.values[0], WithinAbs(3.0, 1e-9));
    CHECK_THAT(eig.values[1], WithinAbs(3.0, 1e-9));
    CHECK_THAT(eig.values[2], WithinAbs(10.0, 1e-9));
    check_orthonormal(eig.vectors);
    check_matrix_close(reconstruct(eig), S);
}

TEST_CASE("eigen_decompose: repeated eigenvalue rotated into a non-diagonal matrix", "[eigen_decomp]") {
    // S = R diag(3,3,10) R^T for a Givens rotation R mixing axes (0,2) by
    // 30 degrees (c=sqrt(3)/2, s=1/2); axis 1 is left alone at its own
    // eigenvalue 3. Worked out by hand from S = sum_k d_k * col_k(R) *
    // col_k(R)^T with col0(R)=(c,0,s), col1(R)=(0,1,0), col2(R)=(-s,0,c):
    //   S(0,0) = 3c^2 + 10s^2, S(2,2) = 3s^2 + 10c^2,
    //   S(0,2) = (3-10)*s*c,   S(1,1) = 3, all other entries 0.
    // Same {3,3,10} spectrum as the diagonal case above, but genuinely
    // non-diagonal -- exercises Jacobi's rotations under a repeated
    // eigenvalue instead of converging immediately because every
    // off-diagonal entry already started at zero.
    using std::sqrt;
    double c = sqrt(3.0) / 2.0, s = 0.5;
    Matrix<double, 3, 3> S;
    S(0, 0) = 3.0 * c * c + 10.0 * s * s;
    S(2, 2) = 3.0 * s * s + 10.0 * c * c;
    S(0, 2) = (3.0 - 10.0) * s * c;
    S(2, 0) = S(0, 2);
    S(1, 1) = 3.0;

    auto eig = eigen_decompose(S);
    CHECK_THAT(eig.values[0], WithinAbs(3.0, 1e-9));
    CHECK_THAT(eig.values[1], WithinAbs(3.0, 1e-9));
    CHECK_THAT(eig.values[2], WithinAbs(10.0, 1e-9));
    check_orthonormal(eig.vectors);
    check_matrix_close(reconstruct(eig), S);
}

TEST_CASE("eigen_decompose: N=5 generic symmetric matrix, reconstruction + orthonormality", "[eigen_decomp]") {
    // No analytically-known eigenbasis here -- this exercises N beyond the
    // N=2,3 that spaces/spd.hpp's closed-form eigen_sym is capped at, via
    // self-consistency (reconstruction + orthonormality) instead.
    double vals[5][5] = {
        {4, 1, 0, 2, 0},
        {1, 3, 1, 0, 1},
        {0, 1, 5, 1, 0},
        {2, 0, 1, 6, 1},
        {0, 1, 0, 1, 2},
    };
    Matrix<double, 5, 5> S;
    for (std::size_t i = 0; i < 5; ++i)
        for (std::size_t j = 0; j < 5; ++j)
            S(i, j) = vals[i][j];

    auto eig = eigen_decompose(S);
    check_ascending(eig.values);
    check_orthonormal(eig.vectors, 1e-8);
    check_matrix_close(reconstruct(eig), S, 1e-8);
}

TEST_CASE("eigen_decompose works with Real50 (generic Scalar path)", "[eigen_decomp][precision]") {
    using T = Real50;
    Matrix<T, 3, 3> S;
    S(0, 0) = T{6.5}; S(1, 1) = T{6.5}; S(2, 2) = T{25};
    S(0, 1) = T{-2.5}; S(1, 0) = T{-2.5};

    auto eig = eigen_decompose(S);
    CHECK(approx_equal(eig.values[0], T{4}));
    CHECK(approx_equal(eig.values[1], T{9}));
    CHECK(approx_equal(eig.values[2], T{25}));

    auto rec = reconstruct(eig);
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            CHECK(approx_equal(rec(i, j), S(i, j)));
}

#if SPATIUM_HAS_EIGEN
#include <spatium/algebra/eigen_interop.hpp>
#include <Eigen/Eigenvalues>
#include <random>

TEST_CASE("eigen_decompose cross-validated against Eigen::SelfAdjointEigenSolver", "[eigen_decomp][eigen]") {
    // Independent-ground-truth cross-check, not just self-consistency --
    // see docs/ROADMAP.md's SE3::exp() bug-catch for why this matters in
    // this codebase: a Jacobi implementation bug that still happens to
    // reconstruct S from a valid-looking orthogonal V could slip past a
    // reconstruction-only check, but not past an independently-computed
    // eigenvalue list.
    constexpr std::size_t N = 6;
    std::mt19937 rng(2026);
    std::uniform_real_distribution<double> d(-5.0, 5.0);

    for (int trial = 0; trial < 20; ++trial) {
        Matrix<double, N, N> S;
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = i; j < N; ++j) {
                double v = d(rng);
                S(i, j) = v;
                S(j, i) = v;
            }

        auto eig = eigen_decompose(S);

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, static_cast<int>(N), static_cast<int>(N)>>
            solver(to_eigen(S));
        auto evals = solver.eigenvalues(); // ascending, matching our convention

        for (std::size_t i = 0; i < N; ++i)
            CHECK_THAT(eig.values[i], WithinAbs(evals(static_cast<int>(i)), 1e-7));
    }
}
#endif
