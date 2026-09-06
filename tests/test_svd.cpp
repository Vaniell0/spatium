#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/algebra/svd.hpp>
#include <spatium/core/precision.hpp>
#include <cmath>

using namespace spatium;
using Catch::Matchers::WithinAbs;

namespace {

template<std::size_t M, std::size_t N>
void check_matrix_close(const Matrix<double, M, N>& a, const Matrix<double, M, N>& b,
                         double tol = 1e-9) {
    for (std::size_t i = 0; i < a.data.size(); ++i)
        CHECK_THAT(a.data[i], WithinAbs(b.data[i], tol));
}

template<std::size_t R, std::size_t K>
void check_orthonormal_columns(const Matrix<double, R, K>& U, double tol = 1e-9) {
    check_matrix_close(U.transpose() * U, Matrix<double, K, K>::identity(), tol);
}

} // namespace

TEST_CASE("svd: diagonal-like 3x2 matrix, hand-verifiable singular values {4,3}", "[svd]") {
    Matrix<double, 3, 2> A;
    A(0, 0) = 3.0;
    A(1, 1) = 4.0;
    // A(2, *) row is all zero -- singular values are exactly the column
    // norms {3,4} since the columns are already orthogonal.

    auto r = svd(A);
    CHECK_THAT(r.singular_values[0], WithinAbs(4.0, 1e-9));
    CHECK_THAT(r.singular_values[1], WithinAbs(3.0, 1e-9));
    check_orthonormal_columns(r.U);
    check_orthonormal_columns(r.V);
    check_matrix_close(reconstruct(r), A);
}

TEST_CASE("svd: 2x2 with a negative diagonal entry -- sign absorbed into U, not the singular value", "[svd]") {
    // Singular values are always >= 0 even though A itself has a negative
    // entry: sigma = {3,2} (not {-3,2}), with the sign folded into U's
    // corresponding column instead.
    Matrix<double, 2, 2> A;
    A(0, 0) = 2.0;
    A(1, 1) = -3.0;

    auto r = svd(A);
    CHECK_THAT(r.singular_values[0], WithinAbs(3.0, 1e-9));
    CHECK_THAT(r.singular_values[1], WithinAbs(2.0, 1e-9));
    check_orthonormal_columns(r.U);
    check_orthonormal_columns(r.V);
    check_matrix_close(reconstruct(r), A);
}

TEST_CASE("svd: wide 2x4 diagonal-like matrix, hand-verifiable singular values {5,3}", "[svd]") {
    Matrix<double, 2, 4> A;
    A(0, 0) = 5.0;
    A(1, 1) = 3.0;

    auto r = svd(A);
    CHECK(SVDResult<double, 2, 4>::K == 2);
    CHECK_THAT(r.singular_values[0], WithinAbs(5.0, 1e-9));
    CHECK_THAT(r.singular_values[1], WithinAbs(3.0, 1e-9));
    check_orthonormal_columns(r.U);
    check_orthonormal_columns(r.V);
    check_matrix_close(reconstruct(r), A);
}

TEST_CASE("svd: symmetric positive-definite square A has U == V and singular values == eigenvalues", "[svd]") {
    // For symmetric PD A, A and A^T A share the same eigenvectors (A^T A =
    // A^2), and A's eigenvalues are already the singular values -- an
    // independently-derivable cross-check of the whole U/sigma/V
    // construction against spaces/spd.hpp's own 45-degree-rotated
    // {4,9,25} test case, using a completely different entry point
    // (svd() instead of eigen_decompose()/detail::eigen_sym()).
    Matrix<double, 3, 3> A;
    A(0, 0) = 6.5; A(1, 1) = 6.5; A(2, 2) = 25.0;
    A(0, 1) = -2.5; A(1, 0) = -2.5;

    auto r = svd(A);
    CHECK_THAT(r.singular_values[0], WithinAbs(25.0, 1e-9));
    CHECK_THAT(r.singular_values[1], WithinAbs(9.0, 1e-9));
    CHECK_THAT(r.singular_values[2], WithinAbs(4.0, 1e-9));
    check_orthonormal_columns(r.U);
    check_orthonormal_columns(r.V);
    check_matrix_close(reconstruct(r), A);

    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t k = 0; k < 3; ++k)
            CHECK_THAT(r.U(i, k), WithinAbs(r.V(i, k), 1e-9));
}

TEST_CASE("svd: tall 5x3 generic matrix, reconstruction + orthonormality", "[svd]") {
    double vals[5][3] = {
        {1, 2, 0},
        {0, 1, 3},
        {4, 0, 1},
        {2, 2, 2},
        {1, 0, 5},
    };
    Matrix<double, 5, 3> A;
    for (std::size_t i = 0; i < 5; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            A(i, j) = vals[i][j];

    auto r = svd(A);
    // Descending order.
    CHECK(r.singular_values[0] >= r.singular_values[1]);
    CHECK(r.singular_values[1] >= r.singular_values[2]);
    check_orthonormal_columns(r.U, 1e-8);
    check_orthonormal_columns(r.V, 1e-8);
    check_matrix_close(reconstruct(r), A, 1e-8);
}

TEST_CASE("svd: rank-deficient 4x3 matrix -- one column a linear combination of the others", "[svd]") {
    // col2 = col0 + col1 exactly, so rank(A) == 2 < min(4,3) == 3: the
    // smallest singular value must come back ~0, and U's corresponding
    // column must still complete an orthonormal basis via the
    // Gram-Schmidt fallback rather than dividing by ~0.
    Matrix<double, 4, 3> A;
    double col0[4] = {1, 2, 3, 4};
    double col1[4] = {0, 1, 0, 1};
    for (std::size_t i = 0; i < 4; ++i) {
        A(i, 0) = col0[i];
        A(i, 1) = col1[i];
        A(i, 2) = col0[i] + col1[i];
    }

    auto r = svd(A);
    CHECK(r.singular_values[0] > 1e-6);
    CHECK(r.singular_values[1] > 1e-6);
    CHECK_THAT(r.singular_values[2], WithinAbs(0.0, 1e-8));

    check_orthonormal_columns(r.U, 1e-8);
    check_orthonormal_columns(r.V, 1e-8);
    check_matrix_close(reconstruct(r), A, 1e-7);
}

TEST_CASE("svd works with Real50 (generic Scalar path)", "[svd][precision]") {
    using T = Real50;
    Matrix<T, 3, 2> A;
    A(0, 0) = T{3};
    A(1, 1) = T{4};

    auto r = svd(A);
    CHECK(approx_equal(r.singular_values[0], T{4}));
    CHECK(approx_equal(r.singular_values[1], T{3}));

    auto rec = reconstruct(r);
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 2; ++j)
            CHECK(approx_equal(rec(i, j), A(i, j)));
}

#if SPATIUM_HAS_EIGEN
#include <spatium/algebra/eigen_interop.hpp>
#include <Eigen/SVD>
#include <random>

TEST_CASE("svd cross-validated against Eigen::JacobiSVD", "[svd][eigen]") {
    // Eigen::JacobiSVD uses two-sided Jacobi/bidiagonalization directly on
    // A, a genuinely different numerical path than this file's
    // eigendecompose(A^T A) construction -- a real independent check, not
    // a restatement of the same computation (see this file's own
    // condition-number caveat in svd.hpp's file comment for why the two
    // approaches aren't expected to be bit-identical, just close).
    constexpr std::size_t M = 6, N = 4;
    std::mt19937 rng(2026);
    std::uniform_real_distribution<double> d(-5.0, 5.0);

    for (int trial = 0; trial < 20; ++trial) {
        Matrix<double, M, N> A;
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t j = 0; j < N; ++j)
                A(i, j) = d(rng);

        auto r = svd(A);

        Eigen::JacobiSVD<Eigen::Matrix<double, static_cast<int>(M), static_cast<int>(N)>> solver(
            to_eigen(A), Eigen::ComputeThinU | Eigen::ComputeThinV);
        auto svals = solver.singularValues(); // descending, matching our convention

        for (std::size_t k = 0; k < SVDResult<double, M, N>::K; ++k)
            CHECK_THAT(r.singular_values[k], WithinAbs(svals(static_cast<int>(k)), 1e-6));
    }
}
#endif
