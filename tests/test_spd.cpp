#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/spaces/spd.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <span>

using namespace spatium;
using Catch::Matchers::WithinAbs;

namespace {

template<std::size_t N>
void check_matrix_close(const Matrix<double, N, N>& a, const Matrix<double, N, N>& b,
                         double tol = 1e-9) {
    for (std::size_t i = 0; i < a.data.size(); ++i)
        CHECK_THAT(a.data[i], WithinAbs(b.data[i], tol));
}

} // namespace

TEST_CASE("SPD<2> satisfies RiemannianManifold, not Surface", "[spd]") {
    static_assert(Manifold<SPD<2>>);
    static_assert(RiemannianManifold<SPD<2>>);
    static_assert(!Surface<SPD<2>>); // no project()/normal() -- see spd.hpp
    SUCCEED();
}

TEST_CASE("SPD<3> satisfies RiemannianManifold, not Surface", "[spd]") {
    static_assert(Manifold<SPD<3>>);
    static_assert(RiemannianManifold<SPD<3>>);
    static_assert(!Surface<SPD<3>>);
    SUCCEED();
}

TEST_CASE("eigen_sym: 2x2 diagonal matrix", "[spd]") {
    Matrix<double, 2, 2> S;
    S(0, 0) = 3.0; S(1, 1) = 7.0;
    auto eig = detail::eigen_sym(S);

    // Values are {3,7} in some order; reconstruction must match regardless.
    auto D = Matrix<double, 2, 2>{};
    D(0, 0) = eig.values[0]; D(1, 1) = eig.values[1];
    auto reconstructed = eig.vectors * D * eig.vectors.transpose();
    check_matrix_close(reconstructed, S);
}

TEST_CASE("eigen_sym: 2x2 non-diagonal matrix, known eigenvalues {4,9}", "[spd]") {
    // S = R diag(4,9) R^T for R = [[c,-s],[s,c]], c=s=1/sqrt(2) (45 degree
    // rotation) -- worked out by hand: S(0,0)=S(1,1)=(4+9)/2=6.5,
    // S(0,1)=S(1,0)=(4-9)/2=-2.5.
    Matrix<double, 2, 2> S;
    S(0, 0) = 6.5; S(1, 1) = 6.5;
    S(0, 1) = -2.5; S(1, 0) = -2.5;

    auto eig = detail::eigen_sym(S);
    double lo = std::min(eig.values[0], eig.values[1]);
    double hi = std::max(eig.values[0], eig.values[1]);
    CHECK_THAT(lo, WithinAbs(4.0, 1e-9));
    CHECK_THAT(hi, WithinAbs(9.0, 1e-9));

    // Eigenvectors orthonormal: U^T U == I.
    auto UtU = eig.vectors.transpose() * eig.vectors;
    check_matrix_close(UtU, Matrix<double, 2, 2>::identity());

    // Reconstruction: U diag(values) U^T == S.
    Matrix<double, 2, 2> D;
    D(0, 0) = eig.values[0]; D(1, 1) = eig.values[1];
    auto reconstructed = eig.vectors * D * eig.vectors.transpose();
    check_matrix_close(reconstructed, S);
}

TEST_CASE("eigen_sym: 3x3 non-diagonal matrix, known eigenvalues {4,9,25}", "[spd]") {
    // Same 45-degree block rotation on axes (0,1), axis 2 untouched.
    Matrix<double, 3, 3> S;
    S(0, 0) = 6.5; S(1, 1) = 6.5; S(2, 2) = 25.0;
    S(0, 1) = -2.5; S(1, 0) = -2.5;

    auto eig = detail::eigen_sym(S);
    std::array<double, 3> values{eig.values[0], eig.values[1], eig.values[2]};
    std::sort(values.begin(), values.end());
    CHECK_THAT(values[0], WithinAbs(4.0, 1e-8));
    CHECK_THAT(values[1], WithinAbs(9.0, 1e-8));
    CHECK_THAT(values[2], WithinAbs(25.0, 1e-8));

    auto UtU = eig.vectors.transpose() * eig.vectors;
    check_matrix_close(UtU, Matrix<double, 3, 3>::identity(), 1e-8);

    Matrix<double, 3, 3> D;
    D(0, 0) = eig.values[0]; D(1, 1) = eig.values[1]; D(2, 2) = eig.values[2];
    auto reconstructed = eig.vectors * D * eig.vectors.transpose();
    check_matrix_close(reconstructed, S, 1e-8);
}

TEST_CASE("SPD<2>::to_spd(from_spd(S)) roundtrip, diagonal and non-diagonal", "[spd]") {
    Matrix<double, 2, 2> S1;
    S1(0, 0) = 3.0; S1(1, 1) = 7.0;

    Matrix<double, 2, 2> S2;
    S2(0, 0) = 6.5; S2(1, 1) = 6.5; S2(0, 1) = -2.5; S2(1, 0) = -2.5;

    for (auto& S : {S1, S2}) {
        auto p = SPD<2>::from_spd(S);
        auto S_back = SPD<2>::to_spd(p);
        check_matrix_close(S_back, S);
    }
}

TEST_CASE("SPD<3>::to_spd(from_spd(S)) roundtrip", "[spd]") {
    Matrix<double, 3, 3> S;
    S(0, 0) = 6.5; S(1, 1) = 6.5; S(2, 2) = 25.0;
    S(0, 1) = -2.5; S(1, 0) = -2.5;

    auto p = SPD<3>::from_spd(S);
    auto S_back = SPD<3>::to_spd(p);
    check_matrix_close(S_back, S, 1e-8);
}

TEST_CASE("SPD<2> exp_map/log_map roundtrip on PointType", "[spd]") {
    // Flat by construction (see spd.hpp) -- this exercises the concept-
    // required signatures, not curvature.
    SPD<2> space;
    SPD<2>::PointType p{0.3, -0.2, 0.5};
    SPD<2>::PointType q{1.1, 0.4, -0.7};
    auto v = space.log_map(p, q);
    auto recovered = space.exp_map(p, v, 1.0);
    for (std::size_t i = 0; i < SPD<2>::M; ++i)
        CHECK_THAT(recovered[i], WithinAbs(q[i], 1e-12));
}

TEST_CASE("frechet_mean: log-Euclidean mean of diagonal SPD matrices is the per-entry geometric mean", "[spd]") {
    // No optimization involved -- see spd.hpp's file comment: log-Euclidean
    // Frechet mean is closed form, average logs then exponentiate.
    Matrix<double, 2, 2> S1; S1(0, 0) = 4.0; S1(1, 1) = 9.0;
    Matrix<double, 2, 2> S2; S2(0, 0) = 16.0; S2(1, 1) = 1.0;
    std::array<Matrix<double, 2, 2>, 2> samples{S1, S2};

    auto mean = frechet_mean<2>(std::span<const Matrix<double, 2, 2>>{samples});

    Matrix<double, 2, 2> expected;
    expected(0, 0) = std::sqrt(4.0 * 16.0); // 8
    expected(1, 1) = std::sqrt(9.0 * 1.0);  // 3
    check_matrix_close(mean, expected);
}

TEST_CASE("frechet_mean: shared-eigenbasis case matches rotating the diagonal geometric mean", "[spd]") {
    // S1, S2 share the same 45-degree eigenbasis as the eigen_sym tests
    // above, with eigenvalues (4,9) and (16,1) respectively. Expected mean
    // eigenvalues: sqrt(4*16)=8, sqrt(9*1)=3, folded back through the same
    // rotation: p=(8+3)/2=5.5, q=(8-3)/2=2.5.
    Matrix<double, 2, 2> S1;
    S1(0, 0) = 6.5; S1(1, 1) = 6.5; S1(0, 1) = -2.5; S1(1, 0) = -2.5;
    Matrix<double, 2, 2> S2;
    S2(0, 0) = 8.5; S2(1, 1) = 8.5; S2(0, 1) = 7.5; S2(1, 0) = 7.5;
    std::array<Matrix<double, 2, 2>, 2> samples{S1, S2};

    auto mean = frechet_mean<2>(std::span<const Matrix<double, 2, 2>>{samples});

    Matrix<double, 2, 2> expected;
    expected(0, 0) = 5.5; expected(1, 1) = 5.5;
    expected(0, 1) = 2.5; expected(1, 0) = 2.5;
    check_matrix_close(mean, expected, 1e-8);
}

TEST_CASE("frechet_mean: SPD<3> mean over three samples matches per-eigenvalue geometric mean", "[spd]") {
    Matrix<double, 3, 3> S1; S1(0, 0) = 4.0; S1(1, 1) = 9.0; S1(2, 2) = 25.0;
    Matrix<double, 3, 3> S2; S2(0, 0) = 16.0; S2(1, 1) = 1.0; S2(2, 2) = 49.0;
    std::array<Matrix<double, 3, 3>, 2> samples{S1, S2};

    auto mean = frechet_mean<3>(std::span<const Matrix<double, 3, 3>>{samples});

    Matrix<double, 3, 3> expected;
    expected(0, 0) = std::sqrt(4.0 * 16.0);  // 8
    expected(1, 1) = std::sqrt(9.0 * 1.0);   // 3
    expected(2, 2) = std::sqrt(25.0 * 49.0); // 35
    check_matrix_close(mean, expected, 1e-8);
}
