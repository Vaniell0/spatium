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

TEST_CASE("eigen_sym: 2x2 scalar multiple of identity (repeated eigenvalue)", "[spd]") {
    // Regression test: S=cI makes every direction an eigenvector, and the
    // old tie-break (match each root to whichever diagonal entry it's
    // closest to) picked the SAME vector for both slots here since both
    // roots are equally close to both (equal) diagonal entries -- producing
    // a rank-1 "eigenvector matrix" that silently reconstructs the wrong
    // matrix despite the eigenvalues themselves being correct.
    Matrix<double, 2, 2> S;
    S(0, 0) = 5.0; S(1, 1) = 5.0;
    auto eig = detail::eigen_sym(S);
    CHECK_THAT(eig.values[0], WithinAbs(5.0, 1e-9));
    CHECK_THAT(eig.values[1], WithinAbs(5.0, 1e-9));

    auto UtU = eig.vectors.transpose() * eig.vectors;
    check_matrix_close(UtU, Matrix<double, 2, 2>::identity()); // not rank-1

    Matrix<double, 2, 2> D;
    D(0, 0) = eig.values[0]; D(1, 1) = eig.values[1];
    check_matrix_close(eig.vectors * D * eig.vectors.transpose(), S);
}

TEST_CASE("eigen_sym: 3x3 scalar multiple of identity (triple repeated eigenvalue)", "[spd]") {
    Matrix<double, 3, 3> S;
    S(0, 0) = 7.0; S(1, 1) = 7.0; S(2, 2) = 7.0;
    auto eig = detail::eigen_sym(S);

    auto UtU = eig.vectors.transpose() * eig.vectors;
    check_matrix_close(UtU, Matrix<double, 3, 3>::identity(), 1e-8);

    Matrix<double, 3, 3> D;
    D(0, 0) = eig.values[0]; D(1, 1) = eig.values[1]; D(2, 2) = eig.values[2];
    check_matrix_close(eig.vectors * D * eig.vectors.transpose(), S, 1e-8);
}

TEST_CASE("eigen_sym: 3x3 axisymmetric case (eigenvalue multiplicity 2, not a scalar multiple of I)", "[spd]") {
    // diag(3,3,10) -- like an axisymmetric rigid body's inertia tensor.
    // Two rows of (S - 3I) are exactly zero and the third is nonzero:
    // rank 1, so the cross product of ANY two rows vanishes regardless of
    // which pair is picked -- the failure mode the largest-norm-row
    // fallback (not just the identity special case) exists for.
    Matrix<double, 3, 3> S;
    S(0, 0) = 3.0; S(1, 1) = 3.0; S(2, 2) = 10.0;
    auto eig = detail::eigen_sym(S);

    // A general cubic solver locates a genuine double root only to about
    // sqrt(machine epsilon), not machine epsilon -- a repeated root's
    // location is a square-root-sensitive function of the polynomial's
    // coefficients (standard numerical-analysis fact, not solve_cubic
    // imprecision), so the two "3.0" roots typically land ~1e-7 apart
    // rather than agreeing to 1e-8.
    std::array<double, 3> values{eig.values[0], eig.values[1], eig.values[2]};
    std::sort(values.begin(), values.end());
    CHECK_THAT(values[0], WithinAbs(3.0, 1e-6));
    CHECK_THAT(values[1], WithinAbs(3.0, 1e-6));
    CHECK_THAT(values[2], WithinAbs(10.0, 1e-8));

    auto UtU = eig.vectors.transpose() * eig.vectors;
    check_matrix_close(UtU, Matrix<double, 3, 3>::identity(), 1e-8); // not rank-deficient

    Matrix<double, 3, 3> D;
    D(0, 0) = eig.values[0]; D(1, 1) = eig.values[1]; D(2, 2) = eig.values[2];
    check_matrix_close(eig.vectors * D * eig.vectors.transpose(), S, 1e-6);
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

TEST_CASE("SPDAffineInvariant<2>/<3> satisfy RiemannianManifold, not Surface", "[spd]") {
    static_assert(Manifold<SPDAffineInvariant<2>>);
    static_assert(RiemannianManifold<SPDAffineInvariant<2>>);
    static_assert(!Surface<SPDAffineInvariant<2>>);
    static_assert(Manifold<SPDAffineInvariant<3>>);
    static_assert(RiemannianManifold<SPDAffineInvariant<3>>);
    static_assert(!Surface<SPDAffineInvariant<3>>);
    SUCCEED();
}

TEST_CASE("SPDAffineInvariant::contains: symmetric positive-definite only", "[spd]") {
    SPDAffineInvariant<2> space;

    Matrix<double, 2, 2> spd;
    spd(0, 0) = 6.5; spd(1, 1) = 6.5; spd(0, 1) = -2.5; spd(1, 0) = -2.5; // eigenvalues 4,9
    CHECK(space.contains(spd));

    Matrix<double, 2, 2> not_symmetric;
    not_symmetric(0, 0) = 1.0; not_symmetric(1, 1) = 1.0;
    not_symmetric(0, 1) = 1.0; not_symmetric(1, 0) = 0.0;
    CHECK_FALSE(space.contains(not_symmetric));

    Matrix<double, 2, 2> not_pd; // symmetric, one negative eigenvalue
    not_pd(0, 0) = 1.0; not_pd(1, 1) = -1.0;
    CHECK_FALSE(space.contains(not_pd));
}

TEST_CASE("SPDAffineInvariant: exp_map/log_map at the identity reduces to plain matrix log/exp", "[spd]") {
    // p = I makes p^{1/2} = p^{-1/2} = I, so the affine-invariant maps
    // collapse to the same eigendecomposition-based matrix log/exp
    // SPD<N,T> uses directly -- a real cross-check between the two classes,
    // not just internal self-consistency.
    SPDAffineInvariant<2> space;
    Matrix<double, 2, 2> I = Matrix<double, 2, 2>::identity();

    Matrix<double, 2, 2> q;
    q(0, 0) = 6.5; q(1, 1) = 6.5; q(0, 1) = -2.5; q(1, 0) = -2.5; // eigenvalues 4,9

    auto v = space.log_map(I, q);
    check_matrix_close(v, SPD<2>::matrix_log_sym(q));

    auto back = space.exp_map(I, v, 1.0);
    check_matrix_close(back, q);
}

TEST_CASE("SPDAffineInvariant<2>: exp_map/log_map roundtrip at a non-identity base point", "[spd]") {
    SPDAffineInvariant<2> space;
    Matrix<double, 2, 2> p;
    p(0, 0) = 3.0; p(1, 1) = 7.0; // diagonal, not identity

    Matrix<double, 2, 2> q;
    q(0, 0) = 6.5; q(1, 1) = 6.5; q(0, 1) = -2.5; q(1, 0) = -2.5; // rotated eigenbasis

    auto v = space.log_map(p, q);
    auto recovered = space.exp_map(p, v, 1.0);
    check_matrix_close(recovered, q, 1e-8);
}

TEST_CASE("SPDAffineInvariant<3>: exp_map/log_map roundtrip at a non-identity base point", "[spd]") {
    SPDAffineInvariant<3> space;
    Matrix<double, 3, 3> p;
    p(0, 0) = 6.5; p(1, 1) = 6.5; p(2, 2) = 25.0; p(0, 1) = -2.5; p(1, 0) = -2.5;

    Matrix<double, 3, 3> q;
    q(0, 0) = 3.0; q(1, 1) = 7.0; q(2, 2) = 12.0;

    auto v = space.log_map(p, q);
    auto recovered = space.exp_map(p, v, 1.0);
    check_matrix_close(recovered, q, 1e-7);
}

TEST_CASE("SPDAffineInvariant::distance is symmetric", "[spd]") {
    SPDAffineInvariant<2> space;
    Matrix<double, 2, 2> p;
    p(0, 0) = 3.0; p(1, 1) = 7.0;
    Matrix<double, 2, 2> q;
    q(0, 0) = 6.5; q(1, 1) = 6.5; q(0, 1) = -2.5; q(1, 0) = -2.5;

    CHECK_THAT(space.distance(p, q), WithinAbs(space.distance(q, p), 1e-8));
    CHECK_THAT(space.distance(p, p), WithinAbs(0.0, 1e-10));
}

TEST_CASE("SPDAffineInvariant::distance matches the generalized-eigenvalue closed form for a shared eigenbasis", "[spd]") {
    // p, q share the standard basis (both diagonal): the affine-invariant
    // distance has an independent, textbook closed form here --
    // sqrt(sum(log(mu_i/lambda_i)^2)) over the per-axis eigenvalue ratios --
    // that log_map()/metric_at() never explicitly compute. Cross-checking
    // against it (not just roundtrip self-consistency) is the same standard
    // this file already holds frechet_mean() to.
    SPDAffineInvariant<2> space;
    Matrix<double, 2, 2> p;
    p(0, 0) = 4.0; p(1, 1) = 9.0;
    Matrix<double, 2, 2> q;
    q(0, 0) = 16.0; q(1, 1) = 1.0;

    double expected = std::sqrt(std::pow(std::log(16.0 / 4.0), 2) + std::pow(std::log(1.0 / 9.0), 2));
    CHECK_THAT(space.distance(p, q), WithinAbs(expected, 1e-8));
}

TEST_CASE("frechet_mean_affine_invariant: shared-eigenbasis case matches per-eigenvalue geometric mean", "[spd]") {
    // When all samples commute (share an eigenbasis), the affine-invariant
    // Karcher mean has to agree with the log-Euclidean closed form: both
    // reduce to the per-eigenvalue geometric mean in that special case. This
    // exercises the actual iteration (no closed form exists here in
    // general, see spd.hpp) against an independently-derivable answer.
    Matrix<double, 2, 2> S1; S1(0, 0) = 4.0; S1(1, 1) = 9.0;
    Matrix<double, 2, 2> S2; S2(0, 0) = 16.0; S2(1, 1) = 1.0;
    std::array<Matrix<double, 2, 2>, 2> samples{S1, S2};

    auto mean = frechet_mean_affine_invariant<2>(std::span<const Matrix<double, 2, 2>>{samples});

    Matrix<double, 2, 2> expected;
    expected(0, 0) = std::sqrt(4.0 * 16.0); // 8
    expected(1, 1) = std::sqrt(9.0 * 1.0);  // 3
    check_matrix_close(mean, expected, 1e-8);
}

TEST_CASE("frechet_mean_affine_invariant: minimizes the sum of squared affine-invariant distances", "[spd]") {
    // Non-commuting samples (rotated relative to each other): no closed
    // form to check against directly, so verify the defining property
    // instead -- the returned mean is a local minimum of the sum of squared
    // distances, i.e. nudging it in any direction can't reduce the sum
    // (the first-order optimality condition sum_i log_map(mean, S_i) ~ 0).
    SPDAffineInvariant<2> space;
    Matrix<double, 2, 2> S1;
    S1(0, 0) = 6.5; S1(1, 1) = 6.5; S1(0, 1) = -2.5; S1(1, 0) = -2.5; // eigenvalues 4,9, rotated
    Matrix<double, 2, 2> S2;
    S2(0, 0) = 3.0; S2(1, 1) = 7.0; // axis-aligned
    std::array<Matrix<double, 2, 2>, 2> samples{S1, S2};

    auto mean = frechet_mean_affine_invariant<2>(std::span<const Matrix<double, 2, 2>>{samples});
    CHECK(space.contains(mean));

    Matrix<double, 2, 2> grad{};
    for (const auto& s : samples) grad = grad + space.log_map(mean, s);
    CHECK_THAT(space.metric_at(mean, grad, grad), WithinAbs(0.0, 1e-16));
}

TEST_CASE("Naive matrix addition can break SPD validity; SPDAffineInvariant::exp_map cannot", "[spd]") {
    // The enforced version of examples/spd_relationship_demo.cpp's claim:
    // repeatedly nudging a relationship matrix in a fixed "strengthening"
    // direction (as a sequence of similar causally-linked events would)
    // breaks positive-definiteness under naive matrix addition within a
    // bounded number of steps -- nothing about plain addition prevents an
    // eigenvalue from crossing zero -- but exp_map's retraction back onto
    // the manifold cannot produce an invalid matrix in exact arithmetic.
    // 14 steps keeps this within double-precision's healthy range: pushed
    // further, the same repeated raw nudge drives the point asymptotically
    // toward the cone's boundary and ordinary floating-point underflow
    // there eventually breaks contains() too -- a real precision
    // phenomenon near the boundary, not a violation of the structural
    // guarantee (see the demo's fuller note), just out of scope here.
    SPDAffineInvariant<2> space;
    Matrix<double, 2, 2> S_naive = Matrix<double, 2, 2>::identity();
    Matrix<double, 2, 2> S_geo = Matrix<double, 2, 2>::identity();
    Matrix<double, 2, 2> nudge;
    nudge(0, 1) = 0.08; nudge(1, 0) = 0.08;

    bool naive_broke = false;
    for (int step = 0; step < 14; ++step) {
        S_naive = S_naive + nudge;
        S_geo = space.exp_map(S_geo, nudge, 1.0);
        if (!space.contains(S_naive)) naive_broke = true;
        CHECK(space.contains(S_geo)); // never fails, by construction
    }
    CHECK(naive_broke); // naive DOES break somewhere in these 14 steps
}
