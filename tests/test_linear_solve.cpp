#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/algebra/linear_solve.hpp>
#include <cmath>

using namespace spatium;
using Catch::Matchers::WithinAbs;

TEST_CASE("solve_direct recovers a known solution exactly", "[linear_solve]") {
    Matrix<double, 3, 3> A;
    A(0, 0) = 4; A(0, 1) = 1; A(0, 2) = 2;
    A(1, 0) = 3; A(1, 1) = 6; A(1, 2) = 1;
    A(2, 0) = 2; A(2, 1) = 1; A(2, 2) = 5;
    Vec<double, 3> x_true{1.0, -2.0, 3.0};
    Vec<double, 3> b = A * x_true;

    auto result = solve_direct(A, b);
    REQUIRE(result.has_value());
    CHECK_THAT((*result)[0], WithinAbs(x_true[0], 1e-9));
    CHECK_THAT((*result)[1], WithinAbs(x_true[1], 1e-9));
    CHECK_THAT((*result)[2], WithinAbs(x_true[2], 1e-9));
}

TEST_CASE("solve_direct reports a singular matrix instead of dividing by zero", "[linear_solve]") {
    Matrix<double, 3, 3> A;
    A(0, 0) = 1; A(0, 1) = 2; A(0, 2) = 3;
    A(1, 0) = 2; A(1, 1) = 4; A(1, 2) = 6;  // row 1 = 2 * row 0 -> singular
    A(2, 0) = 1; A(2, 1) = 0; A(2, 2) = 1;
    Vec<double, 3> b{1.0, 2.0, 3.0};

    auto result = solve_direct(A, b);
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error().code == ErrorCode::SingularMatrix);
}

TEST_CASE("solve_jacobi converges to the exact solution for a strongly diagonally dominant matrix",
          "[linear_solve]") {
    Matrix<double, 4, 4> A;
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 4; ++j)
            A(i, j) = (i == j) ? 10.0 : 0.5;  // |diag|=10 vs. off-diag row sum=1.5 -> ratio ~6.7
    Vec<double, 4> x_true{1.0, 2.0, -1.0, 0.5};
    Vec<double, 4> b = A * x_true;

    Vec<double, 4> x = solve_jacobi(A, b, /*max_iter=*/30);
    for (std::size_t i = 0; i < 4; ++i) CHECK_THAT(x[i], WithinAbs(x_true[i], 1e-6));
}

TEST_CASE("solve_jacobi measurably fails within a small budget once the matrix stops being "
          "diagonally dominant, while solve_direct stays exact",
          "[linear_solve]") {
    // Real, checked claim (not assumed): off-diagonal row sum (3.6) exceeds
    // the diagonal (3.0) on every row -- diagonal_dominance_ratio() < 1,
    // outside Jacobi's guaranteed-convergence regime. This is exactly the
    // two-regime split rsc/include/linear_task.hpp's dispatch signal
    // depends on existing, mirroring rootfind_ops.hpp's Newton-vs-
    // bisection split.
    Matrix<double, 4, 4> A;
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 4; ++j)
            A(i, j) = (i == j) ? 3.0 : 1.2;
    Vec<double, 4> x_true{1.0, 2.0, -1.0, 0.5};
    Vec<double, 4> b = A * x_true;

    auto direct = solve_direct(A, b);
    REQUIRE(direct.has_value());
    for (std::size_t i = 0; i < 4; ++i) CHECK_THAT((*direct)[i], WithinAbs(x_true[i], 1e-9));

    Vec<double, 4> jacobi = solve_jacobi(A, b, /*max_iter=*/20);
    double max_err = 0.0;
    for (std::size_t i = 0; i < 4; ++i) max_err = std::max(max_err, std::abs(jacobi[i] - x_true[i]));
    CHECK(max_err > 0.1);  // did not converge in the budget
}

TEST_CASE("invert recovers A^-1 A = I", "[linear_solve]") {
    Matrix<double, 4, 4> A;
    A(0,0)=4; A(0,1)=1; A(0,2)=0; A(0,3)=2;
    A(1,0)=3; A(1,1)=6; A(1,2)=1; A(1,3)=0;
    A(2,0)=2; A(2,1)=1; A(2,2)=5; A(2,3)=1;
    A(3,0)=0; A(3,1)=2; A(3,2)=1; A(3,3)=7;

    auto result = invert(A);
    REQUIRE(result.has_value());
    Matrix<double, 4, 4> I = A * (*result);
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 4; ++j)
            CHECK_THAT(I(i, j), WithinAbs(i == j ? 1.0 : 0.0, 1e-9));
}

TEST_CASE("invert matches solve_direct column by column", "[linear_solve]") {
    Matrix<double, 3, 3> A;
    A(0,0)=2; A(0,1)=1; A(0,2)=1;
    A(1,0)=1; A(1,1)=3; A(1,2)=2;
    A(2,0)=1; A(2,1)=0; A(2,2)=4;

    auto inv = invert(A);
    REQUIRE(inv.has_value());
    for (std::size_t col = 0; col < 3; ++col) {
        Vec<double, 3> e{};
        e[col] = 1.0;
        auto sol = solve_direct(A, e);
        REQUIRE(sol.has_value());
        for (std::size_t row = 0; row < 3; ++row)
            CHECK_THAT((*inv)(row, col), WithinAbs((*sol)[row], 1e-9));
    }
}

TEST_CASE("invert reports a singular matrix instead of dividing by zero", "[linear_solve]") {
    Matrix<double, 3, 3> A;
    A(0,0)=1; A(0,1)=2; A(0,2)=3;
    A(1,0)=2; A(1,1)=4; A(1,2)=6;  // row 1 = 2 * row 0 -> singular
    A(2,0)=1; A(2,1)=0; A(2,2)=1;

    auto result = invert(A);
    REQUIRE_FALSE(result.has_value());
    CHECK(result.error().code == ErrorCode::SingularMatrix);
}

TEST_CASE("diagonal_dominance_ratio matches hand-computed values", "[linear_solve]") {
    Matrix<double, 3, 3> dominant;
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) dominant(i, j) = (i == j) ? 4.0 : 1.0;
    // Each row: diag=4, off-diag sum=2 -> ratio=2 on every row.
    CHECK_THAT(diagonal_dominance_ratio(dominant), WithinAbs(2.0, 1e-6));

    Matrix<double, 3, 3> weak;
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) weak(i, j) = (i == j) ? 1.0 : 1.0;
    // Each row: diag=1, off-diag sum=2 -> ratio=0.5 on every row.
    CHECK_THAT(diagonal_dominance_ratio(weak), WithinAbs(0.5, 1e-6));
}
