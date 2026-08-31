#pragma once

// General N×N linear system solving: Ax = b.
//
// Matrix<T,R,C> only ever had closed-form inverse() for 2x2/3x3 and
// determinant() for any N via LU -- no general solve existed for N>=4
// before this (a real, surveyed gap in the RSC domain pipeline, not a
// guess: "no general N×N linear algebra (only 2x2/3x3 closed-form
// inverse())"). Two methods, not one, since a genuine method-vs-method
// choice is what every RSC dispatch domain needs:
//   - solve_direct(): Gaussian elimination with partial pivoting. Always
//     gives the exact solution (to floating-point precision) in fixed
//     O(N^3) work, regardless of the matrix's conditioning.
//   - solve_jacobi(): classical Jacobi iteration, fixed iteration budget
//     -- no convergence check, deliberately, same "small enough that a
//     bad case measurably fails" discipline as rootfind_ops.hpp's Newton
//     budget. Cheap per iteration (O(N^2)), but only *guaranteed* to
//     converge when the matrix is strictly diagonally dominant --
//     diagonal_dominance_ratio() below is the real dispatch signal
//     rsc/include/linear_task.hpp builds on.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/matrix.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/core/error.hpp>
#  include <algorithm>
#  include <cmath>
#  include <cstddef>
#  include <limits>
#endif

SPATIUM_EXPORT namespace spatium {
inline namespace algebra {

// Gaussian elimination with partial pivoting, then back-substitution.
// Exact (to floating-point precision) regardless of conditioning; O(N^3).
template<Scalar T, std::size_t N>
Result<Vec<T, N>> solve_direct(Matrix<T, N, N> A, Vec<T, N> b) {
    using std::abs;
    for (std::size_t i = 0; i < N; ++i) {
        std::size_t piv = i;
        for (std::size_t k = i + 1; k < N; ++k)
            if (abs(A(k, i)) > abs(A(piv, i))) piv = k;
        if (abs(A(piv, i)) < epsilon<T>())
            return std::unexpected(Error{ErrorCode::SingularMatrix, "singular matrix in solve_direct"});
        if (piv != i) {
            for (std::size_t c = 0; c < N; ++c) std::swap(A(i, c), A(piv, c));
            std::swap(b[i], b[piv]);
        }
        for (std::size_t k = i + 1; k < N; ++k) {
            T factor = A(k, i) / A(i, i);
            for (std::size_t c = i; c < N; ++c) A(k, c) -= factor * A(i, c);
            b[k] -= factor * b[i];
        }
    }
    Vec<T, N> x{};
    for (std::size_t ii = 0; ii < N; ++ii) {
        std::size_t i = N - 1 - ii;
        T sum = b[i];
        for (std::size_t c = i + 1; c < N; ++c) sum -= A(i, c) * x[c];
        x[i] = sum / A(i, i);
    }
    return x;
}

// General N×N matrix inverse via Gauss-Jordan elimination: one pivoted
// elimination pass carried simultaneously through N right-hand sides
// (the identity's columns), not N independent solve_direct() calls --
// those would redo the same O(N^3) elimination N times over. Added
// alongside solve_direct() specifically for physics/relativity/
// geodesic.hpp's per-step metric inverse, where profiling showed
// 4 solve_direct() calls per Christoffel-symbol evaluation were the
// dominant cost (measured ~3.8x speedup switching to this).
template<Scalar T, std::size_t N>
Result<Matrix<T, N, N>> invert(Matrix<T, N, N> A) {
    using std::abs;
    Matrix<T, N, N> inv = Matrix<T, N, N>::identity();
    for (std::size_t i = 0; i < N; ++i) {
        std::size_t piv = i;
        for (std::size_t k = i + 1; k < N; ++k)
            if (abs(A(k, i)) > abs(A(piv, i))) piv = k;
        if (abs(A(piv, i)) < epsilon<T>())
            return std::unexpected(Error{ErrorCode::SingularMatrix, "singular matrix in invert"});
        if (piv != i) {
            for (std::size_t c = 0; c < N; ++c) {
                std::swap(A(i, c), A(piv, c));
                std::swap(inv(i, c), inv(piv, c));
            }
        }
        T diag = A(i, i);
        for (std::size_t c = 0; c < N; ++c) { A(i, c) /= diag; inv(i, c) /= diag; }
        for (std::size_t k = 0; k < N; ++k) {
            if (k == i) continue;
            T factor = A(k, i);
            if (factor == T{0}) continue;
            for (std::size_t c = 0; c < N; ++c) {
                A(k, c) -= factor * A(i, c);
                inv(k, c) -= factor * inv(i, c);
            }
        }
    }
    return inv;
}

// Classical Jacobi iteration, fixed budget, no convergence check. Only
// *guaranteed* to converge when A is strictly diagonally dominant; off
// that regime this may return an arbitrarily bad (even diverged) answer
// after max_iter steps -- which is the entire point. The dispatcher has
// to recognize which regime it's in from the matrix alone, not from
// watching this function fail after the fact.
template<Scalar T, std::size_t N>
Vec<T, N> solve_jacobi(const Matrix<T, N, N>& A, const Vec<T, N>& b, int max_iter,
                        Vec<T, N> x0 = Vec<T, N>{}) {
    Vec<T, N> x = x0;
    for (int it = 0; it < max_iter; ++it) {
        Vec<T, N> x_new{};
        for (std::size_t i = 0; i < N; ++i) {
            T sum = b[i];
            for (std::size_t j = 0; j < N; ++j)
                if (j != i) sum -= A(i, j) * x[j];
            x_new[i] = sum / A(i, i);
        }
        x = x_new;
    }
    return x;
}

// min over rows of |A_ii| / (sum_{j!=i} |A_ij| + eps) -- the standard
// diagonal-dominance ratio: >1 on every row means A is strictly
// diagonally dominant, the textbook sufficient condition for Jacobi's
// guaranteed convergence. A single scalar summarizing "how well-
// conditioned is this matrix's weakest row" -- the real dispatch signal.
template<Scalar T, std::size_t N>
T diagonal_dominance_ratio(const Matrix<T, N, N>& A) {
    using std::abs;
    T worst = std::numeric_limits<T>::max();
    for (std::size_t i = 0; i < N; ++i) {
        T off_sum{0};
        for (std::size_t j = 0; j < N; ++j)
            if (j != i) off_sum += abs(A(i, j));
        T ratio = abs(A(i, i)) / (off_sum + epsilon<T>());
        if (ratio < worst) worst = ratio;
    }
    return worst;
}

} // namespace algebra
} // namespace spatium
