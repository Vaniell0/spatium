#pragma once

// General M x N Singular Value Decomposition: A = U diag(sigma) V^T.
//
// Approach: eigendecompose A^T A (N x N, symmetric positive-semidefinite)
// via eigen_decomp.hpp -- its eigenvectors are V and its eigenvalues are
// the squared singular values (A^T A v = sigma^2 v follows directly from
// the SVD definition). U's columns then come from A v_i / sigma_i for each
// nonzero sigma_i: (A v_i) is automatically orthogonal to (A v_j) for
// i != j because v_i^T A^T A v_j = sigma_j^2 (v_i . v_j) = 0, so no
// separate orthogonalization step is needed there -- only the near-zero/
// rank-deficient singular values (where dividing by sigma_i would be
// dividing by ~0) need a completion step, via Gram-Schmidt against the
// columns already placed.
//
// This is the standard "normal equations" SVD construction. It is not
// what production numerical libraries use internally (squaring A doubles
// the condition number -- A^T A's condition number is that of A squared,
// so it loses roughly half of double precision's usable digits versus a
// direct bidiagonalization-based SVD for ill-conditioned inputs) but it is
// simple, reuses eigen_decomp.hpp directly, and is entirely adequate for
// well-conditioned inputs and small-to-moderate N -- the same tradeoff
// spaces/spd.hpp's closed-form eigen_sym makes deliberately. A
// condition-number-robust SVD is a real, separate improvement left for if
// it's ever needed (see this file's tests for the accuracy this
// construction actually achieves).
//
// "Thin" SVD, not "full": U is M x K and V is N x K, K = min(M,N), rather
// than the alternative square-orthogonal U (M x M) / V (N x N) padded with
// an arbitrary orthonormal complement of the columns this construction can
// reach. Thin is what every practical consumer (least-squares, low-rank
// approximation, pseudo-inverse) actually wants, and is the form
// Eigen::JacobiSVD::compute() with ComputeThinU/ComputeThinV also produces
// (used directly for this file's cross-validation tests) -- full SVD is a
// straightforward but unneeded extension (pad with any orthonormal
// complement of the K columns already placed) left for a future caller who
// specifically needs it.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/eigen_decomp.hpp>
#  include <spatium/algebra/matrix.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <cmath>
#  include <cstddef>
#endif

SPATIUM_EXPORT namespace spatium {
inline namespace algebra {

// Thin SVD of an M x N matrix: A ~= U * diag(singular_values) * V^T.
// K = min(M,N); singular_values are sorted DESCENDING (the universal SVD
// convention, matching Eigen::JacobiSVD) and are all >= 0. U's and V's
// columns are each orthonormal.
template<Scalar T, std::size_t M, std::size_t N>
struct SVDResult {
    static constexpr std::size_t K = (M < N) ? M : N;
    Vec<T, K> singular_values;
    Matrix<T, M, K> U;
    Matrix<T, N, K> V;
};

// Reconstruct A ~= U diag(singular_values) V^T -- sum of rank-1 terms
// sigma_k * u_k v_k^T, the standard way to verify an SVD without
// re-deriving it independently.
template<Scalar T, std::size_t M, std::size_t N>
Matrix<T, M, N> reconstruct(const SVDResult<T, M, N>& r) {
    Matrix<T, M, N> A;
    for (std::size_t k = 0; k < SVDResult<T, M, N>::K; ++k) {
        T sigma = r.singular_values[k];
        auto u = r.U.col(k);
        auto v = r.V.col(k);
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t j = 0; j < N; ++j)
                A(i, j) += sigma * u[i] * v[j];
    }
    return A;
}

template<Scalar T, std::size_t M, std::size_t N>
SVDResult<T, M, N> svd(const Matrix<T, M, N>& A, std::size_t max_sweeps = 30 * N) {
    using std::sqrt;
    constexpr std::size_t K = SVDResult<T, M, N>::K;

    // A^T A is N x N, symmetric PSD: its eigenpairs give V and sigma^2
    // directly (see file-level comment). eigen_decompose() sorts
    // ascending; SVD wants descending, so column src=N-1-k of eig is
    // placed into slot k below.
    Matrix<T, N, N> AtA = A.transpose() * A;
    auto eig = eigen_decompose(AtA, max_sweeps);

    SVDResult<T, M, N> result;
    T sigma_max{};
    for (std::size_t k = 0; k < K; ++k) {
        std::size_t src = N - 1 - k;
        T lambda = eig.values[src];
        // A^T A is PSD only in exact arithmetic; Jacobi's floating-point
        // result can land a hair below zero for a genuinely-zero
        // eigenvalue (a rank-deficient A) -- clamp rather than feed sqrt()
        // a negative argument.
        T sigma = lambda > T{0} ? sqrt(lambda) : T{0};
        result.singular_values[k] = sigma;
        if (k == 0) sigma_max = sigma;
        for (std::size_t row = 0; row < N; ++row) result.V(row, k) = eig.vectors(row, src);
    }

    // Rank tolerance scaled to the largest singular value, same
    // magnitude-aware idea as eigen_decompose()'s own convergence
    // tolerance -- a singular value can only be trusted relative to the
    // matrix's own scale, not against an absolute constant.
    T tol = relative_epsilon<T>(sigma_max);

    for (std::size_t k = 0; k < K; ++k) {
        T sigma = result.singular_values[k];
        if (sigma > tol) {
            Vec<T, M> u = A * result.V.col(k) / sigma;
            for (std::size_t row = 0; row < M; ++row) result.U(row, k) = u[row];
            continue;
        }

        // Near-zero/rank-deficient singular value: A*v_k/sigma_k would be
        // dividing by ~0 noise, not a usable direction. Complete U's basis
        // instead with a standard basis vector, Gram-Schmidt'd against the
        // U columns already placed -- mirrors spaces/spd.hpp's
        // detail::eigen_sym degenerate-eigenvalue fix: vary the seed axis
        // by column index k so that several near-zero singular values in
        // the same call don't all start Gram-Schmidt from the same
        // candidate (which would leave them numerically parallel until
        // roundoff happens to separate them, rather than genuinely
        // orthogonal).
        Vec<T, M> filled{};
        bool found = false;
        for (std::size_t attempt = 0; attempt < M && !found; ++attempt) {
            std::size_t axis = (k + attempt) % M;
            Vec<T, M> seed = Vec<T, M>::unit(axis);
            for (std::size_t j = 0; j < k; ++j) {
                Vec<T, M> uj = result.U.col(j);
                seed = seed - uj * uj.dot(seed);
            }
            T n = seed.norm();
            if (n > tol) {
                filled = seed / n;
                found = true;
            }
        }
        for (std::size_t row = 0; row < M; ++row) result.U(row, k) = filled[row];
    }

    return result;
}

} // namespace algebra
} // namespace spatium
