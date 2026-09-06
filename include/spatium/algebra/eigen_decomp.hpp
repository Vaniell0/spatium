#pragma once

// General-N eigendecomposition of a real symmetric matrix, via the
// classical cyclic Jacobi eigenvalue algorithm (Golub & Van Loan, "Matrix
// Computations", section 8.4).
//
// spaces/spd.hpp's detail::eigen_sym solves this exactly (root-finding the
// characteristic polynomial via solve_quadratic/solve_cubic) but only for
// N=2,3 -- a general-degree characteristic polynomial has no closed-form
// root formula for N>=5 (Abel-Ruffini), so any N-generic solver has to be
// iterative. Jacobi is the standard choice when N is small-to-moderate and
// simplicity/robustness matter more than the last constant factor: unlike
// shifted-QR (the asymptotically faster method production libraries like
// LAPACK use for large N), it needs no tridiagonalization pass, no shift
// strategy, and no special-casing for clustered eigenvalues -- it just
// repeatedly rotates away the largest sources of asymmetry-from-diagonal
// until none remain, and always converges for a real symmetric input
// (every rotation strictly decreases the sum of squared off-diagonal
// entries). That robustness is why it's still routinely taught and used
// for exactly this regime, not a toy simplification.
//
// Precondition: S is symmetric. Not checked at runtime -- same convention
// as detail::eigen_sym above (see its own file-level comment): this is an
// internal-invariant precondition the caller is responsible for, not a
// boundary function receiving arbitrary untrusted input (see
// docs/conventions.md's Result<T>-vs-assert rule). Passing a non-symmetric
// matrix silently computes garbage rather than diagnosing the misuse.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/matrix.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <algorithm>
#  include <array>
#  include <cmath>
#  include <cstddef>
#endif

SPATIUM_EXPORT namespace spatium {
inline namespace algebra {

// Eigendecomposition of a real symmetric N x N matrix S = V diag(values) V^T.
// values are sorted ASCENDING (matching Eigen::SelfAdjointEigenSolver's
// convention, so cross-validation against it needs no re-sorting); vectors'
// columns are orthonormal, vectors.col(i) is the eigenvector for values[i].
template<Scalar T, std::size_t N>
struct EigenDecomposition {
    Vec<T, N> values;
    Matrix<T, N, N> vectors;
};

// Reconstruct S = V diag(values) V^T -- the standard way to verify an
// eigendecomposition without re-deriving it independently.
template<Scalar T, std::size_t N>
Matrix<T, N, N> reconstruct(const EigenDecomposition<T, N>& eig) {
    Matrix<T, N, N> D;
    for (std::size_t i = 0; i < N; ++i) D(i, i) = eig.values[i];
    return eig.vectors * D * eig.vectors.transpose();
}

// Classical cyclic Jacobi: sweep every (p,q), p<q pair in fixed order
// (rather than always picking the single largest off-diagonal entry, which
// needs an O(N^2) search per rotation -- "classical" cyclic Jacobi visits
// pairs in a fixed pattern instead and still converges, just over a few
// more sweeps), rotating each pair to zero A(p,q) exactly, until the total
// off-diagonal mass is negligible relative to the matrix's own scale.
//
// max_sweeps is a safety cap, not the expected iteration count: cyclic
// Jacobi converges quadratically once the off-diagonal mass is small, and
// typically needs on the order of 5-10 sweeps regardless of N in practice.
// 30*N sweeps (a standard textbook bound -- see Golub & Van Loan) is
// generous headroom against a pathological input rather than a number this
// function expects to spend.
template<Scalar T, std::size_t N>
EigenDecomposition<T, N> eigen_decompose(const Matrix<T, N, N>& S,
                                          std::size_t max_sweeps = 30 * N) {
    using std::abs;
    using std::sqrt;

    Matrix<T, N, N> A = S;
    Matrix<T, N, N> V = Matrix<T, N, N>::identity();

    // Convergence/pivot-skip tolerance, scaled to the matrix's own
    // magnitude via epsilon.hpp's relative_epsilon() -- the same
    // scale-aware-tolerance idea spaces/spd.hpp's dup_tol uses, so a
    // matrix with large entries doesn't get held to an absolute tolerance
    // that's unreachable in floating point, and a matrix with tiny entries
    // doesn't get declared "converged" from the start.
    T frob_sq{};
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < N; ++j)
            frob_sq += A(i, j) * A(i, j);
    T tol = relative_epsilon<T>(sqrt(frob_sq));

    if constexpr (N > 1) {
        for (std::size_t sweep = 0; sweep < max_sweeps; ++sweep) {
            T off_sq{};
            for (std::size_t p = 0; p < N; ++p)
                for (std::size_t q = p + 1; q < N; ++q)
                    off_sq += A(p, q) * A(p, q);
            if (sqrt(off_sq) < tol) break;

            for (std::size_t p = 0; p < N - 1; ++p) {
                for (std::size_t q = p + 1; q < N; ++q) {
                    if (abs(A(p, q)) < tol) continue; // already ~0, rotating would only add roundoff

                    // Standard Jacobi rotation angle (Golub & Van Loan,
                    // algorithm 8.4.1 "sym_schur2"): t = tan(theta/2) of
                    // the rotation that zeroes A(p,q) exactly, taken as the
                    // smaller root of t^2 + 2*cot(2*theta)*t - 1 = 0 (the
                    // sign choice keeps the rotation angle small, which
                    // keeps convergence stable).
                    T theta = (A(q, q) - A(p, p)) / (T{2} * A(p, q));
                    T t = (theta >= T{0} ? T{1} : T{-1})
                        / (abs(theta) + sqrt(theta * theta + T{1}));
                    T c = T{1} / sqrt(t * t + T{1});
                    T s = t * c;

                    T apq = A(p, q);
                    A(p, p) = A(p, p) - t * apq;
                    A(q, q) = A(q, q) + t * apq;
                    A(p, q) = T{0};
                    A(q, p) = T{0};

                    for (std::size_t i = 0; i < N; ++i) {
                        if (i == p || i == q) continue;
                        T aip = A(i, p), aiq = A(i, q);
                        A(i, p) = c * aip - s * aiq;
                        A(p, i) = A(i, p);
                        A(i, q) = s * aip + c * aiq;
                        A(q, i) = A(i, q);
                    }
                    for (std::size_t i = 0; i < N; ++i) {
                        T vip = V(i, p), viq = V(i, q);
                        V(i, p) = c * vip - s * viq;
                        V(i, q) = s * vip + c * viq;
                    }
                }
            }
        }
    }

    // A's diagonal now holds the eigenvalues (up to the tolerance above);
    // sort ascending, permuting V's columns to stay matched to their value.
    std::array<std::size_t, N> order{};
    for (std::size_t i = 0; i < N; ++i) order[i] = i;
    std::sort(order.begin(), order.end(),
              [&A](std::size_t a, std::size_t b) { return A(a, a) < A(b, b); });

    EigenDecomposition<T, N> result;
    for (std::size_t slot = 0; slot < N; ++slot) {
        std::size_t src = order[slot];
        result.values[slot] = A(src, src);
        for (std::size_t r = 0; r < N; ++r) result.vectors(r, slot) = V(r, src);
    }
    return result;
}

} // namespace algebra
} // namespace spatium
