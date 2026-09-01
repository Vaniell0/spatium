#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <spatium/algebra/matrix.hpp>
#  include <spatium/algebra/complex.hpp>
#  include <spatium/algebra/polynomial.hpp>
#  include <cmath>
#  include <cstddef>
#  include <span>
#endif

SPATIUM_EXPORT namespace spatium {

// internal — do not use, no API stability. Closed-form eigendecomposition
// of a real symmetric N×N matrix, for N=2,3 only: the eigenvalues are the
// roots of the characteristic polynomial (degree N), always real for a
// symmetric matrix, so solve_quadratic()/solve_cubic() give them directly
// with no iterative eigensolver. General N needs the not-yet-built native
// SVD/eigendecomposition (ROADMAP → Backlog → Native math); this is
// deliberately scoped to the sizes SPD<N> below actually supports.
namespace detail {

template<Scalar T>
struct EigenSym2 {
    Vec<T, 2> values;
    Matrix<T, 2, 2> vectors; // columns are the (unit) eigenvectors
};

template<Scalar T>
EigenSym2<T> eigen_sym(const Matrix<T, 2, 2>& S) {
    using std::abs;
    T a = S(0, 0), b = S(0, 1), d = S(1, 1);
    auto roots = solve_quadratic(T{1}, -(a + d), a * d - b * b);
    Vec<T, 2> values{roots[0].re, roots[1].re};

    Matrix<T, 2, 2> vecs;
    for (std::size_t i = 0; i < 2; ++i) {
        T lambda = values[i];
        Vec<T, 2> v;
        if (abs(b) > epsilon<T>()) {
            v = Vec<T, 2>{b, lambda - a}.normalized();
        } else {
            // Already diagonal: no rotation to solve for, just match each
            // root back to the diagonal entry it came from.
            v = (abs(lambda - a) <= abs(lambda - d)) ? Vec<T, 2>{T{1}, T{0}}
                                                       : Vec<T, 2>{T{0}, T{1}};
        }
        vecs(0, i) = v[0];
        vecs(1, i) = v[1];
    }
    return {values, vecs};
}

template<Scalar T>
struct EigenSym3 {
    Vec<T, 3> values;
    Matrix<T, 3, 3> vectors;
};

template<Scalar T>
EigenSym3<T> eigen_sym(const Matrix<T, 3, 3>& S) {
    T a = S(0, 0), b = S(0, 1), c = S(0, 2),
                    d = S(1, 1), e = S(1, 2),
                                  f = S(2, 2);
    T trace = a + d + f;
    T q = (d * f - e * e) + (a * f - c * c) + (a * d - b * b); // sum of principal 2x2 minors
    T det = a * (d * f - e * e) - b * (b * f - e * c) + c * (b * e - d * c);

    // lambda^3 - trace*lambda^2 + q*lambda - det = 0
    auto roots = solve_cubic(T{1}, -trace, q, -det);
    Vec<T, 3> values{roots[0].re, roots[1].re, roots[2].re};

    Matrix<T, 3, 3> vecs;
    for (std::size_t i = 0; i < 3; ++i) {
        // Null space of (S - lambda*I): for a generic (rank-2) 3x3 matrix,
        // the cross product of any two independent rows spans it. Picking
        // the pair with the largest cross-product norm keeps this stable
        // near-degenerate cases would otherwise amplify.
        Matrix<T, 3, 3> A = S - Matrix<T, 3, 3>::identity() * values[i];
        Vec<T, 3> r0 = A.row(0), r1 = A.row(1), r2 = A.row(2);
        Vec<T, 3> c01 = r0.cross(r1), c02 = r0.cross(r2), c12 = r1.cross(r2);
        T n01 = c01.norm(), n02 = c02.norm(), n12 = c12.norm();
        Vec<T, 3> v = (n01 >= n02 && n01 >= n12) ? c01 : (n02 >= n12 ? c02 : c12);
        auto vn = v.norm();
        if (vn > epsilon<T>()) v = v / vn;

        // Gram-Schmidt against columns already placed: eigenvectors of a
        // symmetric matrix for distinct eigenvalues are exactly orthogonal
        // already, so this is close to a no-op there -- it earns its keep
        // only when two roots are nearly equal, where the raw cross-product
        // above can land anywhere inside the (near-)degenerate eigenspace
        // instead of orthogonal to the vector already chosen for it.
        for (std::size_t j = 0; j < i; ++j) {
            Vec<T, 3> prev{vecs(0, j), vecs(1, j), vecs(2, j)};
            v = v - prev * prev.dot(v);
        }
        auto gn = v.norm();
        if (gn > epsilon<T>()) v = v / gn;

        vecs(0, i) = v[0]; vecs(1, i) = v[1]; vecs(2, i) = v[2];
    }
    return {values, vecs};
}

} // namespace detail

// SPD(n) — symmetric positive-definite n×n matrices, under the
// log-Euclidean metric (Arsigny, Fillard, Pennec, Ayache, 2006). Real
// applications: SO(3)/SE(3) rotation averaging's natural sibling —
// covariance-matrix descriptors in computer vision, diffusion tensors in
// DTI/medical imaging, spatial-covariance features in BCI/EEG.
//
// The defining move of log-Euclidean: since matrix log is a global
// diffeomorphism SPD(n) <-> Sym(n) (ordinary symmetric matrices, a flat
// vector space), a "point" of this space is represented directly by
// vech(log(S)) -- the independent entries of the LOGARITHM, not of S
// itself. Under that parametrization every Manifold/RiemannianManifold
// operation (distance, exp_map, log_map, metric_at) is the ordinary flat
// Euclidean formula: the curvature of SPD(n) is entirely absorbed into the
// change of coordinates (matrix log/exp), not left for exp_map/log_map to
// deal with the way Sphere/Hyperbolic must. That's the whole reason
// log-Euclidean is called "cheap" in the literature, and it has a sharp,
// checkable consequence used by frechet_mean() below: the Fréchet mean of
// SPD matrices under this metric is a CLOSED FORM (average the logs,
// exponentiate back), not an iterative optimization — riemannian_minimize()
// is genuinely unneeded here, unlike Sphere/Hyperbolic where retraction
// actually has curvature to contend with.
//
// vech() packs a symmetric matrix's independent entries into a
// Vec<T, N*(N+1)/2>, diagonal first then off-diagonals scaled by sqrt(2),
// so the plain dot product on PointType equals the Frobenius inner product
// <A,B>_F = trace(A B) on the matrix itself -- the metric_at() below is the
// ordinary dot product precisely because of this weighting, not despite it.
//
// Deliberately NOT a Surface (no project()/normal()): under log-Euclidean,
// SPD(n) has no meaningful embedding as a hypersurface of some larger
// ambient space the way Sphere/Hyperbolic do -- it's flat and full-
// dimensional in its own log-space coordinates, so there is no normal
// direction to define. Manifold + RiemannianManifold is the honest
// concept-hierarchy fit.
//
// N is restricted to 2 and 3 -- see detail::eigen_sym() above for why.
template<std::size_t N, Scalar T = double>
struct SPD {
    static_assert(N == 2 || N == 3,
        "SPD<N> currently supports N=2,3 only (closed-form eigendecomposition "
        "via solve_quadratic/solve_cubic); general N needs the not-yet-built "
        "native SVD/eigendecomposition, see ROADMAP -> Backlog -> Native math.");

    static constexpr std::size_t M = N * (N + 1) / 2;

    using ScalarType    = T;
    using PointType     = Vec<T, M>; // vech(log(S)) -- log-Euclidean coordinates
    using TangentVector = Vec<T, M>;
    using MatrixType    = Matrix<T, N, N>;

    static constexpr std::size_t dimension = M;
    static constexpr bool is_complete = true;

    // TopologicalSpace: every vech(log(S)) is some SPD matrix's log; the
    // whole M-dimensional log-space is valid.
    constexpr bool contains(const PointType&) const { return true; }

    // MetricSpace / Manifold / RiemannianManifold: flat, by construction --
    // see the file-level comment above for why this isn't approximate.
    ScalarType distance(const PointType& p, const PointType& q) const {
        return (q - p).norm();
    }

    constexpr PointType exp_map(const PointType& p, const TangentVector& v, ScalarType t) const {
        return p + v * t;
    }

    constexpr TangentVector log_map(const PointType& p, const PointType& q) const {
        return q - p;
    }

    constexpr ScalarType metric_at(const PointType&,
                                    const TangentVector& u,
                                    const TangentVector& v) const {
        return u.dot(v);
    }

    // ── Data in/out: actual SPD matrices at the boundary ──────────

    static PointType from_spd(const MatrixType& S) {
        return vech(matrix_log_sym(S));
    }

    static MatrixType to_spd(const PointType& p) {
        return matrix_exp_sym(unvech(p));
    }

    // ── vech <-> symmetric matrix ──────────────────────────────────

    static PointType vech(const MatrixType& L) {
        using std::sqrt; // ADL: lets non-std Scalar T provide its own sqrt
        auto s2 = sqrt(T{2});
        PointType v{};
        for (std::size_t i = 0; i < N; ++i) v[i] = L(i, i);
        std::size_t k = N;
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = i + 1; j < N; ++j)
                v[k++] = L(i, j) * s2;
        return v;
    }

    static MatrixType unvech(const PointType& v) {
        using std::sqrt;
        auto s2 = sqrt(T{2});
        MatrixType L;
        for (std::size_t i = 0; i < N; ++i) L(i, i) = v[i];
        std::size_t k = N;
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = i + 1; j < N; ++j) {
                T off = v[k++] / s2;
                L(i, j) = off;
                L(j, i) = off;
            }
        return L;
    }

    // ── Symmetric matrix log/exp via closed-form eigendecomposition ──
    // Both share detail::eigen_sym(): log(S) needs S's eigenvalues to be
    // positive (SPD's defining property); exp(L) places no such
    // requirement on L (any real symmetric matrix), which is exactly what
    // lets a TangentVector -- an unconstrained point of the flat log-space
    // -- always map back to a genuine SPD matrix via to_spd().

    static MatrixType matrix_log_sym(const MatrixType& S) {
        using std::log;
        auto eig = detail::eigen_sym(S);
        MatrixType D;
        for (std::size_t i = 0; i < N; ++i) D(i, i) = log(eig.values[i]);
        return eig.vectors * D * eig.vectors.transpose();
    }

    static MatrixType matrix_exp_sym(const MatrixType& L) {
        using std::exp;
        auto eig = detail::eigen_sym(L);
        MatrixType D;
        for (std::size_t i = 0; i < N; ++i) D(i, i) = exp(eig.values[i]);
        return eig.vectors * D * eig.vectors.transpose();
    }
};

// Log-Euclidean Fréchet mean: the point minimizing the sum of squared
// log-Euclidean distances to the samples. Closed form (see SPD<N,T>'s
// file-level comment for why) -- average in log-space, exponentiate back.
// No riemannian_minimize() call needed, unlike Sphere/Hyperbolic.
template<std::size_t N, Scalar T = double>
Matrix<T, N, N> frechet_mean(std::span<const Matrix<T, N, N>> matrices) {
    using S = SPD<N, T>;
    typename S::PointType sum{};
    for (const auto& m : matrices) sum += S::from_spd(m);
    typename S::PointType mean_point = sum / T(matrices.size());
    return S::to_spd(mean_point);
}

static_assert(Manifold<SPD<2>>);
static_assert(RiemannianManifold<SPD<2>>);
static_assert(Manifold<SPD<3>>);
static_assert(RiemannianManifold<SPD<3>>);

} // namespace spatium
