#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <spatium/algebra/matrix.hpp>
#  include <spatium/algebra/complex.hpp>
#  include <spatium/algebra/dual.hpp>
#  include <spatium/algebra/polynomial.hpp>
#  include <cmath>
#  include <cstddef>
#  include <limits>
#  include <span>
#endif

SPATIUM_EXPORT namespace spatium {

// internal — do not use, no API stability. Eigendecomposition of a real
// symmetric N×N matrix, for N=2,3 only: the half-angle formulas for 2 by 2
// and cyclic Jacobi rotations for 3 by 3, both orthonormal to rounding
// whatever the gap between the eigenvalues. General N needs the not-yet-built native
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
    // The half-angle formulas, not the characteristic quadratic. Taking the eigenvalues from the
    // quadratic and the vectors as (b, lambda - a) is exact on paper and loses the vectors where
    // the eigenvalues are close: lambda - a is a difference of nearly equal numbers with an error
    // of order sqrt(eps) in the root, so the direction is wrong by that over the gap, and
    // V diag(f(lambda)) V^T missed S by 1e-3 at a distance of 1e-6 from the identity and by 0.15 at
    // 1e-7 (error ~ eps / gap^2; found by a finite difference of sqrt_sym that did not match its
    // derivative). Here c = (a - d)/2 and r = hypot(c, b) are formed without cancellation, the
    // larger eigenvalue is m + r, and its vector is (c + r, b) when c >= 0 and (b, r - c) otherwise,
    // whichever has no subtraction of like signs; the other vector is perpendicular by construction,
    // so the pair is orthonormal to rounding whatever the gap.
    using std::sqrt; using std::abs;
    const T a = S(0, 0), b = S(0, 1), d = S(1, 1);
    const T mean = T((a + d) / T{2});
    const T c = T((a - d) / T{2});
    // Scaled by the larger of |c| and |b|, so that squaring cannot underflow: a float matrix within
    // 1e-23 of a scalar one has c*c + b*b = 0 or a denormal and the normalisation below was a 0/0
    // (the Frechet mean on SPD x S2 in float came back NaN when its iteration converged).
    const T s = (abs(c) >= abs(b)) ? T(abs(c)) : T(abs(b));
    Matrix<T, 2, 2> vecs;
    if (s == T{0}) {                               // a scalar matrix: any orthonormal pair, by slot
        vecs(0, 0) = T{1}; vecs(1, 0) = T{0}; vecs(0, 1) = T{0}; vecs(1, 1) = T{1};
        return {Vec<T, 2>{mean, mean}, vecs};
    }
    const T c1 = T(c / s), b1 = T(b / s);
    const T r1 = T(sqrt(T(c1 * c1 + b1 * b1)));    // in [1, sqrt 2]
    const T r = T(s * r1);
    Vec<T, 2> values{T(mean + r), T(mean - r)};

    Vec<T, 2> up;                                  // eigenvector of the larger eigenvalue
    if (c1 >= T{0}) up = Vec<T, 2>{T(c1 + r1), b1}.normalized();
    else            up = Vec<T, 2>{b1, T(r1 - c1)}.normalized();
    vecs(0, 0) = up[0];
    vecs(1, 0) = up[1];
    vecs(0, 1) = T(-up[1]);
    vecs(1, 1) = up[0];
    return {values, vecs};
}

template<Scalar T>
struct EigenSym3 {
    Vec<T, 3> values;
    Matrix<T, 3, 3> vectors;
};

template<Scalar T>
EigenSym3<T> eigen_sym(const Matrix<T, 3, 3>& S) {
    // Cyclic Jacobi rotations, not the characteristic cubic. The cubic gives the eigenvalues to
    // sqrt(epsilon) near a repeated root, and the vectors taken from them (a cross product of rows of
    // S - lambda I) are wrong where the eigenvalues are close: on I + hM the reconstruction
    // V diag(l) V^T missed by 1.7e-3 at h = 1e-2 and by 0.95, with the vectors no longer orthonormal,
    // at h = 1e-5. A rotation zeroes one off-diagonal entry of a symmetric matrix exactly and keeps
    // V orthogonal to rounding whatever the eigenvalues are, so a repeated or a nearly repeated
    // eigenvalue costs nothing; three entries, a few sweeps (quadratic convergence), 3 by 3 only.
    // Sorted descending, as the 2 by 2 is.
    using std::abs; using std::sqrt;
    Matrix<T, 3, 3> A = S;
    Matrix<T, 3, 3> V = Matrix<T, 3, 3>::identity();
    const T eps = epsilon<T>();
    constexpr std::size_t pairs[3][2] = {{0, 1}, {0, 2}, {1, 2}};
    for (int sweep = 0; sweep < 32; ++sweep) {
        bool rotated = false;
        for (const auto& pq : pairs) {
            const std::size_t p = pq[0], q = pq[1];
            const T apq = A(p, q);
            const T scale = T(abs(A(p, p)) + abs(A(q, q)));
            if (abs(apq) <= eps * scale * T{0.01} || apq == T{0}) { A(p, q) = T{0}; A(q, p) = T{0}; continue; }
            const T tau = T((A(q, q) - A(p, p)) / (T{2} * apq));
            const T t = T((tau >= T{0} ? T{1} : T{-1}) / (abs(tau) + sqrt(T(T{1} + tau * tau))));
            const T c = T(T{1} / sqrt(T(T{1} + t * t)));
            const T s = T(t * c);
            Matrix<T, 3, 3> J = Matrix<T, 3, 3>::identity();
            J(p, p) = c; J(q, q) = c; J(p, q) = s; J(q, p) = T(-s);
            A = J.transpose() * A * J;
            A(p, q) = T{0}; A(q, p) = T{0};            // annihilated by construction; rounding left a speck
            V = V * J;
            rotated = true;
        }
        if (!rotated) break;
    }
    Vec<T, 3> values{A(0, 0), A(1, 1), A(2, 2)};
    Matrix<T, 3, 3> vecs = V;
    for (std::size_t i = 0; i < 2; ++i)                // sort descending, columns with their values
        for (std::size_t j = i + 1; j < 3; ++j)
            if (values[j] > values[i]) {
                const T tmp = values[i]; values[i] = values[j]; values[j] = tmp;
                for (std::size_t k = 0; k < 3; ++k) { const T col = vecs(k, i); vecs(k, i) = vecs(k, j); vecs(k, j) = col; }
            }
    return {values, vecs};
}

// internal — do not use, no API stability. U diag(f(lambda_i)) U^T for a
// real symmetric S = U diag(lambda_i) U^T: the shared reconstruction step
// behind every symmetric matrix function used in this file (log, exp, sqrt,
// inverse-sqrt, inverse) — same eigendecomposition, different scalar map
// over the eigenvalues.
template<std::size_t N, Scalar T, typename F>
Matrix<T, N, N> apply_eigen_sym(const Matrix<T, N, N>& S, F&& f) {
    auto eig = eigen_sym(S);
    Matrix<T, N, N> D;
    for (std::size_t i = 0; i < N; ++i) D(i, i) = f(eig.values[i]);
    return eig.vectors * D * eig.vectors.transpose();
}

// The same spectral function of a matrix of Duals, differentiated by the Daleckii-Krein formula
//     d f(S) = V ( Gamma o (V^T dS V) ) V^T,    Gamma_ij = (f(l_i) - f(l_j)) / (l_i - l_j),  f'(l) when l_i = l_j,
// not by carrying the Dual through the eigendecomposition. Through it, a repeated eigenvalue is a
// point where the derivative of a closed-form root is 0/0 or infinite (sqrt of a discriminant that is
// zero), and the matrix 0 is exactly such a point: exp_map(p, v, t) at t = 0 hands exp the matrix
// t * P^-1/2 V P^-1/2, whose eigenvalues coincide, and its derivative came back NaN, which is why
// the affine-invariant SPD distance could not be differentiated along a geodesic. A function of a
// matrix is analytic where its eigenvalues cross; only the decomposition is not, and this does not
// differentiate the decomposition. f is called once per eigenvalue on a Dual with slope 1, which
// gives f and f' together. The eigenvalues are told apart by sqrt(epsilon) of their scale.
template<std::size_t N, Scalar U, typename F>
Matrix<Dual<U>, N, N> apply_eigen_sym(const Matrix<Dual<U>, N, N>& S, F&& f) {
    using D = Dual<U>;
    Matrix<U, N, N> value, rate;
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < N; ++j) {
            value(i, j) = S(i, j).value;
            rate(i, j) = S(i, j).deriv;
        }
    const auto eig = eigen_sym(value);
    U fv[N], fd[N];
    for (std::size_t i = 0; i < N; ++i) {
        const D r = f(D{eig.values[i], U{1}});
        fv[i] = r.value;
        fd[i] = r.deriv;
    }
    const Matrix<U, N, N> vt = eig.vectors.transpose();
    const Matrix<U, N, N> q = vt * rate * eig.vectors;      // the perturbation in the eigenbasis
    Matrix<U, N, N> lambda, g;
    for (std::size_t i = 0; i < N; ++i) {
        lambda(i, i) = fv[i];
        for (std::size_t j = 0; j < N; ++j) {
            const U gap = eig.values[i] - eig.values[j];
            const double scale = 1.0 + std::abs(primal_double(eig.values[i])) + std::abs(primal_double(eig.values[j]));
            const bool apart = std::abs(primal_double(gap)) > std::sqrt(std::numeric_limits<double>::epsilon()) * scale;
            const U gamma = apart ? U((fv[i] - fv[j]) / gap) : U((fd[i] + fd[j]) / U{2});
            g(i, j) = gamma * q(i, j);
        }
    }
    const Matrix<U, N, N> r = eig.vectors * lambda * vt;
    const Matrix<U, N, N> dr = eig.vectors * g * vt;
    Matrix<D, N, N> out;
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < N; ++j) out(i, j) = D{r(i, j), dr(i, j)};
    return out;
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
        return detail::apply_eigen_sym(S, [](T x) { using std::log; return log(x); });
    }

    static MatrixType matrix_exp_sym(const MatrixType& L) {
        return detail::apply_eigen_sym(L, [](T x) { using std::exp; return exp(x); });
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

// SPD(n) — same matrices as SPD<N,T> above, under the affine-invariant
// metric (Pennec, Fillard, Ayache 2006; Moakher 2005) instead of
// log-Euclidean's flat approximation. Where SPD<N,T> encodes a point as
// vech(log(S)) precisely so every Manifold operation collapses to ordinary
// flat arithmetic, this class keeps a point as the literal SPD matrix and
// lets the real curvature show: the inner product genuinely depends on the
// base point,
//     <U, V>_S = trace(S^-1 U S^-1 V),
// invariant under every congruence S -> A S A^T for invertible A (hence the
// name) -- log-Euclidean's <U,V> = trace(UV) via vech is NOT invariant
// under that action, which is the real mathematical gap this class closes,
// not just an alternate style.
//
// exp_map/log_map are the standard closed forms
//     Exp_S(V) = S^{1/2} exp(S^{-1/2} V S^{-1/2}) S^{1/2}
//     Log_S(Q) = S^{1/2} log(S^{-1/2} Q S^{-1/2}) S^{1/2}
// Both reduce to the same eigendecomposition-based symmetric matrix log/exp
// as SPD<N,T> above (S^{-1/2} V S^{-1/2} is symmetric whenever V is, since
// S^{-1/2} is symmetric) -- reused directly via detail::apply_eigen_sym, no
// new numerical machinery, same N=2,3 restriction for the same reason.
//
// distance() is defined as sqrt(metric_at(p, log_map(p,q), log_map(p,q)))
// rather than the textbook direct formula (generalized eigenvalues of
// S1^-1 S2) -- deliberately: this is the general Riemannian identity
// (distance = norm of the initial geodesic velocity, in the metric at the
// start point), and it keeps distance/log_map/metric_at consistent by
// construction instead of risking the three drifting out of sync under
// independent hand derivations (see the SE3::exp() translation-Jacobian bug
// in ROADMAP.md's SO(3)/SE(3) entry for what that risk actually costs in
// this codebase).
//
// Unlike SPD<N,T>, no from_spd()/to_spd() are needed: PointType already IS
// the SPD matrix, nothing to encode/decode at the boundary.
//
// SPD(n) under this metric is a Hadamard manifold (complete, non-positively
// curved, uniquely geodesic): geodesics genuinely bend away from the
// boundary of the SPD cone (Log_S(Q) diverges as Q approaches a singular
// matrix) instead of being straight lines in disguise. That extra fidelity
// over log-Euclidean's flat approximation is exactly what costs the closed
// form: see frechet_mean_affine_invariant() below.
template<std::size_t N, Scalar T = double>
struct SPDAffineInvariant {
    static_assert(N == 2 || N == 3,
        "SPDAffineInvariant<N> currently supports N=2,3 only, for the same "
        "reason as SPD<N,T> above -- see its static_assert.");

    using ScalarType    = T;
    using MatrixType    = Matrix<T, N, N>;
    using PointType     = MatrixType; // the literal SPD matrix
    using TangentVector = MatrixType; // a symmetric matrix at any base point

    static constexpr std::size_t dimension = N * (N + 1) / 2;
    static constexpr bool is_complete = true; // Hadamard manifold

    // TopologicalSpace: the open SPD cone, not all of Matrix<T,N,N> -- unlike
    // SPD<N,T>'s flat log-space (every vech(log(S)) is valid), a literal
    // matrix point here must actually be symmetric with positive eigenvalues.
    bool contains(const PointType& p) const {
        auto asym = p - p.transpose();
        T off{};
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < N; ++j)
                off += asym(i, j) * asym(i, j);
        if (off > epsilon<T>() * epsilon<T>()) return false;
        auto eig = detail::eigen_sym(p);
        for (std::size_t i = 0; i < N; ++i)
            if (eig.values[i] <= epsilon<T>()) return false;
        return true;
    }

    ScalarType distance(const PointType& p, const PointType& q) const {
        using std::sqrt;
        auto v = log_map(p, q);
        return sqrt(metric_at(p, v, v));
    }

    PointType exp_map(const PointType& p, const TangentVector& v, ScalarType t) const {
        auto p_half = sqrt_sym(p);
        auto p_ih = inv_sqrt_sym(p);
        MatrixType mid = p_ih * (v * t) * p_ih; // symmetric: v and p_ih both are
        MatrixType e = detail::apply_eigen_sym(mid, [](T x) { using std::exp; return exp(x); });
        return p_half * e * p_half;
    }

    TangentVector log_map(const PointType& p, const PointType& q) const {
        auto p_half = sqrt_sym(p);
        auto p_ih = inv_sqrt_sym(p);
        MatrixType mid = p_ih * q * p_ih; // SPD: congruence of an SPD q stays SPD
        MatrixType l = detail::apply_eigen_sym(mid, [](T x) { using std::log; return log(x); });
        return p_half * l * p_half;
    }

    ScalarType metric_at(const PointType& p, const TangentVector& u, const TangentVector& v) const {
        auto p_inv = inv_sym(p);
        return (p_inv * u * p_inv * v).trace();
    }

    // ── Symmetric matrix sqrt/inverse-sqrt/inverse, via the same closed-form
    // eigendecomposition as SPD<N,T>::matrix_log_sym/matrix_exp_sym ──

    static MatrixType sqrt_sym(const MatrixType& S) {
        return detail::apply_eigen_sym(S, [](T x) { using std::sqrt; return sqrt(x); });
    }

    static MatrixType inv_sqrt_sym(const MatrixType& S) {
        return detail::apply_eigen_sym(S, [](T x) { using std::sqrt; return T{1} / sqrt(x); });
    }

    static MatrixType inv_sym(const MatrixType& S) {
        return detail::apply_eigen_sym(S, [](T x) { return T{1} / x; });
    }
};

static_assert(Manifold<SPDAffineInvariant<2>>);
static_assert(RiemannianManifold<SPDAffineInvariant<2>>);
static_assert(Manifold<SPDAffineInvariant<3>>);
static_assert(RiemannianManifold<SPDAffineInvariant<3>>);

// Affine-invariant Fréchet mean (Karcher mean): unlike log-Euclidean's
// frechet_mean() above, this has NO closed form -- the defining tradeoff of
// the fuller metric, per SPDAffineInvariant's file comment. Fixed-point
// iteration (Pennec 2006): repeatedly walk the current estimate along the
// average of log_map() to every sample, until that average tangent vector
// is ~zero. No ambient projection step is needed here (unlike
// algebra::riemannian_minimize() on Sphere/Hyperbolic), since this space's
// tangent space is already the unconstrained space of symmetric matrices --
// hand-rolled rather than calling riemannian_minimize(), which requires
// HasNormal, deliberately absent here (see the class comment above).
template<std::size_t N, Scalar T = double>
Matrix<T, N, N> frechet_mean_affine_invariant(std::span<const Matrix<T, N, N>> matrices,
                                               std::size_t max_iters = 50,
                                               T tol = T{1e-12}) {
    using S = SPDAffineInvariant<N, T>;
    S space;
    typename S::PointType mean = matrices[0]; // any SPD point as the starting estimate
    for (std::size_t iter = 0; iter < max_iters; ++iter) {
        typename S::TangentVector avg{};
        for (const auto& m : matrices) avg = avg + space.log_map(mean, m);
        avg = avg * (T{1} / T(matrices.size()));
        if (space.metric_at(mean, avg, avg) < tol * tol) break;
        mean = space.exp_map(mean, avg, T{1});
    }
    return mean;
}

} // namespace spatium
