#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/concepts.hpp>
#  include <spatium/algebra/dual.hpp>
#  include <spatium/algebra/linear_solve.hpp>
#  include <spatium/algebra/matrix.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <cstddef>
#endif

SPATIUM_EXPORT namespace spatium {
inline namespace algebra {

// An object described as a function, not a hand-derived formula: any plain
// callable F: Domain -> Codomain (lambda, function pointer, functor) — no
// wrapper type required. gradient()/integrate() below work on anything
// matching this shape directly.
template<typename F, typename Domain, typename Codomain>
concept Function = requires(F f, Domain x) {
    { f(x) } -> std::convertible_to<Codomain>;
};

// Gradient of a scalar field f: Vec<T,N> -> T at point p, via forward-mode AD:
// one Dual-seeded evaluation of f per axis (N evaluations for N dimensions).
// f must be generic over its scalar type (e.g. a templated lambda) so it can
// be called with Vec<Dual<T>,N> here and with Vec<T,N> anywhere else — the
// same f serves both plain evaluation and differentiation, unchanged.
// Typical use: surface normal of an implicit surface F(p) = 0 is
// gradient(F, p).normalized().
template<Scalar T, std::size_t N, typename F>
    requires Function<F, Vec<Dual<T>, N>, Dual<T>>
Vec<T, N> gradient(F&& f, const Vec<T, N>& p) {
    Vec<T, N> g{};
    for (std::size_t i = 0; i < N; ++i) {
        Vec<Dual<T>, N> dp;
        for (std::size_t j = 0; j < N; ++j)
            dp[j] = (j == i) ? Dual<T>::variable(p[j]) : Dual<T>::constant(p[j]);
        g[i] = f(dp).deriv;
    }
    return g;
}

namespace calculus_detail {

// One Simpson panel over [a,b] given f at the two endpoints and the midpoint.
template<Scalar T>
constexpr T simpson_panel(T a, T b, T fa, T fm, T fb) {
    return (b - a) / T{6} * (fa + T{4} * fm + fb);
}

// Adaptive refinement: recurse into a half only if its Simpson estimate
// doesn't already agree with the whole-panel estimate to within eps.
// `depth` bounds recursion so a pathological f can't spin forever.
template<Scalar T, typename F>
T adaptive_simpson(F&& f, T a, T b, T fa, T fm, T fb, T whole, T eps, int depth) {
    T m = (a + b) / T{2};
    T lm = (a + m) / T{2};
    T rm = (m + b) / T{2};
    T flm = f(lm);
    T frm = f(rm);
    T left = simpson_panel(a, m, fa, flm, fm);
    T right = simpson_panel(m, b, fm, frm, fb);

    using std::abs;
    if (depth <= 0 || abs(left + right - whole) <= T{15} * eps)
        return left + right + (left + right - whole) / T{15};

    return adaptive_simpson(f, a, m, fa, flm, fm, left, eps / T{2}, depth - 1)
         + adaptive_simpson(f, m, b, fm, frm, fb, right, eps / T{2}, depth - 1);
}

} // namespace calculus_detail

// Definite integral of f: T -> T over [a,b] via adaptive Simpson's rule.
// Scoped to 1-D on purpose — region-aware integration over Spatium's own
// geometry (Box/Sphere/Polygon volumes, intersection overlap) is a separate,
// bigger next step, not folded in here.
template<Scalar T, typename F>
    requires Function<F, T, T>
T integrate(F&& f, T a, T b, T eps = epsilon<T>() * T{1000}) {
    T fa = f(a);
    T fb = f(b);
    T fm = f((a + b) / T{2});
    T whole = calculus_detail::simpson_panel(a, b, fa, fm, fb);
    return calculus_detail::adaptive_simpson(f, a, b, fa, fm, fb, whole, eps, 20);
}

// Calibration, not a trained model: minimizes a scalar loss f: Vec<T,N> -> T
// from an initial guess theta0, via gradient descent with Armijo
// backtracking line search. f must be generic over its scalar type — same
// requirement as gradient() — so it can be evaluated both with Vec<T,N>
// (line search) and Vec<Dual<T>,N> (gradient, via the function above).
//
// Re-solves from scratch on every call, no learning or generalization
// across cases — e.g. fitting a contact solver's (compliance, friction,
// damping) to eliminate penetration in one specific scene. "The model"
// (RL/search-trained) only ever covers which op/family to calibrate here,
// not this fit.
template<Scalar T, std::size_t N, typename F>
    requires Function<F, Vec<Dual<T>, N>, Dual<T>> && Function<F, Vec<T, N>, T>
Vec<T, N> minimize(F&& f, Vec<T, N> theta,
                    T grad_tol = epsilon<T>() * T{1000}, int max_iters = 200) {
    using std::sqrt;
    const T armijo_c = T{1} / T{10000};

    for (int iter = 0; iter < max_iters; ++iter) {
        Vec<T, N> g = gradient(f, theta);
        T g2 = g.dot(g);
        if (sqrt(g2) < grad_tol) break;

        T f0 = f(theta);
        T step = T{1};
        for (int ls = 0; ls < 50; ++ls) {
            Vec<T, N> candidate = theta - g * step;
            if (f(candidate) <= f0 - armijo_c * step * g2) {
                theta = candidate;
                break;
            }
            step = step / T{2};
        }
    }
    return theta;
}

// ── Riemannian optimization ─────────────────────────────────────

// gradient() above returns the ordinary coordinate gradient of f in the
// ambient basis: components d f / d p_i, i.e. a covector, not yet a vector
// in the space's own ambient bilinear form. raise_gradient() converts it
// into the g satisfying metric_at(p, g, e_i) == d f / d p_i for every
// ambient basis vector e_i -- index-raising with the inverse of that form.
// For a Euclidean ambient metric (Euclidean, Sphere) this is the identity;
// for a Minkowski one (Hyperbolic) it isn't -- same reason
// physics/relativity/geodesic.hpp inverts g_{mu nu} via the general
// invert() rather than assuming a diagonal Euclidean shortcut, and needed
// for the same reason here: skipping it silently returns a vector that
// only happens to be right when the ambient metric is Euclidean.
template<typename S>
    requires RiemannianManifold<S>
typename S::TangentVector raise_gradient(const S& space, const typename S::PointType& p,
                                         const typename S::TangentVector& ambient_grad) {
    using T = typename S::ScalarType;
    constexpr auto M = S::PointType::size;

    Matrix<T, M, M> B{};
    for (std::size_t i = 0; i < M; ++i) {
        typename S::TangentVector ei{};
        ei[i] = T{1};
        for (std::size_t j = 0; j < M; ++j) {
            typename S::TangentVector ej{};
            ej[j] = T{1};
            B(i, j) = space.metric_at(p, ei, ej);
        }
    }
    auto inv = invert(B);
    Matrix<T, M, M> Binv = inv ? *inv : Matrix<T, M, M>{};
    return Binv * ambient_grad;
}

// Projects an ambient-space vector v at point p onto the tangent space of
// any Surface (HasNormal) + RiemannianManifold, via the space's own
// metric_at() rather than an assumed-orthonormal Euclidean normal:
//   v_tan = v - n * <n,v> / <n,n>,  n = space.normal(p)
// This is what makes it correct unchanged for both a Euclidean-metric
// space like Sphere and a Minkowski-metric one like Hyperbolic.
template<typename S>
    requires RiemannianManifold<S> && HasNormal<S>
typename S::TangentVector project_tangent(const S& space,
                                          const typename S::PointType& p,
                                          const typename S::TangentVector& v) {
    auto n = space.normal(p);
    auto denom = space.metric_at(p, n, n);
    return v - n * (space.metric_at(p, n, v) / denom);
}

// Riemannian gradient descent: minimizes a scalar field f defined over a
// manifold's ambient embedding, via retraction (exp_map) instead of a flat
// update -- every iterate lands exactly back on the manifold rather than
// drifting off and needing a separate projection step. Reuses the same
// Dual<T>-seeded gradient() as minimize() above; raise_gradient() +
// project_tangent() convert its ambient covector into the actual
// Riemannian gradient before retracting along the manifold.
template<typename S, typename F>
    requires RiemannianManifold<S> && HasNormal<S>
          && Function<F, Vec<Dual<typename S::ScalarType>, S::PointType::size>, Dual<typename S::ScalarType>>
          && Function<F, typename S::PointType, typename S::ScalarType>
typename S::PointType riemannian_minimize(const S& space, F&& f, typename S::PointType theta,
                                          typename S::ScalarType grad_tol
                                              = epsilon<typename S::ScalarType>() * typename S::ScalarType{1000},
                                          int max_iters = 200) {
    using T = typename S::ScalarType;
    using std::sqrt; using std::abs;
    const T armijo_c = T{1} / T{10000};

    for (int iter = 0; iter < max_iters; ++iter) {
        auto g = project_tangent(space, theta, raise_gradient(space, theta, gradient(f, theta)));
        T g2 = abs(space.metric_at(theta, g, g));
        if (sqrt(g2) < grad_tol) break;

        T f0 = f(theta);
        T step = T{1};
        for (int ls = 0; ls < 50; ++ls) {
            auto candidate = space.exp_map(theta, g, -step);
            if (f(candidate) <= f0 - armijo_c * step * g2) {
                theta = candidate;
                break;
            }
            step = step / T{2};
        }
    }
    return theta;
}

} // namespace algebra
} // namespace spatium
