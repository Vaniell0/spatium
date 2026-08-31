#pragma once

// Geometric continuum mechanics scaffolding.
//
// Yavari et al. (2010+) reformulate continuum mechanics in the language of
// Riemannian geometry: the body is a manifold (B, G), space is a manifold
// (S, g), the deformation φ: B → S is a smooth map, and strain emerges as
// the pullback metric C = φ*g compared to the reference G.
//
// This header ships the *types and concepts* needed to start building on
// that picture. Concrete solvers (FEM, XPBD, IPC) are a future addition.
//
// What's here:
//   - `DeformationMap<MFrom, MTo>` — wraps an arbitrary smooth map between
//     two manifolds with auto-differentiation hooks for Jacobian.
//   - `deformation_gradient` — the Jacobian Dφ as a matrix (finite
//     differences on the embedded coordinates).
//   - `pullback_metric` — C = (Dφ)ᵀ · g · (Dφ).
//   - `green_strain` — E = ½(C − I).
//   - `StrainEnergy` concept — anything that maps deformation gradient
//     to a scalar (Saint-Venant-Kirchhoff, Neo-Hookean, custom).
//   - `saint_venant_kirchhoff_energy` — classic linear elasticity in the
//     geometric formulation.
//
// References:
//   - Yavari A. *A Geometric Theory of Thermal Stresses*, J. Math. Phys. 51 (2010).
//   - Goriely A. *The Mathematics and Mechanics of Biological Growth* (Springer 2017).
//   - Marsden & Hughes *Mathematical Foundations of Elasticity* (Dover 1994).

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/algebra/matrix.hpp>
#  include <spatium/core/concepts.hpp>
#  include <utility>
#endif

SPATIUM_EXPORT namespace spatium::physics::mechanics {

// ── DeformationMap: φ : MFrom → MTo ─────────────────────────────
// Wraps a smooth map between two manifolds (typically Euclidean<3> for both
// in classical elasticity). The callable `Phi` is stored by value as a
// template parameter so the optimiser inlines it at every call site —
// avoiding the runtime indirection a `std::function<…>` would impose
// (see docs/concept-driven-physics.md §1.1 for the broader contract).
// Construction is most ergonomic through `make_deformation_map`, which
// captures the callable's deduced type for the user.

template<typename MFrom, typename MTo, typename Phi, Scalar T = double>
struct DeformationMap {
    using FromPoint = typename MFrom::PointType;
    using ToPoint   = typename MTo::PointType;
    static constexpr std::size_t dimFrom = MFrom::dimension;
    static constexpr std::size_t dimTo   = MTo::dimension;

    Phi phi;

    constexpr DeformationMap() = default;
    template<typename F>
    explicit constexpr DeformationMap(F&& f) : phi(std::forward<F>(f)) {}

    ToPoint apply(const FromPoint& X) const { return phi(X); }
};

template<typename MFrom, typename MTo, typename F, Scalar T = double>
constexpr auto make_deformation_map(F&& phi) {
    return DeformationMap<MFrom, MTo, std::decay_t<F>, T>{std::forward<F>(phi)};
}

// Deformation gradient F = Dφ at point X via central finite differences.
// Returns a (dimTo × dimFrom) matrix.
template<typename MFrom, typename MTo, typename Phi, Scalar T = double>
Matrix<T, MTo::PointType::size, MFrom::PointType::size>
deformation_gradient(const DeformationMap<MFrom, MTo, Phi, T>& d,
                     const typename MFrom::PointType& X,
                     T h = T{1e-6})
{
    using FromPoint = typename MFrom::PointType;
    constexpr std::size_t Nf = MFrom::PointType::size;
    constexpr std::size_t Nt = MTo::PointType::size;
    Matrix<T, Nt, Nf> F;

    for (std::size_t j = 0; j < Nf; ++j) {
        FromPoint Xp = X;  Xp[j] += h;
        FromPoint Xm = X;  Xm[j] -= h;
        auto yp = d.phi(Xp);
        auto ym = d.phi(Xm);
        for (std::size_t i = 0; i < Nt; ++i)
            F(i, j) = (yp[i] - ym[i]) / (T{2} * h);
    }
    return F;
}

// Right Cauchy-Green deformation tensor: C = Fᵀ F.
// In Yavari's geometric formulation this IS the pullback of the spatial
// metric g (= identity in flat embedding) via Dφ:  C = (Dφ)ᵀ · g · (Dφ).
template<Scalar T, std::size_t Nt, std::size_t Nf>
Matrix<T, Nf, Nf> right_cauchy_green(const Matrix<T, Nt, Nf>& F) {
    return F.transpose() * F;
}

// Green-Lagrange strain tensor: E = ½(C − I).
// Vanishes for rigid motions, quadratic in displacement gradient. Standard
// finite-strain measure in continuum mechanics.
template<Scalar T, std::size_t N>
Matrix<T, N, N> green_strain(const Matrix<T, N, N>& C) {
    auto I = Matrix<T, N, N>::identity();
    return Matrix<T, N, N>{(C - I) * T{0.5}};
}

// ── StrainEnergy concept ─────────────────────────────────────
// A strain energy density W(F) maps a deformation gradient to a scalar
// energy per unit reference volume. Classical examples: SVK, Neo-Hookean,
// Mooney-Rivlin. Custom hyperelastic models supply only this functional.

template<typename W, typename T, std::size_t N>
concept StrainEnergy = requires(W w, const Matrix<T, N, N>& F) {
    { w(F) } -> std::convertible_to<T>;
};

// Saint-Venant-Kirchhoff energy: W = ½ λ (tr E)² + μ tr(E²).
// Geometric linear elasticity: stress is linear in strain E, but E itself
// is nonlinear in F. Suitable for small-strain large-rotation problems.
template<Scalar T, std::size_t N>
struct SaintVenantKirchhoff {
    T lame_lambda;
    T lame_mu;

    T operator()(const Matrix<T, N, N>& F) const {
        auto C = right_cauchy_green(F);
        auto E = green_strain(C);

        T trE  = T{0};
        T trEE = T{0};
        for (std::size_t i = 0; i < N; ++i) {
            trE += E(i, i);
            for (std::size_t j = 0; j < N; ++j)
                trEE += E(i, j) * E(j, i);
        }
        return T{0.5} * lame_lambda * trE * trE + lame_mu * trEE;
    }
};

static_assert(StrainEnergy<SaintVenantKirchhoff<double, 3>, double, 3>);

} // namespace spatium::physics::mechanics
