#pragma once

// Symplectic manifold abstraction.
//
// A symplectic manifold is a pair (M, ω) where ω is a closed, non-degenerate
// 2-form. The cotangent bundle T*Q of any configuration manifold Q carries
// a canonical such ω = dq ∧ dp, making T*Q the natural "phase space" of
// classical mechanics.
//
// `SymplecticManifold<S>` concept lets integrators dispatch on whatever
// model of phase space the user supplies:
//
//   - flat:   `CotangentBundle<Euclidean<N, T>>` — q ∈ ℝᴺ, p ∈ ℝᴺ.
//   - sphere: `CotangentBundle<Sphere<N, T>>`    — q ∈ Sᴺ ⊂ ℝᴺ⁺¹, p ∈ T*Q.
//   - rigid body: `CotangentBundle<SE3>`          — pose ∈ SE(3), wrench in se(3)*.
//
// This header ships the concept and the cotangent-bundle wrapper. The
// canonical symplectic form is queried via `omega(state)` (returns a function
// of two tangent vectors → scalar), useful for `verify_symplecticity` tests.
// The variational integrators in `variational.hpp` consume this same
// abstraction.
//
// References:
//   - Marsden & Ratiu, *Introduction to Mechanics and Symmetry* (Springer 2nd ed.).
//   - Holm, *Geometric Mechanics* I.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/algebra/vector.hpp>
#endif

SPATIUM_EXPORT namespace spatium::physics::mechanics {

// ── Concept: a symplectic manifold ────────────────────────────
// A type `S` is a `SymplecticManifold` if it provides:
//   - `Configuration`  — element of the base manifold Q (typically a Vec or Lie-group elem).
//   - `Momentum`       — element of T*_q Q (covector); in Vec-based models same shape as Configuration.
//   - `State`          — pair (q, p) on T*Q.
//   - `dimension`      — total phase-space dim (= 2 · dim Q).
//   - `omega(s, dq1, dp1, dq2, dp2)` — symplectic 2-form acting on two tangent vectors.

template<typename S>
concept SymplecticManifold = requires {
    typename S::Configuration;
    typename S::Momentum;
    typename S::State;
    typename S::ScalarType;
    { S::dimension } -> std::convertible_to<std::size_t>;
} && requires(const typename S::Configuration& cfg,
              const typename S::Momentum& mom) {
    // The canonical symplectic 2-form must be reachable as a static
    // method on S: ω((δq₁, δp₁), (δq₂, δp₂)) → Scalar. Without this
    // requirement a struct could pass the concept while leaving the
    // form unimplemented, which would silently break
    // `verify_symplecticity_drift` and any downstream variational
    // integrator that consumes ω.
    { S::omega(cfg, mom, cfg, mom) } -> std::convertible_to<typename S::ScalarType>;
};

// ── Canonical cotangent bundle T*Q ────────────────────────────
// Wraps any configuration manifold M into its phase space with the canonical
// symplectic form ω = dq ∧ dp = Σᵢ dqᵢ ∧ dpᵢ.
//
// Shape conventions (flat vector base):
//   Configuration ≡ M::PointType
//   Momentum      ≡ M::TangentVector  (same shape; identified via metric for now)
//
// On Lie groups one would want `Momentum` to live in the dual algebra g*;
// `lgvi.hpp` already does this directly for SO(3) — extending that
// treatment to SE(3) is a natural future step.

template<typename M>
struct CotangentBundle {
    using Base          = M;
    using Configuration = typename M::PointType;
    using Momentum      = typename M::TangentVector;
    using ScalarType    = typename M::ScalarType;
    static constexpr std::size_t dimension = 2 * M::dimension;

    struct State {
        Configuration q{};
        Momentum      p{};

        constexpr State() = default;
        constexpr State(Configuration qq, Momentum pp)
            : q(std::move(qq)), p(std::move(pp)) {}
    };

    // Canonical symplectic form: ω((δq₁,δp₁), (δq₂,δp₂)) = δq₁·δp₂ − δq₂·δp₁.
    // Independent of the base point on the cotangent bundle (canonical 2-form).
    static constexpr ScalarType omega(
        const Configuration& dq1, const Momentum& dp1,
        const Configuration& dq2, const Momentum& dp2)
    {
        return dq1.dot(dp2) - dq2.dot(dp1);
    }
};

// ── Verify symplecticity of a one-step map numerically ─────────
// Given a step map Φ: T*Q → T*Q (e.g. a Verlet integrator), the discrete
// symplecticity check evaluates dΦ on two tangent perturbations and confirms
// that ω(dΦ·v₁, dΦ·v₂) ≈ ω(v₁, v₂) up to rounding.
//
// Returns the maximum absolute drift over a small set of random perturbations.
// `step_map(state, dt)` should advance the state by one Δt.
//
// This is the empirical test of Verlet/Yoshida4/RKMK4 being symplectic, the
// counterpart of the metric/InnerProduct verifiers in `verify.hpp`.

template<SymplecticManifold S, typename StepMap>
typename S::ScalarType verify_symplecticity_drift(
    const typename S::State& s0, StepMap&& step,
    typename S::ScalarType eps,
    typename S::ScalarType dt)
{
    using T = typename S::ScalarType;
    using Cfg = typename S::Configuration;
    using Mom = typename S::Momentum;

    // The perturbations must be a *conjugate* pair, and the version this
    // replaces was not. It used dq1 = (eps,0,..), dp1 = (0,eps,..),
    // dq2 = (0,eps,..), dp2 = (eps,0,..), for which
    //
    //     omega = dq1·dp2 - dq2·dp1 = eps^2 - eps^2 = 0
    //
    // identically. A 2-form that is zero on the pair being tested measures
    // nothing: preserving zero is what every map does. Measured against
    // deliberately non-symplectic maps -- one inflating momentum by 10% a
    // step, so phase-space volume grows by 1.1 each time, and one pure
    // shear -- it returned 1.7e-12 and 2.3e-12 against Verlet's 1.3e-12.
    // Indistinguishable, and `tests/test_symplectic_integrators.cpp` had
    // been asserting on it.
    //
    // A canonical pair in one coordinate -- perturb q alone, perturb the
    // conjugate p alone -- gives omega = eps^2 instead, which is the
    // largest value the form takes on unit-scale perturbations and so the
    // most sensitive probe available. Every coordinate is swept and the
    // worst plane is reported, because a map can preserve the form in one
    // plane while destroying it in another.
    using std::abs;
    T worst = T{0};

    for (std::size_t axis = 0; axis < Cfg::size; ++axis) {
        Cfg dq1{}; dq1[axis] = eps;      Mom dp1{};
        Cfg dq2{};                       Mom dp2{}; dp2[axis % Mom::size] = eps;

        const auto sP    = step(s0, dt);
        const auto sP_1  = step(typename S::State{Cfg{s0.q + dq1}, Mom{s0.p + dp1}}, dt);
        const auto sP_2  = step(typename S::State{Cfg{s0.q + dq2}, Mom{s0.p + dp2}}, dt);

        const Cfg dq1_after{sP_1.q - sP.q};
        const Mom dp1_after{sP_1.p - sP.p};
        const Cfg dq2_after{sP_2.q - sP.q};
        const Mom dp2_after{sP_2.p - sP.p};

        const T before = S::omega(dq1, dp1, dq2, dp2);
        const T after  = S::omega(dq1_after, dp1_after, dq2_after, dp2_after);

        // Relative to the form's own magnitude, so the answer does not
        // depend on the perturbation size chosen.
        worst = std::max(worst, abs(after - before) / (abs(before) + eps * eps));
    }
    return worst;
}

} // namespace spatium::physics::mechanics
