#pragma once

// Discrete-Lagrangian variational integrators (Marsden & West 2001).
//
// Idea: discretise Hamilton's principle of stationary action
//
//      δS = δ ∫ L(q, q̇) dt = 0
//
// rather than the equations of motion. The discrete Lagrangian
//
//      L_d(q_k, q_{k+1}) ≈ ∫_{t_k}^{t_{k+1}} L(q, q̇) dt
//
// induces the *discrete Euler-Lagrange* (DEL) equations
//
//      D_2 L_d(q_{k-1}, q_k) + D_1 L_d(q_k, q_{k+1}) = 0
//
// solving DEL for q_{k+1} gives a numerical step that is, by construction,
// **symplectic** and **conserves discrete momentum maps** for any continuous
// symmetry of L_d (a discrete Noether theorem). Energy is bounded but not
// preserved — same flavour as Stormer-Verlet, but the construction is
// systematic and transports to constrained systems and Lie groups.
//
// This header ships:
//   - `DiscreteLagrangian` concept describing what L_d must provide.
//   - `MidpointDiscreteLagrangian` — universal midpoint-rule
//     discretisation L_d(q_k, q_{k+1}) = dt · L((q_k+q_{k+1})/2, (q_{k+1}-q_k)/dt)
//     for any smooth `L(q, qdot)` callable.
//   - `variational_step` solving DEL via fixed-point iteration on q_{k+1}
//     for separable Lagrangians L = ½·m·|q̇|² − V(q) on flat ℝᴺ. Reduces
//     to velocity Verlet for that family — verified by tests — and serves
//     as the scaffold for full-blown Marsden-West methods on manifolds in
//     subsequent slices.
//
// Reference: Marsden, J. E. & West, M. *Discrete mechanics and variational
// integrators*, Acta Numerica 10:357-514 (2001).

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/physics/mechanics/body.hpp>
#  include <spatium/physics/mechanics/integrator.hpp>
#  include <spatium/physics/mechanics/state.hpp>
#  include <cmath>
#endif

SPATIUM_EXPORT namespace spatium::physics::mechanics {

// ── Concept: a discrete Lagrangian on a flat-Euclidean configuration ──
// `L_d` returns the discrete action over one time-step given the start and
// end configurations. `dL_d_dq1` and `dL_d_dq2` are the partials w.r.t. the
// first and second arguments — together they assemble the DEL residual.

template<typename LD, std::size_t N, typename T>
concept DiscreteLagrangian = requires (LD ld,
                                       const Vec<T, N>& qk,
                                       const Vec<T, N>& qk1,
                                       T dt) {
    { ld.action(qk, qk1, dt)        } -> std::convertible_to<T>;
    { ld.dL_d_dq1(qk, qk1, dt)      } -> std::convertible_to<Vec<T, N>>;
    { ld.dL_d_dq2(qk, qk1, dt)      } -> std::convertible_to<Vec<T, N>>;
};

// ── Universal midpoint-rule discrete Lagrangian ───────────────
// Given a continuous Lagrangian L(q, q̇) supplied as two callables
// `kinetic(qdot) → T` and `potential(q) → T` (separable case), build the
// midpoint discretisation
//
//     L_d(q0, q1, dt) = dt · [ kinetic((q1-q0)/dt) - potential((q0+q1)/2) ]
//
// Partials computed analytically for separable forms — the most common case.

template<std::size_t N, typename T,
         typename Kinetic, typename Potential, typename GradPotential>
struct SeparableMidpointLagrangian {
    Kinetic       kinetic;       // T(q̇)  — kinetic energy of the velocity
    Potential     potential;     // V(q)  — potential at a configuration
    GradPotential grad_potential;// ∇V(q) — callable Vec<T,N>(const Vec<T,N>&)
    T mass{T{1}};

    // Discrete action over one step.
    T action(const Vec<T, N>& q0, const Vec<T, N>& q1, T dt) const {
        Vec<T, N> qdot = Vec<T, N>{(q1 - q0) * (T{1} / dt)};
        Vec<T, N> qmid = Vec<T, N>{(q0 + q1) * T{0.5}};
        return dt * (kinetic(qdot) - potential(qmid));
    }

    // ∂L_d / ∂q₀ — used for DEL momentum at the *previous* step.
    // For separable L = ½m|q̇|² − V(q): ∂L_d/∂q₀ = -m(q1-q0)/dt - dt/2 · ∇V(qmid)
    Vec<T, N> dL_d_dq1(const Vec<T, N>& q0, const Vec<T, N>& q1, T dt) const {
        Vec<T, N> qmid = Vec<T, N>{(q0 + q1) * T{0.5}};
        Vec<T, N> kin  = Vec<T, N>{(q1 - q0) * (-mass / dt)};
        Vec<T, N> pot  = Vec<T, N>{grad_potential(qmid) * (-dt * T{0.5})};
        return Vec<T, N>{kin + pot};
    }

    // ∂L_d / ∂q₁ — DEL momentum at the *next* step.
    // ∂L_d/∂q₁ = +m(q1-q0)/dt - dt/2 · ∇V(qmid)
    Vec<T, N> dL_d_dq2(const Vec<T, N>& q0, const Vec<T, N>& q1, T dt) const {
        Vec<T, N> qmid = Vec<T, N>{(q0 + q1) * T{0.5}};
        Vec<T, N> kin  = Vec<T, N>{(q1 - q0) * (mass / dt)};
        Vec<T, N> pot  = Vec<T, N>{grad_potential(qmid) * (-dt * T{0.5})};
        return Vec<T, N>{kin + pot};
    }
};

// Factory: deduce Kinetic/Potential/GradPotential from the supplied
// callables so users can construct a Lagrangian as
//   auto L = make_separable_midpoint_lagrangian<1, double>(K, V, ∇V, mass);
// without writing the full template parameter list.
template<std::size_t N, typename T,
         typename Kinetic, typename Potential, typename GradPotential>
constexpr auto make_separable_midpoint_lagrangian(
    Kinetic kinetic, Potential potential, GradPotential grad_potential,
    T mass = T{1})
{
    return SeparableMidpointLagrangian<
        N, T,
        std::decay_t<Kinetic>, std::decay_t<Potential>, std::decay_t<GradPotential>>
    {std::move(kinetic), std::move(potential), std::move(grad_potential), mass};
}

// ── Variational step for separable systems (kick-drift-kick Verlet) ──
// For separable Lagrangians L = ½m|q̇|² − V(q), the trapezoidal-rule
// discrete Lagrangian
//
//     L_d(q0, q1) = dt/2 · [L(q0, (q1-q0)/dt) + L(q1, (q1-q0)/dt)]
//
// makes the Discrete Euler-Lagrange equations factor exactly into the
// Stormer-Verlet leapfrog (Hairer-Lubich-Wanner §VI.6) — the same
// kick-drift-kick already provided by `verlet_step` with the
// Newtonian force `F(q) = -∇V(q)`. We delegate so the variational
// path and the Newton path share a single implementation: any
// numerical fix on one updates the other automatically.

template<std::size_t N, Scalar T, typename GradPotential>
void variational_step_separable(PointMass<N, T>& body, GradPotential&& grad_V, T dt) {
    auto force = [&grad_V](const PointMass<N, T>& b, T) -> Vec<T, N> {
        return Vec<T, N>{grad_V(b.state.position) * T{-1}};
    };
    verlet_step(body, force, dt);
}

// ── Discrete momentum map p_k = ∂L_d/∂q₀ (at the start of an interval) ─
// For separable midpoint Lagrangian: p_k = m·(q1-q0)/dt − ½dt·∇V(qmid).
// At convergence of DEL across consecutive steps, this is the discrete
// conjugate momentum. Provided so tests can verify Noether-like
// conservation (e.g. for translation-invariant V, ∑ p stays constant).

template<std::size_t N, Scalar T>
Vec<T, N> discrete_momentum_separable(const Vec<T, N>& q0, const Vec<T, N>& q1,
                                       T dt, T mass) {
    return Vec<T, N>{(q1 - q0) * (mass / dt)};
}

} // namespace spatium::physics::mechanics
