#pragma once

// Lie-group integrators for systems whose configuration lives on a Lie group
// (rigid body on SO(3)/SE(3), satellite attitude, articulated robotics).
//
// Reference: Celledoni-Çokaj-Leone-Murari-Owren, *"Lie group integrators for
// mechanical systems"*, IJCM (arXiv:2102.12778, 2021); Müller, *"Evaluation
// and implementation of Lie group integration methods for rigid multibody
// systems"*, MSD 2024.
//
// The classic Munthe-Kaas (RKMK) recipe takes any explicit Runge-Kutta method
// for ODEs in the algebra g and lifts it to the group G via `exp`. This
// header ships two methods that exercise the existing SO(3)/SE(3)
// infrastructure without yet pulling in the full BCH machinery needed for
// fourth-order RKMK:
//
//   lie_euler_step      — first-order, single `exp` per step.
//   lie_midpoint_step   — second-order predictor-corrector, two `exp` calls.
//
// The vector-field is supplied as `omega(g) -> AlgebraType` (e.g. for
// torque-free rigid body: `omega(R) = I_inv * (I*ω) × ω` rewritten through R).
//
// `exp_step(g, xi)` defaults to the group's `exp(xi)` left-multiplied onto g
// (so g_{k+1} = g_k · exp(xi · dt)). Callers can override for right-trivialised
// flows by passing a custom `Multiply` policy.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/concepts.hpp>
#  include <spatium/algebra/groups/so3.hpp>
#  include <spatium/algebra/groups/se3.hpp>
#endif

SPATIUM_EXPORT namespace spatium::physics::mechanics {

// Default group action: G_{k+1} = G_k · exp(xi).
struct LeftTrivialised {
    template<typename Group, typename Element, typename Algebra>
    static Element advance(const Group& G, const Element& g, const Algebra& xi) {
        return G.compose(g, G.exp(xi));
    }
};

// First-order Lie-Euler step.
//
//   g_{k+1} = g_k · exp(dt · ω(g_k))
//
template<spatium::algebra::LieGroup G,
         typename Element  = typename G::ElementType,
         typename Algebra  = typename G::AlgebraType,
         typename VectorField,
         typename Scalar   = double,
         typename Multiply = LeftTrivialised>
Element lie_euler_step(const G& group, Element g, Scalar dt, VectorField&& omega,
                       Multiply mult = {})
{
    Algebra xi = omega(g);                  // tangent at g, expressed in g
    // dt * xi: assume Algebra has scalar multiplication. For SO(3)/SE(3)
    // AlgebraType is a Vec, so multiplication is element-wise.
    Algebra step = xi * dt;
    return mult.advance(group, g, step);
}

// Second-order Lie-midpoint (a.k.a. Lie-trapezoidal). Predictor with Euler,
// corrector by averaging slopes at start and predicted endpoint:
//
//   ξ_1 = ω(g_k)
//   ḡ   = g_k · exp(½·dt·ξ_1)
//   ξ_2 = ω(ḡ)
//   g_{k+1} = g_k · exp(dt · ξ_2)
//
// Equivalent to RK2 (midpoint) lifted to the group.
template<spatium::algebra::LieGroup G,
         typename Element  = typename G::ElementType,
         typename Algebra  = typename G::AlgebraType,
         typename VectorField,
         typename Scalar   = double,
         typename Multiply = LeftTrivialised>
Element lie_midpoint_step(const G& group, Element g, Scalar dt,
                          VectorField&& omega, Multiply mult = {})
{
    Algebra xi1   = omega(g);
    Algebra half  = xi1 * (dt * Scalar{0.5});
    Element g_mid = mult.advance(group, g, half);
    Algebra xi2   = omega(g_mid);
    Algebra step  = xi2 * dt;
    return mult.advance(group, g, step);
}

// ── Commutator-free Runge-Kutta-Munthe-Kaas of order 4 ─────────
// Celledoni-Marthinsen-Owren (2003), "Commutator-free Lie group methods".
// Avoids explicit BCH commutators by composing exp's at each stage.
//
//   k_1 = ω(g_n)
//   k_2 = ω(g_n · exp(½·dt·k_1))
//   k_3 = ω(g_n · exp(½·dt·k_2))
//   k_4 = ω(g_n · exp(    dt·k_3))
//   g_{n+1} = g_n · exp(½·dt·(k_1+k_2)/2) · exp(½·dt·(k_3+k_4)/2)
//             — TWO half-step exponentials (not one full step) to keep
//               the order-4 truncation error.
//
// Cost: 4 force-field evaluations + 5 exp calls per step (vs midpoint's
// 2 + 2). Worth it on stiff or chaotic Lie-group ODEs where the ~Δt⁴
// global accuracy lets you take 5-10× larger steps.

template<spatium::algebra::LieGroup G,
         typename Element  = typename G::ElementType,
         typename Algebra  = typename G::AlgebraType,
         typename VectorField,
         typename Scalar   = double,
         typename Multiply = LeftTrivialised>
Element lie_rkmk4_cf_step(const G& group, Element g, Scalar dt,
                          VectorField&& omega, Multiply mult = {})
{
    Algebra k1 = omega(g);

    Algebra h1   = k1 * (dt * Scalar{0.5});
    Element g2   = mult.advance(group, g, h1);
    Algebra k2 = omega(g2);

    Algebra h2   = k2 * (dt * Scalar{0.5});
    Element g3   = mult.advance(group, g, h2);
    Algebra k3 = omega(g3);

    Algebra h3   = k3 * dt;
    Element g4   = mult.advance(group, g, h3);
    Algebra k4 = omega(g4);

    // Two half-step exponentials of averaged stages — the commutator-free
    // way to assemble the order-4 increment without computing brackets.
    Algebra avg12 = Algebra{(k1 + k2) * (dt * Scalar{0.25})};
    Algebra avg34 = Algebra{(k3 + k4) * (dt * Scalar{0.25})};
    Element g_half = mult.advance(group, g, avg12);
    return mult.advance(group, g_half, avg34);
}

} // namespace spatium::physics::mechanics
