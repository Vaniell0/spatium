#pragma once

// Time-step integrators for PointMass-based flat-space mechanics.
//
// Each `step(body, force, dt, t)` advances `body.state` in place by one
// time step. `force` is anything callable as `Vec<T,N>(const PointMass&, T)`.
//
// Methods:
//   euler_step          — explicit Euler. O(Δt). Accumulates energy fast.
//   rk4_step            — classical 4-stage Runge-Kutta. O(Δt⁴). Non-symplectic.
//   semi_implicit_euler — (a.k.a. symplectic Euler). v ← v + a·dt; x ← x + v·dt.
//                         Symplectic — bounded energy error for Hamiltonian systems.
//   verlet_step         — velocity Verlet. O(Δt²) symplectic. Reversible.
//
// Symplectic methods (semi_implicit_euler, verlet) preserve a modified
// Hamiltonian, so total energy oscillates within O(Δt^p) bounds forever
// instead of drifting linearly. This is why they dominate orbital
// integration over long timespans.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/physics/mechanics/body.hpp>
#  include <spatium/physics/mechanics/force.hpp>
#  include <cstddef>
#endif

SPATIUM_EXPORT namespace spatium::physics::mechanics {

// ── Explicit Euler ────────────────────────────────────────────
template<std::size_t N, Scalar T, typename Force>
void euler_step(PointMass<N, T>& body, Force&& force, T dt, T t = T{0}) {
    Vec<T, N> a = Vec<T, N>{force(body, t) / body.mass};
    body.state.position = Vec<T, N>{body.state.position + body.state.velocity * dt};
    body.state.velocity = Vec<T, N>{body.state.velocity + a * dt};
}

// ── Semi-implicit (symplectic) Euler ──────────────────────────
template<std::size_t N, Scalar T, typename Force>
void semi_implicit_euler_step(PointMass<N, T>& body, Force&& force, T dt, T t = T{0}) {
    Vec<T, N> a = Vec<T, N>{force(body, t) / body.mass};
    body.state.velocity = Vec<T, N>{body.state.velocity + a * dt};
    body.state.position = Vec<T, N>{body.state.position + body.state.velocity * dt};
}

// ── Velocity Verlet ───────────────────────────────────────────
// Standard form:
//   x(t+dt) = x(t) + v(t)·dt + ½·a(t)·dt²
//   v(t+dt) = v(t) + ½·(a(t) + a(t+dt))·dt
// Two force evaluations per step but bounded energy drift, time-reversible.
template<std::size_t N, Scalar T, typename Force>
void verlet_step(PointMass<N, T>& body, Force&& force, T dt, T t = T{0}) {
    Vec<T, N> a0 = Vec<T, N>{force(body, t) / body.mass};
    body.state.position = Vec<T, N>{body.state.position
                                    + body.state.velocity * dt
                                    + a0 * (T{0.5} * dt * dt)};
    Vec<T, N> a1 = Vec<T, N>{force(body, t + dt) / body.mass};
    body.state.velocity = Vec<T, N>{body.state.velocity
                                    + (a0 + a1) * (T{0.5} * dt)};
}

// ── Classical RK4 ─────────────────────────────────────────────
// 4-stage Runge-Kutta on the (x, v) phase state. Non-symplectic; great for
// short-time accuracy, slowly accumulates energy on long-time orbits.
template<std::size_t N, Scalar T, typename Force>
void rk4_step(PointMass<N, T>& body, Force&& force, T dt, T t = T{0}) {
    auto eval = [&](const PhaseState<N, T>& s, T tt) -> Derivative<N, T> {
        PointMass<N, T> tmp{body.mass, s.position, s.velocity};
        Vec<T, N> a = Vec<T, N>{force(tmp, tt) / body.mass};
        return Derivative<N, T>{s.velocity, a};
    };

    PhaseState<N, T> s0 = body.state;

    Derivative<N, T> k1 = eval(s0, t);
    Derivative<N, T> k2 = eval(s0 + k1 * (T{0.5} * dt), t + T{0.5} * dt);
    Derivative<N, T> k3 = eval(s0 + k2 * (T{0.5} * dt), t + T{0.5} * dt);
    Derivative<N, T> k4 = eval(s0 + k3 *  dt,           t +          dt);

    PhaseState<N, T> ds = (k1 + (k2 + k3) * T{2} + k4) * (dt / T{6});
    body.state = s0 + ds;
}

// ── Yoshida 4th-order symplectic composition ──────────────────
// Three Verlet sub-steps with magic coefficients that cancel the 2nd-order
// error term, giving overall O(Δt⁴) symplectic accuracy. Same family as
// Forest-Ruth. References: Yoshida, Phys. Lett. A 150 (1990) 262.
//
//   w  = 1 / (2 - 2^(1/3))      ≈ 1.351207...
//   c1 = w · dt
//   c2 = (1 - 2w) · dt           ≈ -1.702415... · dt   (negative middle step)
//
// Three Verlet calls per Yoshida step (versus one for plain Verlet), so
// per-step cost is ~3× but global error scales as Δt⁴ instead of Δt² —
// same accuracy with ~10× larger Δt for typical orbital problems.
template<std::size_t N, Scalar T, typename Force>
void yoshida4_step(PointMass<N, T>& body, Force&& force, T dt, T t = T{0}) {
    constexpr T cube_root_2 = T{1.2599210498948732};   // 2^(1/3)
    const     T w  = T{1} / (T{2} - cube_root_2);
    const     T c1 = w * dt;
    const     T c2 = (T{1} - T{2} * w) * dt;
    verlet_step(body, force, c1, t);
    verlet_step(body, force, c2, t + c1);
    verlet_step(body, force, c1, t + c1 + c2);
}

} // namespace spatium::physics::mechanics
