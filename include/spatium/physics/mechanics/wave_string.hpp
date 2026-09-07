#pragma once

// Explicit finite-difference time-stepping for the 1D wave equation with
// fixed (Dirichlet) boundary conditions -- a string of length `length`,
// clamped at both ends, given an initial displacement (a pluck) or
// velocity and released.
//
//   u_tt(x,t) = c^2 u_xx(x,t),      u(0,t) = u(length,t) = 0
//
// c = sqrt(tension / linear_density) is the wave speed. Deliberately NOT
// routed through algebra/ode.hpp's rk4_step: that solver targets a
// first-order IVP y' = f(t,y) over a small fixed-size Vec<T,N> state, and
// recasting a spatially discretized PDE into that shape would mean
// carrying a 2*node_count stacked (displacement, velocity) state through
// a generic RHS callable for no numerical benefit here. The natural tool
// for u_tt = c^2 u_xx is the standard explicit central-difference
// ("leapfrog") scheme below: second-order accurate in both space and
// time, needs only the previous two time levels, and (the same shape of
// idea as integrator.hpp's verlet_step) updates the whole field with one
// direct local stencil per step instead of multiple RHS evaluations.
// examples/wave_ca_demo.cpp already uses this exact discretization
// inline for a 2D angular-grid/sphere case; this header is the reusable,
// physically-parametrized 1D version, plus the closed-form frequency
// check to verify it against (see tests/test_wave_string.cpp).
//
// Stability: explicit central differences for the wave equation are only
// stable under the CFL condition c*dt/dx <= 1 (Courant, Friedrichs, Lewy
// 1928) -- courant_number() below reports it; violating it blows the
// simulation up rather than silently mis-integrating it, so it is cheap
// to guard against.
//
// Damping: `damping` (gamma, units 1/s) adds the standard velocity-
// proportional loss term to the PDE, u_tt = c^2 u_xx - 2*gamma*u_t,
// discretized centrally (see step()'s derivation in its own comment) --
// the textbook way to add a stand-in for the string's real energy losses
// (air resistance, internal friction, bridge coupling) without modelling
// any of them individually. An early version of this header instead
// multiplied the *entire* newly-computed state by a flat per-step
// (1 - damping) factor; that turned out NOT to be frequency-neutral --
// scaling the "2u^n - u^{n-1}" inertial term along with everything else
// shifts the recursion's characteristic roots by O(damping) in phase,
// not just magnitude (a measurable ~10% fundamental-frequency shift at
// damping values chosen only to fix the *audible* decay time). The
// velocity-proportional form below is the physically correct one: like
// any real damped oscillator its frequency shift is O(gamma^2), i.e.
// negligible for the small gamma a listenable decay needs, which is why
// it doesn't bias the closed-form frequency check in
// tests/test_wave_string.cpp. damping = 0 is physically valid too (an
// ideal string rings forever); it's just not what a listenable demo wants.
//
// Precondition: node_count >= 3 (at least one interior node); not
// checked at runtime -- same convention as algebra/eigen_decomp.hpp's
// symmetric-input precondition (see docs/conventions.md's Result<T>-vs-
// assert rule): an internal-invariant precondition on construction
// parameters, not boundary input from an untrusted source.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <cmath>
#  include <cstddef>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::physics::mechanics {

// ── Closed-form reference ────────────────────────────────────────
// Fundamental frequency of an ideal string under tension, clamped at
// both ends: f1 = c / (2*length) = sqrt(tension/linear_density) / (2*length).
template<Scalar T>
T string_fundamental_frequency(T tension, T linear_density, T length) {
    using std::sqrt;
    T c = sqrt(tension / linear_density);
    return c / (T{2} * length);
}

// nth harmonic (n=1 is the fundamental) of an ideal clamped string.
template<Scalar T>
T string_harmonic_frequency(T tension, T linear_density, T length, int n) {
    return T(n) * string_fundamental_frequency(tension, linear_density, length);
}

// c*dt/dx -- must stay <= 1 for VibratingString::step's explicit scheme
// to be numerically stable (Courant-Friedrichs-Lewy).
template<Scalar T>
T courant_number(T wave_speed, T dt, T dx) {
    return wave_speed * dt / dx;
}

// ── Simulation state ─────────────────────────────────────────────
// A string discretized into `node_count` points spaced dx =
// length/(node_count-1) apart, indices [0, node_count-1]; index 0 and
// node_count-1 are the fixed (Dirichlet) endpoints and never move.
// Build one via make_plucked_string() below rather than filling `current`
// / `previous` by hand -- getting the first step right without a real
// initial-velocity phantom point is easy to get only first-order-accurate.
template<Scalar T = double>
struct VibratingString {
    std::size_t node_count;
    T length;
    T wave_speed;         // c = sqrt(tension / linear_density)
    T damping = T{0};     // gamma in u_tt = c^2 u_xx - 2*gamma*u_t, units 1/s

    std::vector<T> current;   // u^n
    std::vector<T> previous;  // u^{n-1}

    [[nodiscard]] T dx() const { return length / T(node_count - 1); }

    // Advance one time step of size dt. dt should satisfy
    // courant_number(wave_speed, dt, dx()) <= 1.
    //
    // Derivation: central-difference u_tt = c^2 u_xx - 2*gamma*u_t with
    // u_t(t_n) approximated by the centered (u^{n+1}-u^{n-1})/(2dt):
    //
    //   (u^{n+1} - 2u^n + u^{n-1})/dt^2
    //       = c^2*laplacian(u^n)/h^2 - 2*gamma*(u^{n+1}-u^{n-1})/(2dt)
    //
    // Solving for u^{n+1} (the only unknown) gives the update below, which
    // reduces to the plain undamped leapfrog recurrence at damping = 0.
    void step(T dt) {
        const T h    = dx();
        const T r2   = (wave_speed * dt / h) * (wave_speed * dt / h);
        const T gdt  = damping * dt;
        const T denom = T{1} + gdt;

        std::vector<T> next(node_count, T{0});
        for (std::size_t i = 1; i + 1 < node_count; ++i) {
            T laplacian = current[i + 1] - T{2} * current[i] + current[i - 1];
            next[i] = (T{2} * current[i] - (T{1} - gdt) * previous[i] + r2 * laplacian) / denom;
        }

        previous = std::move(current);
        current  = std::move(next);
    }
};

// Builds a plucked string: triangular initial displacement peaking at
// x = pluck_position * length with height pluck_height, zero initial
// velocity -- the standard idealized pluck shape (a string held at one
// point and released, as opposed to a struck string's initial-velocity
// impulse -- contrast wave_membrane.hpp's make_struck_membrane).
//
// Zero initial velocity still needs a second-order-accurate "phantom"
// previous state to start the leapfrog recurrence correctly: Taylor-
// expand u(-dt) = u(0) - dt*u_t(0) + (dt^2/2)*u_tt(0); u_t(0) = 0 (zero
// initial velocity) and u_tt(0) = c^2*u_xx(0) from the PDE itself, giving
// u(-dt) = u(0) + (dt^2/2)*c^2*u_xx(0). Plugging that into VibratingString
// ::step's own (undamped) recurrence at n=0 (2*u^0 - u^{-1} +
// r2*laplacian(u^0)) reproduces exactly the well-known standalone
// first-step formula u^1 = u^0 + (r2/2)*laplacian(u^0) -- so setting
// `previous` this way, rather than just previous = current, makes the
// very first step() call already second-order accurate instead of only
// first-order. This still holds exactly with damping != 0 too: step()'s
// (1-gdt)/(1+gdt) numerator/denominator factors cancel algebraically
// against this same phantom state (u_t(0) = 0 means damping, which acts
// on velocity, has nothing to act on yet at the very first step).
template<Scalar T = double>
VibratingString<T> make_plucked_string(std::size_t node_count, T length,
                                        T wave_speed, T dt,
                                        T pluck_position = T{0.3},
                                        T pluck_height = T{1},
                                        T damping = T{0}) {
    VibratingString<T> s;
    s.node_count  = node_count;
    s.length      = length;
    s.wave_speed  = wave_speed;
    s.damping     = damping;
    s.current.assign(node_count, T{0});
    s.previous.assign(node_count, T{0});

    const T h       = s.dx();
    const T pluck_x = pluck_position * length;
    for (std::size_t i = 0; i < node_count; ++i) {
        T x = T(i) * h;
        if (x <= pluck_x)
            s.current[i] = (pluck_x > T{0}) ? pluck_height * (x / pluck_x) : T{0};
        else
            s.current[i] = pluck_height * (length - x) / (length - pluck_x);
    }
    s.current[0]              = T{0};
    s.current[node_count - 1] = T{0};

    const T r2 = (wave_speed * dt / h) * (wave_speed * dt / h);
    s.previous = s.current;
    for (std::size_t i = 1; i + 1 < node_count; ++i) {
        T laplacian   = s.current[i + 1] - T{2} * s.current[i] + s.current[i - 1];
        s.previous[i] = s.current[i] + T{0.5} * r2 * laplacian;
    }
    return s;
}

} // namespace spatium::physics::mechanics
