#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/concepts.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <cstddef>
#endif

SPATIUM_EXPORT namespace spatium {
inline namespace algebra {

// General first-order IVP ("Cauchy problem"): y'(t) = f(t, y), y(t0) = y0.
// Decoupled from physics/mechanics/integrator.hpp's PointMass-tied
// Euler/semi-implicit-Euler/Verlet on purpose -- those solve Newton's
// second law for a specific state layout; this solves an arbitrary
// first-order vector ODE, the same generality gradient()/integrate()/
// minimize() in calculus.hpp already bring to differentiation/quadrature/
// optimization. A second-order ODE (e.g. the harmonic oscillator) becomes
// first-order here the standard way: stack position and velocity into one
// Vec<T, 2*M> state and let f return their derivatives together.
template<typename F, typename T, std::size_t N>
concept OdeRhs = requires(F f, T t, Vec<T, N> y) {
    { f(t, y) } -> std::convertible_to<Vec<T, N>>;
};

// Explicit Euler: first-order accurate, O(dt) local error. Cheap (one f
// evaluation per step) but its stability region is bounded -- for
// y'=-k*y it only stays bounded while |1-k*dt| <= 1, i.e. k*dt <= 2;
// past that the numerical solution blows up even though the true
// solution decays. That stability boundary is exactly the "recognize
// the delicate case" signal rsc/include/ode_task.hpp dispatches on.
template<Scalar T, std::size_t N, typename F>
    requires OdeRhs<F, T, N>
Vec<T, N> euler_step(F&& f, T t, const Vec<T, N>& y, T dt) {
    return Vec<T, N>{y + f(t, y) * dt};
}

// Classical 4th-order Runge-Kutta: O(dt^4) local error, four f
// evaluations per step, a much larger stability region than Euler's --
// the "expensive but safe" candidate in the same role Real50/bisection
// played for the earlier domains.
template<Scalar T, std::size_t N, typename F>
    requires OdeRhs<F, T, N>
Vec<T, N> rk4_step(F&& f, T t, const Vec<T, N>& y, T dt) {
    Vec<T, N> k1 = f(t, y);
    Vec<T, N> k2 = f(t + dt / T{2}, Vec<T, N>{y + k1 * (dt / T{2})});
    Vec<T, N> k3 = f(t + dt / T{2}, Vec<T, N>{y + k2 * (dt / T{2})});
    Vec<T, N> k4 = f(t + dt, Vec<T, N>{y + k3 * dt});
    return Vec<T, N>{y + (k1 + k2 * T{2} + k3 * T{2} + k4) * (dt / T{6})};
}

// Fixed-step integration from t0 to t0 + steps*dt, stepper = euler_step or
// rk4_step (or any callable with the same (f, t, y, dt) -> Vec shape).
// Adaptive step-size control (Dormand-Prince/RK45) is a deliberate later
// addition, not v1 -- the dispatch domain only needs fixed-step comparison
// to have a real signal (see ode_task.hpp).
template<Scalar T, std::size_t N, typename F, typename Stepper>
Vec<T, N> integrate_fixed(F&& f, T t0, Vec<T, N> y0, T dt, int steps, Stepper&& stepper) {
    Vec<T, N> y = y0;
    T t = t0;
    for (int i = 0; i < steps; ++i) {
        y = stepper(f, t, y, dt);
        t = t + dt;
    }
    return y;
}

} // namespace algebra
} // namespace spatium
