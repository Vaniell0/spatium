#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/concepts.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/core/error.hpp>
#  include <array>
#  include <cmath>
#  include <cstddef>
#  include <optional>
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

// Gragg-Bulirsch-Stoer: the step that reaches the precision of the scalar.
//
// RK4's error is h^4, so a fixed step cannot meet a tolerance that follows
// the arithmetic -- Real50 asks for 1e-33 where double asks for 1e-11, and
// the same step count would have to differ by a factor 10^5. The modified
// midpoint rule of Gragg has an error that is a series in h^2 only, so
// extrapolating a handful of rule evaluations (2, 4, 6, ... substeps) to
// h = 0 gains two orders per column, and a smooth problem converges in a
// few columns at any tolerance. Where it does not -- the step is too long
// for the problem -- the step is halved.
//
// Hence the tolerance is an argument, not a fixed step: a caller passes a
// power of its scalar's epsilon (`machine_epsilon<T>()`), and the same call
// is correct on float, double, Real50 and Dual, each to its own roundoff.
// The convergence test reads the error as a double, `primal_double`, since
// a Dual's derivative carries no magnitude; the derivative parts are
// extrapolated with the values and agree with a finite difference to the
// same order, which the connectivity matrix's derivative cells hold.
//
// Returns NotConverged where it cannot reach `tol` -- a step that halves to
// nothing, or a right-hand side that is not finite, as at a coordinate
// singularity -- never a silently wrong state.

namespace ode_detail {

// Gragg's modified midpoint over [t, t + H] in n substeps, with the
// smoothing step that keeps the error a series in h^2.
template<Scalar T, std::size_t N, typename F>
Vec<T, N> gragg(F& f, T t, const Vec<T, N>& y, T H, int n) {
    const T h = H / T(n);
    Vec<T, N> prev = y;
    Vec<T, N> cur{y + f(t, y) * h};
    for (int m = 1; m < n; ++m) {
        Vec<T, N> next{prev + f(T(t + T(m) * h), cur) * (T{2} * h)};
        prev = cur;
        cur = next;
    }
    return Vec<T, N>{(cur + prev + f(T(t + H), cur) * h) * T{0.5}};
}

// max_i |a_i - b_i| / (1 + |b_i|), as a double.
template<Scalar T, std::size_t N>
double scaled_difference(const Vec<T, N>& a, const Vec<T, N>& b) {
    double worst = 0;
    for (std::size_t i = 0; i < N; ++i) {
        // The difference in the scalar's own arithmetic, then a double: two
        // Real50 values 1e-20 apart are the same double.
        const double bi = primal_double(b[i]);
        const double d = std::abs(primal_double(T(a[i] - b[i]))) / (1.0 + std::abs(bi));
        if (!(d <= worst)) worst = d;      // a NaN becomes the worst, not a skipped term
    }
    return worst;
}

}  // namespace ode_detail

template<Scalar T, std::size_t N, typename F, int MaxColumns = 14>
    requires OdeRhs<F, T, N>
Result<Vec<T, N>> integrate_extrapolated(F&& f, T t0, const Vec<T, N>& y0, T t1, double tol) {
    static_assert(MaxColumns >= 2 && MaxColumns <= 24);
    const double span = std::abs(primal_double(t1) - primal_double(t0));
    if (!std::isfinite(span)) return std::unexpected(Error{ErrorCode::NotConverged, "the interval is not finite"});
    // An empty interval is the identity -- unless the scalar is a Dual, whose
    // end may be 0 in value and still move: d/dt of the flow at t = 0 is the
    // velocity, and returning y0 here dropped it. Such a step runs once, and
    // its tableau is exact in the derivative.
    constexpr bool carries_derivative = requires { t0.deriv; };
    if (span == 0 && !carries_derivative) return y0;

    Vec<T, N> y = y0;
    T t = t0;
    T H = T(t1 - t0);
    for (long steps = 0; steps < 100000; ++steps) {
        const double left = std::abs(primal_double(t1) - primal_double(t));
        const bool last = std::abs(primal_double(H)) >= left * (1.0 - 1e-12);
        const T Hs = last ? T(t1 - t) : H;

        // One attempt at the step: the extrapolation tableau a row at a time
        // -- row j is the rule at 2(j+1) substeps, entry k of it extrapolated
        // k times -- stopped when the last two diagonal entries agree.
        std::array<Vec<T, N>, MaxColumns> above{}, row{};
        std::optional<Vec<T, N>> accepted;
        int used = 0;
        for (int j = 0; j < MaxColumns; ++j) {
            row[0] = ode_detail::gragg(f, t, y, Hs, 2 * (j + 1));
            for (int k = 1; k <= j; ++k) {
                // 1 / ((n_j / n_{j-k})^2 - 1) with n_i = 2(i+1), from the
                // integers: a coefficient rounded to double would cap the
                // extrapolation at double's precision whatever the scalar.
                const long a = j + 1, b = j + 1 - k;
                row[k] = Vec<T, N>{row[k - 1] + (row[k - 1] - above[k - 1]) * (T(b * b) / T(a * a - b * b))};
            }
            used = j + 1;
            if (j >= 1 && ode_detail::scaled_difference(row[j], row[j - 1]) <= tol) {
                accepted = row[j];
                break;
            }
            above = row;
        }

        if (accepted) {
            y = *accepted;
            if (last) return y;
            t = T(t + Hs);
            // Converged in few columns: the step was timid.
            if (used <= 4) H = T(H * T{2});
        } else {
            H = T(Hs * T{0.5});
            if (span == 0 || std::abs(primal_double(H)) < span * 1e-9)
                return std::unexpected(Error{ErrorCode::NotConverged,
                                             "step shrank below 1e-9 of the interval"});
        }
    }
    return std::unexpected(Error{ErrorCode::NotConverged, "more than 1e5 steps"});
}

} // namespace algebra
} // namespace spatium
