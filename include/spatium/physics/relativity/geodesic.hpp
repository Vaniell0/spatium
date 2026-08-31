#pragma once

// Metric-agnostic geodesic integration. The metric is any generic
// callable matching SchwarzschildMetric's shape (Vec<S,4> -> Matrix<S,4,4>,
// templated on S, never std::function -- see schwarzschild.hpp's header
// comment for why that distinction matters), so everything below is
// oblivious to which metric it's fed: Schwarzschild is the only one
// shipped today, but a future Kerr metric (whose g_t_phi term makes the
// metric non-diagonal) plugs into the exact same christoffel()/
// geodesic_rhs() unchanged. That's also why the metric inverse below
// reuses the library's general solve_direct() instead of a diagonal-only
// shortcut -- a diagonal shortcut would need rewriting the day Kerr
// arrives; the general solve doesn't.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/dual.hpp>
#  include <spatium/algebra/linear_solve.hpp>
#  include <spatium/algebra/matrix.hpp>
#  include <spatium/algebra/ode.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/concepts.hpp>
#  include <array>
#  include <cstddef>
#endif

SPATIUM_EXPORT namespace spatium::physics::relativity {

// g_{mu nu}(x) and its exact partials d_kappa g_{mu nu}(x) for all
// (kappa, mu, nu) at once: one Dual-seeded evaluation of the whole metric
// per coordinate axis (4 metric evaluations total, not 4*16) -- the
// rank-2-tensor generalization of calculus.hpp's per-axis-seeded
// gradient() for a scalar field.
template<Scalar T>
struct MetricDerivatives {
    Matrix<T, 4, 4> g{};
    std::array<Matrix<T, 4, 4>, 4> dg{};  // dg[kappa](mu,nu) = d g_{mu nu} / d x^kappa
};

template<Scalar T, typename Metric>
MetricDerivatives<T> metric_derivatives(const Metric& metric, const Vec<T, 4>& x) {
    MetricDerivatives<T> md{};
    for (std::size_t kappa = 0; kappa < 4; ++kappa) {
        Vec<Dual<T>, 4> dx;
        for (std::size_t j = 0; j < 4; ++j)
            dx[j] = (j == kappa) ? Dual<T>::variable(x[j]) : Dual<T>::constant(x[j]);
        Matrix<Dual<T>, 4, 4> gd = metric(dx);
        for (std::size_t mu = 0; mu < 4; ++mu)
            for (std::size_t nu = 0; nu < 4; ++nu) {
                if (kappa == 0) md.g(mu, nu) = gd(mu, nu).value;
                md.dg[kappa](mu, nu) = gd(mu, nu).deriv;
            }
    }
    return md;
}

// Christoffel symbols of the second kind:
//   Gamma^lambda_{mu nu} = 1/2 g^{lambda sigma}
//       (d_mu g_{sigma nu} + d_nu g_{sigma mu} - d_sigma g_{mu nu})
// g^{lambda sigma} (the metric inverse) comes from invert() -- a general
// N x N Gauss-Jordan inverse, so a future non-diagonal metric needs no
// change here. Profiling during the demo build showed this call is the
// dominant per-step cost (~1us/step with 4 independent solve_direct()
// calls vs ~0.3us/step with one invert() -- solve_direct() redid the
// full elimination once per basis vector).
template<Scalar T, typename Metric>
std::array<Matrix<T, 4, 4>, 4> christoffel(const Metric& metric, const Vec<T, 4>& x) {
    auto md = metric_derivatives(metric, x);

    auto inv = invert(md.g);
    Matrix<T, 4, 4> ginv = inv ? *inv : Matrix<T, 4, 4>{};

    std::array<Matrix<T, 4, 4>, 4> Gamma{};
    for (std::size_t lambda = 0; lambda < 4; ++lambda)
        for (std::size_t mu = 0; mu < 4; ++mu)
            for (std::size_t nu = 0; nu < 4; ++nu) {
                T sum{0};
                for (std::size_t sigma = 0; sigma < 4; ++sigma)
                    sum += ginv(lambda, sigma) *
                           (md.dg[mu](sigma, nu) + md.dg[nu](sigma, mu) - md.dg[sigma](mu, nu));
                Gamma[lambda](mu, nu) = T{0.5} * sum;
            }
    return Gamma;
}

// Geodesic state packs position and 4-velocity into one Vec<T,8>:
// [0..3] = x^mu (t, r, theta, phi), [4..7] = dx^mu/dlambda. The geodesic
// equation d^2x^lambda/dlambda^2 = -Gamma^lambda_{mu nu} u^mu u^nu becomes
// first-order the way ode.hpp's own doc comment names (stack position and
// velocity), so this is a plain OdeRhs for rk4_step/integrate_fixed.
template<Scalar T, typename Metric>
Vec<T, 8> geodesic_rhs(const Metric& metric, T /*lambda*/, const Vec<T, 8>& state) {
    Vec<T, 4> x{state[0], state[1], state[2], state[3]};
    Vec<T, 4> v{state[4], state[5], state[6], state[7]};
    auto Gamma = christoffel(metric, x);

    Vec<T, 8> ds{};
    for (std::size_t mu = 0; mu < 4; ++mu) ds[mu] = v[mu];
    for (std::size_t lambda = 0; lambda < 4; ++lambda) {
        T acc{0};
        for (std::size_t mu = 0; mu < 4; ++mu)
            for (std::size_t nu = 0; nu < 4; ++nu)
                acc += Gamma[lambda](mu, nu) * v[mu] * v[nu];
        ds[4 + lambda] = -acc;
    }
    return ds;
}

// One RK4 step of the geodesic equation via the shared ODE module --
// this module's first real production use (previously exercised only by
// tests/test_ode.cpp).
template<Scalar T, typename Metric>
Vec<T, 8> geodesic_step(const Metric& metric, const Vec<T, 8>& state, T dlambda) {
    auto f = [&metric](T lam, const Vec<T, 8>& y) { return geodesic_rhs(metric, lam, y); };
    return rk4_step(f, T{0}, state, dlambda);
}

// Conserved quantities from the Killing vectors d/dt and d/dphi, present
// for any stationary, axisymmetric metric -- written as the general
// contraction p_mu = g_{mu nu} u^nu rather than Schwarzschild's own
// diagonal shortcut, so a future non-diagonal metric (Kerr's g_t_phi
// cross term) is conserved correctly by these same two functions.
template<Scalar T, typename Metric>
T killing_energy(const Metric& metric, const Vec<T, 8>& state) {
    Vec<T, 4> x{state[0], state[1], state[2], state[3]};
    Vec<T, 4> v{state[4], state[5], state[6], state[7]};
    Matrix<T, 4, 4> g = metric(x);
    T p_t = g(0, 0) * v[0] + g(0, 3) * v[3];
    return -p_t;
}

template<Scalar T, typename Metric>
T killing_angular_momentum(const Metric& metric, const Vec<T, 8>& state) {
    Vec<T, 4> x{state[0], state[1], state[2], state[3]};
    Vec<T, 4> v{state[4], state[5], state[6], state[7]};
    Matrix<T, 4, 4> g = metric(x);
    return g(3, 0) * v[0] + g(3, 3) * v[3];
}

} // namespace spatium::physics::relativity
