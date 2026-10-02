#pragma once

// Metric-agnostic geodesic integration in four dimensions. The metric is any
// generic callable matching SchwarzschildMetric's shape (Vec<S,4> ->
// Matrix<S,4,4>, templated on S, never std::function -- see schwarzschild.hpp's
// header comment for why that distinction matters), so everything below is
// oblivious to which metric it's fed: Schwarzschild, Kerr (whose g_t_phi term
// makes the metric non-diagonal), a MetricField, a superposition of holes.
//
// The Christoffel symbols and the geodesic equation are not four-dimensional
// and do not live here: they are spaces/metric_chart.hpp's, for any number of
// coordinates, where the same code derives a space's exp and log from its
// metric alone. This header names them for the relativity code that has
// always used them (`christoffel(metric, x)` with x a Vec<T,4>, state as a
// Vec<T,8> of position and 4-velocity) and keeps what is four-dimensional
// and a spacetime's: the conserved quantities of its Killing vectors.
// Where g does not invert, geodesic_rhs now returns NaN; it used to read
// the symbols as zero and fly straight on, which is a wrong answer.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/dual.hpp>
#  include <spatium/algebra/linear_solve.hpp>
#  include <spatium/algebra/matrix.hpp>
#  include <spatium/algebra/ode.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/concepts.hpp>
#  include <spatium/spaces/metric_chart.hpp>
#  include <array>
#  include <cstddef>
#endif

SPATIUM_EXPORT namespace spatium::physics::relativity {

// g_{mu nu}(x), its exact partials, the Christoffel symbols
//   Gamma^lambda_{mu nu} = 1/2 g^{lambda sigma}
//       (d_mu g_{sigma nu} + d_nu g_{sigma mu} - d_sigma g_{mu nu}),
// the geodesic equation d^2x^lambda/dlambda^2 = -Gamma^lambda_{mu nu} u^mu u^nu
// as a first-order system on [x^mu, dx^mu/dlambda], and one RK4 step of it --
// the generic ones, in four dimensions.
template<Scalar T> using MetricDerivatives = spaces::MetricDerivatives<T, 4>;
using spaces::metric_derivatives;
using spaces::christoffel;
using spaces::geodesic_rhs;
using spaces::geodesic_step;

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
