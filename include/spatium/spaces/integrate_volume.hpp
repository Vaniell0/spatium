#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/quadrature.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/spaces/metric_chart.hpp>
#  include <algorithm>
#  include <cmath>
#  include <cstddef>
#  include <tuple>
#  include <type_traits>
#endif

// The integral of a function over a space given by a metric in coordinates:
//
//   integral of f(x) sqrt|det g(x)| dx_1 ... dx_K
//
// over a product of one-dimensional domains, each of which may be finite, a half
// line or the whole line (algebra/quadrature.hpp's `Finite`, `HalfLine`, `WholeLine`).
// With f = 1 it is the volume -- the area of a surface, the length of a curve -- and
// the metric is whatever `spaces::metric_chart` takes: `pullback_metric` of a map
// (a torus, the sphere through the stereographic map, a hyperboloid through polar
// coordinates), or a metric written down.
//
// It is the door of `integrate` applied once per coordinate, innermost last, so the
// infinite domains of the stereographic chart integrate the whole sphere and nothing
// is sampled at the pole it leaves out. The answer is an `IntegralResult`: its
// status is the worst of any level's, and its error estimate is the outer rule's plus
// the largest estimate any inner integral reported -- an estimate, as everywhere
// else, not a bound, and for a product of domains of measure one it is what the
// inner errors add up to.

SPATIUM_EXPORT namespace spatium::spaces {

namespace volume_detail {

template<class D> struct domain_scalar;
template<class T> struct domain_scalar<Finite<T>> { using type = T; };
template<class T> struct domain_scalar<HalfLine<T>> { using type = T; };
template<class T> struct domain_scalar<WholeLine<T>> { using type = T; };

// Worse statuses win; Divergent and Failed dominate the rest.
constexpr int severity(IntegralStatus s) {
    switch (s) {
        case IntegralStatus::Converged:   return 0;
        case IntegralStatus::Suspicious:  return 1;
        case IntegralStatus::DepthCapped: return 2;
        case IntegralStatus::Divergent:   return 3;
        case IntegralStatus::Failed:      return 4;
    }
    return 4;
}

// What the inner integrals report, kept so that it can be weighed at the end. A status
// is not enough to carry up unweighed: at the farthest nodes of an infinite outer rule
// (x of order 1e39 for the sphere's chart) the inner integral is 1e-117 and its rule
// honestly says it ran out of levels -- immaterial to a total of 12.57, and calling the
// whole answer "depth-capped" for it would be a false alarm. So a failure or a
// divergence (no number, or none that exists) always counts, and a rule that merely
// did not converge counts only if its own error is larger than the outer tolerance
// would let through.
template<class T>
struct Tally {
    IntegralStatus hard = IntegralStatus::Converged;   // Failed or Divergent anywhere
    IntegralStatus soft = IntegralStatus::Converged;   // the worst of Suspicious / DepthCapped
    T soft_error{0};                                   // the largest estimate among those
    T worst_inner_error{0};                            // the largest estimate of any inner integral
    long evaluations = 0;

    void note_inner(const IntegralResult<T>& r) {
        evaluations += r.evaluations;
        if (r.error_estimate > worst_inner_error) worst_inner_error = r.error_estimate;
        if (severity(r.status) >= 3) { if (severity(r.status) > severity(hard)) hard = r.status; }
        else if (severity(r.status) > 0) {
            if (severity(r.status) > severity(soft)) soft = r.status;
            if (r.error_estimate > soft_error) soft_error = r.error_estimate;
        }
    }
};

// One coordinate per recursion level, written as a class rather than a lambda: the
// levels are a compile-time count, and a lambda that carried it would be a type of
// this translation unit exposed in a module's interface.
template<class T, std::size_t K, class Metric, class F, class Domains>
struct Nested {
    const Metric& g;
    F& f;
    const Domains& domains;
    Vec<T, K> x{};
    Tally<T> tally{};

    template<std::size_t I>
    IntegralResult<T> level() {
        const auto inner = [this](T t) -> T {
            x[I] = t;
            if constexpr (I + 1 == K) {
                return T(f(x) * volume_element(g, x));
            } else {
                const IntegralResult<T> r = level<I + 1>();
                tally.note_inner(r);
                return r.value;
            }
        };
        return integrate(inner, std::get<I>(domains));
    }
};

}  // namespace volume_detail

template<class Metric, class F, class D0, class... Ds>
auto integrate_volume(const Metric& g, F&& f, D0 d0, Ds... ds) {
    using T = typename volume_detail::domain_scalar<D0>::type;
    constexpr std::size_t K = 1 + sizeof...(Ds);
    const std::tuple<D0, Ds...> domains{d0, ds...};
    volume_detail::Nested<T, K, Metric, std::remove_reference_t<F>, std::tuple<D0, Ds...>> nested{g, f, domains};
    IntegralResult<T> outer = nested.template level<0>();
    auto& t = nested.tally;
    t.evaluations += outer.evaluations;
    // The outer rule's own tolerance, as it reads it: its default is eps^(3/4) of the size
    // of what it integrates.
    using std::abs;
    const T bar = T(std::pow(machine_epsilon<T>(), 0.75)) * abs(outer.value);
    IntegralStatus status = outer.status;
    if (volume_detail::severity(t.hard) > volume_detail::severity(status)) status = t.hard;
    else if (t.soft_error > bar && volume_detail::severity(t.soft) > volume_detail::severity(status)) status = t.soft;
    outer.error_estimate = T(outer.error_estimate + t.worst_inner_error);
    outer.status = status;
    outer.evaluations = t.evaluations;
    return outer;
}

// The volume itself: f = 1.
template<class Metric, class D0, class... Ds>
auto volume_of(const Metric& g, D0 d0, Ds... ds) {
    using T = typename volume_detail::domain_scalar<D0>::type;
    return integrate_volume(g, [](const auto&) { return T{1}; }, d0, ds...);
}

} // namespace spatium::spaces
