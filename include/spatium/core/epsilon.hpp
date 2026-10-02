#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <algorithm>
#  include <cmath>
#  include <limits>
#  include <type_traits>
#  include <utility>
#endif

SPATIUM_EXPORT namespace spatium {

// Type-aware epsilon for numerical comparisons.
// Uses std::numeric_limits when available, falls back to sensible defaults.
// Scales with type precision: float ~1e-5, double ~1e-10, multiprecision ~1e-40.

template<Scalar T>
constexpr T epsilon() {
    if constexpr (std::numeric_limits<T>::is_specialized) {
        // ~100x machine epsilon: enough for accumulated error in typical algorithms
        return std::numeric_limits<T>::epsilon() * T{128};
    } else {
        // Multiprecision or custom types: conservative default
        return T{1} / T{10000000000}; // 1e-10
    }
}

// The unit roundoff of the arithmetic underneath a scalar, as a double: a
// Dual carries its value's, so Dual<Dual<double>> is double's 2.2e-16 and
// Dual<Real50> is 1e-50. What a tolerance that must follow the scalar --
// an ODE solved to the precision of the numbers it is solved in -- reads;
// `epsilon<T>()` above is a comparison threshold with a fixed fallback for
// types std::numeric_limits does not know, which a Dual is.
template<class T>
constexpr double machine_epsilon() {
    if constexpr (requires(const T& x) { x.value; x.deriv; })
        return machine_epsilon<std::remove_cvref_t<decltype(std::declval<const T&>().value)>>();
    else
        return static_cast<double>(std::numeric_limits<T>::epsilon());
}

// A scalar's value as a double, through any Dual layers: what a decision --
// has it converged, how big is this -- reads, where a Dual's derivative
// carries no magnitude (its ordering compares the value alone).
template<class T>
constexpr double primal_double(const T& x) {
    if constexpr (requires { x.value; x.deriv; }) return primal_double(x.value);
    else return static_cast<double>(x);
}

// Relative epsilon: max(eps, eps * scale)
// Use for comparisons where magnitude matters.
template<Scalar T>
T relative_epsilon(T scale) {
    using std::abs;
    auto eps = epsilon<T>();
    auto s = abs(scale);
    return (s > T{1}) ? eps * s : eps;
}

// Near-zero check
template<Scalar T>
bool near_zero(T value) {
    using std::abs;
    return abs(value) < epsilon<T>();
}

// Approximate equality
template<Scalar T>
bool approx_equal(T a, T b) {
    using std::abs;
    return abs(a - b) < relative_epsilon<T>(std::max(abs(a), abs(b)));
}

} // namespace spatium
