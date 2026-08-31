#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <algorithm>
#  include <cmath>
#  include <limits>
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
