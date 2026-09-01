#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <cmath>
#  include <format>
#endif

SPATIUM_EXPORT namespace spatium {
inline namespace algebra {

// Forward-mode automatic differentiation via dual numbers: value + deriv*eps,
// eps^2 = 0. Dual<T> satisfies Scalar itself, so substituting Dual<T> for T
// in any existing Scalar-templated function (Vec, Matrix, polynomial solvers,
// integrators, ...) computes that function's derivative alongside its value,
// with no change to the function's own code — as long as the function calls
// math via ADL (`using std::sqrt; sqrt(x)`, per this codebase's convention)
// rather than qualified `std::sqrt(x)`, so overload resolution can find the
// spatium:: overloads below.
template<Scalar T = double>
struct Dual {
    T value{};
    T deriv{};

    constexpr Dual() = default;
    constexpr Dual(T v) : value(v), deriv(T{0}) {}
    constexpr Dual(T v, T d) : value(v), deriv(d) {}

    // An independent variable: seed the derivative to 1 so downstream
    // operations accumulate d(result)/d(this variable) via the chain rule.
    static constexpr Dual variable(T v) { return {v, T{1}}; }
    static constexpr Dual constant(T v) { return {v, T{0}}; }

    constexpr Dual operator+(const Dual& o) const { return {value + o.value, deriv + o.deriv}; }
    constexpr Dual operator-(const Dual& o) const { return {value - o.value, deriv - o.deriv}; }
    constexpr Dual operator-() const { return {-value, -deriv}; }

    constexpr Dual operator*(const Dual& o) const {
        return {value * o.value, deriv * o.value + value * o.deriv};
    }

    constexpr Dual operator/(const Dual& o) const {
        return {value / o.value, (deriv * o.value - value * o.deriv) / (o.value * o.value)};
    }

    // Scalar on the left: T op Dual<T>. The member operators above already
    // cover Dual<T> op T via the converting constructor on the right-hand
    // argument; these three cover the other side (e.g. `10.0 * x` in a loss
    // function), matching Complex<T>'s own `friend operator*(T, Complex)`.
    friend constexpr Dual operator+(T s, const Dual& d) { return Dual(s) + d; }
    friend constexpr Dual operator-(T s, const Dual& d) { return Dual(s) - d; }
    friend constexpr Dual operator*(T s, const Dual& d) { return {s * d.value, s * d.deriv}; }

    constexpr Dual& operator+=(const Dual& o) { *this = *this + o; return *this; }
    constexpr Dual& operator-=(const Dual& o) { *this = *this - o; return *this; }
    constexpr Dual& operator*=(const Dual& o) { *this = *this * o; return *this; }
    constexpr Dual& operator/=(const Dual& o) { *this = *this / o; return *this; }

    // Ordering compares the primal value only — the derivative carries no
    // magnitude information, same convention every dual-number AD library uses.
    constexpr bool operator==(const Dual& o) const { return value == o.value; }
    constexpr auto operator<=>(const Dual& o) const { return value <=> o.value; }
};

// ── Chain rule for common transcendentals ─────────────────────
// Not constexpr: same reason Complex<T>'s sqrt/cbrt aren't — std::cmath
// isn't constexpr.

template<Scalar T>
Dual<T> sqrt(const Dual<T>& x) {
    using std::sqrt;
    auto r = sqrt(x.value);
    return {r, x.deriv / (T{2} * r)};
}

template<Scalar T>
Dual<T> sin(const Dual<T>& x) {
    using std::sin, std::cos;
    return {sin(x.value), x.deriv * cos(x.value)};
}

template<Scalar T>
Dual<T> cos(const Dual<T>& x) {
    using std::sin, std::cos;
    return {cos(x.value), -x.deriv * sin(x.value)};
}

template<Scalar T>
Dual<T> exp(const Dual<T>& x) {
    using std::exp;
    auto e = exp(x.value);
    return {e, x.deriv * e};
}

template<Scalar T>
Dual<T> abs(const Dual<T>& x) {
    using std::abs;
    return {abs(x.value), x.value >= T{0} ? x.deriv : -x.deriv};
}

template<Scalar T>
Dual<T> acos(const Dual<T>& x) {
    using std::acos, std::sqrt;
    return {acos(x.value), -x.deriv / sqrt(T{1} - x.value * x.value)};
}

template<Scalar T>
Dual<T> tan(const Dual<T>& x) {
    using std::tan, std::cos;
    auto c = cos(x.value);
    return {tan(x.value), x.deriv / (c * c)};
}

// Real-exponent power; n is a plain constant, not itself a Dual (matches the
// common case — differentiating x^n w.r.t. x, not w.r.t. the exponent).
template<Scalar T>
Dual<T> pow(const Dual<T>& x, T n) {
    using std::pow;
    return {pow(x.value, n), x.deriv * n * pow(x.value, n - T{1})};
}

using Dual64 = Dual<double>;
using Dual32 = Dual<float>;

// derivative(f, x0) == d/dx f(x) at x0, for any f callable as Dual<T> -> Dual<T>.
template<Scalar T, typename F>
T derivative(F&& f, T x0) {
    return f(Dual<T>::variable(x0)).deriv;
}

} // namespace algebra
} // namespace spatium

// std::format support
template<spatium::Scalar T>
struct std::formatter<spatium::Dual<T>> : std::formatter<T> {
    auto format(const spatium::Dual<T>& d, auto& ctx) const {
        return std::format_to(ctx.out(), "({}+{}ε)", d.value, d.deriv);
    }
};
