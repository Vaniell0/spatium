#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <cmath>
#  include <format>
#endif

SPATIUM_EXPORT namespace spatium {
inline namespace algebra {

template<Scalar T = double>
struct Complex {
    T re{};
    T im{};

    constexpr Complex() = default;
    constexpr Complex(T real) : re(real), im(T{0}) {}
    constexpr Complex(T real, T imag) : re(real), im(imag) {}

    static constexpr Complex from_polar(T r, T theta) {
        using std::cos, std::sin; // ADL: lets non-std Scalar T (e.g. Real50) provide these
        return {r * cos(theta), r * sin(theta)};
    }

    // Arithmetic
    constexpr Complex operator+(const Complex& o) const { return {re + o.re, im + o.im}; }
    constexpr Complex operator-(const Complex& o) const { return {re - o.re, im - o.im}; }
    constexpr Complex operator-() const { return {-re, -im}; }

    constexpr Complex operator*(const Complex& o) const {
        return {re * o.re - im * o.im, re * o.im + im * o.re};
    }

    constexpr Complex operator/(const Complex& o) const {
        auto d = o.re * o.re + o.im * o.im;
        return {(re * o.re + im * o.im) / d, (im * o.re - re * o.im) / d};
    }

    // Scalar arithmetic
    constexpr Complex operator*(T s) const { return {re * s, im * s}; }
    constexpr Complex operator/(T s) const { return {re / s, im / s}; }
    friend constexpr Complex operator*(T s, const Complex& c) { return {s * c.re, s * c.im}; }

    // Compound assignment
    constexpr Complex& operator+=(const Complex& o) { re += o.re; im += o.im; return *this; }
    constexpr Complex& operator-=(const Complex& o) { re -= o.re; im -= o.im; return *this; }
    constexpr Complex& operator*=(const Complex& o) { *this = *this * o; return *this; }
    constexpr Complex& operator/=(const Complex& o) { *this = *this / o; return *this; }

    // Properties
    constexpr Complex conjugate() const { return {re, -im}; }
    constexpr T magnitude_sq() const { return re * re + im * im; }
    T magnitude() const { using std::sqrt; return sqrt(magnitude_sq()); }
    T phase() const { using std::atan2; return atan2(im, re); }
    bool is_real(T eps = epsilon<T>()) const { using std::abs; return abs(im) <= eps; }

    constexpr bool operator==(const Complex&) const = default;
};

// Square root of complex number
template<Scalar T>
Complex<T> sqrt(const Complex<T>& z) {
    using std::sqrt; // ADL: lets non-std Scalar T (e.g. Real50) provide its own sqrt
    auto r = z.magnitude();
    auto theta = z.phase();
    return Complex<T>::from_polar(sqrt(r), theta / T{2});
}

// Cube root of complex number (principal root)
template<Scalar T>
Complex<T> cbrt(const Complex<T>& z) {
    using std::cbrt; // ADL: lets non-std Scalar T (e.g. Real50) provide its own cbrt
    auto r = z.magnitude();
    auto theta = z.phase();
    return Complex<T>::from_polar(cbrt(r), theta / T{3});
}

using Complex64 = Complex<double>;
using Complex32 = Complex<float>;

} // namespace algebra
} // namespace spatium

// std::format support
template<spatium::Scalar T>
struct std::formatter<spatium::Complex<T>> : std::formatter<T> {
    auto format(const spatium::Complex<T>& c, auto& ctx) const {
        if (c.im >= T{0})
            return std::format_to(ctx.out(), "({}+{}i)", c.re, c.im);
        else
            return std::format_to(ctx.out(), "({}{}i)", c.re, c.im);
    }
};
