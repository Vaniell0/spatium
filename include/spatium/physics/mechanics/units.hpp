#pragma once

// Compile-time SI dimensional analysis.
//
// Quantity<M, L, T, I, K, N, J, Value> wraps a numeric value and tracks the
// SI dimension exponents at compile time:
//   M  = mass               (kg)
//   L  = length             (m)
//   T  = time               (s)
//   I  = electric current   (A)
//   K  = temperature        (K)
//   N  = amount of substance (mol)
//   J  = luminous intensity (cd)
//
// Multiplication and division add/subtract exponents; addition requires
// matching dimensions and fails to compile otherwise. SI literals expose
// concrete units: 1.0_kg, 9.81_m_per_s2, 5.0_N, etc.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <ratio>
#  include <type_traits>
#endif

SPATIUM_EXPORT namespace spatium::physics::mechanics {

template<typename M, typename L, typename T, typename I, typename K,
         typename N, typename J, typename Value = double>
struct Quantity {
    using mass_dim       = M;
    using length_dim     = L;
    using time_dim       = T;
    using current_dim    = I;
    using temperature_dim= K;
    using amount_dim     = N;
    using luminosity_dim = J;
    using value_type     = Value;

    Value v{};

    constexpr Quantity() = default;
    constexpr explicit Quantity(Value x) : v(x) {}

    [[nodiscard]] constexpr Value value() const noexcept { return v; }

    constexpr Quantity operator+(Quantity rhs) const { return Quantity{v + rhs.v}; }
    constexpr Quantity operator-(Quantity rhs) const { return Quantity{v - rhs.v}; }
    constexpr Quantity operator-() const { return Quantity{-v}; }
    constexpr Quantity& operator+=(Quantity rhs) { v += rhs.v; return *this; }
    constexpr Quantity& operator-=(Quantity rhs) { v -= rhs.v; return *this; }

    constexpr bool operator==(const Quantity&) const = default;
    constexpr auto operator<=>(const Quantity&) const = default;
};

// ── Dimension arithmetic helpers ──────────────────────────────

template<typename A, typename B> using ratio_add = std::ratio_add<A, B>;
template<typename A, typename B> using ratio_sub = std::ratio_subtract<A, B>;

// Multiplication: add exponents.
template<class M1,class L1,class T1,class I1,class K1,class N1,class J1,class V,
                        class M2,class L2,class T2,class I2,class K2,class N2,class J2>
constexpr auto operator*(Quantity<M1,L1,T1,I1,K1,N1,J1,V> a,
                         Quantity<M2,L2,T2,I2,K2,N2,J2,V> b) {
    using R = Quantity<ratio_add<M1,M2>, ratio_add<L1,L2>, ratio_add<T1,T2>,
                       ratio_add<I1,I2>, ratio_add<K1,K2>, ratio_add<N1,N2>,
                       ratio_add<J1,J2>, V>;
    return R{a.value() * b.value()};
}

// Division: subtract exponents.
template<class M1,class L1,class T1,class I1,class K1,class N1,class J1,class V,
                        class M2,class L2,class T2,class I2,class K2,class N2,class J2>
constexpr auto operator/(Quantity<M1,L1,T1,I1,K1,N1,J1,V> a,
                         Quantity<M2,L2,T2,I2,K2,N2,J2,V> b) {
    using R = Quantity<ratio_sub<M1,M2>, ratio_sub<L1,L2>, ratio_sub<T1,T2>,
                       ratio_sub<I1,I2>, ratio_sub<K1,K2>, ratio_sub<N1,N2>,
                       ratio_sub<J1,J2>, V>;
    return R{a.value() / b.value()};
}

// Scalar (dimensionless) multiplication on either side.
template<class M,class L,class T,class I,class K,class N,class J,class V>
constexpr Quantity<M,L,T,I,K,N,J,V> operator*(V scalar, Quantity<M,L,T,I,K,N,J,V> q) {
    return Quantity<M,L,T,I,K,N,J,V>{scalar * q.value()};
}

template<class M,class L,class T,class I,class K,class N,class J,class V>
constexpr Quantity<M,L,T,I,K,N,J,V> operator*(Quantity<M,L,T,I,K,N,J,V> q, V scalar) {
    return Quantity<M,L,T,I,K,N,J,V>{q.value() * scalar};
}

template<class M,class L,class T,class I,class K,class N,class J,class V>
constexpr Quantity<M,L,T,I,K,N,J,V> operator/(Quantity<M,L,T,I,K,N,J,V> q, V scalar) {
    return Quantity<M,L,T,I,K,N,J,V>{q.value() / scalar};
}

// ── Base SI dimensions (single-exponent shortcuts) ────────────

using Z = std::ratio<0>;
using P1 = std::ratio<1>;
using P2 = std::ratio<2>;
using P3 = std::ratio<3>;
using N1 = std::ratio<-1>;
using N2 = std::ratio<-2>;
using N3 = std::ratio<-3>;

using Dimensionless = Quantity<Z, Z, Z, Z, Z, Z, Z, double>;
using Mass          = Quantity<P1, Z,  Z,  Z, Z, Z, Z, double>;  // kg
using Length        = Quantity<Z,  P1, Z,  Z, Z, Z, Z, double>;  // m
using Time          = Quantity<Z,  Z,  P1, Z, Z, Z, Z, double>;  // s
using Current       = Quantity<Z,  Z,  Z,  P1, Z, Z, Z, double>; // A
using Temperature   = Quantity<Z,  Z,  Z,  Z, P1, Z, Z, double>; // K
using Amount        = Quantity<Z,  Z,  Z,  Z, Z, P1, Z, double>; // mol
using Luminosity    = Quantity<Z,  Z,  Z,  Z, Z, Z, P1, double>; // cd

// ── Derived SI quantities ─────────────────────────────────────

using Velocity      = Quantity<Z,  P1, N1, Z, Z, Z, Z, double>;  // m/s
using Acceleration  = Quantity<Z,  P1, N2, Z, Z, Z, Z, double>;  // m/s²
using Force         = Quantity<P1, P1, N2, Z, Z, Z, Z, double>;  // kg·m/s² = N
using Energy        = Quantity<P1, P2, N2, Z, Z, Z, Z, double>;  // kg·m²/s² = J
using Power         = Quantity<P1, P2, N3, Z, Z, Z, Z, double>;  // J/s = W
using Momentum      = Quantity<P1, P1, N1, Z, Z, Z, Z, double>;  // kg·m/s
using AngularVelocity = Quantity<Z, Z, N1, Z, Z, Z, Z, double>;  // rad/s (rad is dimensionless)
using AngularMomentum = Quantity<P1, P2, N1, Z, Z, Z, Z, double>; // kg·m²/s
using Frequency     = Quantity<Z,  Z,  N1, Z, Z, Z, Z, double>;  // 1/s = Hz
using Pressure      = Quantity<P1, N1, N2, Z, Z, Z, Z, double>;  // Pa = N/m²
using Density       = Quantity<P1, N3, Z,  Z, Z, Z, Z, double>;  // kg/m³

// ── Physical constants (SI) ───────────────────────────────────

inline constexpr Acceleration g_earth{9.80665};                  // standard gravity
inline constexpr auto G_newton = Quantity<N1, P3, N2, Z, Z, Z, Z, double>{6.67430e-11};
                                                                 // m³/(kg·s²)
inline constexpr auto c_light  = Velocity{2.99792458e8};         // exact, m/s
inline constexpr auto h_planck = Quantity<P1, P2, N1, Z, Z, Z, Z, double>{6.62607015e-34};
                                                                 // kg·m²/s = J·s
inline constexpr auto k_boltz  = Quantity<P1, P2, N2, Z, N1, Z, Z, double>{1.380649e-23};
                                                                 // J/K

} // namespace spatium::physics::mechanics

// ── User-defined literals (in inline namespace for ADL into top-level scope) ──
SPATIUM_EXPORT namespace spatium::physics::mechanics::literals {

constexpr Mass         operator""_kg(long double v)        { return Mass{static_cast<double>(v)}; }
constexpr Mass         operator""_kg(unsigned long long v) { return Mass{static_cast<double>(v)}; }
constexpr Length       operator""_m (long double v)        { return Length{static_cast<double>(v)}; }
constexpr Length       operator""_m (unsigned long long v) { return Length{static_cast<double>(v)}; }
constexpr Time         operator""_s (long double v)        { return Time{static_cast<double>(v)}; }
constexpr Time         operator""_s (unsigned long long v) { return Time{static_cast<double>(v)}; }
constexpr Force        operator""_N (long double v)        { return Force{static_cast<double>(v)}; }
constexpr Force        operator""_N (unsigned long long v) { return Force{static_cast<double>(v)}; }
constexpr Energy       operator""_J (long double v)        { return Energy{static_cast<double>(v)}; }
constexpr Energy       operator""_J (unsigned long long v) { return Energy{static_cast<double>(v)}; }
constexpr Power        operator""_W (long double v)        { return Power{static_cast<double>(v)}; }
constexpr Power        operator""_W (unsigned long long v) { return Power{static_cast<double>(v)}; }

} // namespace spatium::physics::mechanics::literals
