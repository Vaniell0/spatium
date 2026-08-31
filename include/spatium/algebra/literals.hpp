#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <numbers>
#endif

SPATIUM_EXPORT namespace spatium::literals {

constexpr double operator""_deg(long double d) {
    return static_cast<double>(d) * std::numbers::pi / 180.0;
}

constexpr double operator""_deg(unsigned long long d) {
    return static_cast<double>(d) * std::numbers::pi / 180.0;
}

constexpr double operator""_pi(long double x) {
    return static_cast<double>(x) * std::numbers::pi;
}

constexpr double operator""_pi(unsigned long long x) {
    return static_cast<double>(x) * std::numbers::pi;
}

// ── Vec3 axis literals ────────────────────────────────────────
// Usage: 3.0_x + 2.0_y + 1.0_z → Vec3{3.0, 2.0, 1.0}

constexpr Vec3 operator""_x(long double v) { return Vec3{static_cast<double>(v), 0.0, 0.0}; }
constexpr Vec3 operator""_x(unsigned long long v) { return Vec3{static_cast<double>(v), 0.0, 0.0}; }
constexpr Vec3 operator""_y(long double v) { return Vec3{0.0, static_cast<double>(v), 0.0}; }
constexpr Vec3 operator""_y(unsigned long long v) { return Vec3{0.0, static_cast<double>(v), 0.0}; }
constexpr Vec3 operator""_z(long double v) { return Vec3{0.0, 0.0, static_cast<double>(v)}; }
constexpr Vec3 operator""_z(unsigned long long v) { return Vec3{0.0, 0.0, static_cast<double>(v)}; }

// ── Vec2 axis literals ────────────────────────────────────────
// Usage: 1.0_x2 + 2.0_y2 → Vec2{1.0, 2.0}

constexpr Vec2 operator""_x2(long double v) { return Vec2{static_cast<double>(v), 0.0}; }
constexpr Vec2 operator""_x2(unsigned long long v) { return Vec2{static_cast<double>(v), 0.0}; }
constexpr Vec2 operator""_y2(long double v) { return Vec2{0.0, static_cast<double>(v)}; }
constexpr Vec2 operator""_y2(unsigned long long v) { return Vec2{0.0, static_cast<double>(v)}; }

// ── Vec4 axis literals ────────────────────────────────────────
// Vec4 is its own coherent set so addition stays homogeneous —
// `Vec::operator+` requires same-rank operands, so mixing _x/_y/_z
// (Vec3) with _w (Vec4) does not type-check.  Use the _*4 family
// for full 4D vectors.
// Usage: 1.0_x4 + 2.0_y4 + 3.0_z4 + 4.0_w → Vec4{1.0, 2.0, 3.0, 4.0}

constexpr Vec4 operator""_x4(long double v) { return Vec4{static_cast<double>(v), 0.0, 0.0, 0.0}; }
constexpr Vec4 operator""_x4(unsigned long long v) { return Vec4{static_cast<double>(v), 0.0, 0.0, 0.0}; }
constexpr Vec4 operator""_y4(long double v) { return Vec4{0.0, static_cast<double>(v), 0.0, 0.0}; }
constexpr Vec4 operator""_y4(unsigned long long v) { return Vec4{0.0, static_cast<double>(v), 0.0, 0.0}; }
constexpr Vec4 operator""_z4(long double v) { return Vec4{0.0, 0.0, static_cast<double>(v), 0.0}; }
constexpr Vec4 operator""_z4(unsigned long long v) { return Vec4{0.0, 0.0, static_cast<double>(v), 0.0}; }
constexpr Vec4 operator""_w(long double v) { return Vec4{0.0, 0.0, 0.0, static_cast<double>(v)}; }
constexpr Vec4 operator""_w(unsigned long long v) { return Vec4{0.0, 0.0, 0.0, static_cast<double>(v)}; }

} // namespace spatium::literals
