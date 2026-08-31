#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/core/error.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <spatium/algebra/matrix.hpp>
#  include <cmath>
#  include <format>
#endif

SPATIUM_EXPORT namespace spatium {
inline namespace algebra {

template<Scalar T = double>
struct Quaternion {
    T w{1}, x{0}, y{0}, z{0};  // w + xi + yj + zk

    constexpr Quaternion() = default;
    constexpr Quaternion(T w, T x, T y, T z) : w(w), x(x), y(y), z(z) {}

    // From axis-angle
    static Quaternion from_axis_angle(const Vec<T, 3>& axis, T angle) {
        auto n = axis.normalized();
        T half = angle / T{2};
        using std::sin; using std::cos;
        T s = sin(half);
        return {cos(half), n[0] * s, n[1] * s, n[2] * s};
    }

    // From rotation matrix (3x3)
    static Quaternion from_matrix(const Matrix<T, 3, 3>& m) {
        using std::sqrt;
        T trace = m(0, 0) + m(1, 1) + m(2, 2);
        if (trace > T{0}) {
            T s = T{0.5} / sqrt(trace + T{1});
            return {T{0.25} / s, (m(2, 1) - m(1, 2)) * s,
                    (m(0, 2) - m(2, 0)) * s, (m(1, 0) - m(0, 1)) * s};
        } else if (m(0, 0) > m(1, 1) && m(0, 0) > m(2, 2)) {
            T s = T{2} * sqrt(T{1} + m(0, 0) - m(1, 1) - m(2, 2));
            return {(m(2, 1) - m(1, 2)) / s, T{0.25} * s,
                    (m(0, 1) + m(1, 0)) / s, (m(0, 2) + m(2, 0)) / s};
        } else if (m(1, 1) > m(2, 2)) {
            T s = T{2} * sqrt(T{1} + m(1, 1) - m(0, 0) - m(2, 2));
            return {(m(0, 2) - m(2, 0)) / s, (m(0, 1) + m(1, 0)) / s,
                    T{0.25} * s, (m(1, 2) + m(2, 1)) / s};
        } else {
            T s = T{2} * sqrt(T{1} + m(2, 2) - m(0, 0) - m(1, 1));
            return {(m(1, 0) - m(0, 1)) / s, (m(0, 2) + m(2, 0)) / s,
                    (m(1, 2) + m(2, 1)) / s, T{0.25} * s};
        }
    }

    // Arithmetic
    constexpr Quaternion operator*(const Quaternion& q) const {
        return {
            w * q.w - x * q.x - y * q.y - z * q.z,
            w * q.x + x * q.w + y * q.z - z * q.y,
            w * q.y - x * q.z + y * q.w + z * q.x,
            w * q.z + x * q.y - y * q.x + z * q.w
        };
    }

    constexpr Quaternion operator+(const Quaternion& q) const {
        return {w + q.w, x + q.x, y + q.y, z + q.z};
    }

    constexpr Quaternion operator*(T s) const {
        return {w * s, x * s, y * s, z * s};
    }

    friend constexpr Quaternion operator*(T s, const Quaternion& q) { return q * s; }

    constexpr Quaternion conjugate() const { return {w, -x, -y, -z}; }

    constexpr T norm_squared() const { return w * w + x * x + y * y + z * z; }

    // non-constexpr due to std::sqrt; use norm_squared() in
    // constexpr/static_assert contexts.
    T norm() const { using std::sqrt; return sqrt(norm_squared()); }

    // non-constexpr — wraps norm() for the inverse magnitude.
    Quaternion normalized() const {
        T n = norm();
        if (n < epsilon<T>()) return {};
        return {w / n, x / n, y / n, z / n};
    }

    constexpr Result<Quaternion> inverse() const {
        T n2 = norm_squared();
        if (n2 < epsilon<T>())
            return std::unexpected(Error{ErrorCode::ZeroNorm, "cannot invert zero-norm quaternion"});
        return Quaternion{w / n2, -x / n2, -y / n2, -z / n2};
    }

    // Rotate a 3D vector: q * v * q^(-1)
    constexpr Vec<T, 3> rotate(const Vec<T, 3>& v) const {
        Quaternion p{T{0}, v[0], v[1], v[2]};
        auto r = *this * p * conjugate();
        return {r.x, r.y, r.z};
    }

    // Convert to 3x3 rotation matrix
    constexpr Matrix<T, 3, 3> to_matrix() const {
        T xx = x * x, yy = y * y, zz = z * z;
        T xy = x * y, xz = x * z, yz = y * z;
        T wx = w * x, wy = w * y, wz = w * z;

        Matrix<T, 3, 3> m;
        m(0,0) = T{1} - T{2}*(yy+zz); m(0,1) = T{2}*(xy-wz);       m(0,2) = T{2}*(xz+wy);
        m(1,0) = T{2}*(xy+wz);       m(1,1) = T{1} - T{2}*(xx+zz); m(1,2) = T{2}*(yz-wx);
        m(2,0) = T{2}*(xz-wy);       m(2,1) = T{2}*(yz+wx);       m(2,2) = T{1} - T{2}*(xx+yy);
        return m;
    }

    // Axis and angle extraction
    std::pair<Vec<T, 3>, T> to_axis_angle() const {
        using std::acos; using std::sqrt;
        auto q = normalized();
        T angle = T{2} * acos(std::clamp(q.w, T{-1}, T{1}));
        T s = sqrt(T{1} - q.w * q.w);
        if (s < epsilon<T>())
            return {Vec<T, 3>{T{1}, T{0}, T{0}}, angle};
        return {Vec<T, 3>{q.x / s, q.y / s, q.z / s}, angle};
    }

    // Spherical linear interpolation
    static Quaternion slerp(const Quaternion& a, const Quaternion& b, T t) {
        using std::acos; using std::sin;
        T dot = a.w * b.w + a.x * b.x + a.y * b.y + a.z * b.z;

        Quaternion b2 = b;
        if (dot < T{0}) { b2 = b2 * T{-1}; dot = -dot; }

        if (dot > T{1} - epsilon<T>())
            return (a * (T{1} - t) + b2 * t).normalized();

        T theta = acos(dot);
        T sin_theta = sin(theta);
        return a * (sin((T{1} - t) * theta) / sin_theta)
             + b2 * (sin(t * theta) / sin_theta);
    }

    constexpr bool operator==(const Quaternion&) const = default;
};

using Quat = Quaternion<double>;
using Quatf = Quaternion<float>;

// Free function slerp
template<Scalar T>
Quaternion<T> slerp(const Quaternion<T>& a, const Quaternion<T>& b, T t) {
    return Quaternion<T>::slerp(a, b, t);
}

} // namespace algebra
} // namespace spatium

// std::format support
template<spatium::Scalar T>
struct std::formatter<spatium::Quaternion<T>> {
    constexpr auto parse(auto& ctx) { return ctx.begin(); }
    auto format(const spatium::Quaternion<T>& q, auto& ctx) const {
        return std::format_to(ctx.out(), "({} + {}i + {}j + {}k)", q.w, q.x, q.y, q.z);
    }
};
