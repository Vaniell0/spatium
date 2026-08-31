#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/algebra/vec_expr.hpp>
#  include <spatium/algebra/vec_simd.hpp>
#  include <algorithm>
#  include <array>
#  include <cassert>
#  include <cmath>
#  include <cstddef>
#endif

SPATIUM_EXPORT namespace spatium {
inline namespace algebra {

// SIMD dispatch helpers — true when runtime SIMD is available for this T×N
// internal — do not use, no API stability. The has_simd_* booleans
// are an implementation detail of Vec's SIMD dispatch and may
// disappear or change semantics across releases.
namespace simd_detail {

template<typename T, std::size_t N>
inline constexpr bool has_simd_f4 =
#if defined(__SSE2__)
    std::same_as<T, float> && N == 4;
#else
    false;
#endif

template<typename T, std::size_t N>
inline constexpr bool has_simd_d4 =
#if defined(__AVX2__)
    std::same_as<T, double> && N == 4;
#else
    false;
#endif

template<typename T, std::size_t N>
inline constexpr bool has_simd = has_simd_f4<T, N> || has_simd_d4<T, N>;

} // namespace simd_detail

template<Scalar T, std::size_t N>
struct Vec {
    using scalar_type = T;
    static constexpr std::size_t size = N;
    using is_vec_tag = void;

    std::array<T, N> data{};

    constexpr Vec() = default;

    template<typename... Args>
        requires (sizeof...(Args) == N && (std::convertible_to<Args, T> && ...))
    constexpr Vec(Args... args) : data{static_cast<T>(args)...} {}

    // Converting constructor from any VecLike expression
    template<VecLike E>
        requires std::convertible_to<typename E::scalar_type, T>
              && (E::size == N)
              && (!std::same_as<std::remove_cvref_t<E>, Vec>)
    constexpr Vec(const E& expr) {
        for (std::size_t i = 0; i < N; ++i)
            data[i] = static_cast<T>(expr[i]);
    }

    constexpr T& operator[](std::size_t i) { return data[i]; }
    constexpr const T& operator[](std::size_t i) const { return data[i]; }

    // Arithmetic — return expressions
    constexpr auto operator+(const Vec& rhs) const {
        return VecBinExpr<OpAdd, Vec, Vec>{*this, rhs};
    }

    constexpr auto operator-(const Vec& rhs) const {
        return VecBinExpr<OpSub, Vec, Vec>{*this, rhs};
    }

    constexpr auto operator-() const {
        return VecNegExpr<Vec>{*this};
    }

    constexpr auto operator*(T s) const {
        return VecScalarExpr<OpMul, Vec, T>{*this, s};
    }

    constexpr auto operator/(T s) const {
        assert(s != T{0} && "Vec::operator/: division by zero");
        return VecScalarExpr<OpDiv, Vec, T>{*this, s};
    }

    friend constexpr auto operator*(T s, const Vec& v) {
        return VecScalarExpr<OpMul, Vec, T>{v, s};
    }

    // Compound assignment — accept VecLike, SIMD for Vec+=Vec
    template<VecLike E>
        requires std::convertible_to<typename E::scalar_type, T> && (E::size == N)
    constexpr Vec& operator+=(const E& rhs) {
        if consteval {
            for (std::size_t i = 0; i < N; ++i) data[i] += static_cast<T>(rhs[i]);
        } else {
            if constexpr (simd_detail::has_simd_f4<T, N> && std::same_as<std::remove_cvref_t<E>, Vec>) {
                simd::add_f4(data.data(), rhs.data.data(), data.data());
            } else if constexpr (simd_detail::has_simd_d4<T, N> && std::same_as<std::remove_cvref_t<E>, Vec>) {
                simd::add_d4(data.data(), rhs.data.data(), data.data());
            } else {
                for (std::size_t i = 0; i < N; ++i) data[i] += static_cast<T>(rhs[i]);
            }
        }
        return *this;
    }

    template<VecLike E>
        requires std::convertible_to<typename E::scalar_type, T> && (E::size == N)
    constexpr Vec& operator-=(const E& rhs) {
        if consteval {
            for (std::size_t i = 0; i < N; ++i) data[i] -= static_cast<T>(rhs[i]);
        } else {
            if constexpr (simd_detail::has_simd_f4<T, N> && std::same_as<std::remove_cvref_t<E>, Vec>) {
                simd::sub_f4(data.data(), rhs.data.data(), data.data());
            } else if constexpr (simd_detail::has_simd_d4<T, N> && std::same_as<std::remove_cvref_t<E>, Vec>) {
                simd::sub_d4(data.data(), rhs.data.data(), data.data());
            } else {
                for (std::size_t i = 0; i < N; ++i) data[i] -= static_cast<T>(rhs[i]);
            }
        }
        return *this;
    }

    constexpr Vec& operator*=(T s) {
        if consteval {
            for (std::size_t i = 0; i < N; ++i) data[i] *= s;
        } else {
            if constexpr (simd_detail::has_simd_f4<T, N>) {
                simd::mul_scalar_f4(data.data(), s, data.data());
            } else if constexpr (simd_detail::has_simd_d4<T, N>) {
                simd::mul_scalar_d4(data.data(), s, data.data());
            } else {
                for (std::size_t i = 0; i < N; ++i) data[i] *= s;
            }
        }
        return *this;
    }

    constexpr Vec& operator/=(T s) {
        assert(s != T{0} && "Vec::operator/=: division by zero");
        if consteval {
            for (std::size_t i = 0; i < N; ++i) data[i] /= s;
        } else {
            if constexpr (simd_detail::has_simd_f4<T, N>) {
                simd::div_scalar_f4(data.data(), s, data.data());
            } else if constexpr (simd_detail::has_simd_d4<T, N>) {
                simd::div_scalar_d4(data.data(), s, data.data());
            } else {
                for (std::size_t i = 0; i < N; ++i) data[i] /= s;
            }
        }
        return *this;
    }

    constexpr bool operator==(const Vec&) const = default;

    // Dot product — SIMD for Vec.dot(Vec)
    template<VecLike E>
        requires std::convertible_to<typename E::scalar_type, T> && (E::size == N)
    constexpr T dot(const E& rhs) const {
        if consteval {
            T sum{0};
            for (std::size_t i = 0; i < N; ++i)
                sum += data[i] * static_cast<T>(rhs[i]);
            return sum;
        } else {
            if constexpr (simd_detail::has_simd_f4<T, N> && std::same_as<std::remove_cvref_t<E>, Vec>) {
                return simd::dot_f4(data.data(), rhs.data.data());
            } else if constexpr (simd_detail::has_simd_d4<T, N> && std::same_as<std::remove_cvref_t<E>, Vec>) {
                return simd::dot_d4(data.data(), rhs.data.data());
            } else {
                T sum{0};
                for (std::size_t i = 0; i < N; ++i)
                    sum += data[i] * static_cast<T>(rhs[i]);
                return sum;
            }
        }
    }

    constexpr T norm_squared() const { return dot(*this); }

    // non-constexpr due to std::sqrt; use norm_squared() in
    // constexpr/static_assert contexts.
    T norm() const {
        using std::sqrt;
        return sqrt(norm_squared());
    }

    // non-constexpr — wraps norm() for the inverse magnitude.
    Vec normalized() const {
        auto n = norm();
        if (n < epsilon<T>()) return Vec{};
        return Vec{*this / n};
    }

    // Linear interpolation: this + t * (other - this)
    constexpr Vec lerp(const Vec& other, T t) const {
        Vec result;
        for (std::size_t i = 0; i < N; ++i)
            result[i] = data[i] + t * (other[i] - data[i]);
        return result;
    }

    // Distance to another vector
    T distance_to(const Vec& other) const {
        return Vec{*this - other}.norm();
    }

    // Reflect across a unit normal: v - 2*(v.n)*n
    constexpr Vec reflect(const Vec& unit_normal) const {
        return Vec{*this - unit_normal * (T{2} * dot(unit_normal))};
    }

    // Static factories
    static constexpr Vec zero() { return Vec{}; }

    static constexpr Vec unit(std::size_t axis) {
        Vec v{};
        v[axis] = T{1};
        return v;
    }

    // Component-wise operations
    constexpr Vec abs() const {
        Vec result;
        for (std::size_t i = 0; i < N; ++i)
            result[i] = data[i] < T{0} ? -data[i] : data[i];
        return result;
    }

    static constexpr Vec min(const Vec& a, const Vec& b) {
        Vec result;
        for (std::size_t i = 0; i < N; ++i)
            result[i] = a[i] < b[i] ? a[i] : b[i];
        return result;
    }

    static constexpr Vec max(const Vec& a, const Vec& b) {
        Vec result;
        for (std::size_t i = 0; i < N; ++i)
            result[i] = a[i] > b[i] ? a[i] : b[i];
        return result;
    }

    constexpr Vec clamp(const Vec& lo, const Vec& hi) const {
        return max(lo, min(*this, hi));
    }

    // Cross product (3D only)
    template<VecLike E>
        requires (N == 3)
              && std::convertible_to<typename E::scalar_type, T>
              && (E::size == 3)
    constexpr Vec cross(const E& rhs) const {
        return Vec{
            data[1] * static_cast<T>(rhs[2]) - data[2] * static_cast<T>(rhs[1]),
            data[2] * static_cast<T>(rhs[0]) - data[0] * static_cast<T>(rhs[2]),
            data[0] * static_cast<T>(rhs[1]) - data[1] * static_cast<T>(rhs[0])
        };
    }
};

using Vec2 = Vec<double, 2>;
using Vec3 = Vec<double, 3>;
using Vec4 = Vec<double, 4>;

using Vec2f = Vec<float, 2>;
using Vec3f = Vec<float, 3>;
using Vec4f = Vec<float, 4>;

} // namespace algebra
} // namespace spatium
