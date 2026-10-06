#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/core/error.hpp>
#  include <spatium/algebra/quaternion.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <cmath>
#  include <cstddef>
#  include <format>
#endif

SPATIUM_EXPORT namespace spatium {
inline namespace algebra {

// The octonions: the normed division algebra of dimension 8 (Hurwitz, 1898),
// the last of the four. Not associative -- (ab)c is not always a(bc) -- but
// alternative: the subalgebra any two elements generate is associative, so
// a(ab) = (aa)b, and |ab| = |a||b| holds. That identity is why they are here:
// the 7-dimensional cross product, a global tangent frame on S^7, and the
// row of docs/dimensions.md that said "allowed by the theorem, not written".
//
// Built by the Cayley-Dickson doubling of the quaternions the library already
// has, rather than from a table of structure constants:
//
//   (a, b)(c, d) = (a c - d* b,  d a + b c*)         a, b, c, d quaternions
//
// so an octonion is a + b l with l^2 = -1, and the first four components are
// a `Quaternion` multiplied by `Quaternion`'s own product. Components e0..e7
// are a.w a.x a.y a.z b.w b.x b.y b.z; e1, e2, e3 are i, j, k and e4 is l.
//
// The doubling is the only definition here, and it is held by properties
// that determine the octonions up to isomorphism rather than by a second
// copy of the table: composition (|ab| = |a||b|), alternativity, the Moufang
// identities, non-associativity, and its embeddings of Quaternion and Complex
// (tests/test_octonion.cpp).

namespace octonion_detail {
template<Scalar T>
constexpr Quaternion<T> sub(const Quaternion<T>& p, const Quaternion<T>& q) {
    return {p.w - q.w, p.x - q.x, p.y - q.y, p.z - q.z};
}
}  // namespace octonion_detail

template<Scalar T = double>
struct Octonion {
    Quaternion<T> a{};                                      // e0 e1 e2 e3: the real part is a.w
    Quaternion<T> b{T{0}, T{0}, T{0}, T{0}};                // e4 e5 e6 e7: a + b l

    constexpr Octonion() = default;                          // 1
    constexpr Octonion(const Quaternion<T>& lo, const Quaternion<T>& hi) : a(lo), b(hi) {}
    constexpr Octonion(T e0, T e1, T e2, T e3, T e4, T e5, T e6, T e7)
        : a(e0, e1, e2, e3), b(e4, e5, e6, e7) {}

    // The unit e_i: e0 = 1, e1..e7 imaginary.
    static constexpr Octonion basis(std::size_t i) {
        Octonion o{T{0}, T{0}, T{0}, T{0}, T{0}, T{0}, T{0}, T{0}};
        o[i] = T{1};
        return o;
    }

    static constexpr Octonion from_vec8(const Vec<T, 8>& v) {
        return {v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7]};
    }
    // The imaginary octonion with the seven components v: what a 7-vector is.
    static constexpr Octonion from_imag(const Vec<T, 7>& v) {
        return {T{0}, v[0], v[1], v[2], v[3], v[4], v[5], v[6]};
    }

    constexpr T& operator[](std::size_t i) {
        switch (i) {
            case 0: return a.w;  case 1: return a.x;  case 2: return a.y;  case 3: return a.z;
            case 4: return b.w;  case 5: return b.x;  case 6: return b.y;  default: return b.z;
        }
    }
    constexpr const T& operator[](std::size_t i) const {
        switch (i) {
            case 0: return a.w;  case 1: return a.x;  case 2: return a.y;  case 3: return a.z;
            case 4: return b.w;  case 5: return b.x;  case 6: return b.y;  default: return b.z;
        }
    }

    constexpr Vec<T, 8> to_vec8() const {
        return Vec<T, 8>{a.w, a.x, a.y, a.z, b.w, b.x, b.y, b.z};
    }
    constexpr T real() const { return a.w; }
    constexpr Vec<T, 7> imag() const { return Vec<T, 7>{a.x, a.y, a.z, b.w, b.x, b.y, b.z}; }

    constexpr Octonion operator+(const Octonion& o) const { return {a + o.a, b + o.b}; }
    constexpr Octonion operator-(const Octonion& o) const {
        return {octonion_detail::sub(a, o.a), octonion_detail::sub(b, o.b)};
    }
    constexpr Octonion operator-() const { return {a * T{-1}, b * T{-1}}; }
    constexpr Octonion operator*(T s) const { return {a * s, b * s}; }
    friend constexpr Octonion operator*(T s, const Octonion& o) { return o * s; }

    // (a, b)(c, d) = (a c - d* b, d a + b c*)
    constexpr Octonion operator*(const Octonion& o) const {
        const Quaternion<T>& c = o.a;
        const Quaternion<T>& d = o.b;
        return {octonion_detail::sub(a * c, d.conjugate() * b), d * a + b * c.conjugate()};
    }

    constexpr Octonion conjugate() const { return {a.conjugate(), b * T{-1}}; }
    constexpr T norm_squared() const { return a.norm_squared() + b.norm_squared(); }
    T norm() const { using std::sqrt; return sqrt(norm_squared()); }
    constexpr T dot(const Octonion& o) const {
        return a.w * o.a.w + a.x * o.a.x + a.y * o.a.y + a.z * o.a.z +
               b.w * o.b.w + b.x * o.b.x + b.y * o.b.y + b.z * o.b.z;
    }

    // a^-1 = a* / |a|^2 -- two-sided, because the octonions are alternative.
    constexpr Result<Octonion> inverse() const {
        const T n2 = norm_squared();
        if (n2 < epsilon<T>())
            return std::unexpected(Error{ErrorCode::ZeroNorm, "cannot invert zero-norm octonion"});
        return conjugate() * (T{1} / n2);
    }

    constexpr bool operator==(const Octonion&) const = default;
};

using Oct = Octonion<double>;

// The cross product of R^7: the imaginary part of the product of two
// imaginary octonions, x y = -x.y + x cross y, so that
//   x cross y is orthogonal to x and to y, and
//   |x cross y|^2 = |x|^2 |y|^2 - (x.y)^2
// -- the identities that make a binary cross product, which exist at
// dimensions 3 and 7 only (Brown and Gray, 1967). Unlike the 3-dimensional
// one it does not satisfy the Jacobi identity.
template<Scalar T>
Vec<T, 7> cross7(const Vec<T, 7>& x, const Vec<T, 7>& y) {
    return (Octonion<T>::from_imag(x) * Octonion<T>::from_imag(y)).imag();
}

} // namespace algebra
} // namespace spatium

// std::format support
template<spatium::Scalar T>
struct std::formatter<spatium::Octonion<T>> {
    constexpr auto parse(auto& ctx) { return ctx.begin(); }
    auto format(const spatium::Octonion<T>& o, auto& ctx) const {
        return std::format_to(ctx.out(), "(e0 {}, e1 {}, e2 {}, e3 {}, e4 {}, e5 {}, e6 {}, e7 {})",
                              o[0], o[1], o[2], o[3], o[4], o[5], o[6], o[7]);
    }
};
