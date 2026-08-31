#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <cassert>
#  include <cmath>
#  include <cstddef>
#  include <type_traits>
#endif

SPATIUM_EXPORT namespace spatium {
inline namespace algebra {

// Forward-declare Vec for storage specialization
template<Scalar T, std::size_t N> struct Vec;

// ── VecLike concept ───────────────────────────────────────────

template<typename E>
concept VecLike = requires(const E& e, std::size_t i) {
    typename E::scalar_type;
    { E::size } -> std::convertible_to<std::size_t>;
    { e[i] } -> std::convertible_to<typename E::scalar_type>;
};

// ── Storage trait ─────────────────────────────────────────────

// internal — do not use, no API stability. Names inside `detail`
// (expr_storage_t, vec_dot/vec_norm/vec_cross helpers) back the
// SPATIUM_VEC_EXPR_EAGER_METHODS macro and may change without
// notice. Public users should call Vec::dot/.norm/.cross instead.
namespace expr_detail {

// All expressions stored by value. Vec is small (array<T,N>),
// and value semantics preserve constexpr compatibility.
template<typename E>
using expr_storage_t = E;

// ── Shared implementations for expr methods ───────────────────

template<VecLike A, VecLike B>
    requires std::same_as<typename A::scalar_type, typename B::scalar_type>
          && (A::size == B::size)
constexpr typename A::scalar_type vec_dot(const A& a, const B& b) {
    using T = typename A::scalar_type;
    T sum{0};
    for (std::size_t i = 0; i < A::size; ++i)
        sum += static_cast<T>(a[i]) * static_cast<T>(b[i]);
    return sum;
}

template<VecLike E>
constexpr typename E::scalar_type vec_norm_squared(const E& e) {
    return vec_dot(e, e);
}

template<VecLike E>
typename E::scalar_type vec_norm(const E& e) {
    using std::sqrt;
    return sqrt(vec_norm_squared(e));
}

template<VecLike E>
Vec<typename E::scalar_type, E::size> vec_normalized(const E& e) {
    using T = typename E::scalar_type;
    constexpr auto N = E::size;
    auto n = vec_norm(e);
    if (n < epsilon<T>()) return Vec<T, N>{};
    Vec<T, N> result{e};
    for (std::size_t i = 0; i < N; ++i) result[i] /= n;
    return result;
}

template<VecLike A, VecLike B>
    requires (A::size == 3) && (B::size == 3)
          && std::same_as<typename A::scalar_type, typename B::scalar_type>
constexpr Vec<typename A::scalar_type, 3> vec_cross(const A& a, const B& b) {
    using T = typename A::scalar_type;
    return Vec<T, 3>{
        static_cast<T>(a[1]) * static_cast<T>(b[2]) - static_cast<T>(a[2]) * static_cast<T>(b[1]),
        static_cast<T>(a[2]) * static_cast<T>(b[0]) - static_cast<T>(a[0]) * static_cast<T>(b[2]),
        static_cast<T>(a[0]) * static_cast<T>(b[1]) - static_cast<T>(a[1]) * static_cast<T>(b[0])
    };
}

} // namespace expr_detail

// ── Operations ────────────────────────────────────────────────

struct OpAdd { template<typename T> constexpr T operator()(T a, T b) const { return a + b; } };
struct OpSub { template<typename T> constexpr T operator()(T a, T b) const { return a - b; } };
struct OpMul { template<typename T> constexpr T operator()(T a, T b) const { return a * b; } };
struct OpDiv { template<typename T> constexpr T operator()(T a, T b) const { return a / b; } };

// Eager method set — added to each expression type (preserves aggregate)
#define SPATIUM_VEC_EXPR_EAGER_METHODS                                        \
    template<VecLike E>                                                       \
        requires std::convertible_to<typename E::scalar_type, scalar_type>    \
              && (E::size == size)                                             \
    constexpr scalar_type dot(const E& rhs) const {                           \
        return expr_detail::vec_dot(*this, rhs);                                   \
    }                                                                         \
    constexpr scalar_type norm_squared() const {                              \
        return expr_detail::vec_norm_squared(*this);                               \
    }                                                                         \
    scalar_type norm() const {                                                \
        return expr_detail::vec_norm(*this);                                       \
    }                                                                         \
    Vec<scalar_type, size> normalized() const {                               \
        return expr_detail::vec_normalized(*this);                                 \
    }                                                                         \
    template<VecLike E>                                                       \
        requires (size == 3) && (E::size == 3)                                \
              && std::convertible_to<typename E::scalar_type, scalar_type>    \
    constexpr Vec<scalar_type, 3> cross(const E& rhs) const {                \
        return expr_detail::vec_cross(*this, rhs);                                 \
    }

// ── Binary expression ─────────────────────────────────────────

template<typename Op, VecLike LHS, VecLike RHS>
    requires std::same_as<typename LHS::scalar_type, typename RHS::scalar_type>
          && (LHS::size == RHS::size)
struct VecBinExpr {
    using scalar_type = typename LHS::scalar_type;
    static constexpr std::size_t size = LHS::size;

    expr_detail::expr_storage_t<LHS> lhs;
    expr_detail::expr_storage_t<RHS> rhs;

    constexpr scalar_type operator[](std::size_t i) const {
        return Op{}(lhs[i], rhs[i]);
    }

    SPATIUM_VEC_EXPR_EAGER_METHODS
};

// ── Scalar broadcast expression ───────────────────────────────

template<typename Op, VecLike Expr, Scalar S>
    requires std::convertible_to<S, typename Expr::scalar_type>
struct VecScalarExpr {
    using scalar_type = typename Expr::scalar_type;
    static constexpr std::size_t size = Expr::size;

    expr_detail::expr_storage_t<Expr> expr;
    scalar_type scalar;

    constexpr scalar_type operator[](std::size_t i) const {
        return Op{}(expr[i], scalar);
    }

    SPATIUM_VEC_EXPR_EAGER_METHODS
};

// ── Unary negate expression ───────────────────────────────────

template<VecLike Expr>
struct VecNegExpr {
    using scalar_type = typename Expr::scalar_type;
    static constexpr std::size_t size = Expr::size;

    expr_detail::expr_storage_t<Expr> expr;

    constexpr scalar_type operator[](std::size_t i) const {
        return -expr[i];
    }

    SPATIUM_VEC_EXPR_EAGER_METHODS
};

#undef SPATIUM_VEC_EXPR_EAGER_METHODS

// ── Free operators on VecLike (non-Vec) ───────────────────────

template<VecLike L, VecLike R>
    requires std::same_as<typename L::scalar_type, typename R::scalar_type>
          && (L::size == R::size)
          && (!requires { typename L::is_vec_tag; } || !requires { typename R::is_vec_tag; })
constexpr auto operator+(const L& l, const R& r) {
    return VecBinExpr<OpAdd, L, R>{l, r};
}

template<VecLike L, VecLike R>
    requires std::same_as<typename L::scalar_type, typename R::scalar_type>
          && (L::size == R::size)
          && (!requires { typename L::is_vec_tag; } || !requires { typename R::is_vec_tag; })
constexpr auto operator-(const L& l, const R& r) {
    return VecBinExpr<OpSub, L, R>{l, r};
}

template<VecLike E>
    requires (!requires { typename E::is_vec_tag; })
constexpr auto operator-(const E& e) {
    return VecNegExpr<E>{e};
}

template<VecLike E>
    requires (!requires { typename E::is_vec_tag; })
constexpr auto operator*(const E& e, typename E::scalar_type s) {
    return VecScalarExpr<OpMul, E, typename E::scalar_type>{e, s};
}

template<VecLike E>
    requires (!requires { typename E::is_vec_tag; })
constexpr auto operator*(typename E::scalar_type s, const E& e) {
    return VecScalarExpr<OpMul, E, typename E::scalar_type>{e, s};
}

template<VecLike E>
    requires (!requires { typename E::is_vec_tag; })
constexpr auto operator/(const E& e, typename E::scalar_type s) {
    return VecScalarExpr<OpDiv, E, typename E::scalar_type>{e, s};
}

} // namespace algebra
} // namespace spatium
