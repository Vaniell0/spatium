#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/concepts.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <span>
#endif

SPATIUM_EXPORT namespace spatium {
inline namespace algebra {

// Free function wrappers for common vector operations.
// Mathematical convention: dot(a, b) instead of a.dot(b).

template<Scalar T, std::size_t N>
constexpr T dot(const Vec<T, N>& a, const Vec<T, N>& b) { return a.dot(b); }

template<Scalar T>
constexpr Vec<T, 3> cross(const Vec<T, 3>& a, const Vec<T, 3>& b) { return a.cross(b); }

template<Scalar T, std::size_t N>
Vec<T, N> normalize(const Vec<T, N>& v) { return v.normalized(); }

template<Scalar T, std::size_t N>
constexpr Vec<T, N> lerp(const Vec<T, N>& a, const Vec<T, N>& b, T t) { return a.lerp(b, t); }

template<Scalar T, std::size_t N>
T distance(const Vec<T, N>& a, const Vec<T, N>& b) { return a.distance_to(b); }

// ── Generic algebraic functions ───────────────────────────────

// Group: repeated composition (power) via square-and-multiply.
// a^n for n>0, uses inverse for n<0, identity for n==0.
template<Group G>
auto power(const G& g, const typename G::ElementType& a, int n) -> typename G::ElementType {
    if (n == 0) return g.identity();
    auto base = (n > 0) ? a : g.inverse(a);
    auto abs_n = (n > 0) ? n : -n;
    auto result = g.identity();
    auto curr = base;
    while (abs_n > 0) {
        if (abs_n & 1) result = g.compose(result, curr);
        curr = g.compose(curr, curr);
        abs_n >>= 1;
    }
    return result;
}

// Group: commutator [a, b] = a · b · a⁻¹ · b⁻¹
template<Group G>
auto commutator(const G& g, const typename G::ElementType& a,
                const typename G::ElementType& b) -> typename G::ElementType {
    return g.compose(g.compose(a, b), g.compose(g.inverse(a), g.inverse(b)));
}

// LieGroup: adjoint action Ad_g(v) ≈ log(g · exp(εv) · g⁻¹) / ε
template<LieGroup G>
auto adjoint(const G& g, const typename G::ElementType& elem,
             const typename G::AlgebraType& v,
             double eps = 1e-6) -> typename G::AlgebraType {
    auto gv = g.compose(elem, g.exp(v * eps));
    auto conj = g.compose(gv, g.inverse(elem));
    return g.log(conj) / eps;
}

// Ring: polynomial evaluation via Horner's method.
// coeffs = {a₀, a₁, ..., aₙ} → a₀ + a₁x + ... + aₙxⁿ
template<Ring R>
auto poly_eval(const R& ring, std::span<const typename R::ElementType> coeffs,
               const typename R::ElementType& x) -> typename R::ElementType {
    if (coeffs.empty()) return ring.zero();
    auto result = coeffs.back();
    for (std::size_t i = coeffs.size() - 1; i > 0; --i)
        result = ring.add(ring.multiply(result, x), coeffs[i - 1]);
    return result;
}

} // namespace algebra
} // namespace spatium
