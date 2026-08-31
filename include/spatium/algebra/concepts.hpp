#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <concepts>
#endif

SPATIUM_EXPORT namespace spatium::inline algebra {

// ── Algebraic structure concepts ───────────────────────────────
// Model the algebraic branch of the mathematical spaces hierarchy.
// Each level adds axioms over the previous.

// Magma: closed binary operation
template<typename G>
concept Magma = requires {
    typename G::ElementType;
} && requires(const G& g, const typename G::ElementType& a, const typename G::ElementType& b) {
    { g.compose(a, b) } -> std::convertible_to<typename G::ElementType>;
};

// Semigroup: associative Magma
// (associativity checked via verify, not compile-time)
template<typename G>
concept Semigroup = Magma<G>;

// Monoid: Semigroup with identity
template<typename G>
concept Monoid = Semigroup<G>
    && requires(const G& g) {
        { g.identity() } -> std::convertible_to<typename G::ElementType>;
    };

// Group: Monoid with inverses
template<typename G>
concept Group = Monoid<G>
    && requires(const G& g, const typename G::ElementType& a) {
        { g.inverse(a) } -> std::convertible_to<typename G::ElementType>;
    };

// AbelianGroup: commutative Group
// (commutativity checked via verify, not compile-time)
template<typename G>
concept AbelianGroup = Group<G>;

// Ring: two operations (+, *) where + is AbelianGroup, * is Monoid, * distributes over +
template<typename R>
concept Ring = requires {
    typename R::ElementType;
} && requires(const R& r, const typename R::ElementType& a, const typename R::ElementType& b) {
    { r.add(a, b) }      -> std::convertible_to<typename R::ElementType>;
    { r.negate(a) }       -> std::convertible_to<typename R::ElementType>;
    { r.zero() }          -> std::convertible_to<typename R::ElementType>;
    { r.multiply(a, b) }  -> std::convertible_to<typename R::ElementType>;
    { r.one() }           -> std::convertible_to<typename R::ElementType>;
};

// Field: Ring where nonzero elements have multiplicative inverse
template<typename F>
concept Field = Ring<F>
    && requires(const F& f, const typename F::ElementType& a) {
        { f.reciprocal(a) } -> std::convertible_to<typename F::ElementType>;
    };

// ── Lie group / Lie algebra ────────────────────────────────────

// LieGroup: Group + smooth manifold structure (exp/log on the algebra)
template<typename G>
concept LieGroup = Group<G>
    && requires {
        typename G::AlgebraType; // tangent at identity = Lie algebra element
    }
    && requires(const G& g,
                const typename G::AlgebraType& v,
                const typename G::ElementType& a) {
        { g.exp(v) }       -> std::convertible_to<typename G::ElementType>;
        { g.log(a) }       -> std::convertible_to<typename G::AlgebraType>;
    };

// LieAlgebra: vector space with bracket [x, y]
template<typename A>
concept LieAlgebra = requires {
    typename A::ElementType;
    typename A::ScalarType;
} && requires(const A& a,
              const typename A::ElementType& x,
              const typename A::ElementType& y,
              typename A::ScalarType s) {
    { a.bracket(x, y) }  -> std::convertible_to<typename A::ElementType>;
    { a.add(x, y) }      -> std::convertible_to<typename A::ElementType>;
    { a.scale(s, x) }    -> std::convertible_to<typename A::ElementType>;
};

} // namespace spatium::algebra
