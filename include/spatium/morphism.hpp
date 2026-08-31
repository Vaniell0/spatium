#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/error.hpp>
#  include <spatium/point.hpp>
#  include <functional>
#  include <optional>
#  include <utility>
#endif

SPATIUM_EXPORT namespace spatium {

// ── Morphism ───────────────────────────────────────────────────
// Runtime-polymorphic map between spaces. Uses std::function for flexibility.

template<Set From, Set To>
struct Morphism {
    using DomainPoint = typename From::PointType;
    using CodomainPoint = typename To::PointType;
    using MapFn = std::function<CodomainPoint(const DomainPoint&)>;

    MapFn forward;
    std::optional<std::function<DomainPoint(const CodomainPoint&)>> inverse;

    // Apply: morphism(point)
    CodomainPoint operator()(const DomainPoint& p) const {
        return forward(p);
    }

    // Apply to typed Point
    Point<To> operator()(const Point<From>& p) const {
        return Point<To>{forward(p.raw())};
    }

    bool has_inverse() const { return inverse.has_value(); }

    // Return the inverse morphism (To → From)
    Morphism<To, From> invert() const {
        return {
            .forward = *inverse,
            .inverse = forward,
        };
    }
};

// ── Composition ────────────────────────────────────────────────
// g * f = g ∘ f  (apply f first, then g)

template<Set A, Set B, Set C>
Morphism<A, C> operator*(const Morphism<B, C>& g, const Morphism<A, B>& f) {
    auto inv = (f.inverse && g.inverse)
        ? std::optional<std::function<typename A::PointType(const typename C::PointType&)>>(
            [f_inv = *f.inverse, g_inv = *g.inverse](const typename C::PointType& p) {
                return f_inv(g_inv(p));
            })
        : std::nullopt;
    return {
        .forward = [f_fn = f.forward, g_fn = g.forward](const typename A::PointType& p) {
            return g_fn(f_fn(p));
        },
        .inverse = std::move(inv),
    };
}

// ── Pipe operator ──────────────────────────────────────────────
// point | morphism  (apply morphism to point)

template<Set From, Set To>
Point<To> operator|(const Point<From>& p, const Morphism<From, To>& m) {
    return m(p);
}

// Chain morphisms: f | g = g ∘ f  (pipe order: apply f, then g)
template<Set A, Set B, Set C>
Morphism<A, C> operator|(const Morphism<A, B>& f, const Morphism<B, C>& g) {
    return g * f;
}

// Pipe through Result: Result<Point> | morphism  (auto-unwrap)
template<Set From, Set To>
Result<Point<To>> operator|(const Result<Point<From>>& p, const Morphism<From, To>& m) {
    if (!p) return std::unexpected(p.error());
    return m(*p);
}

// ── Factory helpers ────────────────────────────────────────────

// Create a morphism from a lambda
template<Set From, Set To, typename F>
    requires std::invocable<F, const typename From::PointType&>
Morphism<From, To> morph(F&& f) {
    return {.forward = std::forward<F>(f), .inverse = std::nullopt};
}

// Create a morphism with inverse
template<Set From, Set To, typename F, typename G>
    requires std::invocable<F, const typename From::PointType&>
          && std::invocable<G, const typename To::PointType&>
Morphism<From, To> morph(F&& f, G&& inv) {
    return {
        .forward = std::forward<F>(f),
        .inverse = std::forward<G>(inv),
    };
}

// ── Identity morphism ──────────────────────────────────────────

template<Set S>
Morphism<S, S> identity() {
    return {.forward = [](const typename S::PointType& p) { return p; },
            .inverse = [](const typename S::PointType& p) { return p; }};
}

} // namespace spatium
