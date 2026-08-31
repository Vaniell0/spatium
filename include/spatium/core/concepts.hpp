#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <concepts>
#  include <cstddef>
#endif

SPATIUM_EXPORT namespace spatium {

// Sentinel for runtime-determined dimension
inline constexpr std::size_t kDynamic = static_cast<std::size_t>(-1);

// ── Scalar ─────────────────────────────────────────────────────

template<typename T>
concept Scalar = std::regular<T>
    && std::totally_ordered<T>
    && requires(T a, T b) {
        { a + b } -> std::convertible_to<T>;
        { a - b } -> std::convertible_to<T>;
        { a * b } -> std::convertible_to<T>;
        { a / b } -> std::convertible_to<T>;
        { -a }    -> std::convertible_to<T>;
        { T{0} };
        { T{1} };
    };

// ── Set ────────────────────────────────────────────────────────
// Most basic: a type with points, scalars, and dimension.

template<typename S>
concept Set = requires {
    typename S::PointType;
    typename S::ScalarType;
    requires Scalar<typename S::ScalarType>;
    { S::dimension } -> std::convertible_to<std::size_t>;
} && std::equality_comparable<typename S::PointType>;

// ── TopologicalSpace ───────────────────────────────────────────
// Adds containment check (membership in the space).

template<typename S>
concept TopologicalSpace = Set<S>
    && requires(const S& space, const typename S::PointType& p) {
        { space.contains(p) } -> std::convertible_to<bool>;
    };

// ── MetricSpace ────────────────────────────────────────────────
// Adds distance function. Axioms (symmetry, triangle inequality,
// d(x,y)=0 iff x=y) enforced via tests, not compile-time.

template<typename S>
concept MetricSpace = TopologicalSpace<S>
    && requires(const S& space,
                const typename S::PointType& p,
                const typename S::PointType& q) {
        { space.distance(p, q) } -> std::convertible_to<typename S::ScalarType>;
    };

// ── VectorSpace ────────────────────────────────────────────────
// Points support linear operations via VectorType.

template<typename S>
concept VectorSpace = Set<S>
    && requires {
        typename S::VectorType;
    }
    && requires(const typename S::PointType& p,
                const typename S::VectorType& v,
                const typename S::VectorType& w,
                typename S::ScalarType a) {
        { p + v } -> std::convertible_to<typename S::PointType>;
        { v + w } -> std::convertible_to<typename S::VectorType>;
        { a * v } -> std::convertible_to<typename S::VectorType>;
        { -v }    -> std::convertible_to<typename S::VectorType>;
    };

// ── NormedSpace ────────────────────────────────────────────────
// VectorSpace + MetricSpace with a norm on vectors.

template<typename S>
concept NormedSpace = VectorSpace<S> && MetricSpace<S>
    && requires(const S& space, const typename S::VectorType& v) {
        { space.norm(v) } -> std::convertible_to<typename S::ScalarType>;
    };

// ── InnerProductSpace ──────────────────────────────────────────
// NormedSpace with an inner product inducing the norm.

template<typename S>
concept InnerProductSpace = NormedSpace<S>
    && requires(const S& space,
                const typename S::VectorType& u,
                const typename S::VectorType& v) {
        { space.inner(u, v) } -> std::convertible_to<typename S::ScalarType>;
    };

// ── Completeness (opt-in tag) ──────────────────────────────────
// Cauchy sequences converge. Can't verify at compile-time;
// spaces declare it via static constexpr bool.

template<typename S>
concept Complete = requires {
    { S::is_complete } -> std::convertible_to<bool>;
    requires S::is_complete;
};

// ── Composed classifications ───────────────────────────────────

template<typename S>
concept BanachSpace = NormedSpace<S> && Complete<S>;

template<typename S>
concept HilbertSpace = InnerProductSpace<S> && Complete<S>;

template<typename S>
concept EuclideanSpace = HilbertSpace<S>
    && requires {
        requires (S::dimension != kDynamic);
    };

// ── Manifold sub-requirements ──────────────────────────────────
// Each named concept below checks exactly one structural requirement.
// Composing them via && in the higher-level Manifold/Surface/...
// concepts gives precise diagnostics: when a user type fails, the
// compiler reports the specific HasXxx<MyType> that did not hold,
// not an opaque "associated constraints not satisfied".

template<typename S>
concept HasTangentVector = requires {
    typename S::TangentVector;
};

template<typename S>
concept HasExpMap = HasTangentVector<S>
    && requires {
        typename S::PointType;
        typename S::ScalarType;
    }
    && requires(const S& space,
                const typename S::PointType& p,
                const typename S::TangentVector& v,
                typename S::ScalarType t) {
        { space.exp_map(p, v, t) } -> std::convertible_to<typename S::PointType>;
    };

template<typename S>
concept HasLogMap = HasTangentVector<S>
    && requires { typename S::PointType; }
    && requires(const S& space,
                const typename S::PointType& p,
                const typename S::PointType& q) {
        { space.log_map(p, q) } -> std::convertible_to<typename S::TangentVector>;
    };

template<typename S>
concept HasRiemannianMetric = HasTangentVector<S>
    && requires { typename S::PointType; typename S::ScalarType; }
    && requires(const S& space,
                const typename S::PointType& p,
                const typename S::TangentVector& u,
                const typename S::TangentVector& v) {
        { space.metric_at(p, u, v) } -> std::convertible_to<typename S::ScalarType>;
    };

template<typename S>
concept HasProject = requires { typename S::PointType; }
    && requires(const S& space, const typename S::PointType& p) {
        { space.project(p) } -> std::convertible_to<typename S::PointType>;
    };

template<typename S>
concept HasNormal = HasTangentVector<S>
    && requires { typename S::PointType; }
    && requires(const S& space, const typename S::PointType& p) {
        { space.normal(p) } -> std::convertible_to<typename S::TangentVector>;
    };

// ── Manifold ───────────────────────────────────────────────────
// Locally Euclidean topological space with exp/log maps.

template<typename S>
concept Manifold = TopologicalSpace<S> && HasExpMap<S> && HasLogMap<S>;

// ── RiemannianManifold ─────────────────────────────────────────
// Manifold + MetricSpace + pointwise inner product on tangent space.

template<typename S>
concept RiemannianManifold =
    Manifold<S> && MetricSpace<S> && HasRiemannianMetric<S>;

// ── Surface ────────────────────────────────────────────────────
// Manifold embedded in ambient space with projection and normals.
// Required by the mesh/subdivision system.

template<typename S>
concept Surface = Manifold<S> && HasProject<S> && HasNormal<S>;

} // namespace spatium
