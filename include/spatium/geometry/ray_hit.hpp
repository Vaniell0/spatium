#pragma once

// Concept-driven ray-hit dispatch.
//
// `ray_hit(ray, shape)` returns an optional unified `RayHit3<T>`
// regardless of whether the shape is a triangle (Möller-Trumbore),
// a quadric (closed-form), a torus (quartic), or a user-defined
// type that supplies its own overload. This is the dispatch
// channel BVH and any future ray-tracing front-end can use without
// branching on `is_same_v<Shape, Triangle<3, T>>`.
//
// Adding a new ray-hittable type means: write `std::optional<
// RayHit3<T>> ray_hit(const Ray<3, T>&, const MyShape&)`. The
// `RayHittable<S>` concept then matches it automatically.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/geometry/intersection.hpp>
#  include <spatium/geometry/line.hpp>
#  include <spatium/geometry/ray_surface.hpp>
#  include <spatium/geometry/triangle.hpp>
#  include <optional>
#endif

SPATIUM_EXPORT namespace spatium::geometry {

template<Scalar T>
struct RayHit3 {
    T t;                       // ray parameter; point = ray.origin + t·direction
    Vec<T, 3> point;
    Vec<T, 3> normal{};        // unit-length surface normal (zero if undefined)
    T u{};                     // barycentric weight on vertex 1 (Triangle only)
    T v{};                     // barycentric weight on vertex 2 (Triangle only)
};

// Triangle<3, T> — Möller-Trumbore, returns full barycentric + normal.
template<Scalar T>
inline std::optional<RayHit3<T>> ray_hit(const Ray<3, T>& ray,
                                         const Triangle<3, T>& tri) {
    auto h = ray_triangle(ray, tri);
    if (!h) return std::nullopt;
    return RayHit3<T>{h->t, h->point, h->normal, h->u, h->v};
}

// Quadric<T> — first forward (t ≥ 0) closed-form hit, if any.
template<Scalar T>
inline std::optional<RayHit3<T>> ray_hit(const Ray<3, T>& ray,
                                         const Quadric<T>& q) {
    auto hits = ray_quadric(ray, q);
    if (hits.empty()) return std::nullopt;
    auto& h = hits.front();
    return RayHit3<T>{h.t, h.point, h.normal, T{0}, T{0}};
}

// BoundedQuadric<T> — nearest forward hit that also lies inside the clip.
// Not simply the nearest hit: a ray crossing a truncated cylinder beyond
// its end passes through the infinite surface twice, and both of those
// solutions are real. Rejecting the ones outside the box is what makes
// the truncation mean something rather than just labelling the shape
// finite.
//
// Deliberately does not go through ray_quadric(): that returns a
// std::vector, so it heap-allocates on every call that hits and then
// sorts a list of at most two. Harmless for a direct call, but this
// overload is a BVH leaf test invoked many times per ray, where the
// allocation dominates the arithmetic it is wrapping. A quadratic has
// at most two roots, so the nearest valid one can be picked in place.
template<Scalar T>
inline std::optional<RayHit3<T>> ray_hit(const Ray<3, T>& ray,
                                         const BoundedQuadric<T>& bq) {
    auto [a, b, c] = detail::quadric_coeffs(ray, bq.surface);
    std::optional<RayHit3<T>> best;
    for (const auto& root : solve_quadratic(a, b, c)) {
        if (!root.is_real() || root.re < T{0}) continue;
        if (best && root.re >= best->t) continue;
        Vec<T, 3> pt{ray.origin + ray.direction * root.re};
        if (!bq.within_clip(pt)) continue;
        best = RayHit3<T>{root.re, pt, bq.surface.normal(pt), T{0}, T{0}};
    }
    return best;
}

// Torus<T> — first forward quartic hit, if any.
template<Scalar T>
inline std::optional<RayHit3<T>> ray_hit(const Ray<3, T>& ray,
                                         const Torus<T>& t) {
    auto hits = ray_torus(ray, t);
    if (hits.empty()) return std::nullopt;
    auto& h = hits.front();
    return RayHit3<T>{h.t, h.point, h.normal, T{0}, T{0}};
}

// Concept: any type with a matching `ray_hit` overload qualifies.
// The default template scalar is double; the BVH derives its T
// from the shape and instantiates the requirement with that.
template<typename S, typename T = double>
concept RayHittable = requires(const Ray<3, T>& ray, const S& s) {
    { ray_hit(ray, s) } -> std::convertible_to<std::optional<RayHit3<T>>>;
};

} // namespace spatium::geometry
