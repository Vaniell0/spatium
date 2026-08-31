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
