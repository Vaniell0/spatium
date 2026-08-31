#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/geometry/box.hpp>
#  include <spatium/geometry/circle.hpp>
#  include <spatium/geometry/hyperplane.hpp>
#  include <spatium/geometry/line.hpp>
#  include <spatium/geometry/polygon.hpp>
#  include <spatium/geometry/triangle.hpp>
#  include <initializer_list>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::geometry {

// ── Triangle: tri(a, b, c) ─────────────────────────────────────

template<std::size_t N, Scalar T>
constexpr Triangle<N, T> tri(Vec<T, N> a, Vec<T, N> b, Vec<T, N> c) {
    return Triangle<N, T>(a, b, c);
}

// ── Segment: seg(a, b) ─────────────────────────────────────────

template<std::size_t N, Scalar T>
constexpr Segment<N, T> seg(Vec<T, N> a, Vec<T, N> b) {
    return {a, b};
}

// ── Box: box(min, max) or box(center, half) ────────────────────

template<std::size_t N, Scalar T>
constexpr Box<N, T> box(Vec<T, N> min_corner, Vec<T, N> max_corner) {
    return {min_corner, max_corner};
}

// ── Line: line(origin, direction) → Result ─────────────────────

template<std::size_t N, Scalar T>
Result<Line<N, T>> line(Vec<T, N> origin, Vec<T, N> direction) {
    return Line<N, T>::from(origin, direction);
}

// ── Ray: ray(origin, direction) → Result ───────────────────────

template<std::size_t N, Scalar T>
Result<Ray<N, T>> ray(Vec<T, N> origin, Vec<T, N> direction) {
    return Ray<N, T>::from(origin, direction);
}

// ── Plane: plane(normal, point) → Result ───────────────────────

template<std::size_t N, Scalar T>
Result<Hyperplane<N, T>> plane(Vec<T, N> normal, Vec<T, N> point) {
    return Hyperplane<N, T>::from_normal_and_point(normal, point);
}

// From 3 points (3D)
template<Scalar T>
Result<Hyperplane<3, T>> plane(Vec<T, 3> a, Vec<T, 3> b, Vec<T, 3> c) {
    return Hyperplane<3, T>::from_points(a, b, c);
}

// ── Circle: circle(center, radius, normal) ─────────────────────

template<std::size_t N, Scalar T>
Circle<N, T> circle(Vec<T, N> center, T radius, Vec<T, N> normal = {}) {
    return {center, radius, normal};
}

// ── Disk: disk(center, radius, normal) ─────────────────────────

template<std::size_t N, Scalar T>
Disk<N, T> disk(Vec<T, N> center, T radius, Vec<T, N> normal = {}) {
    return {{center, radius, normal}};
}

// ── Polygon: poly({v0, v1, v2, ...}) ───────────────────────────

template<std::size_t N, Scalar T>
Polygon<N, T> poly(std::initializer_list<Vec<T, N>> verts) {
    return {std::vector<Vec<T, N>>(verts)};
}

template<std::size_t N, Scalar T>
Polygon<N, T> poly(std::vector<Vec<T, N>> verts) {
    return {std::move(verts)};
}

// ── operator| for intersection / pipe ─────────────────────────
// a | b  is syntactic sugar for intersect(a, b).
// Operator convention:  | = intersect/pipe,  & = boolean intersect,
//                       + = union,  - = difference.

template<typename A, typename B>
    requires requires(const A& a, const B& b) { intersect(a, b); }
auto operator|(const A& a, const B& b) {
    return intersect(a, b);
}

// Pipe through Result: ray(o,d) | triangle  (no need to dereference)
template<typename A, typename B>
    requires requires(const A& a, const B& b) { intersect(a, b); }
auto operator|(const Result<A>& a, const B& b)
    -> decltype(intersect(std::declval<const A&>(), b))
{
    if (!a) return std::unexpected(a.error());
    return intersect(*a, b);
}

} // namespace spatium::geometry
