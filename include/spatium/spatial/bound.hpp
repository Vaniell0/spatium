#pragma once
// What a node of a bounding-volume hierarchy needs from its bound, stated
// once so one tree serves several: an axis-aligned box, a ball, and later
// a geodesic ball on a manifold. A box and a ball differ only here --
// traversal, the stack, the leaf tests and the build all stay the same.
//
// A bound B provides, found by ADL:
//
//   merge(a, b)                 a bound containing both
//   lower_distance(a, p)        no point inside `a` is closer to `p` than
//                               this; 0 when `p` is inside
//   sah_measure(a)              how likely a random ray is to cross it --
//                               the surface measure, not the volume, by
//                               Cauchy-Crofton
//   ray_interval(ray, a, tmax)  the parameter interval the ray spends in
//                               `a`, if it enters before `tmax`
//   center(a)                   a point to sort by when building
//   overlaps(a, box)            whether it meets an axis-aligned box
//   bound_of<B>(shape)          the bound of one shape
//
// Every one of them may be loose on the safe side and none may be loose on
// the other: a lower distance that is too high, or an interval that misses
// part of the ray, prunes a node that holds the answer.
//
// The box's operations are the functions the tree always called
// (`union_with`, `distance`, `surface_measure`, `intersect_parameters`), so
// a tree over boxes builds and walks exactly as it did before this existed.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/concepts.hpp>
#  include <spatium/geometry/box.hpp>
#  include <spatium/geometry/intersection.hpp>
#  include <spatium/geometry/line.hpp>
#  include <cmath>
#  include <concepts>
#  include <limits>
#  include <numbers>
#  include <optional>
#  include <utility>
#endif

SPATIUM_EXPORT namespace spatium::spatial {

template<typename B>
concept Bound = requires(const B& a, const B& b,
                         const Vec<typename B::ScalarType, B::ambient_dimension>& p,
                         const geometry::Ray<B::ambient_dimension, typename B::ScalarType>& ray,
                         const geometry::Box<B::ambient_dimension, typename B::ScalarType>& box,
                         typename B::ScalarType tmax) {
    { merge(a, b) } -> std::same_as<B>;
    { lower_distance(a, p) } -> std::same_as<typename B::ScalarType>;
    { sah_measure(a) } -> std::same_as<typename B::ScalarType>;
    { ray_interval(ray, a, tmax) }
        -> std::same_as<std::optional<std::pair<typename B::ScalarType, typename B::ScalarType>>>;
    { center(a) } -> std::same_as<Vec<typename B::ScalarType, B::ambient_dimension>>;
    { overlaps(a, box) } -> std::same_as<bool>;
};

// The bound of one shape. A box is what every Bounded shape already gives;
// another bound is built from it unless the shape offers something tighter.
template<typename B, typename Shape>
    requires std::same_as<B, geometry::Box<Shape::ambient_dimension, typename Shape::ScalarType>>
B bound_of(const Shape& s) {
    return s.bounding_box();
}

} // namespace spatium::spatial

// The box's side of the concept, beside Box itself so ADL finds it.
SPATIUM_EXPORT namespace spatium::geometry {

template<std::size_t N, Scalar T>
Box<N, T> merge(const Box<N, T>& a, const Box<N, T>& b) { return a.union_with(b); }

template<std::size_t N, Scalar T>
T lower_distance(const Box<N, T>& a, const Vec<T, N>& p) { return a.distance(p); }

template<std::size_t N, Scalar T>
T sah_measure(const Box<N, T>& a) { return a.surface_measure(); }

template<std::size_t N, Scalar T>
std::optional<std::pair<T, T>> ray_interval(const Ray<N, T>& ray, const Box<N, T>& a, T tmax) {
    auto hit = intersect_parameters(ray, a);
    if (!hit || hit->first > tmax) return std::nullopt;
    return *hit;
}

template<std::size_t N, Scalar T>
Vec<T, N> center(const Box<N, T>& a) { return a.centroid(); }

template<std::size_t N, Scalar T>
bool overlaps(const Box<N, T>& a, const Box<N, T>& box) { return a.intersects(box); }

} // namespace spatium::geometry

SPATIUM_EXPORT namespace spatium::spatial {

// A ball in the ambient space: the bound a point cloud, a swept point or
// the cells of a chart search want, where a box would be loose on every
// diagonal. It is also the Euclidean case of the geodesic ball a tree on a
// manifold is built from, which is why it is here rather than in geometry/.
template<std::size_t N, Scalar T = double>
struct Ball {
    using ScalarType = T;
    using PointType = Vec<T, N>;
    static constexpr std::size_t ambient_dimension = N;

    PointType c{};
    T r{};
};

// Rounding in the arithmetic below is a few ulps of the coordinates; the
// pad keeps a merged ball containing both of its children, which the
// pruning relies on.
namespace detail {
template<std::size_t N, Scalar T>
T ball_pad(const Ball<N, T>& b) {
    using std::abs;
    T scale = b.r;
    for (std::size_t i = 0; i < N; ++i) scale = std::max(scale, abs(b.c[i]));
    return scale * std::numeric_limits<T>::epsilon() * T{8};
}
}  // namespace detail

template<std::size_t N, Scalar T>
Ball<N, T> merge(const Ball<N, T>& a, const Ball<N, T>& b) {
    const Vec<T, N> ab{b.c - a.c};
    const T d = ab.norm();
    if (d + b.r <= a.r) return a;
    if (d + a.r <= b.r) return b;
    const T r = (d + a.r + b.r) / T{2};
    Ball<N, T> out{Vec<T, N>{a.c + ab * ((r - a.r) / d)}, r};
    out.r += detail::ball_pad(out);
    return out;
}

template<std::size_t N, Scalar T>
T lower_distance(const Ball<N, T>& a, const Vec<T, N>& p) {
    const T d = Vec<T, N>{p - a.c}.norm() - a.r;
    return d > T{0} ? d : T{0};
}

// The surface measure of the sphere bounding it -- 2 pi r in the plane,
// 4 pi r^2 in space -- which is what a ray's chance of crossing it scales
// with. Only ratios of it are ever used, within one tree.
template<std::size_t N, Scalar T>
T sah_measure(const Ball<N, T>& a) {
    constexpr T pi = std::numbers::pi_v<T>;
    if constexpr (N == 2) return T{2} * pi * a.r;
    else if constexpr (N == 3) return T{4} * pi * a.r * a.r;
    else {
        T m = T{1};
        for (std::size_t i = 1; i < N; ++i) m *= a.r;
        return m;
    }
}

template<std::size_t N, Scalar T>
std::optional<std::pair<T, T>> ray_interval(const geometry::Ray<N, T>& ray, const Ball<N, T>& a, T tmax) {
    using std::sqrt;
    const Vec<T, N> oc{ray.origin - a.c};
    const T dd = ray.direction.dot(ray.direction);
    const T b = ray.direction.dot(oc);
    const T cc = oc.dot(oc) - a.r * a.r;
    const T disc = b * b - dd * cc;
    if (disc < T{0}) return std::nullopt;
    const T s = sqrt(disc);
    const T t1 = (-b + s) / dd;
    if (t1 < T{0}) return std::nullopt;
    T t0 = (-b - s) / dd;
    if (t0 < T{0}) t0 = T{0};
    if (t0 > tmax) return std::nullopt;
    return std::pair{t0, t1};
}

template<std::size_t N, Scalar T>
Vec<T, N> center(const Ball<N, T>& a) { return a.c; }

template<std::size_t N, Scalar T>
bool overlaps(const Ball<N, T>& a, const geometry::Box<N, T>& box) {
    return box.distance(a.c) <= a.r;
}

// A shape's ball, circumscribing its box unless the shape gives its own.
template<typename B, typename Shape>
    requires std::same_as<B, Ball<Shape::ambient_dimension, typename Shape::ScalarType>>
B bound_of(const Shape& s) {
    using T = typename Shape::ScalarType;
    if constexpr (requires { { s.bounding_ball() } -> std::same_as<B>; }) {
        return s.bounding_ball();
    } else {
        const auto box = s.bounding_box();
        B out{box.centroid(), box.extents().norm() / T{2}};
        out.r += detail::ball_pad(out);
        return out;
    }
}

} // namespace spatium::spatial
