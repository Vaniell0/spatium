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
#  include <algorithm>
#  include <optional>
#  include <limits>
#  include <span>
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
        // Phrased as "keep only t >= 0" rather than "skip t < 0" on
        // purpose: solve_quadratic divides by 2a with no guard, so a ray
        // whose direction lies in the quadric's null direction (down a
        // cylinder's axis, along a cone's generator) comes back NaN. Every
        // comparison against NaN is false, so the two phrasings disagree
        // exactly there -- the negated form admits it and reports a hit on
        // a surface the ray never touches.
        if (!root.is_real() || !(root.re >= T{0})) continue;
        if (best && !(root.re < best->t)) continue;
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

// ── One shape, placed: the instance leaf ─────────────────────────
//
// A shape that is *somewhere else*. It points at geometry shared with
// every other instance of the same shape and carries the placement that
// puts it there, so a BVH holding these by value is holding placements,
// not meshes.
//
// **Measured on 19 800 specks**, against the flattened form a scene
// builds when every object carries its own triangles:
//
//   leaf                  leaves     build       cast      memory
//   flattened triangles   237 600   319.3 ms   336.8 ms   17.1 MB
//   Instanced<Box>         19 800    19.8 ms   310.4 ms    0.79 MB
//
// 7 740 hits either way. 16.1x on build, 1.08x on traversal, 21.6x on
// memory.
//
// **The leaf must be ONE shape, and that took a wrong turn to learn.**
// The first version held a span of the shape's triangles and scanned them
// linearly, on the reasoning that a dozen triangles beat a tree. It lost
// traversal 2.9x. The comparison was never a 12-triangle scan against a
// 12-triangle tree -- it was a scan against a top-level BVH that had
// already narrowed the ray to *one* triangle. Arithmetic settles it
// without a second measurement: flattened is depth ~18 plus one
// Möller-Trumbore, instanced-with-12 is depth ~15 plus twelve, so three
// slab tests are saved and eleven triangle tests added; at a triangle
// test costing roughly three slabs that predicts 2.4x worse against a
// measured 2.9x. Bounding-box overlap was never involved.
//
// So the leaf holds a shape and does one `ray_hit`. For a speck that
// never rotates, that shape is a `Box` -- and a box is not an
// approximation of a cube, it *is* the cube, tested in six slab
// comparisons rather than twelve Möller-Trumbore.
//
// **The placement is a translation and a uniform scale, deliberately.**
// That is what `VecField::is_placement()` can extract, and it is the case
// where transforming the ray is exact and free: with
// `local = (world - translation) / scale`, the ray parameter `t` is
// unchanged and the normal is unchanged, so a hit needs no correction on
// the way back out. A rotation would join them easily. A *non-uniform*
// scale would not -- normals would want the inverse transpose -- and
// shear or a projective map are a different conversation again. Named
// here rather than discovered by someone whose normals quietly go wrong.
//
// Nothing is wired into a renderer yet.

// A Box is a shape a ray can hit, not only a bound. Six slab comparisons,
// cheaper than a single triangle test -- which is what makes it the right
// leaf for anything whose geometry really is a box.
template<Scalar T>
inline std::optional<RayHit3<T>> ray_hit(const Ray<3, T>& ray,
                                         const Box<3, T>& box) {
    auto span = intersect_parameters(ray, box);
    if (!span) return std::nullopt;
    T t = span->first;
    if (!(t >= T{0})) t = span->second;       // origin inside: take the exit
    if (!(t >= T{0})) return std::nullopt;

    Vec<T, 3> p{ray.origin + ray.direction * t};
    // The face the hit landed on: whichever coordinate sits on its slab.
    Vec<T, 3> n{};
    T closest = std::numeric_limits<T>::max();
    for (std::size_t i = 0; i < 3; ++i) {
        using std::abs;
        T d_lo = abs(p[i] - box.min_corner[i]);
        T d_hi = abs(p[i] - box.max_corner[i]);
        if (d_lo < closest) { closest = d_lo; n = Vec<T, 3>{}; n[i] = T{-1}; }
        if (d_hi < closest) { closest = d_hi; n = Vec<T, 3>{}; n[i] = T{1}; }
    }
    return RayHit3<T>{t, p, n, T{0}, T{0}};
}

template<typename S>
struct Instanced {
    using ScalarType = typename S::ScalarType;
    using PointType = typename S::PointType;
    static constexpr std::size_t ambient_dimension = S::ambient_dimension;

    // A pointer, not a copy. This is the whole point: a thousand instances
    // of one shape cost a thousand placements and one geometry, and the
    // BVH holding them by value is holding placements, not meshes.
    const S* shape = nullptr;
    Vec<ScalarType, 3> translation{};
    ScalarType scale = ScalarType{1};

    PointType centroid() const {
        if (!shape) return PointType{translation};
        return PointType{shape->centroid() * scale + translation};
    }

    Box<3, ScalarType> bounding_box() const {
        if (!shape) return Box<3, ScalarType>{translation, translation};
        auto b = shape->bounding_box();
        return Box<3, ScalarType>{
            Vec<ScalarType, 3>{b.min_corner * scale + translation},
            Vec<ScalarType, 3>{b.max_corner * scale + translation}};
    }
};

// The ray goes into the object's own space; the hit comes back out.
//
// `t` survives the round trip untouched, which is the reason to restrict
// the placement to translation plus uniform scale: dividing both origin
// and direction by the same scale leaves `origin + t * direction` meaning
// the same point, so no hit needs rescaling and no normal needs fixing.
template<typename S>
inline std::optional<RayHit3<typename S::ScalarType>> ray_hit(
    const Ray<3, typename S::ScalarType>& ray, const Instanced<S>& inst) {
    using T = typename S::ScalarType;
    if (!inst.shape || inst.scale == T{0}) return std::nullopt;

    const Ray<3, T> local{Vec<T, 3>{(ray.origin - inst.translation) / inst.scale},
                          Vec<T, 3>{ray.direction / inst.scale}};

    using geometry::ray_hit;
    auto h = ray_hit(local, *inst.shape);
    if (!h) return std::nullopt;
    h->point = Vec<T, 3>{h->point * inst.scale + inst.translation};
    return h;
}

// Concept: any type with a matching `ray_hit` overload qualifies.
// The default template scalar is double; the BVH derives its T
// from the shape and instantiates the requirement with that.
template<typename S, typename T = double>
concept RayHittable = requires(const Ray<3, T>& ray, const S& s) {
    { ray_hit(ray, s) } -> std::convertible_to<std::optional<RayHit3<T>>>;
};

} // namespace spatium::geometry
