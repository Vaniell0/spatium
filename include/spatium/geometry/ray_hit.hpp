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
// A shape that is *somewhere else*. It does not own its geometry -- it
// points at triangles shared with every other instance of the same shape
// -- and carries the placement that puts it there.
//
// This is what lets a BVH hold one leaf per object instead of one leaf
// per triangle. The donut's dust is 19 800 copies of the same 12-triangle
// cube: 396 000 triangle leaves today, 19 800 instance leaves with this.
// Measured on that scene, the dust is 56% of the frame and scales
// linearly with the particle count, and leaf count is what both the build
// and the traversal depth follow.
//
// **The placement is a translation and a uniform scale, deliberately.**
// That is what `VecField::is_placement()` can extract, and it is the case
// where transforming the ray is exact and free: with
// `local = (world - translation) / scale`, the ray parameter `t` is
// unchanged and the surface normal is unchanged, so a hit needs no
// correction on the way back out. A rotation would join them easily. A
// *non-uniform* scale would not: normals would need the inverse transpose
// and `t` would need rescaling, and shear or a projective map are a
// different conversation again. The boundary is named here rather than
// discovered by someone whose normals quietly go wrong.
//
// **Measured, and it does not pay yet -- nothing uses this.** On 19 800
// specks of 12 triangles each, against the flattened 237 600-triangle
// form the demo builds today:
//
//   BVH build   304.2 ms -> 20.1 ms   (15x faster, 237 600 -> 19 800 leaves)
//   geometry      6.5 ms -> 0.2 ms
//   200 000 rays  330.6 ms -> 959.6 ms  (2.9x SLOWER; 7 740 hits either way)
//
// The traversal number is the one that decides it, and the cause is this
// comment's original claim, which was wrong: "at a dozen triangles a
// linear scan beats a tree". The comparison is not a 12-triangle scan
// against a 12-triangle tree. It is a 12-triangle scan against a top-level
// BVH that had already narrowed the ray down to *one* triangle. Putting a
// linear scan at the leaf reintroduces exactly what the tree exists to
// remove, and 15x on build does not buy back 2.9x on traversal when a
// frame is a million rays.
//
// **The cause, settled by arithmetic and then confirmed.** Flattened:
// 237 600 leaves, depth ~18, one triangle per leaf -- 18 AABB tests plus
// one Möller-Trumbore. Instanced: 19 800 leaves, depth ~15, twelve
// triangles per leaf -- 15 AABB plus twelve Möller-Trumbore. Three AABB
// tests saved, eleven triangle tests added; at a triangle test costing
// roughly three slab tests that predicts ~2.4x worse, and the measurement
// said 2.9x. The overlap of instance bounding boxes was never the
// problem.
//
// So the fix is not a nested acceleration structure -- that would be
// 19 800 subtrees to build. **The leaf has to be one shape.** For a speck
// that never rotates, that shape is a `Box`: axis-aligned in local space,
// so the slab test is exact rather than an approximation, and it is
// cheaper than a single triangle test. Measured on the same scene:
//
//   leaf                  leaves     build       cast
//   flattened triangles   237 600   314.7 ms   318.0 ms
//   instance of 12 tris    19 800    19.9 ms   982.9 ms
//   instance of one box    19 800    19.8 ms   293.5 ms
//
// 7 740 hits in all three. One-shape leaves win on both axes: 15.9x on
// build and 8% on traversal.
//
// This type stays as written -- an instance over a span of triangles is
// the general case, and it is the right general case for a shape that
// really is a mesh. What the numbers say is that the *dust* should not be
// a mesh at all; a cube is a Box. Nothing is wired into a renderer yet.
template<Scalar T>
struct Instanced {
    using ScalarType = T;
    using PointType = Vec<T, 3>;
    static constexpr std::size_t ambient_dimension = 3;

    std::span<const Triangle<3, T>> triangles{};  // shared, not owned
    Vec<T, 3> translation{};
    T scale = T{1};

    PointType centroid() const {
        Vec<T, 3> c{};
        if (triangles.empty()) return PointType{translation};
        for (const auto& tri : triangles) c = Vec<T, 3>{c + tri.centroid()};
        return PointType{c * (scale / static_cast<T>(triangles.size())) + translation};
    }

    Box<3, T> bounding_box() const {
        if (triangles.empty()) return Box<3, T>{translation, translation};
        auto b = triangles.front().bounding_box();
        Vec<T, 3> lo = b.min_corner, hi = b.max_corner;
        for (const auto& tri : triangles) {
            auto tb = tri.bounding_box();
            for (std::size_t i = 0; i < 3; ++i) {
                lo[i] = std::min(lo[i], tb.min_corner[i]);
                hi[i] = std::max(hi[i], tb.max_corner[i]);
            }
        }
        return Box<3, T>{Vec<T, 3>{lo * scale + translation},
                         Vec<T, 3>{hi * scale + translation}};
    }
};

// The ray goes into the object's own space; the hit comes back out.
//
// `t` survives the round trip untouched, which is the reason to restrict
// the placement to translation plus uniform scale: dividing both origin
// and direction by the same scale leaves `origin + t * direction` meaning
// the same point, so no hit needs rescaling and no normal needs fixing.
template<Scalar T>
inline std::optional<RayHit3<T>> ray_hit(const Ray<3, T>& ray,
                                         const Instanced<T>& inst) {
    if (inst.triangles.empty() || inst.scale == T{0}) return std::nullopt;

    const Ray<3, T> local{Vec<T, 3>{(ray.origin - inst.translation) / inst.scale},
                          Vec<T, 3>{ray.direction / inst.scale}};

    std::optional<RayHit3<T>> best;
    for (const auto& tri : inst.triangles) {
        auto h = ray_hit(local, tri);
        if (!h) continue;
        if (best && !(h->t < best->t)) continue;
        best = h;
    }
    if (!best) return std::nullopt;

    best->point = Vec<T, 3>{best->point * inst.scale + inst.translation};
    return best;
}

// Concept: any type with a matching `ray_hit` overload qualifies.
// The default template scalar is double; the BVH derives its T
// from the shape and instantiates the requirement with that.
template<typename S, typename T = double>
concept RayHittable = requires(const Ray<3, T>& ray, const S& s) {
    { ray_hit(ray, s) } -> std::convertible_to<std::optional<RayHit3<T>>>;
};

} // namespace spatium::geometry
