#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/io/scene.hpp>
#  include <spatium/mesh/mesh.hpp>
#  include <spatium/mesh/primitives.hpp>
#  include <spatium/spaces/chart.hpp>
#  include <spatium/spaces/euclidean.hpp>
#  include <spatium/spaces/offset.hpp>
#  include <spatium/spaces/parametric.hpp>
#  include <spatium/spaces/sample.hpp>
#  include <spatium/algebra/quaternion.hpp>
#  include <spatium/geometry/ray_surface.hpp>
#  include <spatium/io/field.hpp>
#  include <any>
#  include <cstdint>
#  include <functional>
#  include <initializer_list>
#  include <optional>
#  include <stdexcept>
#  include <unordered_map>
#  include <string>
#  include <utility>
#  include <vector>
#endif

// Declarative scene builder, space-first: a Trace is a flat, inspectable
// record of operations (Space / Offset / Scatter / Compose) over real
// spatium spaces (ParametricSurface -- torus, cylinder, any u,v->R^3
// formula), addressed by index rather than closures. This IS the trace:
// walk trace.node(i) for any i, dump it, or write a second walker over
// the same data later (e.g. codegen) without touching this one.
//
// Analytic until the last mile: Space and Offset nodes resolve to real
// ParametricSurface function compositions (offset_surface()) -- no mesh
// exists anywhere until materialize() tessellates, purely for display.
// Scatter placement is analytic too (spaces/sample.hpp's area-weighted
// rejection sampling over the first fundamental form) -- no intermediate
// mesh/Voronoi graph.
//
// Every node carries one motion/mutation slot: `Vec3(Vec3 point, T t) ->
// Vec3`, defaulting to identity. A constant shift, a Morphism pipe
// (`[m](p,T){ return m(p); }`), or genuine time-dependence all fit the
// same slot -- not three different mechanisms.

SPATIUM_EXPORT namespace spatium::io::build {

// ── Field slots ───────────────────────────────────────────────────
//
// A node's motion and thickness are `Field`s -- expressions whose leaves
// may be opaque callables -- not bare erased functions. See io/field.hpp
// for why the leaf rather than the whole field is the unit of opacity.
//
// The motion slot used to be a `std::move_only_function`, chosen so a
// hook could own move-only state and so the signature could be const.
// Both of those were about capability, and the second one still holds
// (the field's call operator is const). The first turned out to aim at a
// problem nobody had: nothing in the tree owns move-only state, while the
// problem that *was* measured -- 19 800 closures each copying a 512-byte
// noise table -- is fixed by sharing the table, which the demo already
// does. Meanwhile move-only cost something real: it made a `TraceNode`
// uncopyable, so a `Trace` could not be copied, so `cook()` could only
// consume a trace rather than copy it, and `std::move_only_function`
// exposes no `target_type()`, so a motion could not even be counted by
// type. A field that owns its pool gives back all three.
//
// A lambda of the historical `(point, t)` shape still converts
// implicitly, so every scene reads unchanged. What it gains is that the
// conversion is now visible: an identity motion is a `Point` node and a
// constant shift is `Point + Const`, both structural and both reported as
// such, where they used to be opaque one-line lambdas indistinguishable
// from a noise field.
template<Scalar T = double>
using PointField = VecField<T>;

template<Scalar T = double>
using ScalarField = Field<T>;

enum class Kind { Space, Offset, Scatter, Compose, Literal };

inline const char* kind_name(Kind k) {
    switch (k) {
        case Kind::Space:   return "Space";
        case Kind::Offset:  return "Offset";
        case Kind::Scatter: return "Scatter";
        case Kind::Compose: return "Compose";
        case Kind::Literal: return "Literal";
    }
    return "?";
}

// What happens at the rim when a shell is built over a surface that has
// one. Only a shell needs this: offsetting a closed surface produces a
// closed surface with no rim to rule on, which is why `offset()` and
// `offset_shell()` are separate operations rather than one with a flag.
//
// `ZeroThickness` is the rim the donut's icing already had, now named
// rather than improvised: the thickness field is expected to fall to
// zero before the edge, so the shell meets its base there and closes
// against it. It is the only rule that carries no data, which is exactly
// why it is the only one here today -- a round cap carries a radius, and
// extending to another surface carries a reference to another node.
// Those want to be fields, and fields are the next piece of work.
enum class EdgeRule { ZeroThickness };

// Which of a scattered item's own axes points along the target's normal.
//
// `Z` is the default and the only behaviour there used to be: an item is
// modelled standing up, and scattering stands it on the surface. That is
// right for a bristle, a blade of grass, a tree -- anything whose
// interesting direction *is* the normal.
//
// It is wrong for anything that lies *across* a surface, and the donut
// demo is the proof: a sprinkle is a small cylinder lying down, and with
// only `Z` available the only way to express that was to permute the
// item's mesh vertices by hand -- which a closed form cannot follow,
// since `cylinder()`'s exact `BoundedQuadric` is a cylinder about z and
// stays one. So the sprinkle had to be a `Literal`, and a `Literal` has
// no closed form, and an object with no closed form cannot be instanced.
// One missing enum cost 11 000 objects their instancing.
//
// The frame's columns are cyclically permuted rather than swapped, which
// keeps the basis right-handed: a cyclic permutation leaves a
// determinant alone, a swap negates it, and a left-handed frame would
// mirror every item placed through it.
enum class SeatAxis { X, Y, Z };

// How a renderer should turn this node into ray hits. Three levels that
// answer the same question at three prices, measured on the shape
// donut_demo actually builds (benchmarks/bench_raycast.cpp):
//
//   Exact         closed form -- the real surface, no approximation.
//                 A torus is 20.9 ns/ray as one BVH leaf.
//   Tessellated   triangles in a BVH. 40.8 ns/ray, plus 10.8 ms per
//                 frame to build the 25 600 of them first.
//   Newton        Newton's method on the (u,v) map. 23 612 ns/ray --
//                 three orders of magnitude, and the only way to hit a
//                 general parametric surface exactly.
//
// The gap is why this is a choice rather than an implementation detail.
// But note what the ratio does and does not say: the exact path's
// advantage grows with the scene, because it replaces *leaves*, not
// because a leaf test is faster. One torus against 25 600 triangles is
// 2x; a hundred tori against 2.5M triangles is a difference in traversal
// depth. Reading 2x as the ceiling is the mistake to avoid.
//
// Tessellated is the default even for a node that has a (u,v) map,
// because Newton is three orders of magnitude dearer and nobody should
// pay that by accident. Exactness on a general parametric surface is
// something you ask for.
enum class RenderLevel { Exact, Tessellated, Newton };

inline const char* render_level_name(RenderLevel l) {
    switch (l) {
        case RenderLevel::Exact:       return "Exact";
        case RenderLevel::Tessellated: return "Tessellated";
        case RenderLevel::Newton:      return "Newton";
    }
    return "?";
}

template<Scalar T = double>
struct TraceNode {
    Kind kind{};

    // Space
    std::optional<ParametricSurface<T>> surface;
    std::size_t u_steps = 48, v_steps = 24;

    // The exact analytic form, when this node has one -- a Torus, a
    // BoundedQuadric, or any user type that is Bounded and RayHittable.
    // A renderer that can consume it never tessellates this node at all;
    // one that cannot ignores the slot entirely.
    //
    // std::any rather than a closed variant, so adding a shape kind does
    // not mean editing this header. And deliberately NOT the shape of
    // io/scene.hpp's ResolvedShape, which erases to a std::function
    // returning hits: that buys the same openness at the price of an
    // indirect call at every leaf test. The renderer needs the concrete
    // type *back* in order to bucket it into a monomorphic BVH, and
    // std::any keeps the type tag that makes that possible. Erasure
    // belongs to describing a scene; evaluating one stays monomorphic.
    //
    // Known constraint, same family as ScalarField's above: std::any
    // requires the stored type to be copy-constructible, so a shape that
    // owns move-only state (a device handle, a voxel grid behind a
    // unique_ptr) cannot go here. Shapes are value types today so
    // nothing hits it, but the failure would be a std::any error message
    // that says nothing about this design. Tracked in ROADMAP with the
    // other copy-constructibility leaks rather than left as a comment.
    //
    // This is an invariant, not a field. It must agree with `surface`,
    // and any operation that can move them apart is required to clear it
    // -- see Handle::moving(). A stale exact form renders a picture that
    // is correct for the shape and wrong for the scene, silently, which
    // is the same class of defect as a NaN walking through a filter.
    std::any exact;

    // Unset means "infer from what this node is" -- see
    // Placed::render_level(). Set by .rendered_as(), which refuses on the
    // spot a level this node cannot actually serve.
    std::optional<RenderLevel> level;

    // Literal -- a precomputed mesh, for shapes with no natural single
    // (u,v)->R^3 formula (a cube, most obviously). The escape hatch out
    // of the analytic-first path, not the default.
    std::optional<mesh::Mesh<Euclidean<3, T>>> literal_mesh;

    // Offset -- `edge` is meaningful only for a shell, i.e. when the
    // base is open. A closed base has no edge and never consults it.
    std::size_t base = 0;
    ScalarField<T> thickness;
    EdgeRule edge = EdgeRule::ZeroThickness;

    // Scatter
    std::size_t item = 0, target = 0;
    std::size_t count = 0;
    std::uint32_t seed = 42;
    SeatAxis seat_axis = SeatAxis::Z;

    // How deep a scattered item sits. 1 rests it on the surface (its
    // lowest point touching), 0 puts its own origin there -- half sunk,
    // for an item modelled around its centre. Anything between is a
    // press into a soft surface.
    //
    // It has to be a named knob rather than an offset baked into the
    // item's mesh, and the reason is not style. `scatter_lift()` derives
    // the rise from the item's own lowest point, so shifting every vertex
    // down by d lowers min_z by d and raises the lift by d: the two
    // cancel exactly and the nudge does nothing at all. The donut demo
    // had such a nudge, with a comment explaining what it was for, and it
    // had silently stopped doing anything.
    T seat = T{1};

    // Compose
    std::vector<std::size_t> children;

    Material<T> material{};
    PointField<T> transform;

    // Optional: color as a function of (a representative point, t) --
    // evaluated once per node per materialize() call (at the
    // materialized mesh's own centroid), same "evaluate at t" spirit as
    // `transform`, just producing a color instead of a position. When
    // set, overrides `material.base_color` for that materialize() call.
    PointField<T> color_fn;

    // Emission, same shape and same reason as color_fn: a field, because
    // a node that glows only sometimes cannot say so with a constant.
    // Unset means the identity VecField, which is falsy -- so "is there
    // one" is a question about the expression rather than a null check.
    PointField<T> emissive_fn;
};

template<Scalar T>
class Trace;

// One object of a materialized scene.
//
// Deliberately NOT a mesh. The whole argument for describing a scene as
// spaces rather than triangles is that every operation in the library --
// geodesics, contact, physics, Riemannian optimisation -- works on a
// scene object for free, because the object *is* a space. Handing back
// `{mesh, material}` threw that away one step before the renderer: the
// trace stayed analytic right up to materialize() and then became
// triangles, which left "why not just load an OBJ" a fair question.
//
// So this is a view onto the node instead: it says what the object is,
// and derives the rest on request. `surface()` is present whenever the
// node really is a space, with the node's own motion composed into the
// map, so a caller that can consume an exact surface never pays for a
// tessellation. `mesh()` builds the triangle view for callers that need
// one.
//
// Lifetime follows `Handle`: this points into the trace and is valid as
// long as the trace is. `mesh()` builds on every call rather than
// caching, so bind it once (`auto m = obj.mesh();`) if you need it more
// than once -- explicit, rather than a hidden cost inside an innocent
// looking member access.
template<Scalar T = double>
struct Placed {
    const Trace<T>* trace;
    std::size_t index;
    T t;

    // The triangle view, built on demand.
    mesh::Mesh<Euclidean<3, T>> mesh() const;

    // The shape before this object's own motion moved it. Two copies of
    // one speck differ only in where they were put, so this is the form
    // they can share; `mesh()` is this with the placement already applied
    // and is therefore unique per object by construction.
    mesh::Mesh<Euclidean<3, T>> rest_mesh() const {
        return materialize_mesh(*trace, index, t, /*placed=*/false);
    }

    // The transform this object's motion is, if it is one. Present when
    // the motion is affine in the point -- a placement, which moves the
    // object without reshaping it -- and absent when it deforms per
    // vertex, because then there is no transform to name and nothing to
    // share. A renderer uses this to decide whether an object can be an
    // instance of a shared shape or has to be its own triangles.
    std::optional<typename VecField<T>::Placement> placement() const {
        const auto& n = trace->node(index);
        if (!n.transform.is_placement()) return std::nullopt;
        return n.transform.placement_at(MotionEnv<T>{Vec<T, 3>{}, t});
    }

    // The analytic description, when there is one. Present for Space and
    // Offset nodes; absent for Literal (a precomputed mesh with no
    // (u,v) map), Scatter and Compose (many objects, not one surface).
    std::optional<ParametricSurface<T>> surface() const;

    bool is_analytic() const;

    // Whether this object has an exact closed form a renderer can hit
    // directly, rather than only a (u,v) map to tessellate or run Newton
    // against. Visible on purpose: losing the exact form is a silent
    // downgrade from a closed-form hit to a tessellation, and a caller
    // that cares should be able to see it rather than infer it from a
    // frame time. `exact_type()` says which shape it is, for a renderer
    // bucketing nodes into one monomorphic tree per type.
    bool is_exact() const;
    const std::type_info& exact_type() const;

    // The level a renderer should use for this object: whatever
    // .rendered_as() asked for, or the inference below when nothing
    // asked. Resolved rather than raw, so a renderer never has to
    // reimplement the default and drift from it.
    RenderLevel render_level() const;

    // The exact form itself, when it is the type asked for. Returns null
    // when the node has no exact form or has a different one, so the
    // bucketing loop is a cast attempt rather than a type interrogation
    // followed by a cast.
    template<class Shape>
    const Shape* exact_as() const;

    // Resolved material. A node's color_fn needs a representative point,
    // which today means the mesh centroid, so a node that has one pays
    // for its mesh here; a node with a plain Material does not.
    Material<T> material() const;
};

// Orthonormal basis {t1, t2} spanning the plane perpendicular to a unit
// normal n (Duff, Burgess, Christensen, Hery, Kensler, Liani, Villemin
// 2017) -- orients a scattered item so its local z maps to the target's
// normal at the site it was placed.
template<Scalar T>
std::pair<Vec<T, 3>, Vec<T, 3>> basis_from_normal(const Vec<T, 3>& n) {
    T sign = n[2] >= T{0} ? T{1} : T{-1};
    T a = T{-1} / (sign + n[2]);
    T b = n[0] * n[1] * a;
    Vec<T, 3> t1{T{1} + sign * n[0] * n[0] * a, sign * b, -sign * n[0]};
    Vec<T, 3> t2{b, sign + n[1] * n[1] * a, -n[1]};
    return {t1, t2};
}

// The full orientation of one scattered item: the site's tangent frame,
// turned in-plane by an angle that is a pure function of (seed, index).
//
// Without the in-plane turn every item at every site shares the same
// {t1, t2} up to the smooth drift `basis_from_normal` has with the
// surface's own curvature, which reads as "all facing the same way" -- a
// real complaint about the donut's sprinkles, and the reason this angle
// exists at all.
//
// It is a named function rather than a loop body because there are two
// consumers that must agree to the last bit: `materialize_mesh()` bakes
// this frame into vertices, and `cook()` hands it to an instance as a
// transform. They are two renderings of the same object and any drift
// between them is a disagreement no test would phrase as such -- it would
// surface as a picture that changes when a node's render level changes.
// Deterministic per (seed, index): same trace, same t, same result.
//
// The columns are {r1, r2, normal}, so `frame * v` maps an item's local
// (x, y, z) onto (in-plane, in-plane, along the normal) exactly as
// scattering means it.
// How far a scattered item has to rise along the site normal so that it
// rests ON the surface instead of halfway inside it: the depth of its own
// geometry below its local origin.
//
// This was never decided; it was inherited from a bug. `offset()` used to
// ignore its base's tessellation and render at 48x24, and a coarse mesh
// of a convex surface sits *below* the analytic surface it approximates,
// which lifted every item placed on that analytic surface into view by
// accident. Making offset() honour the request removed the accident and
// the sprinkles sank into the icing -- the same picture the ROADMAP entry
// about items "lying fully flush, not slightly proud" describes from the
// other side. Either way the answer is to place them deliberately.
//
// Computed from the item's own rest geometry rather than taken as a
// parameter: the caller already said how big the item is by building it,
// and asking twice is how the two answers drift apart.
template<Scalar T>
T scatter_lift(const mesh::Mesh<Euclidean<3, T>>& item, SeatAxis axis = SeatAxis::Z) {
    // Measured along whichever axis this item stands on -- the same one
    // `scatter_frame` sends to the normal. Measuring z regardless would
    // seat a lying cylinder by its length instead of its radius.
    const std::size_t a = axis == SeatAxis::X ? 0 : axis == SeatAxis::Y ? 1 : 2;
    T lowest{};
    for (const auto& v : item.vertices) lowest = std::min(lowest, v[a]);
    return -lowest;
}

template<Scalar T>
Matrix<T, 3, 3> scatter_frame(const Vec<T, 3>& normal, std::uint32_t seed,
                              std::size_t index, SeatAxis axis = SeatAxis::Z) {
    auto [t1, t2] = basis_from_normal(normal);

    std::uint32_t h = seed * 2654435761u + static_cast<std::uint32_t>(index) * 40503u;
    h ^= h >> 13; h *= 0x85ebca6bu; h ^= h >> 16;
    T theta = (static_cast<T>(h % 360000) / T{100000}) * T{6.283185307179586};

    using std::cos, std::sin;
    Vec<T, 3> r1{Vec<T, 3>{t1 * cos(theta) + t2 * sin(theta)}};
    Vec<T, 3> r2{Vec<T, 3>{t2 * cos(theta) - t1 * sin(theta)}};

    // Cyclic, not a swap. [r1 r2 n] is right-handed, and a cyclic
    // permutation of columns leaves the determinant alone while a swap
    // negates it -- a left-handed frame would mirror every item placed
    // through it, which is the kind of thing that looks like a modelling
    // mistake rather than a sign error.
    Vec<T, 3> c0 = r1, c1 = r2, c2 = normal;
    if (axis == SeatAxis::X) { c0 = normal; c1 = r1; c2 = r2; }   // local x -> n
    else if (axis == SeatAxis::Y) { c0 = r2; c1 = normal; c2 = r1; }   // local y -> n

    Matrix<T, 3, 3> f{};
    for (std::size_t i = 0; i < 3; ++i) {
        f(i, 0) = c0[i];
        f(i, 1) = c1[i];
        f(i, 2) = c2[i];
    }
    return f;
}

template<Scalar T = double>
class Trace;

template<Scalar T = double>
struct Handle {
    Trace<T>* trace;
    std::size_t index;

    Handle colored(Material<T> m) const;
    Handle rendered_as(RenderLevel level) const;
    Handle colored(PointField<T> color_fn) const;

    // Emission, as a field for the same reason colour is one: the donut's
    // dust has to *catch* light as it nears its letterform and lose it
    // again as it drifts off, and a constant on the node could only say
    // "always" or "never".
    Handle glowing(PointField<T> emissive_fn) const;
    Handle glowing(Vec<T, 3> emissive) const;
    Handle moving(PointField<T> f) const;
};

template<Scalar T>
ParametricSurface<T> resolve_surface(const Trace<T>& trace, std::size_t idx);

template<Scalar T>
mesh::Mesh<Euclidean<3, T>> materialize_mesh(const Trace<T>& trace, std::size_t idx, T t = T{0}, bool placed = true);

template<Scalar T>
std::vector<Placed<T>> materialize(const Trace<T>& trace, std::size_t idx, T t = T{0});

template<Scalar T>
class Trace {
public:
    const TraceNode<T>& node(std::size_t i) const { return nodes_[i]; }
    TraceNode<T>& node(std::size_t i) { return nodes_[i]; }
    std::size_t size() const { return nodes_.size(); }

    Handle<T> space(ParametricSurface<T> surface, std::size_t u_steps = 48, std::size_t v_steps = 24) {
        TraceNode<T> n{};
        n.kind = Kind::Space;
        n.surface = std::move(surface);
        n.u_steps = u_steps;
        n.v_steps = v_steps;
        return push(std::move(n));
    }

    // The same door, for a space that is not itself a chart but has one.
    // `chart_of` is found by ADL, so a space defined in a user's own
    // namespace enters the DSL by having an overload written beside it --
    // no edit to this header, no registry, no indirect call: the chart is
    // built once here, at trace-build time, and the node holds the result.
    //
    // Note what this does *not* open. `Kind` stays closed, and adding a
    // space does not go near it: Space/Offset/Scatter/Compose/Literal are
    // operations, not shapes, and a sphere is a new chart on the existing
    // Space node rather than a new kind of node. See ROADMAP for the one
    // thing that would open `Kind`, and why it is not this.
    template<typename S>
        requires Chartable<S, T> && (!std::same_as<std::remove_cvref_t<S>, ParametricSurface<T>>)
    Handle<T> space(const S& s, std::size_t u_steps = 48, std::size_t v_steps = 24) {
        return space(chart_of(s), u_steps, v_steps);
    }

    // The factories that know their shape exactly record it alongside the
    // (u,v) map, so a renderer can hit the real surface instead of a
    // tessellation of it. space() deliberately does not: an arbitrary
    // (u,v)->R^3 formula has no closed form to record.
    Handle<T> torus(T major_r, T minor_r, std::size_t u_steps = 48, std::size_t v_steps = 24) {
        auto h = space(make_torus<T>(major_r, minor_r), u_steps, v_steps);
        node(h.index).exact =
            geometry::Torus<T>{.major_radius = major_r, .minor_radius = minor_r};
        return h;
    }

    Handle<T> cylinder(T radius, T height, std::size_t u_steps = 12, std::size_t v_steps = 4) {
        auto h = space(make_cylinder<T>(radius, height), u_steps, v_steps);
        // make_cylinder puts v in [0, height], so the clip matches the map.
        node(h.index).exact =
            geometry::BoundedQuadric<T>::cylinder_z(radius, T{0}, height);
        return h;
    }

    // A sphere, entering through its own chart (chart_of(Sphere<2,T>))
    // and recording the exact form a renderer can hit directly.
    //
    // Worth having as a factory rather than as `space(chart_of(...))`,
    // for the same reason `torus()` and `cylinder()` are: `space()` cannot
    // know what it was handed. A chart is a (u,v) map and nothing more,
    // so a sphere that arrives that way is a surface with no closed form
    // recorded, and every renderer downstream has to tessellate it. The
    // factory knows, so it says.
    Handle<T> sphere(T radius, std::size_t u_steps = 24, std::size_t v_steps = 12) {
        auto h = space(chart_of(Sphere<2, T>{radius}), u_steps, v_steps);
        node(h.index).exact = geometry::BoundedQuadric<T>::sphere(radius);
        return h;
    }

    // A flake: a shallow spherical cap, clipped to a thin slab. What a
    // speck of dust or ash actually is -- a curved sheet, not a ball --
    // and the shape to reach for when there will be millions of them,
    // because the clip box *is* the speck's size and the box is what a
    // tree traverses.
    //
    // The chart and the exact form are built from the same three numbers
    // here, which is the point of it being a factory: they cannot drift.
    // `half` is the flake's own extent; `bulge` is how much wider the
    // sphere is than the slab it is cut by, and so how curved the sheet
    // is. Keep it near 1 — a sphere much larger than its clip is a nearly
    // flat patch wearing a huge bounding box.
    Handle<T> flake(const Vec<T, 3>& half, T bulge = T{1.35},
                    std::size_t u_steps = 8, std::size_t v_steps = 3) {
        using std::acos, std::sqrt, std::sin, std::cos, std::min, std::max;
        const T r = bulge * max(half[0], max(half[1], half[2]));
        const Vec<T, 3> centre{T{0}, T{0}, half[2] - r};
        // How far down the sphere the slab still admits: the polar angle
        // at which the cap's radius reaches the slab's half-width.
        const T rim = min(half[0], half[1]);
        const T v_max = acos(max(T{-1}, min(T{1}, sqrt(max(T{0}, r * r - rim * rim)) / r)));

        auto h = space(ParametricSurface<T>(
                           [r, centre](T u, T v) -> Vec<T, 3> {
                               return {centre[0] + r * sin(v) * cos(u),
                                       centre[1] + r * sin(v) * sin(u),
                                       centre[2] + r * cos(v)};
                           },
                           {T{0}, T{2} * acos(T{-1}), T{0}, v_max},
                           /*periodic_u=*/true, /*periodic_v=*/false),
                       u_steps, v_steps);
        node(h.index).exact = geometry::BoundedQuadric<T>::flake(half, bulge);
        return h;
    }

    Handle<T> literal(mesh::Mesh<Euclidean<3, T>> m) {
        TraceNode<T> n{};
        n.kind = Kind::Literal;
        n.literal_mesh = std::move(m);
        return push(std::move(n));
    }

    Handle<T> cube(Vec<T, 3> half_extents = {T{0.5}, T{0.5}, T{0.5}}) {
        return literal(mesh::box_mesh<T>(half_extents));
    }

    // The parallel surface: every point of a closed base pushed along
    // its own normal. Closed in, closed out -- no rim anywhere, nothing
    // to decide at an edge, and the result is a surface in exactly the
    // sense the base was one.
    //
    // A base that is open, or that is not a surface at all, is refused
    // here, at the call that makes the mistake. That is the whole reason
    // this is a separate operation from offset_shell(): "offset" used to
    // mean both this and building a shell over a band, and the second
    // meaning quietly required an edge rule that nobody was stating.
    Handle<T> offset(Handle<T> base, ScalarField<T> thickness) {
        require_surface(base.index, "offset");
        if (!is_closed(resolve_surface(*this, base.index)))
            throw std::invalid_argument(
                "offset: the base surface has an edge, so its offset is a shell and "
                "needs a rule for what happens at that edge -- use offset_shell(base, "
                "thickness, EdgeRule). offset() is the parallel surface of a closed "
                "base and produces a closed surface.");
        TraceNode<T> n{};
        n.kind = Kind::Offset;
        n.base = base.index;
        n.thickness = std::move(thickness);
        // The base's tessellation, not this node's defaults. An offset
        // *is* its base pushed along its own normals, so the two are the
        // same grid by construction and a different one here would be a
        // second answer to a question the base already answered.
        //
        // It used to default to 48x24 regardless, which meant
        // `torus(2.0, 1.0, 160, 80)` followed by `offset()` silently
        // rendered at 48x24: the caller's request was accepted and
        // discarded, with nothing saying so. The donut demo had asked for
        // 160x80 since it was written and had been drawing 48x24. Same
        // family as every other defect this tree has caught -- a value
        // wearing the costume of an answer -- and the one place it shows
        // is a surface detail that never appears no matter how fine the
        // field describing it gets.
        n.u_steps = node(base.index).u_steps;
        n.v_steps = node(base.index).v_steps;
        return push(std::move(n));
    }

    Handle<T> offset(Handle<T> base, T thickness) {
        return offset(base, ScalarField<T>{[thickness](T, T) { return thickness; }});
    }

    // A shell over any base, open or closed, with the rim rule stated.
    // Same construction as offset() -- the base's own surface pushed
    // along its own normal -- but the caller has said what the edge
    // means instead of leaving it to whatever the renderer happened to
    // do with it.
    Handle<T> offset_shell(Handle<T> base, ScalarField<T> thickness, EdgeRule edge) {
        require_surface(base.index, "offset_shell");
        TraceNode<T> n{};
        n.kind = Kind::Offset;
        n.base = base.index;
        n.thickness = std::move(thickness);
        n.edge = edge;
        n.u_steps = node(base.index).u_steps;   // see offset(), same reason
        n.v_steps = node(base.index).v_steps;
        return push(std::move(n));
    }

    Handle<T> offset_shell(Handle<T> base, T thickness, EdgeRule edge) {
        return offset_shell(base, ScalarField<T>{[thickness](T, T) { return thickness; }}, edge);
    }

    // The target must be a surface, for the same reason and with the
    // same timing as offset()'s base: scattering is placement *on* a
    // space, and a Literal mesh or a group is not one. Unlike offset(),
    // an open target is fine -- sampling a band by its own area element
    // is well posed, and the rim never comes up.
    Handle<T> scatter(Handle<T> item, Handle<T> target, std::size_t count,
                      std::uint32_t seed = 42, T seat = T{1},
                      SeatAxis axis = SeatAxis::Z) {
        require_surface(target.index, "scatter");
        TraceNode<T> n{};
        n.kind = Kind::Scatter;
        n.item = item.index;
        n.target = target.index;
        n.count = count;
        n.seed = seed;
        n.seat = seat;
        n.seat_axis = axis;
        return push(std::move(n));
    }

    Handle<T> compose(std::initializer_list<Handle<T>> children) {
        return compose(std::vector<Handle<T>>(children));
    }

    // For a group whose size isn't known until runtime (e.g. one scatter
    // per color in a palette) -- same node, just not initializer_list-shaped.
    Handle<T> compose(const std::vector<Handle<T>>& children) {
        TraceNode<T> n{};
        n.kind = Kind::Compose;
        for (auto h : children) n.children.push_back(h.index);
        return push(std::move(n));
    }

private:
    // Both offset flavours and scatter need their base or target to be a
    // space. It always was required -- resolve_surface() threw otherwise
    // -- but it threw from inside materialize(), a long way from the
    // line that got it wrong, and only if the scene was ever
    // materialized at all.
    void require_surface(std::size_t idx, const char* op) const {
        Kind k = nodes_[idx].kind;
        if (k == Kind::Space || k == Kind::Offset) return;
        throw std::invalid_argument(
            std::string(op) + ": node " + std::to_string(idx) + " is a " + kind_name(k) +
            ", which is not a surface. Only Space and Offset nodes carry a (u,v) map; a "
            "Literal is a precomputed mesh and Scatter and Compose are many objects rather "
            "than one surface.");
    }

    Handle<T> push(TraceNode<T> n) {
        std::size_t idx = nodes_.size();
        nodes_.push_back(std::move(n));
        return Handle<T>{this, idx};
    }

    std::vector<TraceNode<T>> nodes_;
};

template<Scalar T>
Handle<T> Handle<T>::colored(Material<T> m) const {
    trace->node(index).material = std::move(m);
    return *this;
}

// Refuses on the spot, in the house pattern offset() established: a
// level this node cannot serve is a mistake in the line that asked for
// it, not a surprise for materialize() to raise later or -- worse -- for
// a renderer to paper over by silently picking something else. Silently
// picking something else is exactly how a caller ends up paying Newton's
// three orders of magnitude without ever having asked.
template<Scalar T>
Handle<T> Handle<T>::rendered_as(RenderLevel level) const {
    const auto& n = trace->node(index);
    auto refuse = [&](const char* why) {
        throw std::invalid_argument(std::string("rendered_as(") + render_level_name(level) +
                                     "): node " + std::to_string(index) + " is a " +
                                     kind_name(n.kind) + " and " + why);
    };
    if (level == RenderLevel::Exact && !n.exact.has_value())
        refuse("has no exact closed form. Only factories that know their shape record "
               "one -- torus() and cylinder() do, space() cannot, and .moving() drops it.");
    if (level == RenderLevel::Newton && n.kind != Kind::Space && n.kind != Kind::Offset)
        refuse("has no (u,v) map for Newton to iterate on.");
    if (level == RenderLevel::Tessellated && n.kind == Kind::Compose)
        refuse("is a group of objects rather than one, so it has no single mesh.");
    trace->node(index).level = level;
    return *this;
}

// Last call wins, deliberately: a color is a value a point maps to, and
// two such maps have no meaningful composition (unlike motion below).
template<Scalar T>
Handle<T> Handle<T>::colored(PointField<T> color_fn) const {
    trace->node(index).color_fn = std::move(color_fn);
    return *this;
}

template<Scalar T>
Handle<T> Handle<T>::glowing(PointField<T> emissive_fn) const {
    trace->node(index).emissive_fn = std::move(emissive_fn);
    return *this;
}

template<Scalar T>
Handle<T> Handle<T>::glowing(Vec<T, 3> emissive) const {
    trace->node(index).material.emissive = emissive;
    return *this;
}

// Composes, rather than replacing: `.moving(f).moving(g)` applies f then
// g, i.e. g(f(p, t), t). This used to be a plain assignment, so the
// earlier motion was silently dropped and the chaining syntax quietly
// meant something other than what it reads as -- with no diagnostic.
// Composition is also what makes the eventual structural form work: two
// expression trees compose by substituting one into the other's point
// slot, and that has to agree with what the callable form does here.
template<Scalar T>
Handle<T> Handle<T>::moving(PointField<T> f) const {
    // Refused on a Compose, at the call site, because until now it was
    // accepted and silently did nothing: materialize() never builds a
    // Placed for a Compose node -- it flattens to the children's -- so
    // the motion had nowhere to be read from. A call that returns a
    // handle and changes nothing is the same defect as a filter that
    // skips nothing; see conventions.md.
    //
    // Refused rather than implemented, and that is the design decision
    // rather than the cheap way out. Giving a Compose its own transform
    // would make a node's effective transform the product of its
    // ancestors' -- a scene graph, walked by following parents, which is
    // exactly what addressing a flat array by index exists to avoid.
    // Group motion belongs at the moment operations are expanded into
    // objects, where it folds into each object's own transform once and
    // the array stays flat. Until that expansion exists, move the
    // children.
    if (trace->node(index).kind == Kind::Compose)
        throw std::logic_error(
            "moving() on a Compose node: a Compose is grouping, not an object -- it "
            "materializes to its children, so a motion here would have nothing to read "
            "it. Apply .moving() to each child, or compose the moved children.");

    // The exact form survives a *placement* and not a deformation, and
    // this is the question that used to be unanswerable.
    //
    // The rule was "clear it, always", with the reason written here: an
    // arbitrary point map does not in general send a torus to a torus, so
    // a recorded shape that no longer agrees with the map is worse than
    // no recorded shape -- the render would be right about the shape and
    // wrong about the scene. And recognising an isometry meant asking an
    // opaque callable what it does, which a callable cannot answer. The
    // note ended: *it becomes answerable once motion has a structural
    // form, and that is where it belongs.*
    //
    // It has one. `is_placement()` is derived by walking, not asked of a
    // closure: the motion is affine in the point, so it translates and
    // uniformly scales and nothing else. That sends a sphere to a sphere
    // and a torus to a torus. The recorded form stays the *rest* shape --
    // where the object is belongs to the placement, which a renderer reads
    // separately -- so the two cannot drift apart.
    //
    // A deformation still clears it, for the original reason, unchanged.
    if (!f.is_placement()) trace->node(index).exact.reset();

    auto& slot = trace->node(index).transform;
    if (!slot) {
        slot = std::move(f);
        return *this;
    }
    slot = PointField<T>{
        [prev = std::move(slot), next = std::move(f)](const Vec<T, 3>& p, T t) -> Vec<T, 3> {
            return next(prev(p, t), t);
        }};
    return *this;
}

// ── Resolution: Space/Offset nodes only, stays analytic ──────────

template<Scalar T>
ParametricSurface<T> resolve_surface(const Trace<T>& trace, std::size_t idx) {
    const auto& n = trace.node(idx);
    if (n.kind == Kind::Space) return *n.surface;
    // T cannot be deduced from a Field through offset_surface's
    // std::function parameter, so the erasure is spelled here. This is
    // the boundary where the field's structure stops travelling and only
    // its value continues: the node keeps the Field, so the report and
    // any lowering pass still see the expression.
    if (n.kind == Kind::Offset)
        return offset_surface<T>(resolve_surface(trace, n.base),
                                 std::function<T(T, T)>{n.thickness});
    throw std::logic_error("resolve_surface: node is not a Space or Offset");
}

// ── Materialization: the one place a mesh gets built, for display ──

template<Scalar T>
mesh::Mesh<Euclidean<3, T>> materialize_mesh(const Trace<T>& trace, std::size_t idx, T t, bool placed) {
    const auto& n = trace.node(idx);
    mesh::Mesh<Euclidean<3, T>> out;
    std::vector<Vec<T, 3>> scatter_seats;   // one per site, for the per-instance origin below

    if (n.kind == Kind::Space || n.kind == Kind::Offset) {
        auto surface = resolve_surface(trace, idx);
        auto pm = mesh::parametric_mesh(surface, n.u_steps, n.v_steps);
        out.vertices.assign(pm.vertices.begin(), pm.vertices.end());
        out.faces = pm.faces;
    } else if (n.kind == Kind::Literal) {
        out = *n.literal_mesh;
    } else if (n.kind == Kind::Scatter) {
        auto target_surface = resolve_surface(trace, n.target);
        auto sites = sample_surface_uniform(target_surface, n.count, n.seed);
        auto item_mesh = materialize_mesh(trace, n.item, t, /*placed=*/false);
        const T lift = scatter_lift<T>(item_mesh, n.seat_axis) * n.seat;
        scatter_seats.reserve(sites.size());

        out.vertices.reserve(item_mesh.vertex_count() * sites.size());
        out.faces.reserve(item_mesh.face_count() * sites.size());
        for (std::size_t i = 0; i < sites.size(); ++i) {
            const auto& site = sites[i];

            // The same frame cook() gives the instance, from the same
            // function, so the baked picture and the instanced one cannot
            // drift apart. See scatter_frame().
            auto frame = scatter_frame<T>(site.normal, n.seed, i, n.seat_axis);

            uint32_t base_idx = static_cast<uint32_t>(out.vertices.size());
            Vec<T, 3> seat{site.position + site.normal * lift};
            scatter_seats.push_back(seat);
            for (const auto& iv : item_mesh.vertices)
                out.vertices.push_back(Vec<T, 3>{seat + frame * iv});
            for (const auto& f : item_mesh.faces)
                out.faces.push_back({f[0] + base_idx, f[1] + base_idx, f[2] + base_idx});
        }
    } else {
        throw std::logic_error("materialize_mesh: Compose has no single mesh, use materialize()");
    }

    // A Scatter applies its motion per site, with that site's seat as the
    // instance origin -- so a field reading `origin` can send every
    // instance somewhere different while staying affine in the point.
    // Vertices were appended per site with a fixed stride above, so the
    // split is arithmetic rather than bookkeeping.
    //
    // This has to agree, site for site, with what cook() computes for the
    // same node; the vertex-for-vertex test between the two paths is what
    // says it does.
    if (placed && n.transform && n.kind == Kind::Scatter && !scatter_seats.empty()) {
        const std::size_t per = out.vertices.size() / scatter_seats.size();
        for (std::size_t i = 0; i < scatter_seats.size(); ++i)
            for (std::size_t k = 0; k < per; ++k) {
                auto& v = out.vertices[i * per + k];
                v = n.transform(MotionEnv<T>{v, t, scatter_seats[i]});
            }
        return out;
    }

    // `placed = false` gives the *rest* geometry -- the shape before this
    // node's own motion moves it. That is what cook() has to hash and
    // share: two dust specks are the same little cube and become two
    // different meshes only because each has been moved somewhere else.
    // Hashing the placed form would make every instance unique by
    // construction, which is exactly the way the first version of cook()
    // silently found 19 801 shapes for 19 801 objects.
    if (placed && n.transform)
        for (auto& v : out.vertices) v = n.transform(v, t);
    return out;
}

// ── Placed: what the object is, not what it tessellates into ──

template<Scalar T>
mesh::Mesh<Euclidean<3, T>> Placed<T>::mesh() const {
    return materialize_mesh(*trace, index, t);
}

template<Scalar T>
bool Placed<T>::is_exact() const {
    return trace->node(index).exact.has_value();
}

template<Scalar T>
const std::type_info& Placed<T>::exact_type() const {
    return trace->node(index).exact.type();
}

template<Scalar T>
template<class Shape>
const Shape* Placed<T>::exact_as() const {
    return std::any_cast<Shape>(&trace->node(index).exact);
}

// Inference, in one place so a renderer cannot drift from it: an exact
// form is used when there is one, and everything else tessellates.
// Newton is never inferred -- at 23 612 ns/ray it is something a caller
// opts into, never something a default hands them.
template<Scalar T>
RenderLevel Placed<T>::render_level() const {
    const auto& n = trace->node(index);
    if (n.level) return *n.level;
    return n.exact.has_value() ? RenderLevel::Exact : RenderLevel::Tessellated;
}

template<Scalar T>
bool Placed<T>::is_analytic() const {
    auto k = trace->node(index).kind;
    return k == Kind::Space || k == Kind::Offset;
}

template<Scalar T>
std::optional<ParametricSurface<T>> Placed<T>::surface() const {
    if (!is_analytic()) return std::nullopt;
    auto base = resolve_surface(*trace, index);
    if (!trace->node(index).transform) return base;

    // The node's motion composed into the map rather than applied to
    // vertices afterwards, so the result is still a surface and not a
    // deformed mesh. ParametricSurface derives normals by finite
    // differences of this map, so a deforming motion is accounted for
    // without anyone hand-deriving a Jacobian.
    //
    // Captures the trace and the index rather than the node itself: the
    // node's motion slot is a move_only_function, and ParametricSurface
    // stores its map in a std::function, which demands a copyable
    // callable. Same lifetime contract as the rest of this type.
    const Trace<T>* tr = trace;
    std::size_t i = index;
    T tt = t;
    return ParametricSurface<T>(
        [tr, i, tt, base](T u, T v) -> Vec<T, 3> {
            auto p = base.evaluate(u, v);
            const auto& node = tr->node(i);
            return node.transform ? node.transform(p, tt) : p;
        },
        base.domain(), base.periodic_u(), base.periodic_v());
}

// A node's material, with its colour and emission fields evaluated at one
// representative point.
//
// A free function rather than a method because there are two callers that
// must not disagree: `Placed::material()`, which answers for a whole node,
// and `cook()`, which answers per object. They differ only in *which*
// point is representative -- a node's placed centroid, or one instance's
// -- and keeping the rest in one place is what stops that difference from
// quietly growing into two different materials.
template<Scalar T>
Material<T> resolve_material(const TraceNode<T>& n, const Vec<T, 3>& at, T t) {
    Material<T> mat = n.material;
    if (n.color_fn)    mat.base_color = n.color_fn(at, t);
    if (n.emissive_fn) mat.emissive   = n.emissive_fn(at, t);
    return mat;
}

template<Scalar T>
Material<T> Placed<T>::material() const {
    const auto& n = trace->node(index);
    // The early out is not a micro-optimisation: mesh() rebuilds on every
    // call, so a node with no fields would otherwise tessellate itself to
    // answer a question whose answer does not depend on the geometry.
    if (!n.color_fn && !n.emissive_fn) return n.material;
    return resolve_material(n, Vec<T, 3>{mesh().centroid()}, t);
}

// Walks a Compose node into its leaves. Each leaf comes back as a view
// onto its own node -- nothing is tessellated here, and a caller that
// only wants the analytic form never causes a tessellation at all.
template<Scalar T>
std::vector<Placed<T>> materialize(const Trace<T>& trace, std::size_t idx, T t) {
    const auto& n = trace.node(idx);
    if (n.kind == Kind::Compose) {
        std::vector<Placed<T>> out;
        for (auto c : n.children) {
            auto sub = materialize(trace, c, t);
            out.insert(out.end(), sub.begin(), sub.end());
        }
        return out;
    }
    return {Placed<T>{&trace, idx, t}};
}

// ── The field report ─────────────────────────────────────────────
//
// What the trace's fields actually are, walked once. Cheap enough to call
// per frame -- measured at 2.8-5.1 ns per field, so 112-204 us across the
// donut demo's ~39 600 against an 80 ms frame -- which is why there is no
// cache here and no need for one. Note that this is a property of the
// trace being immutable once built: if a scene ever mutates between
// frames, computing it once at that boundary is the obvious move.
//
// Reports recognized and unknown separately, and prints unknown even when
// it is zero. A single "distinct types" number cannot distinguish "they
// really are all one closure" from "they all fell into one bucket nothing
// could classify", and that is the defect class conventions.md names.
template<Scalar T>
FieldStats field_report(const Trace<T>& trace) {
    FieldStats stats{};
    std::vector<std::type_index> seen;
    for (std::size_t i = 0; i < trace.size(); ++i) {
        const auto& n = trace.node(i);
        if (n.kind == Kind::Offset) accumulate(stats, n.thickness, seen);
        accumulate(stats, n.transform, seen);
        accumulate(stats, n.color_fn, seen);
    }
    return stats;
}

// ── cook(): operations become objects, once ──────────────────────
//
// A `Trace` is mutated while it is built and read-only while it is
// rendered. The DSL had that boundary in practice and never said so.
// `cook()` says it.
//
// **Takes the trace by const reference, not by value.** `Cooked` holds
// objects, shapes and the report -- it never keeps the trace -- so
// consuming one bought nothing and cost the caller the scene they still
// need to render. Expressing "the build phase is over" is the job of
// `Cooked` having no builders, not of destroying the input.
//
// **A type, not a flag.** A `frozen_` bool would be a value that looks
// like a guarantee and holds none: nothing stops `trace.torus()` after it
// is set, and no compiler notices. `Cooked<T>` holds the guarantee
// structurally -- it has no `torus()`, no `offset()`, no `scatter()`.
// Mutation after cooking is not discouraged; it does not compile.
//
// **Five things that were being tracked separately are one piece of
// work**, because each of them needs the whole scene to be known:
// expanding operations into objects, deduplicating shapes, folding a
// group's transform into its members, computing the field report once,
// and (later) compaction and AoS->SoA. One well-placed boundary closing
// several problems is the sign that they were one problem.
//
// Deduplication is *part of the expansion*, not an optimisation on top of
// it: without it there is no "one shape, N transforms" and so no
// instancing at all. Compaction and layout are the optimisations, and
// they are a different kind of thing.
//
// **`materialize()` is deliberately left alone.** It answers "what is at
// this node" and returns one `Placed` per node, which is what the demos,
// the benchmarks and a dozen tests rely on -- a `Scatter` is one `Placed`
// holding a merged mesh. Cooking answers a different question, "what
// objects are in this scene", and a `Scatter` is N of them. Two
// questions, two answers, neither pretending to be the other.

// One instance in a cooked scene: which shape, where, and what it looks
// like. The shape is an index into the shape table, so N instances of one
// speck cost N transforms and one geometry.
template<Scalar T = double>
struct Object {
    std::size_t shape = 0;        // index into Cooked::shapes()
    std::size_t source_node = 0;  // the trace node it came from, for reporting
    Vec<T, 3> translation{};      // folded from ancestors and from the node
    T scale = T{1};               // uniform, from the node's placement

    // Orientation, and it is not decoration. `materialize_mesh()` has
    // oriented scattered items since the sprinkles complaint -- each site
    // gets the {t1, t2, normal} frame plus a deterministic in-plane turn
    // -- while cooking kept only `site.position` and dropped the frame on
    // the floor. So the two answers to "where is this item" disagreed for
    // every scattered object, and the disagreement was invisible because
    // no renderer consumed the cooked one yet.
    //
    // Orthonormal, and paired with a *uniform* scale rather than folded
    // into one general 3x3, for the reason `VecField::Placement` records:
    // that pair is what lets a ray test go into local space without
    // rescaling `t` or repairing a normal.
    Matrix<T, 3, 3> rotation = Matrix<T, 3, 3>::identity();

    // The same rotation, compactly -- 32 bytes against the matrix's 72,
    // which is 80 MB at two million objects.
    //
    // **The matrix is authoritative while both are here.** Everything
    // downstream reads `rotation`; this is carried and compared, and
    // nothing decides anything by it. Two live sources of truth for one
    // quantity is how they drift apart in silence, so one of them is
    // explicitly not a source.
    //
    // That is also what makes this step checkable. The renderer's output
    // cannot move, because the thing it draws with has not changed -- so a
    // frame that is not byte-identical means the plumbing broke, and
    // nothing subtler. The round trip's own error arrives in the *next*
    // step, when the matrix goes and this becomes the source; measured at
    // about 8 ulp, with zero of four thousand rotations surviving
    // bit-exact (tests/test_build_dsl.cpp).
    Quaternion<T> rotation_q{};

    bool instanceable = true;     // false when the motion deforms per vertex
    Material<T> material{};
};

// The shape table: one entry per distinct geometry in the scene.
template<Scalar T = double>
struct Shape {
    mesh::Mesh<Euclidean<3, T>> geometry;
    std::size_t instances = 0;    // how many objects point here

    // The closed form, when the node this shape came from had one --
    // carried across the deduplication rather than dropped at it.
    //
    // Dropping it had a consequence worth stating, because it made the
    // instancing story incomplete in exactly the place the DSL says
    // "many": a node that is its own object (the donut's dust, one node
    // per particle) kept its exact form and could be instanced, while a
    // `Scatter` -- the one operation whose entire meaning is "N of
    // these" -- came out of cook() as N meshes of one mesh, with the
    // BoundedQuadric that `cylinder()` had carefully recorded nowhere to
    // be found. The dedup key is `content_hash` of the geometry node, so
    // two objects sharing a shape share its exact form too, by the same
    // argument that lets them share the mesh.
    std::any exact;
};

template<Scalar T = double>
class Cooked {
public:
    const std::vector<Object<T>>& objects() const { return objects_; }
    const std::vector<Shape<T>>& shapes() const { return shapes_; }
    const FieldStats& fields() const { return fields_; }

    std::size_t object_count() const { return objects_.size(); }
    std::size_t shape_count() const { return shapes_.size(); }

    // Two counters, because one would lie in the same way a single
    // "distinct types" number lied: an object that could not be instanced
    // because its motion is opaque, and an object that was never a
    // candidate because its shape is unique, are different facts. Rolled
    // into one "not instanced" they read identically, and the honest
    // reading -- "the mechanism works, the scene simply has nothing to
    // share" -- becomes indistinguishable from "the mechanism is
    // refusing everything".
    //
    // They do not sum to object_count(), and that is correct: an object
    // can be unique *and* structural, which is neither.
    std::size_t shared_objects() const { return shared_; }
    std::size_t opaque_refused() const { return refused_; }

    // The nodes that were refused, so the report can point at them rather
    // than state a number nobody can act on.
    const std::vector<std::size_t>& refused_nodes() const { return refused_nodes_; }

    // The vertices a renderer would build if every object carried its own
    // copy, against what the shape table actually holds. The ratio is the
    // whole point of the pass, and it is reported rather than assumed.
    std::size_t vertices_without_instancing() const {
        std::size_t n = 0;
        for (const auto& o : objects_) n += shapes_[o.shape].geometry.vertex_count();
        return n;
    }
    std::size_t vertices_stored() const {
        std::size_t n = 0;
        for (const auto& s : shapes_) n += s.geometry.vertex_count();
        return n;
    }

private:
    template<Scalar U> friend Cooked<U> cook(const Trace<U>&, std::size_t, U);
    std::vector<Object<T>> objects_;
    std::vector<Shape<T>> shapes_;
    FieldStats fields_{};
    std::size_t shared_ = 0, refused_ = 0;
    std::vector<std::size_t> refused_nodes_;
};

// The key two nodes must share to be the same shape. Node identity is not
// enough: the donut's dust is 19 800 separate Literal nodes all holding
// the result of the same `dust_speck(0.014)` call -- equal without being
// the same object.
//
// Three kinds of node, three different keys, and they cost different
// things:
//
//   Literal   the mesh's own content. O(vertices) per node, which is the
//             expensive one -- ~158 000 vertices hashed across the
//             donut's dust. Acceptable because cook() runs once per
//             scene; if it ever runs per frame, this is the line to look
//             at first.
//   Space     the exact analytic form's *type* plus the chart sampled at
//             fixed (u,v), plus the tessellation steps. Sampling the
//             chart rather than reading parameters back out of the node
//             is deliberate: it hashes what the geometry will actually
//             be, so it cannot drift away from what gets tessellated the
//             way a reconstructed parameter list could. The exact form's
//             type is folded in so two nodes with the same chart but
//             different closed forms never merge.
//   Offset    the base's key, plus the thickness field sampled the same
//             way. An offset of the same base by the same thickness is
//             the same shape.
template<Scalar T>
inline std::size_t content_hash(const Trace<T>& trace, std::size_t idx) {
    std::size_t h = 1469598103934665603ull;
    auto mix = [&h](std::size_t v) { h = (h ^ v) * 1099511628211ull; };
    auto mix_scalar = [&](T v) { mix(std::hash<double>{}(static_cast<double>(v))); };

    const auto& n = trace.node(idx);
    mix(static_cast<std::size_t>(n.kind));
    mix(n.u_steps);
    mix(n.v_steps);

    if (n.kind == Kind::Literal) {
        const auto& m = *n.literal_mesh;
        mix(m.vertex_count());
        mix(m.face_count());
        for (const auto& v : m.vertices)
            for (std::size_t i = 0; i < 3; ++i) mix_scalar(v[i]);
        return h;
    }

    if (n.kind == Kind::Offset) {
        mix(content_hash(trace, n.base));
        for (int i = 0; i <= 4; ++i)
            for (int j = 0; j <= 4; ++j)
                mix_scalar(n.thickness(T(i) * T{0.25}, T(j) * T{0.25}));
        return h;
    }

    if (n.kind == Kind::Space) {
        if (n.exact.has_value()) mix(n.exact.type().hash_code());
        const auto& s = *n.surface;
        auto [u0, u1, v0, v1] = s.domain();
        mix_scalar(u0); mix_scalar(u1); mix_scalar(v0); mix_scalar(v1);
        for (int i = 0; i <= 3; ++i)
            for (int j = 0; j <= 3; ++j) {
                auto p = s.evaluate(u0 + (u1 - u0) * T(i) / T{3}, v0 + (v1 - v0) * T(j) / T{3});
                for (std::size_t k = 0; k < 3; ++k) mix_scalar(p[k]);
            }
        return h;
    }

    // Scatter: keyed by its item, since that is the geometry being
    // repeated; the placements are what differ and they live on objects.
    mix(content_hash(trace, n.item));
    mix(n.count);
    mix(n.seed);
    return h;
}

template<Scalar T = double>
Cooked<T> cook(const Trace<T>& trace, std::size_t root, T t = T{0}) {
    Cooked<T> out;
    std::unordered_map<std::size_t, std::size_t> shape_of_key;

    // Where one instance of an operation sits and which way it faces. The
    // frame is carried alongside the position rather than recomputed
    // later, because it belongs to the site -- dropping it is exactly the
    // bug this pair fixes.
    //
    // Declared here rather than inside the walk below: a local class
    // written in a *generic* lambda's body resolves `T` against the
    // lambda's own invented parameter list under clang, which reports
    // `Matrix<T, 3, 3>` with T deduced as the lambda type itself. gcc
    // accepts it. The construct buys nothing, so it is simply not used.
    struct Spot {
        Vec<T, 3> position{};
        Matrix<T, 3, 3> frame = Matrix<T, 3, 3>::identity();
    };

    // A group's transform is folded into its members here, which is why
    // `Compose` never needs one of its own at render time and why the
    // result stays a flat array: the accumulator lives in the walk, not
    // in the output. Nothing downstream has to follow a parent pointer.
    auto walk = [&](auto&& self, std::size_t idx, const Vec<T, 3>& carried) -> void {
        const auto& n = trace.node(idx);

        if (n.kind == Kind::Compose) {
            for (auto c : n.children) self(self, c, carried);
            return;
        }

        // How many instances this operation is -- see Spot above.
        std::vector<Spot> placements;
        std::size_t geometry_node = idx;

        if (n.kind == Kind::Scatter) {
            auto target = resolve_surface(trace, n.target);
            auto sites = sample_surface_uniform(target, n.count, n.seed);
            const T lift =
                scatter_lift<T>(materialize_mesh(trace, n.item, t, /*placed=*/false), n.seat_axis) *
                n.seat;
            placements.reserve(sites.size());
            for (std::size_t i = 0; i < sites.size(); ++i)
                placements.push_back(
                    Spot{Vec<T, 3>{sites[i].position + sites[i].normal * lift},
                         scatter_frame<T>(sites[i].normal, n.seed, i, n.seat_axis)});
            geometry_node = n.item;
        } else {
            placements.push_back(Spot{});
        }

        // A renderer can instance this object -- draw one geometry many
        // times, each with its own transform -- exactly when the node's
        // motion is a *placement*: affine in the point, so it moves the
        // object without reshaping it.
        //
        // Note what the test is not. It is not "is the motion structural".
        // The donut's dust is Perlin noise and will always be opaque, yet
        // it is a placement, because the noise decides *where the particle
        // is* and never reads the vertex. Opacity and point-dependence are
        // different questions, and only the second one costs anything:
        // a term that ignores the point is evaluated once per object, a
        // term that reads it once per vertex.
        const bool placeable = n.transform.is_placement();

        // ── The refused path, which used to hand back unusable geometry ──
        //
        // A deformation cannot be decomposed into per-instance transforms
        // -- that is what makes it a deformation -- so there is nothing to
        // instance and nothing to share. The geometry *is* the answer, and
        // it has to be the **placed** geometry: storing the rest shape here
        // left a consumer with a mesh and no way to recover the motion that
        // moved it, which would have drawn the donut demo's exploding cube
        // unexploded. Nothing noticed for as long as nothing consumed a
        // cooked scene.
        //
        // No deduplication either, and that is not a shortcut: a placed
        // shape is unique by construction, which the comment on the shared
        // path below has always said. Hashing it would find no matches and
        // cost a full pass over the vertices to find none.
        //
        // And one object, not N: a refused `Scatter` is a single merged
        // mesh, exactly the answer `materialize()` gives for it. The
        // refusal *count* still adds the instances, because the number
        // answers "how much sharing was lost", which is a question about
        // instances rather than about objects.
        if (!placeable) {
            const std::size_t shape_index = out.shapes_.size();
            auto placed_mesh = materialize_mesh(trace, idx, t, /*placed=*/true);
            auto at = placed_mesh.centroid();
            out.shapes_.push_back(Shape<T>{std::move(placed_mesh), 1, std::any{}});
            out.objects_.push_back(Object<T>{.shape = shape_index,
                                             .source_node = idx,
                                             .translation = carried,
                                             .scale = T{1},
                                             .rotation = Matrix<T, 3, 3>::identity(),
                                             .rotation_q = Quaternion<T>{},
                                             .instanceable = false,
                                             .material = resolve_material(n, at, t)});
            out.refused_ += placements.size();
            out.refused_nodes_.push_back(idx);
            return;
        }

        const auto key = content_hash(trace, geometry_node);
        auto it = shape_of_key.find(key);
        std::size_t shape_index;
        if (it == shape_of_key.end()) {
            shape_index = out.shapes_.size();
            // Rest geometry, not placed: the shape before this node's own
            // motion moved it. Sharing is only possible between shapes
            // that have not yet been put anywhere.
            out.shapes_.push_back(
                Shape<T>{materialize_mesh(trace, geometry_node, t, /*placed=*/false), 0,
                         trace.node(geometry_node).exact});
            shape_of_key.emplace(key, shape_index);
        } else {
            shape_index = it->second;
        }
        const auto rest_centroid = Vec<T, 3>{out.shapes_[shape_index].geometry.centroid()};

        // Composing the two transforms, written out because getting the
        // order wrong is silent. A vertex `v` of the rest item becomes
        //
        //     carried + b + R*s*(site + F*v)
        //             = [carried + b + R*(s*site)] + (R*F)*(s*v)
        //
        // which is where each line below comes from. Note the `s*site`:
        // the old form added the site position *outside* the node's
        // scale, so a scattered node whose motion also scaled placed its
        // items at the wrong distance -- invisible until now only because
        // nothing consumed the cooked scene and every scattering scene so
        // far had scale 1.
        for (const auto& p : placements) {
            // Resolved per instance, with that instance's seat as its
            // origin -- which is what lets one motion field send every
            // scattered item somewhere different while staying a
            // placement. For a node that is its own object the seat is
            // the origin and this reduces to what it was.
            const auto pl = n.transform.placement_at(MotionEnv<T>{Vec<T, 3>{}, t, p.position});

            Object<T> o{.shape = shape_index,
                        .source_node = idx,
                        .translation = Vec<T, 3>{carried + pl.translation +
                                                 pl.rotation * Vec<T, 3>{p.position * pl.scale}},
                        .scale = pl.scale,
                        .rotation = pl.rotation * p.frame,
                        .rotation_q = {},
                        .instanceable = true,
                        .material = {}};
            o.rotation_q = Quaternion<T>::from_matrix(o.rotation);
            // Resolved where this object actually ends up, and computed
            // from the rest centroid rather than by building the placed
            // mesh -- which is the same point, and is the point
            // `Placed::material()` uses, so the two paths cannot answer
            // differently.
            o.material = resolve_material(
                n, Vec<T, 3>{o.rotation * Vec<T, 3>{rest_centroid * o.scale} + o.translation}, t);
            out.objects_.push_back(std::move(o));
            ++out.shapes_[shape_index].instances;
        }
    };

    walk(walk, root, Vec<T, 3>{});

    // Counted after the walk, because "shared" is a property of the
    // finished grouping: an object is shared if the shape it points at
    // ended up with more than one instance, which is not knowable while
    // the group is still being filled.
    for (const auto& o : out.objects_)
        if (out.shapes_[o.shape].instances > 1) ++out.shared_;

    out.fields_ = field_report(trace);
    return out;
}

} // namespace spatium::io::build
