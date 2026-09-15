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
#  include <spatium/geometry/ray_surface.hpp>
#  include <any>
#  include <cstdint>
#  include <functional>
#  include <initializer_list>
#  include <optional>
#  include <stdexcept>
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
// A node's motion/thickness/color hooks are std::move_only_function,
// not std::function, for two reasons that are about capability rather
// than speed (the indirect call itself measures ~3.3 ns/vertex, ~4% of
// a realistic motion callee -- see benchmarks/bench_trace.cpp):
//
//   - std::function requires its callable to be copy-constructible, so
//     anything a hook wants to own has to be copyable too. That makes
//     capturing heavy state by value the path of least resistance:
//     donut_demo captured a 512-byte PerlinNoise into each of 19 800
//     closures, ~9.7 MB of identical tables, because sharing it would
//     have needed a hand-rolled shared_ptr dance. move_only_function
//     accepts move-only state directly.
//   - The signatures are const-qualified, so a hook is callable through
//     the `const Trace&` materialize() actually holds. std::function's
//     operator() is const but happily calls a non-const callable, a
//     known hole this type closes.
template<Scalar T = double>
using PointField = std::move_only_function<Vec<T, 3>(const Vec<T, 3>&, T) const>;

// Deliberately still std::function, unlike PointField: a thickness flows
// through offset_surface() into a ParametricSurface, whose ParamFn is
// itself a std::function and therefore requires a copy-constructible
// callable. The constraint belongs to ParametricSurface, not to the DSL,
// and moving this slot needs that type to change first.
template<Scalar T = double>
using ScalarField = std::function<T(T, T)>;

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

template<Scalar T = double>
class Trace;

template<Scalar T = double>
struct Handle {
    Trace<T>* trace;
    std::size_t index;

    Handle colored(Material<T> m) const;
    Handle rendered_as(RenderLevel level) const;
    Handle colored(PointField<T> color_fn) const;
    Handle moving(PointField<T> f) const;
};

template<Scalar T>
ParametricSurface<T> resolve_surface(const Trace<T>& trace, std::size_t idx);

template<Scalar T>
mesh::Mesh<Euclidean<3, T>> materialize_mesh(const Trace<T>& trace, std::size_t idx, T t = T{0});

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
    Handle<T> scatter(Handle<T> item, Handle<T> target, std::size_t count, std::uint32_t seed = 42) {
        require_surface(target.index, "scatter");
        TraceNode<T> n{};
        n.kind = Kind::Scatter;
        n.item = item.index;
        n.target = target.index;
        n.count = count;
        n.seed = seed;
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

// Composes, rather than replacing: `.moving(f).moving(g)` applies f then
// g, i.e. g(f(p, t), t). This used to be a plain assignment, so the
// earlier motion was silently dropped and the chaining syntax quietly
// meant something other than what it reads as -- with no diagnostic.
// Composition is also what makes the eventual structural form work: two
// expression trees compose by substituting one into the other's point
// slot, and that has to agree with what the callable form does here.
template<Scalar T>
Handle<T> Handle<T>::moving(PointField<T> f) const {
    // The exact form goes, always. `f` is an arbitrary point map, so in
    // general it does not send a torus to a torus, and a recorded shape
    // that no longer agrees with the map is worse than no recorded shape
    // at all: the render would be correct for the shape and wrong for
    // the scene, with nothing to notice.
    //
    // Always, not "unless f is an isometry", on purpose. Recognising an
    // isometry means asking an opaque callable what it does, which is
    // exactly the question a callable cannot answer -- it becomes
    // answerable once motion has a structural form, and that is where it
    // belongs. Until then this costs the exact path on a moving node and
    // keeps the invariant true, which is the cheaper of the two mistakes.
    trace->node(index).exact.reset();

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
    if (n.kind == Kind::Offset) return offset_surface(resolve_surface(trace, n.base), n.thickness);
    throw std::logic_error("resolve_surface: node is not a Space or Offset");
}

// ── Materialization: the one place a mesh gets built, for display ──

template<Scalar T>
mesh::Mesh<Euclidean<3, T>> materialize_mesh(const Trace<T>& trace, std::size_t idx, T t) {
    const auto& n = trace.node(idx);
    mesh::Mesh<Euclidean<3, T>> out;

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
        auto item_mesh = materialize_mesh(trace, n.item, t);

        out.vertices.reserve(item_mesh.vertex_count() * sites.size());
        out.faces.reserve(item_mesh.face_count() * sites.size());
        for (std::size_t i = 0; i < sites.size(); ++i) {
            const auto& site = sites[i];
            auto [t1, t2] = basis_from_normal(site.normal);

            // Per-instance in-plane rotation about the normal -- without
            // it every item at every site shares the same {t1,t2} up to
            // the smooth drift basis_from_normal has with the surface's
            // own curvature, which reads as "all facing the same way"
            // (a real complaint on the donut's sprinkles). Deterministic
            // per (seed, site index) -- same trace, same t, same result.
            std::uint32_t h = n.seed * 2654435761u + static_cast<std::uint32_t>(i) * 40503u;
            h ^= h >> 13; h *= 0x85ebca6bu; h ^= h >> 16;
            T theta = (static_cast<T>(h % 360000) / T{100000}) * T{6.283185307179586};
            using std::cos, std::sin;
            Vec<T, 3> r1{Vec<T, 3>{t1 * cos(theta) + t2 * sin(theta)}};
            Vec<T, 3> r2{Vec<T, 3>{t2 * cos(theta) - t1 * sin(theta)}};

            uint32_t base_idx = static_cast<uint32_t>(out.vertices.size());
            for (const auto& iv : item_mesh.vertices)
                out.vertices.push_back(Vec<T, 3>{
                    site.position + r1 * iv[0] + r2 * iv[1] + site.normal * iv[2]});
            for (const auto& f : item_mesh.faces)
                out.faces.push_back({f[0] + base_idx, f[1] + base_idx, f[2] + base_idx});
        }
    } else {
        throw std::logic_error("materialize_mesh: Compose has no single mesh, use materialize()");
    }

    if (n.transform)
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

template<Scalar T>
Material<T> Placed<T>::material() const {
    const auto& n = trace->node(index);
    Material<T> mat = n.material;
    if (n.color_fn) mat.base_color = n.color_fn(mesh().centroid(), t);
    return mat;
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

} // namespace spatium::io::build
