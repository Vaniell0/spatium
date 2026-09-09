#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/io/scene.hpp>
#  include <spatium/mesh/mesh.hpp>
#  include <spatium/mesh/primitives.hpp>
#  include <spatium/spaces/euclidean.hpp>
#  include <spatium/spaces/offset.hpp>
#  include <spatium/spaces/parametric.hpp>
#  include <spatium/spaces/sample.hpp>
#  include <cstdint>
#  include <functional>
#  include <initializer_list>
#  include <optional>
#  include <stdexcept>
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

enum class Kind { Space, Offset, Scatter, Compose, Literal };

template<Scalar T = double>
struct TraceNode {
    Kind kind{};

    // Space
    std::optional<ParametricSurface<T>> surface;
    std::size_t u_steps = 48, v_steps = 24;

    // Literal -- a precomputed mesh, for shapes with no natural single
    // (u,v)->R^3 formula (a cube, most obviously). The escape hatch out
    // of the analytic-first path, not the default.
    std::optional<mesh::Mesh<Euclidean<3, T>>> literal_mesh;

    // Offset
    std::size_t base = 0;
    std::function<T(T, T)> thickness;

    // Scatter
    std::size_t item = 0, target = 0;
    std::size_t count = 0;
    std::uint32_t seed = 42;

    // Compose
    std::vector<std::size_t> children;

    Material<T> material{};
    std::function<Vec<T, 3>(const Vec<T, 3>&, T)> transform;

    // Optional: color as a function of (a representative point, t) --
    // evaluated once per node per materialize() call (at the
    // materialized mesh's own centroid), same "evaluate at t" spirit as
    // `transform`, just producing a color instead of a position. When
    // set, overrides `material.base_color` for that materialize() call.
    std::function<Vec<T, 3>(const Vec<T, 3>&, T)> color_fn;
};

template<Scalar T = double>
struct Placed {
    mesh::Mesh<Euclidean<3, T>> mesh;
    Material<T> material;
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
    Handle colored(std::function<Vec<T, 3>(const Vec<T, 3>&, T)> color_fn) const;
    Handle moving(std::function<Vec<T, 3>(const Vec<T, 3>&, T)> f) const;
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

    Handle<T> torus(T major_r, T minor_r, std::size_t u_steps = 48, std::size_t v_steps = 24) {
        return space(make_torus<T>(major_r, minor_r), u_steps, v_steps);
    }

    Handle<T> cylinder(T radius, T height, std::size_t u_steps = 12, std::size_t v_steps = 4) {
        return space(make_cylinder<T>(radius, height), u_steps, v_steps);
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

    Handle<T> offset(Handle<T> base, std::function<T(T, T)> thickness) {
        TraceNode<T> n{};
        n.kind = Kind::Offset;
        n.base = base.index;
        n.thickness = std::move(thickness);
        return push(std::move(n));
    }

    Handle<T> offset(Handle<T> base, T thickness) {
        return offset(base, std::function<T(T, T)>{[thickness](T, T) { return thickness; }});
    }

    Handle<T> scatter(Handle<T> item, Handle<T> target, std::size_t count, std::uint32_t seed = 42) {
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

template<Scalar T>
Handle<T> Handle<T>::colored(std::function<Vec<T, 3>(const Vec<T, 3>&, T)> color_fn) const {
    trace->node(index).color_fn = std::move(color_fn);
    return *this;
}

template<Scalar T>
Handle<T> Handle<T>::moving(std::function<Vec<T, 3>(const Vec<T, 3>&, T)> f) const {
    trace->node(index).transform = std::move(f);
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
    auto m = materialize_mesh(trace, idx, t);
    Material<T> mat = n.material;
    if (n.color_fn) mat.base_color = n.color_fn(m.centroid(), t);
    return {Placed<T>{std::move(m), mat}};
}

} // namespace spatium::io::build
