#pragma once
// A cooked scene laid out for a ray tracer: world-space triangles in one
// tree, instanced closed forms in another.
//
// This lived inside `donut_demo.cpp` until a second renderer needed it.
// The GPU path has to see exactly the scene the CPU path sees, and two
// copies of "how a Cooked becomes triangles" would disagree the first
// time one of them changed -- the same way `materialize_mesh()` and
// `cook()` once disagreed about every scattered object's orientation
// because nothing compared them.
//
// Not part of the `spatium.render` module, for the reason `io/build.hpp`
// is not part of `spatium.io`: it sits on the scene DSL, which is
// header-only outside the module build.
#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/geometry/ray_hit.hpp>
#  include <spatium/geometry/ray_surface.hpp>
#  include <spatium/geometry/triangle.hpp>
#  include <spatium/io/build.hpp>
#  include <spatium/mesh/mesh.hpp>
#  include <spatium/spaces/euclidean.hpp>
#  include <spatium/spatial/bvh.hpp>
#  include <any>
#  include <array>
#  include <cstddef>
#  include <cstdint>
#  include <utility>
#  include <vector>
#endif

namespace spatium::render {

// Smooth normals per face *corner*, not per vertex, so an edge can stay
// sharp. Flat per-triangle normals on a coarse tessellation read as
// faceted plastic, and averaging every face around a vertex -- what this
// did until the cube's diagonal bands were noticed -- smooths across edges
// that are meant to be there: a cube's corner normal came out pointing
// along the diagonal, and each face shaded as a gradient between three.
//
// So a corner averages only the faces around its vertex that lie within
// `crease_cos` of its own face (60 degrees by default). A smooth surface's
// neighbours are all well inside that, and for it the result is the
// area-weighted average it always was, summed in the same face order and
// therefore to the same bits; a cube's are at 90 degrees and drop out.
template<Scalar T = double>
std::vector<std::array<Vec<T, 3>, 3>> corner_normals(const mesh::Mesh<Euclidean<3, T>>& m,
                                                     T crease_cos = T{0.5}) {
    const std::size_t F = m.faces.size(), V = m.vertex_count();
    std::vector<Vec<T, 3>> area_n(F), unit_n(F);
    for (std::size_t f = 0; f < F; ++f) {
        const auto& a = m.vertices[m.faces[f][0]];
        const auto& b = m.vertices[m.faces[f][1]];
        const auto& c = m.vertices[m.faces[f][2]];
        area_n[f] = Vec<T, 3>{Vec<T, 3>{b - a}.cross(Vec<T, 3>{c - a})};  // area-weighted
        const T len = area_n[f].norm();
        unit_n[f] = len > T{0} ? Vec<T, 3>{area_n[f] / len} : Vec<T, 3>{};
    }

    // Faces around each vertex, as one flat array with offsets, filled in
    // face order so the sums below run in the order the old per-vertex
    // accumulation did.
    std::vector<std::uint32_t> start(V + 1, 0), around;
    for (const auto& face : m.faces)
        for (auto v : face) ++start[v + 1];
    for (std::size_t v = 0; v < V; ++v) start[v + 1] += start[v];
    around.resize(start[V]);
    std::vector<std::uint32_t> fill(start.begin(), start.end() - 1);
    for (std::size_t f = 0; f < F; ++f)
        for (auto v : m.faces[f]) around[fill[v]++] = static_cast<std::uint32_t>(f);

    std::vector<std::array<Vec<T, 3>, 3>> out(F);
    for (std::size_t f = 0; f < F; ++f)
        for (std::size_t k = 0; k < 3; ++k) {
            const auto v = m.faces[f][k];
            Vec<T, 3> sum{0, 0, 0};
            for (auto i = start[v]; i < start[v + 1]; ++i) {
                const auto g = around[i];
                if (unit_n[g].dot(unit_n[f]) >= crease_cos) sum = Vec<T, 3>{sum + area_n[g]};
            }
            const T len = sum.norm();
            out[f][k] = len > T{1e-12} ? Vec<T, 3>{sum / len} : Vec<T, 3>{0, 0, 1};
        }
    return out;
}

// A triangle as the shader sees it. No opacity, and that is the current
// behaviour rather than an omission made here: only instances are ever
// see-through.
template<Scalar T = double>
struct ShadedTriangle {
    geometry::Triangle<3, T> tri;
    std::array<Vec<T, 3>, 3> vertex_normals;
    Vec<T, 3> color;
    T roughness;
    Vec<T, 3> emissive;
};

// What an instance looks like, parallel to the instance tree's input.
template<Scalar T = double>
struct InstanceLook {
    Vec<T, 3> color;
    T roughness;
    Vec<T, 3> emissive;
    T opacity;
};

// Two trees, because the scene has two kinds of object in it.
//
// Most of it is ordinary geometry, each piece different, and a triangle
// tree is the right home for it. The rest is many copies of one closed
// form, and an instance holds a *reference* to the shape plus where this
// copy sits, so the tree holds one leaf per copy rather than one per
// triangle of every copy.
//
// What makes an object eligible is not that it looks small. It is that
// its motion is a *placement* -- affine in the point, so it moves the
// object without reshaping it (see VecField::is_placement) -- and that
// its shape carries a closed form.
//
// Move-only: every instance points into `quadrics`. Moving a vector keeps
// its buffer, copying does not, so a copy would leave the instance tree
// pointing into the original.
template<Scalar T = double>
struct CookedScene {
    std::vector<ShadedTriangle<T>> prims;          // parallel to the triangle tree's input
    std::vector<geometry::BoundedQuadric<T>> quadrics;  // one slot per shape, never grown
    std::vector<InstanceLook<T>> looks;            // parallel to the instance tree's input
    spatial::BVH<geometry::Triangle<3, T>> triangles;
    spatial::BVH<geometry::Instanced<geometry::BoundedQuadric<T>>> instances;

    CookedScene() = default;
    CookedScene(CookedScene&&) noexcept = default;
    CookedScene& operator=(CookedScene&&) noexcept = default;
    CookedScene(const CookedScene&) = delete;
    CookedScene& operator=(const CookedScene&) = delete;
};

// Lay a cooked scene out for tracing. `extra` is geometry that is drawn
// but is not in the scene -- the build-up frames' gizmos are the editor,
// not the subject, and arrive as Placed views onto their own Trace.
//
// `skip_nodes`, when given, is indexed by trace node: objects from a node
// marked in it are left out -- the instances a device kernel moves itself,
// which would otherwise be laid out twice.
template<Scalar T = double>
CookedScene<T> lay_out(const io::build::Cooked<T>& cooked,
                       const std::vector<io::build::Placed<T>>& extra = {},
                       const std::vector<char>* skip_nodes = nullptr) {
    using Tri = geometry::Triangle<3, T>;
    using Quadric = geometry::BoundedQuadric<T>;

    CookedScene<T> out;
    std::vector<Tri> tris;
    std::vector<geometry::Instanced<Quadric>> insts;

    // One quadric per distinct *shape*, not per object -- which is the
    // whole point of reading a cooked scene rather than a list of
    // materialized nodes. Sized up front and never grown, because every
    // instance below holds a pointer into it.
    out.quadrics.resize(cooked.shape_count());
    std::vector<char> has_quadric(cooked.shape_count(), 0);
    std::vector<std::vector<std::array<Vec<T, 3>, 3>>> shape_normals(cooked.shape_count());
    for (std::size_t i = 0; i < cooked.shape_count(); ++i) {
        const auto& sh = cooked.shapes()[i];
        if (const auto* q = std::any_cast<Quadric>(&sh.exact)) {
            out.quadrics[i] = *q;
            has_quadric[i] = 1;
        }
    }

    // A shape's geometry is its *rest* form when the object placing it is
    // instanceable, and its *placed* form when the object deforms -- in
    // which case cook() hands back an identity transform. So one formula
    // covers both, and the deforming case is not a special case here.
    // Takes the rotation already expanded, rather than reaching into the
    // object for it: this runs per *vertex*, and a quaternion unpacked
    // here would spend the arithmetic the compact storage was supposed to
    // be paying for. The expansion happens once per object, below.
    auto to_world = [](const Vec<T, 3>& v, const Matrix<T, 3, 3>& R,
                       const io::build::Object<T>& o) {
        return Vec<T, 3>{R * Vec<T, 3>{v * o.scale} + o.translation};
    };

    auto emit_triangles = [&](const mesh::Mesh<Euclidean<3, T>>& m,
                              const std::vector<std::array<Vec<T, 3>, 3>>& cn,
                              const io::build::Object<T>& o, const Matrix<T, 3, 3>& R,
                              const io::Material<T>& mat) {
        for (std::size_t i = 0; i < m.faces.size(); ++i) {
            const auto& f = m.faces[i];
            Tri t(to_world(m.vertices[f[0]], R, o), to_world(m.vertices[f[1]], R, o),
                  to_world(m.vertices[f[2]], R, o));
            tris.push_back(t);
            // Normals turn by the rotation alone: the scale is uniform, so
            // it divides out of the inverse transpose.
            out.prims.push_back(ShadedTriangle<T>{t,
                                                  {Vec<T, 3>{R * cn[i][0]},
                                                   Vec<T, 3>{R * cn[i][1]},
                                                   Vec<T, 3>{R * cn[i][2]}},
                                                  mat.base_color, mat.roughness, mat.emissive});
        }
    };

    for (const auto& obj : cooked.objects()) {
        if (skip_nodes && obj.source_node < skip_nodes->size() && (*skip_nodes)[obj.source_node])
            continue;
        const auto& mat = obj.material;
        // Once per object. Everything below uses this, including the
        // instance, whose rotation stays a matrix precisely so the
        // traversal never has to unpack anything.
        const Matrix<T, 3, 3> R = obj.rotation_q.to_matrix();

        // Scaled to nothing is not there. The donut's dust dissolves by
        // scaling to zero, so at t=3.9 35 200 of its 46 192 objects were
        // still in the instance tree as points: never hit -- `ray_hit`
        // refuses a zero scale -- and still walked, 31 nodes a pixel
        // against 4.6 for a scene without them.
        if (obj.scale == T{0}) continue;

        // Two conditions, and neither is a guess about what the object
        // looks like: the shape must carry a closed form, and the object's
        // motion must be a placement so copies differ only by where they
        // are.
        if (obj.instanceable && has_quadric[obj.shape]) {
            insts.push_back({&out.quadrics[obj.shape], obj.translation, obj.scale, R});
            out.looks.push_back({mat.base_color, mat.roughness, mat.emissive, mat.opacity});
            continue;
        }

        // Smooth normals once per shape rather than once per object: a
        // thousand instances of one sprinkle used to recompute the same
        // normals a thousand times, through a merged mesh that had to be
        // built first.
        auto& vn = shape_normals[obj.shape];
        if (vn.empty()) vn = corner_normals(cooked.shapes()[obj.shape].geometry);
        emit_triangles(cooked.shapes()[obj.shape].geometry, vn, obj, R, mat);
    }

    for (const auto& g : extra) {
        auto m = g.mesh();
        auto cn = corner_normals(m);
        auto mat = g.material();
        for (std::size_t i = 0; i < m.faces.size(); ++i) {
            const auto& f = m.faces[i];
            Tri t(m.vertices[f[0]], m.vertices[f[1]], m.vertices[f[2]]);
            tris.push_back(t);
            out.prims.push_back(ShadedTriangle<T>{t, cn[i], mat.base_color, mat.roughness,
                                                  mat.emissive});
        }
    }

    // Moved, not copied. `build` takes its shapes by value and keeps
    // them, so handing it an lvalue leaves two copies of the array alive
    // for the rest of the frame -- 112 bytes per instance, which is 214 MB
    // at two million and was simply being spent.
    out.triangles = spatial::BVH<Tri>::build(std::move(tris));
    out.instances = spatial::BVH<geometry::Instanced<Quadric>>::build(std::move(insts));
    return out;
}

}  // namespace spatium::render
