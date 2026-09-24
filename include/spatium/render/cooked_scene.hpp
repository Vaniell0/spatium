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
#  include <utility>
#  include <vector>
#endif

namespace spatium::render {

// Smooth per-vertex normals (area-weighted face-normal average) -- flat
// per-triangle normals on a coarse tessellation read as faceted plastic;
// this is the standard fix, needs nothing beyond the mesh's own topology.
template<Scalar T = double>
std::vector<Vec<T, 3>> smooth_normals(const mesh::Mesh<Euclidean<3, T>>& m) {
    std::vector<Vec<T, 3>> n(m.vertex_count(), Vec<T, 3>{0, 0, 0});
    for (const auto& f : m.faces) {
        auto& a = m.vertices[f[0]];
        auto& b = m.vertices[f[1]];
        auto& c = m.vertices[f[2]];
        Vec<T, 3> face_n{Vec<T, 3>{b - a}.cross(Vec<T, 3>{c - a})}; // area-weighted (unnormalized)
        n[f[0]] = Vec<T, 3>{n[f[0]] + face_n};
        n[f[1]] = Vec<T, 3>{n[f[1]] + face_n};
        n[f[2]] = Vec<T, 3>{n[f[2]] + face_n};
    }
    for (auto& v : n) {
        T len = v.norm();
        v = len > T{1e-12} ? Vec<T, 3>{v / len} : Vec<T, 3>{0, 0, 1};
    }
    return n;
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
template<Scalar T = double>
CookedScene<T> lay_out(const io::build::Cooked<T>& cooked,
                       const std::vector<io::build::Placed<T>>& extra = {}) {
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
    std::vector<std::vector<Vec<T, 3>>> shape_normals(cooked.shape_count());
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
                              const std::vector<Vec<T, 3>>& vn,
                              const io::build::Object<T>& o, const Matrix<T, 3, 3>& R,
                              const io::Material<T>& mat) {
        for (const auto& f : m.faces) {
            Tri t(to_world(m.vertices[f[0]], R, o), to_world(m.vertices[f[1]], R, o),
                  to_world(m.vertices[f[2]], R, o));
            tris.push_back(t);
            // Normals turn by the rotation alone: the scale is uniform, so
            // it divides out of the inverse transpose.
            out.prims.push_back(ShadedTriangle<T>{t,
                                                  {Vec<T, 3>{R * vn[f[0]]},
                                                   Vec<T, 3>{R * vn[f[1]]},
                                                   Vec<T, 3>{R * vn[f[2]]}},
                                                  mat.base_color, mat.roughness, mat.emissive});
        }
    };

    for (const auto& obj : cooked.objects()) {
        const auto& mat = obj.material;
        // Once per object. Everything below uses this, including the
        // instance, whose rotation stays a matrix precisely so the
        // traversal never has to unpack anything.
        const Matrix<T, 3, 3> R = obj.rotation_q.to_matrix();

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
        if (vn.empty()) vn = smooth_normals(cooked.shapes()[obj.shape].geometry);
        emit_triangles(cooked.shapes()[obj.shape].geometry, vn, obj, R, mat);
    }

    for (const auto& g : extra) {
        auto m = g.mesh();
        auto vn = smooth_normals(m);
        auto mat = g.material();
        for (const auto& f : m.faces) {
            Tri t(m.vertices[f[0]], m.vertices[f[1]], m.vertices[f[2]]);
            tris.push_back(t);
            out.prims.push_back(ShadedTriangle<T>{t, {vn[f[0]], vn[f[1]], vn[f[2]]},
                                                  mat.base_color, mat.roughness, mat.emissive});
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
