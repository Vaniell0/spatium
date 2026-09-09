#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/mesh/mesh.hpp>
#  include <spatium/spaces/euclidean.hpp>
#  include <functional>
#endif

SPATIUM_EXPORT namespace spatium::mesh {

// conform_to_surface: drape a guide mesh onto a target Surface.
//
// Every guide vertex is projected onto `target` (target.project()) and
// pushed out along the target's local normal by `thickness` (a per-point
// field or a constant) -- topology (faces) carried over unchanged. This is
// the "icing on the donut" primitive: build the guide shape in whatever
// space is convenient (its own parametric patch, a flat grid, ...), then
// wrap it onto any Surface -- torus, sphere, implicit blob -- with one
// call. Output lives in ambient R^3 (Euclidean<3,T>), since the draped
// result in general no longer satisfies the guide's own Surface concept.

template<Surface GuideS, class Target>
Mesh<Euclidean<3, typename GuideS::ScalarType>>
conform_to_surface(const Mesh<GuideS>& guide, const Target& target,
                    std::function<typename GuideS::ScalarType(const Vec<typename GuideS::ScalarType, 3>&)> thickness)
{
    using T = typename GuideS::ScalarType;
    Mesh<Euclidean<3, T>> out;
    out.vertices.reserve(guide.vertices.size());
    for (const auto& v : guide.vertices) {
        Vec<T, 3> p = target.project(Vec<T, 3>{v});
        Vec<T, 3> n = target.normal(p);
        out.vertices.push_back(Vec<T, 3>{p + n * thickness(p)});
    }
    out.faces = guide.faces;
    return out;
}

template<Surface GuideS, class Target>
Mesh<Euclidean<3, typename GuideS::ScalarType>>
conform_to_surface(const Mesh<GuideS>& guide, const Target& target,
                    typename GuideS::ScalarType thickness = typename GuideS::ScalarType{0})
{
    return conform_to_surface(guide, target,
        std::function<typename GuideS::ScalarType(const Vec<typename GuideS::ScalarType, 3>&)>{
            [thickness](const Vec<typename GuideS::ScalarType, 3>&) { return thickness; }});
}

} // namespace spatium::mesh
