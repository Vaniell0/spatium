#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/mesh/mesh.hpp>
#  include <spatium/mesh/geodesic.hpp>
#endif

SPATIUM_EXPORT namespace spatium::mesh {

// Parallel transport of a tangent vector along a geodesic path
// using Schild's ladder. Works on any Manifold (uses only exp_map + log_map).

template<Surface S>
    requires Manifold<S>
typename S::TangentVector parallel_transport(
    const S& space,
    const Mesh<S>& mesh,
    const GeodesicPath<S>& path,
    const typename S::TangentVector& vector)
{
    using T = typename S::ScalarType;

    if (path.vertices.size() < 2) return vector;

    auto v = vector;

    // Zero vector stays zero — skip the ladder
    auto eps = spatium::epsilon<T>();
    if (v.norm_squared() < eps * eps)
        return v;

    for (std::size_t i = 0; i + 1 < path.vertices.size(); ++i) {
        auto& p = mesh.vertices[path.vertices[i]];
        auto& q = mesh.vertices[path.vertices[i + 1]];

        // Schild's ladder:
        // 1. x = exp(p, v) — tip of vector at p
        auto x = space.exp_map(p, v, T{1});

        // 2. c = midpoint of geodesic from x to q
        auto c = space.exp_map(x, space.log_map(x, q), T{0.5});

        // 3. x' = reflect p through c: exp(p, 2 * log(p, c))
        auto xprime = space.exp_map(p, space.log_map(p, c), T{2});

        // 4. v' = log(q, x') — transported vector at q
        v = space.log_map(q, xprime);
    }

    return v;
}

} // namespace spatium::mesh
