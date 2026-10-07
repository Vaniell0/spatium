#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/mesh/mesh.hpp>
#  include <spatium/mesh/geodesic.hpp>
#  include <cmath>
#endif

SPATIUM_EXPORT namespace spatium::mesh {

// Parallel transport of a tangent vector along a geodesic path
// using Schild's ladder. Works on any Manifold (uses only exp_map + log_map).
//
// Schild's ladder is parallel transport in the limit of a small vector: its error
// grows with the size of the vector (a unit vector laddered over a quarter of a great
// circle came out 46% short, a vector of length 2 not even pointing the right way --
// findings 16). Transport is linear, so the vector is scaled to a small size h,
// laddered, and scaled back -- which leaves a relative error of order h, the ladder's
// h^2 truncation divided by the scale. A second ladder at h/2 and the extrapolation
// 2 T(h/2) - T(h) remove that term, leaving h^2. The other side is the rounding of the
// exp and log the ladder is built from: a log between points h apart that goes through
// acos (a sphere's does) is good to eps / h in the angle, eps / h^2 relative to the
// vector. The two balance at h = eps^(1/4), about 1e-8 relative in double -- the floor
// of this method until those logs are stable for close points (stage 6, precision).

template<Surface S>
    requires Manifold<S>
typename S::TangentVector parallel_transport(
    const S& space,
    const Mesh<S>& mesh,
    const GeodesicPath<S>& path,
    const typename S::TangentVector& vector)
{
    using T = typename S::ScalarType;
    using Tangent = typename S::TangentVector;

    if (path.vertices.size() < 2) return vector;

    // Zero vector stays zero — skip the ladder
    auto eps = spatium::epsilon<T>();
    if (vector.norm_squared() < eps * eps)
        return vector;

    // The ladder along the whole path for a vector of length `size` in the direction of
    // `vector`, scaled back to the length of `vector`.
    using std::sqrt;
    const T length = T(sqrt(vector.norm_squared()));
    const auto ladder = [&](T size) -> Tangent {
        const T scale = T(size / length);
        Tangent v = vector * scale;       // not auto: that would be the expression, not the vector

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
        return Tangent{v * T(T{1} / scale)};
    };

    const T h = T(std::pow(machine_epsilon<T>(), 1.0 / 4.0));
    const Tangent coarse = ladder(h);
    const Tangent fine = ladder(T(h / T{2}));
    return Tangent{fine * T{2} - coarse};
}

} // namespace spatium::mesh
