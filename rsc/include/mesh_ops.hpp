#pragma once

// Domain 5: mesh-strategy dispatch (uniform UV step vs. an
// anisotropy-adapted UV step) for a parametrized-surface patch -- the
// general principle behind the Klein-bottle seam mesh-twist bug (see
// examples/primitives_demo.cpp's KleinBottle comment): a UV-uniform step
// maps to wildly different R^3 edge lengths wherever the parametrization's
// local anisotropy is high, producing long, thin ("sliver") triangles.
// spatium::ParametricSurface::parametrization_anisotropy() (first
// fundamental form eigenvalue ratio) is the cheap, single-evaluation
// feature that predicts this without ever building a mesh -- same
// "recognize the delicate case before paying for it" shape as
// rootfind_ops.hpp's |f'(x0)| and ode_ops.hpp's freq*dt.
//
// Real, measured correction made while building this: the first design
// used face-*normal* angular deviation (uniform-mesh normal vs. the
// analytic normal) as the ground-truth error, matching how
// subdivide_adaptive() already judges curvature. Measured directly (a
// standalone probe, not assumed): that error came out essentially
// IDENTICAL across torus aspect ratios from 1.2 to 30 -- a real,
// surprising fact, not a bug in the measurement. Reason: sweeping u
// around a torus rotates the tangent plane by exactly du radians per
// unit du regardless of major_r (it's a plain circle in the xy-plane),
// so normal-angle error is driven by *curvature*, not by the fu/fv
// *magnitude ratio* anisotropy actually measures -- two genuinely
// different quantities. Switched the ground-truth metric to ambient
// *edge-length aspect ratio* of the UV-uniform quad instead
// (spatium::mesh::face_aspect_ratio()'s own hi/lo convention, just
// computed directly from two surface evaluations) -- this one really is
// governed by |fu|/|fv|, confirmed by the same kind of standalone probe
// before committing to it (uniform_aspect empirically tracks
// parametrization_anisotropy's own range at every tested level).
//
// Two test families, chosen for two different real causes of high
// anisotropy, not the same knob twice:
//   - Torus: aspect ratio R/r -- a thin-ring torus has fu (around the big
//     loop) growing much faster than fv (around the tube) almost
//     everywhere.
//   - Cone: distance from the apex -- fu (around the rim) shrinks toward
//     zero near the apex while fv (up the slant) stays roughly constant,
//     a genuine degenerate point, not just a global ratio.

#include <registry.hpp>
#include <spatium/spaces/parametric.hpp>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

namespace rsc {

enum class MeshFamily { Torus = 0, Cone = 1 };

namespace mesh_detail {

using Surf = spatium::ParametricSurface<double>;

inline Surf make_surface(MeshFamily family, double shape_param) {
    switch (family) {
    case MeshFamily::Torus: return spatium::make_torus(shape_param, 1.0);
    case MeshFamily::Cone: return spatium::make_cone(1.0, shape_param);
    }
    return spatium::make_torus();
}

struct QuadProbe {
    double u0, v0, u1, v1;
};

// Ambient-space aspect ratio (>= 1, longer edge / shorter edge, same
// convention spatium::mesh::face_aspect_ratio() uses) of the two edges
// leaving (u0,v0) for a UV-uniform quad -- exactly what a plain,
// unadapted tessellation at this UV resolution produces.
inline double uniform_aspect(const Surf& surf, const QuadProbe& q) {
    auto p00 = surf.evaluate(q.u0, q.v0);
    auto p10 = surf.evaluate(q.u1, q.v0);
    auto p01 = surf.evaluate(q.u0, q.v1);
    double eu = (p10 - p00).norm();
    double ev = (p01 - p00).norm();
    double lo = std::min(eu, ev), hi = std::max(eu, ev);
    return lo > 1e-15 ? hi / lo : std::numeric_limits<double>::max();
}

// Anisotropy-adapted strategy: a genuinely different triangulation
// choice, not just "more compute" on the same one -- probe the local
// speed |fu|, |fv| once (one extra pair of surface evaluations, the real
// small added cost this candidate pays for), then pick per-axis step
// sizes du', dv' that (a) equalize the resulting ambient edge lengths
// (du'*speed_u == dv'*speed_v) and (b) keep the same total UV-area
// budget as the uniform quad (du'*dv' == (u1-u0)*(v1-v0)) -- two
// equations, two unknowns, closed form below.
inline double adapted_aspect(const Surf& surf, const QuadProbe& q) {
    constexpr double h = 1e-4;
    auto p0 = surf.evaluate(q.u0, q.v0);
    double speed_u = (surf.evaluate(q.u0 + h, q.v0) - p0).norm() / h;
    double speed_v = (surf.evaluate(q.u0, q.v0 + h) - p0).norm() / h;
    if (speed_u < 1e-12 || speed_v < 1e-12) return uniform_aspect(surf, q);

    double footprint = (q.u1 - q.u0) * (q.v1 - q.v0);
    double du_p = std::sqrt(footprint * speed_v / speed_u);
    double dv_p = std::sqrt(footprint * speed_u / speed_v);
    if (!std::isfinite(du_p) || !std::isfinite(dv_p) || du_p <= 0.0 || dv_p <= 0.0)
        return uniform_aspect(surf, q);

    return uniform_aspect(surf, QuadProbe{q.u0, q.v0, q.u0 + du_p, q.v0 + dv_p});
}

} // namespace mesh_detail

// Registered for real usable dispatch. in = [family_id, shape_param, u0,
// v0, u1, v1], out = [aspect_ratio].
inline Registry build_mesh_registry() {
    Registry reg;

    auto make_fn = [](auto strategy) {
        return [strategy](std::span<const double> in, std::span<double> out) {
            auto family = static_cast<MeshFamily>(static_cast<int>(in[0]));
            auto surf = mesh_detail::make_surface(family, in[1]);
            mesh_detail::QuadProbe q{in[2], in[3], in[4], in[5]};
            out[0] = strategy(surf, q);
        };
    };

    reg.add({.name = "mesh_uniform",
             .tier = Tier::General,
             .in_size = 6,
             .out_size = 1,
             .input_names = {"family", "shape_param", "u0", "v0", "u1", "v1"},
             .output_names = {"aspect_ratio"}},
            make_fn([](const mesh_detail::Surf& s, const mesh_detail::QuadProbe& q) {
                return mesh_detail::uniform_aspect(s, q);
            }));

    reg.add({.name = "mesh_adapted",
             .tier = Tier::General,
             .in_size = 6,
             .out_size = 1,
             .input_names = {"family", "shape_param", "u0", "v0", "u1", "v1"},
             .output_names = {"aspect_ratio"}},
            make_fn([](const mesh_detail::Surf& s, const mesh_detail::QuadProbe& q) {
                return mesh_detail::adapted_aspect(s, q);
            }));

    return reg;
}

} // namespace rsc
