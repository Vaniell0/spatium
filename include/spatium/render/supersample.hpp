#pragma once

// CPU-raytracer supersampling (antialiasing).
//
// Moved into the engine (2026-08-26) after three independent example
// renderers (blackhole_demo, wormhole_demo, parametric_analytical_demo)
// each fired exactly one ray per pixel and each separately grew visible
// aliasing symptoms from it -- jagged horizon/silhouette edges, a
// fragmented photon-ring, aliased spiral-arm structure -- that were
// getting patched one at a time as they were noticed. The actual missing
// piece was general: nothing averaged multiple samples per pixel. That
// belongs here, not duplicated per demo, since any future CPU raytracer
// built on Spatium's ray-surface primitives (ray_parametric, ray_quadric,
// ray_torus, or a caller's own integrator) needs the same fix for the
// same reason.
//
// Ordinary NxN jittered-grid supersampling (SSAA): the caller's own
// per-ray color function is called aa*aa times per output pixel, on a
// jittered sub-pixel grid, and averaged in linear space before a single
// clamp -- silhouette edges, disk/ring boundaries, and thin high-
// frequency detail all blend smoothly instead of hard-aliasing,
// automatically, for any scene.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <algorithm>
#  include <cstdint>
#  include <utility>
#endif

SPATIUM_EXPORT namespace spatium::render {

inline constexpr int kDefaultAA = 2;  // 2x2 = 4 rays/pixel

// `ray_color(sx, sy)` receives normalized screen coordinates in the
// tan(fov/2)-scaled NDC convention a pinhole-camera CPU raytracer
// computes per pixel:
//   sx = (2*(x+jitter)/width  - 1) * tan_fov * aspect
//   sy = -(2*(y+jitter)/height - 1) * tan_fov
// and returns a linear color in [0, 255] scale *before* clamping --
// clamping happens once here, after averaging every subsample, not per
// sample.
// HDR variant: the same jittered-grid average, WITHOUT the final clamp
// -- for a caller that needs extended range before its own tone-mapping
// (e.g. a bloom pass). Clamping before bloom is a real bug, not a
// simplification: every sufficiently bright pixel flattens to the same
// clamped value first, so bloom can only blur an already-flat block
// instead of a true HDR highlight (found in blackhole_gr_demo.cpp,
// 2026-08-26, as the cause of a solid-white blowout with no internal
// structure at close camera range).
template<class RayColorFn>
Vec<double, 3> supersample_pixel_hdr(int x, int y, int width, int height, double tan_fov,
                                      double aspect, RayColorFn&& ray_color, int aa = kDefaultAA) {
    Vec<double, 3> accum{0.0, 0.0, 0.0};
    for (int j = 0; j < aa; ++j) {
        double jy = (j + 0.5) / aa;
        double sy = -(2.0 * (y + jy) / height - 1.0) * tan_fov;
        for (int i = 0; i < aa; ++i) {
            double jx = (i + 0.5) / aa;
            double sx = (2.0 * (x + jx) / width - 1.0) * tan_fov * aspect;
            accum = Vec<double, 3>{accum + ray_color(sx, sy)};
        }
    }
    return Vec<double, 3>{accum * (1.0 / (aa * aa))};
}

template<class RayColorFn>
void supersample_pixel(int x, int y, int width, int height, double tan_fov, double aspect,
                        RayColorFn&& ray_color, std::uint8_t px[3], int aa = kDefaultAA) {
    Vec<double, 3> c =
        supersample_pixel_hdr(x, y, width, height, tan_fov, aspect,
                               std::forward<RayColorFn>(ray_color), aa);
    px[0] = static_cast<std::uint8_t>(std::clamp(c[0], 0.0, 255.0));
    px[1] = static_cast<std::uint8_t>(std::clamp(c[1], 0.0, 255.0));
    px[2] = static_cast<std::uint8_t>(std::clamp(c[2], 0.0, 255.0));
}

}  // namespace spatium::render
