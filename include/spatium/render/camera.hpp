#pragma once

// Pinhole camera for CPU raytracers.
//
// Every CPU raytracer example in the tree (parametric_analytical_demo.cpp,
// collatz_demo.cpp, blackhole_demo.cpp, wormhole_demo.cpp) independently
// reimplemented the same position/target/up/fov_deg -> orthonormal basis
// -> per-pixel ray direction math -- found duplicated 4 separate times
// during the 2026-08-26 architecture audit (`project_spatium_
// architecture_audit` memory). None of it depends on anything example-
// specific, so it belongs here, not copy-pasted per demo.
//
// Two screen-coordinate conventions were in use and are actually
// identical: `ny = (1 - 2t) * tan_half` (parametric_analytical_demo.cpp)
// and `sy = -(2t - 1) * tan_fov` (blackhole_demo.cpp/wormhole_demo.cpp)
// -- `1-2t == -(2t-1)` algebraically, so `camera_pixel_dir` below serves
// both call shapes with one implementation.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/concepts.hpp>
#  include <cmath>
#  include <numbers>
#endif

SPATIUM_EXPORT namespace spatium::render {

template<Scalar T = double>
struct Camera {
    Vec<T, 3> position;
    Vec<T, 3> target;
    Vec<T, 3> up{0, 0, 1};
    T fov_deg = T{38};
};

// Precomputed orthonormal basis + half-FOV tangent -- built once per
// frame, not per ray.
template<Scalar T>
struct CameraBasis {
    Vec<T, 3> fwd, right, up;
    T tan_half;
};

template<Scalar T>
CameraBasis<T> make_camera_basis(const Camera<T>& cam) {
    Vec<T, 3> fwd = Vec<T, 3>{(cam.target - cam.position).normalized()};
    Vec<T, 3> right = Vec<T, 3>{fwd.cross(cam.up).normalized()};
    Vec<T, 3> up = Vec<T, 3>{right.cross(fwd).normalized()};
    T tan_half = std::tan(cam.fov_deg * std::numbers::pi_v<T> / T{360});
    return CameraBasis<T>{fwd, right, up, tan_half};
}

// Ray direction for normalized screen coordinates (nx, ny) -- the
// tan(fov/2)-scaled NDC convention every Spatium CPU raytracer uses.
template<Scalar T>
Vec<T, 3> camera_ray_dir(const CameraBasis<T>& basis, T nx, T ny) {
    return Vec<T, 3>{(basis.fwd + basis.right * nx + basis.up * ny).normalized()};
}

// Ray direction for pixel (px, py) of a width x height image, with
// optional sub-pixel jitter (jx, jy in [0,1); default 0.5 = pixel
// center) -- the same jitter hook `render::supersample_pixel()` needs.
template<Scalar T>
Vec<T, 3> camera_pixel_dir(const Camera<T>& cam, const CameraBasis<T>& basis, int px, int py,
                           int width, int height, T jx = T{0.5}, T jy = T{0.5}) {
    T aspect = static_cast<T>(width) / static_cast<T>(height);
    T nx = (T{2} * (static_cast<T>(px) + jx) / static_cast<T>(width) - T{1}) * aspect *
           basis.tan_half;
    T ny = (T{1} - T{2} * (static_cast<T>(py) + jy) / static_cast<T>(height)) * basis.tan_half;
    return camera_ray_dir(basis, nx, ny);
}

}  // namespace spatium::render
