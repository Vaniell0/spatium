#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/spaces/parametric.hpp>
#  include <algorithm>
#  include <cstdint>
#  include <random>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium {

template<Scalar T = double>
struct SurfaceSample {
    T u, v;
    Vec<T, 3> position;
    Vec<T, 3> normal;
};

// sample_surface_uniform: place `count` points evenly by *surface area*
// (not by parameter (u,v)) across an analytic ParametricSurface, via
// rejection sampling weighted by the first fundamental form's area
// element sqrt(EG-F^2) -- the "sprinkles" primitive, fully analytic: no
// mesh, no geodesic-Voronoi graph. Sampling (u,v) uniformly instead would
// bunch points wherever the parametrization compresses ambient space
// (e.g. a thin-ring torus's inner rim); this doesn't.
template<Scalar T = double>
std::vector<SurfaceSample<T>> sample_surface_uniform(
    const ParametricSurface<T>& surface, std::size_t count, std::uint32_t seed = 42)
{
    auto [u_min, u_max, v_min, v_max] = surface.domain();

    // Coarse grid to bound the rejection envelope.
    constexpr int grid = 24;
    T max_area{0};
    for (int j = 0; j <= grid; ++j) {
        T v = v_min + (v_max - v_min) * static_cast<T>(j) / static_cast<T>(grid);
        for (int i = 0; i <= grid; ++i) {
            T u = u_min + (u_max - u_min) * static_cast<T>(i) / static_cast<T>(grid);
            max_area = std::max(max_area, surface.area_element(u, v));
        }
    }
    if (max_area <= T{0}) max_area = T{1};

    std::mt19937 rng(seed);
    std::uniform_real_distribution<T> du(u_min, u_max), dv(v_min, v_max), unit(T{0}, T{1});

    std::vector<SurfaceSample<T>> out;
    out.reserve(count);
    std::size_t guard = 0, guard_max = count * 200 + 10000;
    while (out.size() < count && guard < guard_max) {
        ++guard;
        T u = du(rng), v = dv(rng);
        if (unit(rng) * max_area <= surface.area_element(u, v))
            out.push_back({u, v, surface.evaluate(u, v), surface.normal_at(u, v)});
    }
    return out;
}

} // namespace spatium
