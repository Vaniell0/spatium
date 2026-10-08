#pragma once

// What more than one pair uses: the rays through the torus's box (the closed form against fifty
// digits, and against the chart), and a point from a draw of numbers.

#include "pair.hpp"

#include <spatium/core/precision.hpp>
#include <spatium/geometry/ray_surface.hpp>

#include <vector>

namespace symmetry {

using spatium::Vec;

struct TorusRay { Vec<double, 3> o, d; };

inline TorusRay torus_ray(std::size_t i) {
    const auto u = [&](std::size_t k) { return halton(i, k); };
    Vec<double, 3> o{-6 + 12 * u(0), -6 + 12 * u(1), -2 + 4 * u(2)};
    const Vec<double, 3> target{-2.6 + 5.2 * u(3), -2.6 + 5.2 * u(4), -0.6 + 1.2 * u(5)};
    Vec<double, 3> d = target - o;
    d = d / d.norm();
    return {o, d};
}

inline std::vector<double> torus_hits_double(const TorusRay& s) {
    using namespace spatium::geometry;
    Torus<double> torus{.major_radius = 2.0, .minor_radius = 0.5};
    std::vector<double> t;
    for (const auto& h : ray_torus(Ray<3, double>{s.o, s.d}, torus)) t.push_back(h.t);
    return t;
}

inline std::vector<double> torus_hits_real50(const TorusRay& s) {
    using namespace spatium::geometry;
    using R = spatium::Real50;
    Torus<R> torus{.major_radius = R(2.0), .minor_radius = R(0.5)};
    const Vec<R, 3> o{R(s.o[0]), R(s.o[1]), R(s.o[2])};
    const Vec<R, 3> d{R(s.d[0]), R(s.d[1]), R(s.d[2])};
    std::vector<double> t;
    for (const auto& h : ray_torus(Ray<3, R>{o, d}, torus)) t.push_back(spatium::primal_double(h.t));
    return t;
}

template<class T>
Vec<T, 3> geometry_point(const double* a) { return Vec<T, 3>{T(a[0]), T(a[1]), T(a[2])}; }

}  // namespace symmetry
