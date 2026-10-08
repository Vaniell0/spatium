#include "../pair.hpp"
#include "../common.hpp"

#include <spatium/geometry/make.hpp>
#include <spatium/geometry/ray_parametric.hpp>
#include <spatium/geometry/ray_surface.hpp>
#include <spatium/spaces/parametric.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <type_traits>
#include <vector>

namespace symmetry {
namespace {

using spatium::Vec;

// ── ray_torus: the quartic against the general path through the chart ──
//
// The closed form against the same torus as a ParametricSurface, hit by ray_parametric, a
// Newton search over the chart that knows nothing about quartics. Two paths that share no
// code, on the rays of the first pair. The route's own promise is its Newton tolerance, 1e-6;
// over the 3 000 rays the worst distance measured is 3.5e-7 and no hit is gained or lost.
inline Entry ray_torus_against_chart(std::size_t n = 3000) {
    return {"ray_torus: the quartic against the parametric search", n, [n] {
        const auto surface = spatium::make_torus<double>(2.0, 0.5);
        return check<TorusRay>(
            "ray_torus: the quartic against the parametric search", Kind::Function, n, torus_ray,
            torus_hits_double,
            [&](const TorusRay& s) {
                std::vector<double> t;
                for (const auto& h : spatium::geometry::ray_parametric(spatium::geometry::Ray<3, double>{s.o, s.d}, surface))
                    t.push_back(h.t);
                return t;
            },
            [](const std::vector<double>& a, const std::vector<double>& b) {
                if (a.size() != b.size()) return std::numeric_limits<double>::infinity();
                double worst = 0.0;
                for (std::size_t k = 0; k < a.size(); ++k) worst = std::max(worst, std::abs(a[k] - b[k]));
                return worst;
            },
            [](const TorusRay& s) { return 1e-6 * (1.0 + s.o.norm()); });
    }};
}

}  // namespace

const Registrar registered_ray_torus_chart{"ray_torus_against_chart", 3000, [](std::size_t n) { return ray_torus_against_chart(n); }};

}  // namespace symmetry
