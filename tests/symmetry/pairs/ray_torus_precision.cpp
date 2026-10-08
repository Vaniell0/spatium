#include "../pair.hpp"
#include "../common.hpp"

#include <spatium/core/precision.hpp>
#include <spatium/geometry/ray_surface.hpp>

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

// ── ray_torus: the closed form on double against the same code on Real50 ──
//
// A ray from a point of the box around the torus (R = 2, r = 0.5), aimed at a
// point inside the torus's own box, so most rays hit and some graze. The
// reference is the same function on fifty digits; the route is `double`. The
// bound is 1e-9 of the scale of the ray. A hit with two coinciding roots could
// only promise sqrt(eps) = 1.5e-8, but none of the 6 000 rays is that close, the
// worst distance is 8e-13, and a bound nearer the observation is what lets the
// pair see a fault of 1e-6.
//
// First run (20 000 rays, 15 846 with hits): `ray_torus<Real50>` did not compile
// (a qualified std::abs in torus_basis), and once it did, one ray (the 5 037th)
// had both hits off by 4.3e-4 while the other roots were 3.8 away. The cause was
// the closed forms' absolute branch tests; the roots are now polished by Newton
// steps on the polynomial itself, and the worst error over the 20 000 is 2e-11.

inline Entry ray_torus_precision(std::size_t n = 6000) {
    return {"ray_torus: double against Real50", n, [n] {
        return check<TorusRay>(
            "ray_torus: double against Real50", Kind::Function, n, torus_ray,
            torus_hits_real50, torus_hits_double,
            [](const std::vector<double>& a, const std::vector<double>& b) {
                if (a.size() != b.size()) return std::numeric_limits<double>::infinity();          // a hit gained or lost is not a rounding
                double worst = 0.0;
                for (std::size_t k = 0; k < a.size(); ++k) worst = std::max(worst, std::abs(a[k] - b[k]));
                return worst;
            },
            [](const TorusRay& s) { return 1e-9 * (1.0 + s.o.norm()); });
    }};
}

}  // namespace

const Registrar registered_ray_torus_precision{"ray_torus_precision", 6000, [](std::size_t n) { return ray_torus_precision(n); }};

}  // namespace symmetry
