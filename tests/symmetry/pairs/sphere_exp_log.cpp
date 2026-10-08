#include "../pair.hpp"

#include <spatium/core/access.hpp>
#include <spatium/spaces/sphere.hpp>

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

// ── exp then log on the sphere: a relation, not a second path ──
//
// log_p(exp_p(v)) = v for |v| inside the injectivity radius (pi r). Checked on tangent
// vectors of length up to 2.75 r, away from the cut locus where acos loses its footing.
// The bound is 1e-11 r: acos near +-1 is what limits it and the vectors are kept off both,
// and the worst distance measured over 2 000 draws is 2e-14.
struct SphereVector { Vec<double, 3> p, v; };

inline SphereVector sphere_vector(std::size_t i, double radius) {
    const double w = 2 * halton(i, 0) - 1, phi = 2 * std::numbers::pi * halton(i, 1);
    const double rho = std::sqrt(std::max(0.0, 1 - w * w));
    const Vec<double, 3> n{rho * std::cos(phi), rho * std::sin(phi), w};
    const Vec<double, 3> a{2 * halton(i, 2) - 1, 2 * halton(i, 3) - 1, 2 * halton(i, 4) - 1};
    Vec<double, 3> t = a - n * a.dot(n);
    const double length = radius * (0.05 + 2.7 * halton(i, 5));
    t = t * (length / t.norm());
    return {n * radius, t};
}

inline Entry sphere_exp_log(std::size_t n = 2000) {
    return {"Sphere: log after exp is the identity", n, [n] {
        const spatium::Sphere<2> s{.radius = 2.0};
        return check<SphereVector>(
            "Sphere: log after exp is the identity", Kind::Relation, n,
            [](std::size_t i) { return sphere_vector(i, 2.0); },
            [](const SphereVector& x) { return x.v; },
            [&](const SphereVector& x) {
                const auto q = spatium::spaces::exp_map(s, x.p, x.v, 1.0);
                return Vec<double, 3>{spatium::spaces::log_map(s, x.p, q)};
            },
            [](const Vec<double, 3>& a, const Vec<double, 3>& b) { return (a - b).norm(); },
            [](const SphereVector&) { return 1e-11 * 2.0; });
    }};
}

}  // namespace

const Registrar registered_sphere_exp_log{"sphere_exp_log", 2000, [](std::size_t n) { return sphere_exp_log(n); }};

}  // namespace symmetry
