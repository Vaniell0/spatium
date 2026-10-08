#include "../pair.hpp"
#include "../common.hpp"

#include <spatium/algebra/dual.hpp>
#include <spatium/geometry/geometry.hpp>

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

// ── The derivative of a distance to a primitive: Dual against central differences ──
//
// The gradient of the distance from a point to a segment, a triangle and a box with respect to
// the point, taken through `Dual` and by central differences with a step of 1e-6. Away from the
// places where the nearest feature switches the distance is smooth, and the two must agree to
// what a difference allows: 1e-8. At a switch the distance is C^1 but its second derivative jumps,
// and a difference whose step reaches across it is off by the jump times the step; and near a
// feature, at perpendicular distance a, the third derivative is of order 1/a^2 and the truncation
// grows with it. A sweep of 300 000 draws found two of each kind (9e-8 at worst). The difference
// quotient knows when it is unreliable: the bound adds three times the distance between the
// quotients at step h and h/2, so it opens exactly where the quotient is not to be trusted and
// stays 1e-8 everywhere else.
template<class T>
std::vector<double> distance_gradients(const std::vector<double>& u, const Vec<double, 3>& p0, int axis, double shift) {
    using namespace spatium::geometry;
    const double* d = u.data();
    Vec<T, 3> p{T(p0[0]), T(p0[1]), T(p0[2])};
    if constexpr (std::is_same_v<T, spatium::Dual<double>>) {
        p[axis] = T::variable(p0[axis]);
    } else {
        p[axis] = T(p0[axis] + shift);
    }
    const Segment<3, T> s{geometry_point<T>(d + 3), geometry_point<T>(d + 6)};
    const Triangle<3, T> t(geometry_point<T>(d + 15), geometry_point<T>(d + 18), geometry_point<T>(d + 21));
    const Vec<T, 3> lo = geometry_point<T>(d + 33);
    const Box<3, T> box{lo, lo + Vec<T, 3>{T(1), T(1), T(1)}};
    const T results[3] = {distance(p, s), distance(p, t), distance(p, box)};
    std::vector<double> out;
    for (const T& r : results) {
        if constexpr (std::is_same_v<T, spatium::Dual<double>>) out.push_back(r.deriv);
        else out.push_back(r);
    }
    return out;
}

inline std::vector<double> distance_difference_quotients(const std::vector<double>& u, double h) {
    std::vector<double> g;
    const Vec<double, 3> p{u[0], u[1], u[2]};
    for (int axis = 0; axis < 3; ++axis) {
        const auto hi = distance_gradients<double>(u, p, axis, +h);
        const auto lo = distance_gradients<double>(u, p, axis, -h);
        for (std::size_t k = 0; k < hi.size(); ++k) g.push_back((hi[k] - lo[k]) / (2 * h));
    }
    return g;
}

inline Entry distance_derivatives(std::size_t n = 2000) {
    return {"geometry: the derivative of a distance, Dual against differences", n, [n] {
        return check<std::vector<double>>(
            "geometry: the derivative of a distance, Dual against differences", Kind::Function, n,
            [](std::size_t i) {
                std::vector<double> u(40);
                for (std::size_t k = 0; k < 40; ++k) u[k] = -2 + 4 * halton(k < 32 ? i : i + 104729, k % 32);
                return u;
            },
            [](const std::vector<double>& u) {
                std::vector<double> g;
                const Vec<double, 3> p{u[0], u[1], u[2]};
                for (int axis = 0; axis < 3; ++axis) {
                    const auto v = distance_gradients<spatium::Dual<double>>(u, p, axis, 0.0);
                    g.insert(g.end(), v.begin(), v.end());
                }
                return g;
            },
            [](const std::vector<double>& u) { return distance_difference_quotients(u, 1e-6); },
            [](const std::vector<double>& a, const std::vector<double>& b) {
                double worst = 0.0;
                for (std::size_t k = 0; k < a.size(); ++k) worst = std::max(worst, std::abs(a[k] - b[k]));
                return worst;
            },
            [](const std::vector<double>& u) {
                const auto a = distance_difference_quotients(u, 1e-6);
                const auto b = distance_difference_quotients(u, 0.5e-6);
                double apart = 0.0;
                for (std::size_t k = 0; k < a.size(); ++k) apart = std::max(apart, std::abs(a[k] - b[k]));
                return 1e-8 + 3.0 * apart;
            });
    }};
}

}  // namespace

const Registrar registered_distance_derivatives{"distance_derivatives", 2000, [](std::size_t n) { return distance_derivatives(n); }};

}  // namespace symmetry
