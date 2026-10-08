#include "../pair.hpp"

#include <spatium/core/precision.hpp>
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

// ── Polygons: hull, areas of a boolean, clipping and distances, double against Real50 ──
//
// Two convex hulls of six points each in [-1, 1]^2, their areas, the area of their intersection,
// of their difference and of their symmetric difference, a point clipped to a hull and to a box,
// and the distance from the point to a hull and between the hulls. A result that is an error in
// one scalar and a value in the other changes the length of the answer and fails the pair.
//
// First run: `Polygon::measure` called `std::abs` qualified, so none of this compiled on Real50.
template<class T>
std::vector<double> polygon_answers(const std::vector<double>& u) {
    using namespace spatium::geometry;
    std::vector<double> out;
    const auto V = [&](std::size_t k) { return Vec<T, 2>{T(u[2 * k]), T(u[2 * k + 1])}; };
    std::vector<Vec<T, 2>> pa, pb;
    for (std::size_t k = 0; k < 6; ++k) pa.push_back(V(k));
    for (std::size_t k = 6; k < 12; ++k) pb.push_back(V(k));
    const auto ha = convex_hull(pa);
    const auto hb = convex_hull(pb);
    if (!ha || !hb) return out;
    const auto push = [&](const T& x) { out.push_back(spatium::primal_double(x)); };
    push(ha->area());
    push(hb->area());
    const auto ir = intersection_region(*ha, *hb);
    if (ir) push(ir->area()); else out.push_back(-1.0);
    const auto da = difference_area<2, T>(*ha, *hb);
    if (da) push(*da); else out.push_back(-1.0);
    const auto sa = symmetric_difference_area<2, T>(*ha, *hb);
    if (sa) push(*sa); else out.push_back(-1.0);
    const Vec<T, 2> p = V(12);
    const auto c = clip(p, *ha);
    if (c) { push((*c)[0]); push((*c)[1]); } else out.push_back(-1.0);
    const Box<2, T> box{Vec<T, 2>{T(-0.5), T(-0.5)}, Vec<T, 2>{T(0.5), T(0.5)}};
    const auto cb = clip(p, box);
    if (cb) { push((*cb)[0]); push((*cb)[1]); } else out.push_back(-1.0);
    push(distance(p, *ha));
    push(distance(*ha, *hb));
    return out;
}

inline Entry polygon_precision(std::size_t n = 3000) {
    return {"geometry: polygons, double against Real50", n, [n] {
        return check<std::vector<double>>(
            "geometry: polygons, double against Real50", Kind::Function, n,
            [](std::size_t i) {
                std::vector<double> u(26);
                for (std::size_t k = 0; k < 26; ++k) u[k] = -1 + 2 * halton(i, k);
                return u;
            },
            polygon_answers<spatium::Real50>, polygon_answers<double>,
            [](const std::vector<double>& a, const std::vector<double>& b) {
                if (a.size() != b.size()) return std::numeric_limits<double>::infinity();
                double worst = 0.0;
                for (std::size_t k = 0; k < a.size(); ++k) worst = std::max(worst, std::abs(a[k] - b[k]));
                return worst;
            },
            [](const std::vector<double>&) { return 1e-9; });
    }};
}

}  // namespace

const Registrar registered_polygon_precision{"polygon_precision", 3000, [](std::size_t n) { return polygon_precision(n); }};

}  // namespace symmetry
