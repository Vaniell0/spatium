#include "../pair.hpp"
#include "../common.hpp"

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

// ── Distances and hits between primitives: double against Real50 ──
//
// Ten answers from one draw of 40 numbers in [-2, 2]: a point against a segment, a triangle
// and a box; two segments, a segment and a triangle, two triangles, two lines; a ray against a
// triangle and against a box. The same function templates on fifty digits are the reference.
// A hit gained or lost is a failure (the flag is part of the answer); the bound is 1e-9, the
// worst measured distance over the draws is below it by orders of magnitude.
//
// First run: none of this compiled on Real50. Qualified std::abs / std::clamp / std::min /
// std::max / std::swap in the segment-segment and line-line distances, in the ray-triangle and
// ray-box intersections and in the barycentric coordinates of a triangle blocked the scalar's
// own functions, and expression templates (`auto t1 = (a - b) * c`) broke std::swap and
// std::max. A pair that exercises a function is what makes it worth fixing.

template<class T>
std::vector<double> geometry_answers(const std::vector<double>& u) {
    using namespace spatium::geometry;
    const double* d = u.data();
    std::vector<double> out;
    const auto push = [&](const T& x) { out.push_back(spatium::primal_double(x)); };
    const Vec<T, 3> p = geometry_point<T>(d);
    const Segment<3, T> s1{geometry_point<T>(d + 3), geometry_point<T>(d + 6)};
    const Segment<3, T> s2{geometry_point<T>(d + 9), geometry_point<T>(d + 12)};
    const Triangle<3, T> t1(geometry_point<T>(d + 15), geometry_point<T>(d + 18), geometry_point<T>(d + 21));
    const Triangle<3, T> t2(geometry_point<T>(d + 24), geometry_point<T>(d + 27), geometry_point<T>(d + 30));
    const Vec<T, 3> lo = geometry_point<T>(d + 33);
    const Box<3, T> box{lo, lo + Vec<T, 3>{T(1), T(1), T(1)}};
    push(distance(p, s1));
    push(distance(p, t1));
    push(distance(p, box));
    push(distance(s1, s2));
    push(distance(s1, t1));
    push(distance(t1, t2));
    const auto l1 = Line<3, T>::from(geometry_point<T>(d + 3), geometry_point<T>(d + 6));
    const auto l2 = Line<3, T>::from(geometry_point<T>(d + 9), geometry_point<T>(d + 12));
    if (l1 && l2) push(distance(*l1, *l2)); else out.push_back(-1.0);
    const auto ray = Ray<3, T>::from(p, geometry_point<T>(d + 36) - p);
    if (ray) {
        const auto h = intersect(*ray, t1);
        out.push_back(h ? 1.0 : 0.0);
        if (h) push((*h)[0]);
        const auto hb = intersect(*ray, box);
        out.push_back(hb ? 1.0 : 0.0);
        if (hb) push((*hb)[0]);
    }
    return out;
}

inline Entry geometry_precision(std::size_t n = 3000) {
    return {"geometry: distances and hits, double against Real50", n, [n] {
        return check<std::vector<double>>(
            "geometry: distances and hits, double against Real50", Kind::Function, n,
            [](std::size_t i) {
                std::vector<double> u(40);
                for (std::size_t k = 0; k < 40; ++k) u[k] = -2 + 4 * halton(k < 32 ? i : i + 104729, k % 32);
                return u;
            },
            geometry_answers<spatium::Real50>, geometry_answers<double>,
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

const Registrar registered_geometry_precision{"geometry_precision", 3000, [](std::size_t n) { return geometry_precision(n); }};

}  // namespace symmetry
