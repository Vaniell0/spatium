#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/error.hpp>
#  include <spatium/geometry/polygon.hpp>
#  include <algorithm>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::geometry {

// 2D convex hull via Andrew's monotone chain algorithm.
// Returns a Polygon<2, T> with vertices in CCW order.
// Fewer than three input points cannot define a polygon, so the
// function returns `unexpected(DegenerateInput)`. Three or more
// collinear points still yield a valid (degenerate) Polygon — the
// monotone chain handles that case downstream.

template<Scalar T>
Result<Polygon<2, T>> convex_hull(std::vector<Vec<T, 2>> points) {
    auto n = points.size();
    if (n < 3)
        return std::unexpected(Error{ErrorCode::DegenerateInput,
                                      "convex_hull requires at least 3 input points"});

    // Sort by x, then by y
    std::sort(points.begin(), points.end(), [](const auto& a, const auto& b) {
        return a[0] < b[0] || (a[0] == b[0] && a[1] < b[1]);
    });

    auto cross = [](const Vec<T, 2>& O, const Vec<T, 2>& A, const Vec<T, 2>& B) -> T {
        return (A[0] - O[0]) * (B[1] - O[1]) - (A[1] - O[1]) * (B[0] - O[0]);
    };

    std::vector<Vec<T, 2>> hull(2 * n);
    std::size_t k = 0;

    // Lower hull
    for (std::size_t i = 0; i < n; ++i) {
        while (k >= 2 && cross(hull[k - 2], hull[k - 1], points[i]) <= T{0})
            --k;
        hull[k++] = points[i];
    }

    // Upper hull
    auto lower_size = k + 1;
    for (std::size_t i = n - 1; i-- > 0;) {
        while (k >= lower_size && cross(hull[k - 2], hull[k - 1], points[i]) <= T{0})
            --k;
        hull[k++] = points[i];
    }

    hull.resize(k - 1); // remove last point (duplicate of first)
    return Polygon<2, T>{std::move(hull)};
}

} // namespace spatium::geometry
