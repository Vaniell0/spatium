#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/core/error.hpp>
#  include <spatium/geometry/box.hpp>
#  include <spatium/geometry/concepts.hpp>
#  include <spatium/geometry/line.hpp>
#  include <spatium/geometry/triangle.hpp>
#  include <cmath>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::geometry {

template<std::size_t N, Scalar T = double>
struct Polygon {
    using ScalarType = T;
    using PointType = Vec<T, N>;
    static constexpr std::size_t ambient_dimension = N;

    std::vector<PointType> vertices; // ordered, CCW orientation assumed

    std::size_t size() const { return vertices.size(); }

    Segment<N, T> edge(std::size_t i) const {
        return {vertices[i], vertices[(i + 1) % vertices.size()]};
    }

    // Measure (2D: shoelace formula)
    T measure() const requires (N == 2) {
        T sum{0};
        auto n = vertices.size();
        for (std::size_t i = 0; i < n; ++i) {
            auto j = (i + 1) % n;
            sum += vertices[i][0] * vertices[j][1];
            sum -= vertices[j][0] * vertices[i][1];
        }
        return std::abs(sum) * T{0.5};
    }

    // Measure (3D+: sum of triangle fan areas from centroid)
    T measure() const requires (N >= 3) {
        if (vertices.size() < 3) return T{0};
        T total{0};
        auto c = centroid();
        for (std::size_t i = 0; i < vertices.size(); ++i) {
            auto j = (i + 1) % vertices.size();
            Triangle<N, T> tri(c, vertices[i], vertices[j]);
            total += tri.measure();
        }
        return total;
    }

    T area() const { return measure(); }

    T perimeter() const {
        T total{0};
        for (std::size_t i = 0; i < vertices.size(); ++i)
            total += edge(i).length();
        return total;
    }

    PointType centroid() const {
        PointType sum{};
        for (const auto& v : vertices)
            sum = sum + v;
        return sum / T(vertices.size());
    }

    // Normal (3D, from first two edges)
    Vec<T, N> normal() const requires (N == 3) {
        if (vertices.size() < 3) return {};
        auto e1 = vertices[1] - vertices[0];
        auto e2 = vertices[2] - vertices[0];
        return e1.cross(e2).normalized();
    }

    // Point containment (2D: winding number algorithm)
    bool contains(const PointType& p) const requires (N == 2) {
        int winding = 0;
        auto n = vertices.size();
        for (std::size_t i = 0; i < n; ++i) {
            auto j = (i + 1) % n;
            if (vertices[i][1] <= p[1]) {
                if (vertices[j][1] > p[1]) {
                    if (cross_2d(vertices[j] - vertices[i], p - vertices[i]) > T{0})
                        ++winding;
                }
            } else {
                if (vertices[j][1] <= p[1]) {
                    if (cross_2d(vertices[j] - vertices[i], p - vertices[i]) < T{0})
                        --winding;
                }
            }
        }
        return winding != 0;
    }

    // Point containment (3D+: project to polygon plane, test in 2D)
    bool contains(const PointType& p) const requires (N >= 3) {
        // Check if point is on the polygon's plane
        if (vertices.size() < 3) return false;
        auto n = normal();
        auto d = n.dot(p - vertices[0]);
        if (std::abs(d) > epsilon<T>()) return false;

        // Use barycentric test via triangulation
        auto tris = triangulate();
        for (const auto& tri : tris) {
            if (tri.contains(p))
                return true;
        }
        return false;
    }

    T distance(const PointType& p) const {
        if (contains(p)) return T{0};
        T min_dist = std::numeric_limits<T>::max();
        for (std::size_t i = 0; i < vertices.size(); ++i) {
            auto d = edge(i).distance(p);
            min_dist = std::min(min_dist, d);
        }
        return min_dist;
    }

    PointType project(const PointType& p) const {
        T min_dist_sq = std::numeric_limits<T>::max();
        PointType best{};
        for (std::size_t i = 0; i < vertices.size(); ++i) {
            auto proj = edge(i).project(p);
            auto dist_sq = (p - proj).norm_squared();
            if (dist_sq < min_dist_sq) {
                min_dist_sq = dist_sq;
                best = proj;
            }
        }
        return best;
    }

    // Ear-clipping triangulation (2D)
    std::vector<Triangle<N, T>> triangulate() const requires (N == 2) {
        std::vector<Triangle<N, T>> result;
        if (vertices.size() < 3) return result;

        std::vector<std::size_t> indices(vertices.size());
        for (std::size_t i = 0; i < vertices.size(); ++i)
            indices[i] = i;

        // Determine winding sign: positive cross = convex vertex for CCW polygon
        T winding_sign = signed_area_2d() > T{0} ? T{1} : T{-1};

        while (indices.size() > 2) {
            bool found = false;
            auto n = indices.size();
            for (std::size_t i = 0; i < n; ++i) {
                auto pi = (i + n - 1) % n;
                auto ni = (i + 1) % n;
                const auto& prev = vertices[indices[pi]];
                const auto& curr = vertices[indices[i]];
                const auto& next = vertices[indices[ni]];

                // Convex vertex test: cross(curr-prev, next-curr) has same sign as winding
                auto cross = cross_2d(curr - prev, next - curr);
                if (cross * winding_sign < T{0}) continue; // reflex vertex, skip

                // Check no other vertex inside this ear triangle
                Triangle<N, T> tri(prev, curr, next);
                bool has_point_inside = false;
                for (std::size_t j = 0; j < n; ++j) {
                    if (j == pi || j == i || j == ni) continue;
                    if (tri.contains(vertices[indices[j]])) {
                        has_point_inside = true;
                        break;
                    }
                }
                if (has_point_inside) continue;

                result.push_back(tri);
                indices.erase(indices.begin() + static_cast<std::ptrdiff_t>(i));
                found = true;
                break;
            }
            if (!found) break; // degenerate polygon
        }
        return result;
    }

    // Fan triangulation for 3D+ (simpler, assumes convex or uses centroid fan)
    std::vector<Triangle<N, T>> triangulate() const requires (N >= 3) {
        std::vector<Triangle<N, T>> result;
        if (vertices.size() < 3) return result;
        for (std::size_t i = 1; i + 1 < vertices.size(); ++i) {
            result.push_back(Triangle<N, T>(vertices[0], vertices[i], vertices[i + 1]));
        }
        return result;
    }

    Box<N, T> bounding_box() const {
        if (vertices.empty()) return {};
        return Box<N, T>::from_points(vertices);
    }

    // Subspace: the plane this polygon lies on (3D only).
    Result<Hyperplane<N, T>> subspace() const requires (N == 3) {
        if (vertices.size() < 3)
            return std::unexpected(Error{ErrorCode::DegenerateInput, "need 3+ vertices"});
        return Hyperplane<N, T>::from_points(vertices[0], vertices[1], vertices[2]);
    }

private:
    static constexpr T cross_2d(const Vec<T, 2>& a, const Vec<T, 2>& b) requires (N == 2) {
        return a[0] * b[1] - a[1] * b[0];
    }

    T signed_area_2d() const requires (N == 2) {
        T sum{0};
        auto n = vertices.size();
        for (std::size_t i = 0; i < n; ++i) {
            auto j = (i + 1) % n;
            sum += vertices[i][0] * vertices[j][1];
            sum -= vertices[j][0] * vertices[i][1];
        }
        return sum * T{0.5};
    }
};

// Concept checks
static_assert(Shape<Polygon<2>>);
static_assert(ClosedShape<Polygon<2>>);
static_assert(Measurable<Polygon<2>>);
static_assert(Bounded<Polygon<2>>);

static_assert(Shape<Polygon<3>>);
static_assert(ClosedShape<Polygon<3>>);
static_assert(Measurable<Polygon<3>>);

// Aliases
using Polygon2 = Polygon<2>;
using Polygon3 = Polygon<3>;

} // namespace spatium::geometry
