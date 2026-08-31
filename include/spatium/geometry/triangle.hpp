#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/core/error.hpp>
#  include <spatium/geometry/box.hpp>
#  include <spatium/geometry/concepts.hpp>
#  include <spatium/geometry/hyperplane.hpp>
#  include <spatium/geometry/line.hpp>
#  include <array>
#  include <cmath>
#endif

SPATIUM_EXPORT namespace spatium::geometry {

template<std::size_t N, Scalar T = double>
struct Triangle {
    using ScalarType = T;
    using PointType = Vec<T, N>;
    static constexpr std::size_t ambient_dimension = N;

    std::array<PointType, 3> vertices;

    constexpr Triangle() = default;
    constexpr Triangle(PointType a, PointType b, PointType c)
        : vertices{a, b, c} {}

    constexpr const PointType& operator[](std::size_t i) const { return vertices[i]; }
    constexpr PointType& operator[](std::size_t i) { return vertices[i]; }

    // Edges
    constexpr Segment<N, T> edge(std::size_t i) const {
        // edge i connects vertex (i+1)%3 to (i+2)%3 (opposite to vertex i)
        return {vertices[(i + 1) % 3], vertices[(i + 2) % 3]};
    }

    // Edge vectors from vertex 0
    constexpr PointType e01() const { return vertices[1] - vertices[0]; }
    constexpr PointType e02() const { return vertices[2] - vertices[0]; }

    // Measure (area) via cross product (3D) or 2D formula
    T measure() const requires (N == 3) {
        return e01().cross(e02()).norm() * T{0.5};
    }

    T measure() const requires (N == 2) {
        auto a = e01();
        auto b = e02();
        return std::abs(a[0] * b[1] - a[1] * b[0]) * T{0.5};
    }

    // General N-D measure via Gram determinant: sqrt(det(G))/2
    // where G = [[e01·e01, e01·e02], [e01·e02, e02·e02]]
    T measure() const requires (N > 3) {
        auto a = e01();
        auto b = e02();
        auto aa = a.dot(a);
        auto ab = a.dot(b);
        auto bb = b.dot(b);
        return std::sqrt(aa * bb - ab * ab) * T{0.5};
    }

    T area() const { return measure(); }

    T perimeter() const {
        return edge(0).length() + edge(1).length() + edge(2).length();
    }

    // Normal (3D only)
    Vec<T, N> normal() const requires (N == 3) {
        return e01().cross(e02()).normalized();
    }

    // Supporting plane (3D)
    Result<Hyperplane<N, T>> supporting_plane() const requires (N == 3) {
        return Hyperplane<N, T>::from_points(vertices[0], vertices[1], vertices[2]);
    }

    // Subspace: the plane this triangle is a bounded region of.
    Result<Hyperplane<N, T>> subspace() const requires (N == 3) {
        return supporting_plane();
    }

    constexpr PointType centroid() const {
        return (vertices[0] + vertices[1] + vertices[2]) / T{3};
    }

    // Barycentric coordinates of point p with respect to this triangle.
    // Returns (u, v, w) where p ≈ u*v0 + v*v1 + w*v2, u+v+w=1.
    // Works by projecting onto the triangle's plane in N-D.
    Vec<T, 3> barycentric(const PointType& p) const {
        auto v0 = e01();
        auto v1 = e02();
        auto v2 = p - vertices[0];

        auto d00 = v0.dot(v0);
        auto d01 = v0.dot(v1);
        auto d11 = v1.dot(v1);
        auto d20 = v2.dot(v0);
        auto d21 = v2.dot(v1);

        auto denom = d00 * d11 - d01 * d01;
        if (std::abs(denom) < epsilon<T>() * epsilon<T>())
            return Vec<T, 3>{T{1} / T{3}, T{1} / T{3}, T{1} / T{3}}; // degenerate

        auto v = (d11 * d20 - d01 * d21) / denom;
        auto w = (d00 * d21 - d01 * d20) / denom;
        auto u = T{1} - v - w;
        return Vec<T, 3>{u, v, w};
    }

    bool contains(const PointType& p, T eps = spatium::epsilon<T>()) const {
        auto bary = barycentric(p);
        return bary[0] >= -eps && bary[1] >= -eps && bary[2] >= -eps;
    }

    // Closest point on triangle to p
    PointType project(const PointType& p) const {
        auto bary = barycentric(p);

        // If inside, project onto the triangle's affine plane
        if (bary[0] >= T{0} && bary[1] >= T{0} && bary[2] >= T{0}) {
            return vertices[0] * bary[0] + vertices[1] * bary[1] + vertices[2] * bary[2];
        }

        // Outside: closest point is on the nearest edge
        T best_dist_sq = std::numeric_limits<T>::max();
        PointType best;
        for (std::size_t i = 0; i < 3; ++i) {
            auto proj = edge(i).project(p);
            auto dist_sq = (p - proj).norm_squared();
            if (dist_sq < best_dist_sq) {
                best_dist_sq = dist_sq;
                best = proj;
            }
        }
        return best;
    }

    T distance(const PointType& p) const {
        return (p - project(p)).norm();
    }

    // Midpoint subdivision: split each edge at midpoint → 4 triangles
    constexpr std::array<Triangle, 4> subdivide() const {
        auto m01 = (vertices[0] + vertices[1]) * T{0.5};
        auto m12 = (vertices[1] + vertices[2]) * T{0.5};
        auto m20 = (vertices[2] + vertices[0]) * T{0.5};
        return {{
            Triangle(vertices[0], m01, m20),
            Triangle(m01, vertices[1], m12),
            Triangle(m20, m12, vertices[2]),
            Triangle(m01, m12, m20),
        }};
    }

    Box<N, T> bounding_box() const {
        Box<N, T> result{vertices[0], vertices[0]};
        for (std::size_t i = 1; i < 3; ++i) {
            for (std::size_t j = 0; j < N; ++j) {
                result.min_corner[j] = std::min(result.min_corner[j], vertices[i][j]);
                result.max_corner[j] = std::max(result.max_corner[j], vertices[i][j]);
            }
        }
        return result;
    }

    constexpr bool operator==(const Triangle&) const = default;
};

// Concept checks
static_assert(Shape<Triangle<3>>);
static_assert(ClosedShape<Triangle<3>>);
static_assert(Measurable<Triangle<3>>);
static_assert(Bounded<Triangle<3>>);
static_assert(DistanceQueryable<Triangle<3>>);
static_assert(BoundedRegion<Triangle<3>>);
static_assert(Shape<Triangle<2>>);
static_assert(Measurable<Triangle<2>>);

// Aliases
using Triangle2 = Triangle<2>;
using Triangle3 = Triangle<3>;

} // namespace spatium::geometry
