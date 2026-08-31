#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/geometry/box.hpp>
#  include <spatium/geometry/concepts.hpp>
#  include <array>
#  include <cmath>
#endif

SPATIUM_EXPORT namespace spatium::geometry {

// K-simplex in N-dimensional space.
// K=0: point, K=1: segment, K=2: triangle, K=3: tetrahedron, ...
// Has K+1 vertices. Requires K <= N.

template<std::size_t N, std::size_t K, Scalar T = double>
    requires (K <= N)
struct Simplex {
    using ScalarType = T;
    using PointType = Vec<T, N>;
    static constexpr std::size_t ambient_dimension = N;
    static constexpr std::size_t intrinsic_dimension = K;

    std::array<PointType, K + 1> vertices;

    constexpr Simplex() = default;

    // Variadic constructor: Simplex(v0, v1, ..., vK)
    template<typename... Args>
        requires (sizeof...(Args) == K + 1 && (std::convertible_to<Args, PointType> && ...))
    constexpr Simplex(Args&&... args) : vertices{static_cast<PointType>(std::forward<Args>(args))...} {}

    constexpr const PointType& operator[](std::size_t i) const { return vertices[i]; }
    constexpr PointType& operator[](std::size_t i) { return vertices[i]; }

    // i-th face: (K-1)-simplex formed by removing vertex i
    constexpr auto face(std::size_t i) const requires (K > 0) {
        Simplex<N, K - 1, T> result;
        std::size_t idx = 0;
        for (std::size_t j = 0; j <= K; ++j) {
            if (j != i)
                result.vertices[idx++] = vertices[j];
        }
        return result;
    }

    constexpr PointType centroid() const {
        PointType sum = vertices[0];
        for (std::size_t i = 1; i <= K; ++i)
            sum = sum + vertices[i];
        return sum / T(K + 1);
    }

    // Measure: length (K=1), area (K=2), volume (K=3), etc.
    // Uses the Gram matrix determinant: measure = sqrt(det(G)) / K!
    // where G[i][j] = edge_i · edge_j, edges from vertex 0.
    T measure() const {
        if constexpr (K == 0) {
            return T{1}; // point has measure 1 (counting measure)
        } else if constexpr (K == 1) {
            return (vertices[1] - vertices[0]).norm();
        } else {
            // Build Gram matrix
            std::array<PointType, K> edges;
            for (std::size_t i = 0; i < K; ++i)
                edges[i] = vertices[i + 1] - vertices[0];

            // Gram matrix G[i][j] = edges[i] · edges[j]
            std::array<std::array<T, K>, K> G;
            for (std::size_t i = 0; i < K; ++i)
                for (std::size_t j = 0; j < K; ++j)
                    G[i][j] = edges[i].dot(edges[j]);

            // Determinant of KxK matrix (small K: direct formulas)
            T det = gram_det(G);

            // K! factorial
            T fact{1};
            for (std::size_t i = 2; i <= K; ++i)
                fact *= T(i);

            return std::sqrt(std::abs(det)) / fact;
        }
    }

    // Barycentric coordinates: K+1 weights summing to 1
    Vec<T, K + 1> barycentric(const PointType& p) const {
        // Solve via Gram matrix: project p-v0 onto edge basis
        std::array<PointType, K> edges;
        for (std::size_t i = 0; i < K; ++i)
            edges[i] = vertices[i + 1] - vertices[0];
        auto dp = p - vertices[0];

        // G * lambda = b, where b[i] = dp · edges[i]
        std::array<std::array<T, K>, K> G;
        std::array<T, K> b;
        for (std::size_t i = 0; i < K; ++i) {
            b[i] = dp.dot(edges[i]);
            for (std::size_t j = 0; j < K; ++j)
                G[i][j] = edges[i].dot(edges[j]);
        }

        // Solve (Gaussian elimination for small K)
        auto lambda = solve_linear(G, b);

        Vec<T, K + 1> result;
        T sum{0};
        for (std::size_t i = 0; i < K; ++i) {
            result[i + 1] = lambda[i];
            sum += lambda[i];
        }
        result[0] = T{1} - sum;
        return result;
    }

    bool contains(const PointType& p, T eps = spatium::epsilon<T>()) const {
        auto bary = barycentric(p);
        for (std::size_t i = 0; i <= K; ++i)
            if (bary[i] < -eps)
                return false;
        return true;
    }

    Box<N, T> bounding_box() const {
        Box<N, T> result{vertices[0], vertices[0]};
        for (std::size_t i = 1; i <= K; ++i)
            for (std::size_t j = 0; j < N; ++j) {
                result.min_corner[j] = std::min(result.min_corner[j], vertices[i][j]);
                result.max_corner[j] = std::max(result.max_corner[j], vertices[i][j]);
            }
        return result;
    }

    constexpr bool operator==(const Simplex&) const = default;

private:
    // Determinant of a KxK matrix (up to 4x4, general case via cofactor expansion)
    static T gram_det(const std::array<std::array<T, K>, K>& M) {
        if constexpr (K == 1) {
            return M[0][0];
        } else if constexpr (K == 2) {
            return M[0][0] * M[1][1] - M[0][1] * M[1][0];
        } else if constexpr (K == 3) {
            return M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
                 - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
                 + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
        } else {
            // Cofactor expansion along first row
            T det{0};
            for (std::size_t j = 0; j < K; ++j) {
                std::array<std::array<T, K - 1>, K - 1> sub;
                for (std::size_t r = 1; r < K; ++r) {
                    std::size_t col = 0;
                    for (std::size_t c = 0; c < K; ++c) {
                        if (c != j)
                            sub[r - 1][col++] = M[r][c];
                    }
                }
                T sign = (j % 2 == 0) ? T{1} : T{-1};
                // Recursive call would need Simplex<N, K-1> gram_det...
                // For K>3, use LU decomposition in practice
                // This compiles only for K<=4 in typical use
                det += sign * M[0][j] * gram_det_sub<K - 1>(sub);
            }
            return det;
        }
    }

    template<std::size_t M>
    static T gram_det_sub(const std::array<std::array<T, M>, M>& mat) {
        if constexpr (M == 1) {
            return mat[0][0];
        } else if constexpr (M == 2) {
            return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
        } else if constexpr (M == 3) {
            return mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1])
                 - mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0])
                 + mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
        } else {
            // For simplicity, not supporting K > 4 in this version
            static_assert(M <= 3, "Simplex measure for K > 4 not yet supported");
            return T{0};
        }
    }

    // Solve KxK linear system via Gaussian elimination (small K)
    static std::array<T, K> solve_linear(
        std::array<std::array<T, K>, K> A,
        std::array<T, K> b)
    {
        // Forward elimination
        for (std::size_t col = 0; col < K; ++col) {
            // Partial pivoting
            std::size_t max_row = col;
            for (std::size_t row = col + 1; row < K; ++row)
                if (std::abs(A[row][col]) > std::abs(A[max_row][col]))
                    max_row = row;
            std::swap(A[col], A[max_row]);
            std::swap(b[col], b[max_row]);

            if (std::abs(A[col][col]) < epsilon<T>() * epsilon<T>()) {
                // Singular: return zeros
                std::array<T, K> result{};
                return result;
            }

            for (std::size_t row = col + 1; row < K; ++row) {
                T factor = A[row][col] / A[col][col];
                for (std::size_t j = col; j < K; ++j)
                    A[row][j] -= factor * A[col][j];
                b[row] -= factor * b[col];
            }
        }

        // Back substitution
        std::array<T, K> x{};
        for (std::size_t i = K; i-- > 0;) {
            x[i] = b[i];
            for (std::size_t j = i + 1; j < K; ++j)
                x[i] -= A[i][j] * x[j];
            x[i] /= A[i][i];
        }
        return x;
    }
};

// Concept checks
static_assert(Shape<Simplex<3, 2>>);
static_assert(ClosedShape<Simplex<3, 2>>);
static_assert(Measurable<Simplex<3, 2>>);
static_assert(Bounded<Simplex<3, 2>>);

} // namespace spatium::geometry
