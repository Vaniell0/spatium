#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/matrix.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/error.hpp>
#  include <spatium/morphism.hpp>
#  include <spatium/spaces/euclidean.hpp>
#  include <cmath>
#  include <numbers>
#endif

SPATIUM_EXPORT namespace spatium::geometry {

// Affine transform in N-D Euclidean space.
// Stores (N+1)x(N+1) homogeneous matrix. Acts as Morphism<Euclidean<N>, Euclidean<N>>.

template<std::size_t N, Scalar T = double>
struct AffineTransform {
    Matrix<T, N + 1, N + 1> matrix;

    AffineTransform() : matrix(Matrix<T, N + 1, N + 1>::identity()) {}
    explicit AffineTransform(Matrix<T, N + 1, N + 1> m) : matrix(m) {}

    // Apply to point
    Vec<T, N> operator()(const Vec<T, N>& p) const {
        Vec<T, N + 1> h;
        for (std::size_t i = 0; i < N; ++i) h[i] = p[i];
        h[N] = T{1};
        auto r = matrix * h;
        Vec<T, N> result;
        for (std::size_t i = 0; i < N; ++i) result[i] = r[i];
        return result;
    }

    // Apply to typed Point
    Point<Euclidean<N, T>> operator()(const Point<Euclidean<N, T>>& p) const {
        return Point<Euclidean<N, T>>{(*this)(p.raw())};
    }

    // Compose: this * other (apply other first, then this)
    AffineTransform operator*(const AffineTransform& other) const {
        return AffineTransform{matrix * other.matrix};
    }

    // Compute inverse transform (via linear sub-matrix + translation)
    Result<AffineTransform> inverse() const {
        // Extract NxN linear part and translation
        Matrix<T, N, N> lin;
        for (std::size_t r = 0; r < N; ++r)
            for (std::size_t c = 0; c < N; ++c)
                lin(r, c) = matrix(r, c);

        auto lin_inv = lin.inverse();
        if (!lin_inv) return std::unexpected(lin_inv.error());

        // inv_translation = -inv_linear * translation
        Vec<T, N> trans;
        for (std::size_t i = 0; i < N; ++i) trans[i] = matrix(i, N);

        Vec<T, N> neg_trans = Vec<T, N>{trans * T{-1}};
        auto inv_trans = (*lin_inv) * neg_trans;

        auto m = Matrix<T, N + 1, N + 1>::identity();
        for (std::size_t r = 0; r < N; ++r)
            for (std::size_t c = 0; c < N; ++c)
                m(r, c) = (*lin_inv)(r, c);
        for (std::size_t i = 0; i < N; ++i) m(i, N) = inv_trans[i];
        return AffineTransform{m};
    }

    // Convert to Morphism for pipe syntax (with inverse when possible)
    operator Morphism<Euclidean<N, T>, Euclidean<N, T>>() const {
        auto apply = [](const Matrix<T, N+1, N+1>& mat, const Vec<T, N>& p) {
            Vec<T, N + 1> h;
            for (std::size_t i = 0; i < N; ++i) h[i] = p[i];
            h[N] = T{1};
            auto r = mat * h;
            Vec<T, N> result;
            for (std::size_t i = 0; i < N; ++i) result[i] = r[i];
            return result;
        };
        auto fwd_mat = matrix;
        auto inv = inverse();
        if (inv) {
            auto inv_mat = inv->matrix;
            return {
                .forward = [apply, fwd_mat](const Vec<T, N>& p) { return apply(fwd_mat, p); },
                .inverse = [apply, inv_mat](const Vec<T, N>& p) { return apply(inv_mat, p); },
            };
        }
        return {
            .forward = [apply, fwd_mat](const Vec<T, N>& p) { return apply(fwd_mat, p); },
            .inverse = std::nullopt,
        };
    }

    // ── Factories ──────────────────────────────────────────────

    static AffineTransform translation(Vec<T, N> offset) {
        auto m = Matrix<T, N + 1, N + 1>::identity();
        for (std::size_t i = 0; i < N; ++i) m(i, N) = offset[i];
        return AffineTransform{m};
    }

    static AffineTransform scaling(Vec<T, N> factors) {
        auto m = Matrix<T, N + 1, N + 1>::identity();
        for (std::size_t i = 0; i < N; ++i) m(i, i) = factors[i];
        return AffineTransform{m};
    }

    static AffineTransform uniform_scaling(T factor) {
        auto m = Matrix<T, N + 1, N + 1>::identity();
        for (std::size_t i = 0; i < N; ++i) m(i, i) = factor;
        return AffineTransform{m};
    }

    // 3D rotation around axis (angle in radians)
    static AffineTransform rotation(Vec<T, 3> axis, T angle) requires (N == 3) {
        using std::cos; using std::sin;
        auto u = axis.normalized();
        auto c = cos(angle);
        auto s = sin(angle);
        auto t = T{1} - c;

        auto m = Matrix<T, 4, 4>::identity();
        m(0, 0) = t * u[0] * u[0] + c;
        m(0, 1) = t * u[0] * u[1] - s * u[2];
        m(0, 2) = t * u[0] * u[2] + s * u[1];
        m(1, 0) = t * u[0] * u[1] + s * u[2];
        m(1, 1) = t * u[1] * u[1] + c;
        m(1, 2) = t * u[1] * u[2] - s * u[0];
        m(2, 0) = t * u[0] * u[2] - s * u[1];
        m(2, 1) = t * u[1] * u[2] + s * u[0];
        m(2, 2) = t * u[2] * u[2] + c;
        return AffineTransform{m};
    }

    // 2D rotation (angle in radians)
    static AffineTransform rotation(T angle) requires (N == 2) {
        using std::cos; using std::sin;
        auto m = Matrix<T, 3, 3>::identity();
        m(0, 0) = cos(angle);  m(0, 1) = -sin(angle);
        m(1, 0) = sin(angle);  m(1, 1) = cos(angle);
        return AffineTransform{m};
    }
};

// Pipe: point | transform
template<std::size_t N, Scalar T>
Point<Euclidean<N, T>> operator|(const Point<Euclidean<N, T>>& p, const AffineTransform<N, T>& t) {
    return t(p);
}

// Pipe through Result<Vec>: intersect result | transform  (auto-unwrap)
template<std::size_t N, Scalar T>
Result<Vec<T, N>> operator|(const Result<Vec<T, N>>& v, const AffineTransform<N, T>& t) {
    if (!v) return std::unexpected(v.error());
    return t(*v);
}

// Pipe through Result<Point>: Result<Point> | transform  (auto-unwrap)
template<std::size_t N, Scalar T>
Result<Point<Euclidean<N, T>>> operator|(const Result<Point<Euclidean<N, T>>>& p,
                                         const AffineTransform<N, T>& t) {
    if (!p) return std::unexpected(p.error());
    return t(*p);
}

using Transform2 = AffineTransform<2>;
using Transform3 = AffineTransform<3>;

// ── Shorthand factories ───────────────────────────────────────

inline Transform3 translate(Vec3 offset) { return Transform3::translation(offset); }
inline Transform3 translate(double x, double y, double z) { return translate(Vec3{x, y, z}); }
inline Transform3 scale(double s) { return Transform3::uniform_scaling(s); }
inline Transform3 scale(Vec3 factors) { return Transform3::scaling(factors); }
inline Transform3 scale(double x, double y, double z) { return Transform3::scaling(Vec3{x, y, z}); }
inline Transform3 rotate_x(double angle) { return Transform3::rotation(Vec3{1, 0, 0}, angle); }
inline Transform3 rotate_y(double angle) { return Transform3::rotation(Vec3{0, 1, 0}, angle); }
inline Transform3 rotate_z(double angle) { return Transform3::rotation(Vec3{0, 0, 1}, angle); }

inline Transform2 translate(Vec2 offset) { return Transform2::translation(offset); }
inline Transform2 rotate(double angle) { return Transform2::rotation(angle); }

// ── Lazy transform expression templates ──────────────────────
// Usage: lazy(translate(1,0,0)) * lazy(rotate_z(0.5)) * lazy(scale(2.0))
// Defers matrix multiplication until .collapse() or .apply().

template<typename T>
concept TransformLike = requires(const T& t) {
    typename T::scalar_type;
    { T::dim } -> std::convertible_to<std::size_t>;
};

template<std::size_t N, Scalar T>
struct TransformLeaf {
    using scalar_type = T;
    static constexpr std::size_t dim = N;

    AffineTransform<N, T> t;

    Vec<T, N> apply(const Vec<T, N>& p) const { return t(p); }

    Point<Euclidean<N, T>> apply(const Point<Euclidean<N, T>>& p) const { return t(p); }

    AffineTransform<N, T> collapse() const { return t; }
};

template<typename L, typename R>
    requires TransformLike<L> && TransformLike<R>
          && (L::dim == R::dim)
struct TransformExpr {
    using scalar_type = typename L::scalar_type;
    static constexpr std::size_t dim = L::dim;

    L lhs;  // applied second
    R rhs;  // applied first

    Vec<scalar_type, dim> apply(const Vec<scalar_type, dim>& p) const {
        return lhs.apply(rhs.apply(p));
    }

    Point<Euclidean<dim, scalar_type>> apply(const Point<Euclidean<dim, scalar_type>>& p) const {
        return lhs.apply(rhs.apply(p));
    }

    AffineTransform<dim, scalar_type> collapse() const {
        return lhs.collapse() * rhs.collapse();
    }
};

// lazy() wrapper
template<std::size_t N, Scalar T>
TransformLeaf<N, T> lazy(const AffineTransform<N, T>& t) {
    return {t};
}

// operator* for lazy composition
template<TransformLike L, TransformLike R>
    requires (L::dim == R::dim)
TransformExpr<L, R> operator*(const L& lhs, const R& rhs) {
    return {lhs, rhs};
}

// Pipe: point | lazy_chain
template<TransformLike E, std::size_t N, Scalar T>
    requires (E::dim == N)
Point<Euclidean<N, T>> operator|(const Point<Euclidean<N, T>>& p, const E& expr) {
    return expr.apply(p);
}

// Pipe: Result<Vec> | lazy_chain
template<TransformLike E>
Result<Vec<typename E::scalar_type, E::dim>> operator|(
    const Result<Vec<typename E::scalar_type, E::dim>>& v, const E& expr) {
    if (!v) return std::unexpected(v.error());
    return expr.apply(*v);
}

} // namespace spatium::geometry
