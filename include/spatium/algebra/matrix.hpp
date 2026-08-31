#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/core/error.hpp>
#  include <array>
#  include <cmath>
#  include <cstddef>
#  include <format>
#  include <ostream>
#endif

SPATIUM_EXPORT namespace spatium {
inline namespace algebra {

// Column-major RxC matrix.
// M[col][row] in storage, but operator()(row, col) for math convention.

template<Scalar T, std::size_t R, std::size_t C>
struct Matrix {
    std::array<T, R * C> data{};

    constexpr Matrix() = default;

    // Access: (row, col)
    constexpr T& operator()(std::size_t row, std::size_t col) {
        return data[col * R + row];
    }
    constexpr const T& operator()(std::size_t row, std::size_t col) const {
        return data[col * R + row];
    }

    // Identity (square only)
    static constexpr Matrix identity() requires (R == C) {
        Matrix m;
        for (std::size_t i = 0; i < R; ++i) m(i, i) = T{1};
        return m;
    }

    // From row vectors
    static constexpr Matrix from_rows(const std::array<Vec<T, C>, R>& rows) {
        Matrix m;
        for (std::size_t r = 0; r < R; ++r)
            for (std::size_t c = 0; c < C; ++c)
                m(r, c) = rows[r][c];
        return m;
    }

    // Arithmetic
    constexpr Matrix operator+(const Matrix& rhs) const {
        Matrix result;
        for (std::size_t i = 0; i < R * C; ++i)
            result.data[i] = data[i] + rhs.data[i];
        return result;
    }

    constexpr Matrix operator-(const Matrix& rhs) const {
        Matrix result;
        for (std::size_t i = 0; i < R * C; ++i)
            result.data[i] = data[i] - rhs.data[i];
        return result;
    }

    constexpr Matrix operator*(T s) const {
        Matrix result;
        for (std::size_t i = 0; i < R * C; ++i)
            result.data[i] = data[i] * s;
        return result;
    }

    friend constexpr Matrix operator*(T s, const Matrix& m) { return m * s; }

    // Matrix-matrix multiply: (R x C) * (C x C2) = (R x C2)
    template<std::size_t C2>
    constexpr Matrix<T, R, C2> operator*(const Matrix<T, C, C2>& rhs) const {
        Matrix<T, R, C2> result;
        for (std::size_t c = 0; c < C2; ++c)
            for (std::size_t r = 0; r < R; ++r)
                for (std::size_t k = 0; k < C; ++k)
                    result(r, c) += (*this)(r, k) * rhs(k, c);
        return result;
    }

    // Matrix-vector multiply: (R x C) * Vec<C> = Vec<R>
    constexpr Vec<T, R> operator*(const Vec<T, C>& v) const {
        Vec<T, R> result;
        for (std::size_t r = 0; r < R; ++r)
            for (std::size_t c = 0; c < C; ++c)
                result[r] += (*this)(r, c) * v[c];
        return result;
    }

    constexpr Matrix<T, C, R> transpose() const {
        Matrix<T, C, R> result;
        for (std::size_t r = 0; r < R; ++r)
            for (std::size_t c = 0; c < C; ++c)
                result(c, r) = (*this)(r, c);
        return result;
    }

    constexpr bool operator==(const Matrix&) const = default;

    // Trace (square only)
    constexpr T trace() const requires (R == C) {
        T sum{};
        for (std::size_t i = 0; i < R; ++i) sum += (*this)(i, i);
        return sum;
    }

    // Determinant (square only)
    constexpr T determinant() const requires (R == C && R == 1) {
        return data[0];
    }

    constexpr T determinant() const requires (R == C && R == 2) {
        return (*this)(0, 0) * (*this)(1, 1) - (*this)(0, 1) * (*this)(1, 0);
    }

    constexpr T determinant() const requires (R == C && R == 3) {
        return (*this)(0, 0) * ((*this)(1, 1) * (*this)(2, 2) - (*this)(1, 2) * (*this)(2, 1))
             - (*this)(0, 1) * ((*this)(1, 0) * (*this)(2, 2) - (*this)(1, 2) * (*this)(2, 0))
             + (*this)(0, 2) * ((*this)(1, 0) * (*this)(2, 1) - (*this)(1, 1) * (*this)(2, 0));
    }

    // General determinant via cofactor expansion (N >= 4)
    T determinant() const requires (R == C && R >= 4) {
        // LU decomposition
        Matrix<T, R, C> L, U;
        T det{1};
        auto copy = *this;
        for (std::size_t i = 0; i < R; ++i) {
            // Partial pivot
            std::size_t max_row = i;
            using std::abs;
            for (std::size_t k = i + 1; k < R; ++k)
                if (abs(copy(k, i)) > abs(copy(max_row, i)))
                    max_row = k;
            if (max_row != i) {
                for (std::size_t c = 0; c < C; ++c)
                    std::swap(copy(i, c), copy(max_row, c));
                det = -det;
            }
            if (abs(copy(i, i)) < epsilon<T>()) return T{0};
            det *= copy(i, i);
            for (std::size_t k = i + 1; k < R; ++k) {
                auto factor = copy(k, i) / copy(i, i);
                for (std::size_t j = i + 1; j < C; ++j)
                    copy(k, j) -= factor * copy(i, j);
            }
        }
        return det;
    }

    // Inverse (square, small N direct formulas)
    constexpr Result<Matrix> inverse() const requires (R == C && R == 2) {
        auto det = determinant();
        using std::abs;
        if (abs(det) < epsilon<T>())
            return std::unexpected(Error{ErrorCode::SingularMatrix, "singular 2x2 matrix"});
        auto inv_det = T{1} / det;
        Matrix result;
        result(0, 0) = (*this)(1, 1) * inv_det;
        result(0, 1) = -(*this)(0, 1) * inv_det;
        result(1, 0) = -(*this)(1, 0) * inv_det;
        result(1, 1) = (*this)(0, 0) * inv_det;
        return result;
    }

    constexpr Result<Matrix> inverse() const requires (R == C && R == 3) {
        auto det = determinant();
        using std::abs;
        if (abs(det) < epsilon<T>())
            return std::unexpected(Error{ErrorCode::SingularMatrix, "singular 3x3 matrix"});
        auto inv_det = T{1} / det;
        Matrix result;
        result(0, 0) = ((*this)(1, 1) * (*this)(2, 2) - (*this)(1, 2) * (*this)(2, 1)) * inv_det;
        result(0, 1) = ((*this)(0, 2) * (*this)(2, 1) - (*this)(0, 1) * (*this)(2, 2)) * inv_det;
        result(0, 2) = ((*this)(0, 1) * (*this)(1, 2) - (*this)(0, 2) * (*this)(1, 1)) * inv_det;
        result(1, 0) = ((*this)(1, 2) * (*this)(2, 0) - (*this)(1, 0) * (*this)(2, 2)) * inv_det;
        result(1, 1) = ((*this)(0, 0) * (*this)(2, 2) - (*this)(0, 2) * (*this)(2, 0)) * inv_det;
        result(1, 2) = ((*this)(0, 2) * (*this)(1, 0) - (*this)(0, 0) * (*this)(1, 2)) * inv_det;
        result(2, 0) = ((*this)(1, 0) * (*this)(2, 1) - (*this)(1, 1) * (*this)(2, 0)) * inv_det;
        result(2, 1) = ((*this)(0, 1) * (*this)(2, 0) - (*this)(0, 0) * (*this)(2, 1)) * inv_det;
        result(2, 2) = ((*this)(0, 0) * (*this)(1, 1) - (*this)(0, 1) * (*this)(1, 0)) * inv_det;
        return result;
    }

    // Column/row access
    constexpr Vec<T, R> col(std::size_t c) const {
        Vec<T, R> result;
        for (std::size_t r = 0; r < R; ++r) result[r] = (*this)(r, c);
        return result;
    }

    constexpr Vec<T, C> row(std::size_t r) const {
        Vec<T, C> result;
        for (std::size_t c = 0; c < C; ++c) result[c] = (*this)(r, c);
        return result;
    }

    // Stream output
    friend std::ostream& operator<<(std::ostream& os, const Matrix& m) {
        os << '[';
        for (std::size_t r = 0; r < R; ++r) {
            if (r > 0) os << "; ";
            for (std::size_t c = 0; c < C; ++c) {
                if (c > 0) os << ", ";
                os << m(r, c);
            }
        }
        return os << ']';
    }
};

// Aliases
using Mat2 = Matrix<double, 2, 2>;
using Mat3 = Matrix<double, 3, 3>;
using Mat4 = Matrix<double, 4, 4>;

using Mat2f = Matrix<float, 2, 2>;
using Mat3f = Matrix<float, 3, 3>;
using Mat4f = Matrix<float, 4, 4>;

} // namespace algebra
} // namespace spatium

// std::format
template<spatium::Scalar T, std::size_t R, std::size_t C>
struct std::formatter<spatium::Matrix<T, R, C>> {
    constexpr auto parse(auto& ctx) { return ctx.begin(); }
    auto format(const spatium::Matrix<T, R, C>& m, auto& ctx) const {
        auto out = ctx.out();
        *out++ = '[';
        for (std::size_t r = 0; r < R; ++r) {
            if (r > 0) { *out++ = ';'; *out++ = ' '; }
            for (std::size_t c = 0; c < C; ++c) {
                if (c > 0) { *out++ = ','; *out++ = ' '; }
                ctx.advance_to(out);
                out = std::formatter<T>{}.format(m(r, c), ctx);
            }
        }
        *out++ = ']';
        return out;
    }
};
