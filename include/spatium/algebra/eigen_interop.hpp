#pragma once

// Eigen interop for spatium::Vec and spatium::Matrix.
// Only available when SPATIUM_HAS_EIGEN is defined and nonzero.
// Not included in the umbrella header (spatium.hpp) — opt-in only.

#include <spatium/_export_macro.hpp>

#if defined(SPATIUM_HAS_EIGEN) && SPATIUM_HAS_EIGEN

#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/algebra/matrix.hpp>
#  include <Eigen/Core>
#  include <climits>
#endif

SPATIUM_EXPORT namespace spatium {
inline namespace algebra {

// ── Vec <-> Eigen::Vector ────────────────────────────────────

template<Scalar T, std::size_t N>
Eigen::Matrix<T, static_cast<int>(N), 1> to_eigen(const Vec<T, N>& v) {
    static_assert(N <= INT_MAX, "Vec dimension too large for Eigen");
    Eigen::Matrix<T, static_cast<int>(N), 1> ev;
    for (std::size_t i = 0; i < N; ++i)
        ev(static_cast<int>(i)) = v[i];
    return ev;
}

template<Scalar T, std::size_t N>
Vec<T, N> from_eigen(const Eigen::Matrix<T, static_cast<int>(N), 1>& ev) {
    static_assert(N <= INT_MAX, "Vec dimension too large for Eigen");
    Vec<T, N> v;
    for (std::size_t i = 0; i < N; ++i)
        v[i] = ev(static_cast<int>(i));
    return v;
}

// ── Matrix <-> Eigen::Matrix ─────────────────────────────────

template<Scalar T, std::size_t R, std::size_t C>
Eigen::Matrix<T, static_cast<int>(R), static_cast<int>(C)> to_eigen(const Matrix<T, R, C>& m) {
    static_assert(R <= INT_MAX && C <= INT_MAX, "Matrix dimensions too large for Eigen");
    Eigen::Matrix<T, static_cast<int>(R), static_cast<int>(C)> em;
    for (std::size_t r = 0; r < R; ++r)
        for (std::size_t c = 0; c < C; ++c)
            em(static_cast<int>(r), static_cast<int>(c)) = m(r, c);
    return em;
}

template<Scalar T, std::size_t R, std::size_t C>
Matrix<T, R, C> from_eigen(const Eigen::Matrix<T, static_cast<int>(R), static_cast<int>(C)>& em) {
    static_assert(R <= INT_MAX && C <= INT_MAX, "Matrix dimensions too large for Eigen");
    Matrix<T, R, C> m;
    for (std::size_t r = 0; r < R; ++r)
        for (std::size_t c = 0; c < C; ++c)
            m(r, c) = em(static_cast<int>(r), static_cast<int>(c));
    return m;
}

// ── Zero-copy Eigen::Map views ───────────────────────────────
// Caller must ensure the source outlives the Map.
// Safe because Vec<T,N> stores std::array<T,N> (contiguous),
// and Matrix<T,R,C> stores std::array<T,R*C> in column-major order,
// matching Eigen's default ColMajor layout.

template<Scalar T, std::size_t N>
Eigen::Map<const Eigen::Matrix<T, static_cast<int>(N), 1>> eigen_view(const Vec<T, N>& v) {
    static_assert(N <= INT_MAX, "Vec dimension too large for Eigen");
    return Eigen::Map<const Eigen::Matrix<T, static_cast<int>(N), 1>>(v.data.data());
}

template<Scalar T, std::size_t N>
Eigen::Map<Eigen::Matrix<T, static_cast<int>(N), 1>> eigen_view(Vec<T, N>& v) {
    static_assert(N <= INT_MAX, "Vec dimension too large for Eigen");
    return Eigen::Map<Eigen::Matrix<T, static_cast<int>(N), 1>>(v.data.data());
}

template<Scalar T, std::size_t R, std::size_t C>
Eigen::Map<const Eigen::Matrix<T, static_cast<int>(R), static_cast<int>(C)>>
eigen_view(const Matrix<T, R, C>& m) {
    static_assert(R <= INT_MAX && C <= INT_MAX, "Matrix dimensions too large for Eigen");
    return Eigen::Map<const Eigen::Matrix<T, static_cast<int>(R), static_cast<int>(C)>>(m.data.data());
}

template<Scalar T, std::size_t R, std::size_t C>
Eigen::Map<Eigen::Matrix<T, static_cast<int>(R), static_cast<int>(C)>>
eigen_view(Matrix<T, R, C>& m) {
    static_assert(R <= INT_MAX && C <= INT_MAX, "Matrix dimensions too large for Eigen");
    return Eigen::Map<Eigen::Matrix<T, static_cast<int>(R), static_cast<int>(C)>>(m.data.data());
}

} // namespace algebra
} // namespace spatium

#endif // SPATIUM_HAS_EIGEN
