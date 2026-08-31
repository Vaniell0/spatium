#pragma once

// Triangle face metrics — single source of truth for the area /
// cross-product / unit-normal arithmetic used by mesh modules.
//
// Before this header six places (`mesh.hpp`, `differential.hpp`,
// `dec.hpp`, `quality.hpp`, `operations.hpp`, plus the central
// triangle helpers in `geometry/`) each open-coded `½‖(b−a)×(c−a)‖`
// or its 2D / Gram-determinant variants. Splitting them apart lets
// the optimiser inline a single body and keeps any future change
// (precision tweak, alternative formulation, SIMD path) localised.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/concepts.hpp>
#  include <cmath>
#endif

SPATIUM_EXPORT namespace spatium::mesh {

// Edge cross-product e1 × e2 for a 3D triangle. Magnitude equals
// twice the unsigned area; direction is the outward face normal.
// Not normalised — area-weighted on purpose so callers can sum
// face contributions without re-multiplying by 2A.
template<Scalar T>
inline Vec<T, 3> face_cross_3d(const Vec<T, 3>& a,
                               const Vec<T, 3>& b,
                               const Vec<T, 3>& c) {
    return (b - a).cross(c - a);
}

// Unsigned area of a 3D triangle.
template<Scalar T>
inline T face_area_3d(const Vec<T, 3>& a,
                      const Vec<T, 3>& b,
                      const Vec<T, 3>& c) {
    return T{0.5} * face_cross_3d(a, b, c).norm();
}

// Unsigned area of a 2D triangle, computed from the planar cross-z.
template<Scalar T>
inline T face_area_2d(const Vec<T, 2>& a,
                      const Vec<T, 2>& b,
                      const Vec<T, 2>& c) {
    auto e1 = b - a;
    auto e2 = c - a;
    T cross_z = e1[0] * e2[1] - e1[1] * e2[0];
    using std::abs;
    return T{0.5} * abs(cross_z);
}

// Unsigned area for an arbitrary-dimension triangle via the Gram
// determinant ½√(|e1|²|e2|² − (e1·e2)²). Reduces to ½‖e1×e2‖ in 3D
// and to ½|e1[0]e2[1] − e1[1]e2[0]| in 2D, but stays scalar in any
// embedding dimension. Used by `Mesh::area` for N-D coverage.
template<Scalar T, std::size_t N>
inline T face_area_gram(const Vec<T, N>& a,
                        const Vec<T, N>& b,
                        const Vec<T, N>& c) {
    auto ab = b - a;
    auto ac = c - a;
    auto aa = ab.dot(ab);
    auto bb = ac.dot(ac);
    auto ab_dot = ab.dot(ac);
    using std::abs;
    using std::sqrt;
    return T{0.5} * sqrt(abs(aa * bb - ab_dot * ab_dot));
}

} // namespace spatium::mesh
