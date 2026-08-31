#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/concepts.hpp>
#  include <spatium/algebra/matrix.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <cmath>
#  include <numbers>
#endif

SPATIUM_EXPORT namespace spatium::inline algebra {

// SO(3) — Special Orthogonal Group in 3D.
// Elements: 3x3 rotation matrices (det=1, R^T R = I).
// Lie algebra so(3): 3-vectors (angular velocity), mapped to skew-symmetric matrices.

struct SO3 {
    using ElementType = Matrix<double, 3, 3>;
    using AlgebraType = Vec3; // angular velocity vector (compact form of skew-symmetric)
    using ScalarType = double;

    // Group operations
    ElementType identity() const { return ElementType::identity(); }

    ElementType compose(const ElementType& a, const ElementType& b) const {
        return a * b;
    }

    ElementType inverse(const ElementType& a) const {
        return a.transpose(); // orthogonal matrix: inverse = transpose
    }

    // Lie group: exponential map (Rodrigues formula)
    // v = axis * angle (axis-angle representation)
    ElementType exp(const AlgebraType& v) const {
        using std::sin; using std::cos; using std::sqrt;
        auto angle = v.norm();
        if (angle < epsilon<double>()) return identity();

        auto axis = v / angle;
        auto K = skew(axis);
        auto s = sin(angle);
        auto c = cos(angle);

        // R = I + sin(θ)K + (1-cos(θ))K²
        return identity() + K * s + (K * K) * (1.0 - c);
    }

    // Logarithmic map: rotation matrix → axis-angle vector
    AlgebraType log(const ElementType& R) const {
        using std::acos; using std::sqrt; using std::abs;

        auto trace = R(0, 0) + R(1, 1) + R(2, 2);
        auto cos_angle = (trace - 1.0) * 0.5;
        cos_angle = std::clamp(cos_angle, -1.0, 1.0);
        auto angle = acos(cos_angle);

        if (angle < epsilon<double>()) return AlgebraType{};

        if (abs(angle - std::numbers::pi) < epsilon<double>()) {
            // angle ≈ π: extract axis from R + I
            // Find column of R + I with largest norm
            auto RpI = R + identity();
            int best = 0;
            double best_norm = 0;
            for (int i = 0; i < 3; ++i) {
                auto col = RpI.col(i);
                auto n = col.norm();
                if (n > best_norm) { best_norm = n; best = i; }
            }
            auto axis = RpI.col(best).normalized();
            return axis * angle;
        }

        // General case: axis from skew-symmetric part
        auto s = sin(angle);
        AlgebraType axis{
            (R(2, 1) - R(1, 2)) / (2.0 * s),
            (R(0, 2) - R(2, 0)) / (2.0 * s),
            (R(1, 0) - R(0, 1)) / (2.0 * s),
        };
        return axis * angle;
    }

    // Action: rotate a point
    Vec3 act(const ElementType& R, const Vec3& p) const {
        return R * p;
    }

    // ── Convenience factories ──────────────────────────────────

    // Rotation around axis by angle (radians)
    ElementType from_axis_angle(const Vec3& axis, double angle) const {
        return exp(axis.normalized() * angle);
    }

    // Rotation around X/Y/Z axes
    ElementType rx(double angle) const { return from_axis_angle(Vec3{1, 0, 0}, angle); }
    ElementType ry(double angle) const { return from_axis_angle(Vec3{0, 1, 0}, angle); }
    ElementType rz(double angle) const { return from_axis_angle(Vec3{0, 0, 1}, angle); }

private:
    // Skew-symmetric matrix from 3-vector: [v]_× such that [v]_× w = v × w
    static ElementType skew(const Vec3& v) {
        ElementType K;
        K(0, 1) = -v[2]; K(0, 2) =  v[1];
        K(1, 0) =  v[2]; K(1, 2) = -v[0];
        K(2, 0) = -v[1]; K(2, 1) =  v[0];
        return K;
    }
};

// Concept checks
static_assert(Group<SO3>);
static_assert(LieGroup<SO3>);

} // namespace spatium::algebra
