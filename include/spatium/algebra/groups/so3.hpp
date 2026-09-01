#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/concepts.hpp>
#  include <spatium/algebra/dual.hpp>
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
//
// Templated on Scalar so Dual<T> substitutes for T like everywhere else in
// the tree: differentiating through a chain of rotations (pose-graph /
// SLAM-style optimization) then needs no separate analytical Jacobian, the
// same way calculus.hpp's gradient() differentiates any plain function.
template<Scalar T = double>
struct SO3 {
    using ElementType = Matrix<T, 3, 3>;
    using AlgebraType = Vec<T, 3>; // angular velocity vector (compact form of skew-symmetric)
    using ScalarType = T;
    using Vec3 = Vec<T, 3>;

    // Group operations
    ElementType identity() const { return ElementType::identity(); }

    ElementType compose(const ElementType& a, const ElementType& b) const {
        return a * b;
    }

    ElementType inverse(const ElementType& a) const {
        return a.transpose(); // orthogonal matrix: inverse = transpose
    }

    // Lie group: exponential map (Rodrigues formula).
    // R = I + a(θ)K + b(θ)K², K = skew(v) (UNNORMALIZED — K = θ·skew(axis)),
    // a(θ) = sin(θ)/θ, b(θ) = (1-cos θ)/θ². Built from θ² = v·v (a smooth
    // polynomial in v) rather than branching on the normalized axis v/θ:
    // that division, and the old `if (angle < eps) return identity();`
    // short-circuit, were both exact at θ=0 in VALUE but silently returned
    // a v-INDEPENDENT constant there — under Dual<T> this zeroed the
    // gradient of exp() at the identity, the single most common
    // optimization starting point, even though exp() is analytically
    // smooth there (its own Taylor series in v has no singularity; only
    // the angle=sqrt(θ²) and axis=v/θ intermediates do).
    ElementType exp(const AlgebraType& v) const {
        using std::sin; using std::cos; using std::sqrt;
        T theta2 = v.dot(v);
        auto angle = sqrt(theta2);
        auto K = skew(v);

        T a, b;
        if (angle < epsilon<T>()) {
            // Taylor series in θ² (smooth at v=0): sin(θ)/θ = 1 - θ²/6 + ...,
            // (1-cos θ)/θ² = 1/2 - θ²/24 + ...
            a = T{1} - theta2 / T{6};
            b = T{0.5} - theta2 / T{24};
        } else {
            a = sin(angle) / angle;
            b = (T{1} - cos(angle)) / theta2;
        }

        return identity() + K * a + (K * K) * b;
    }

    // Logarithmic map: rotation matrix → axis-angle vector.
    // v = raw · θ/(2 sin θ), raw = vee(R - Rᵀ) (linear in R, always smooth).
    AlgebraType log(const ElementType& R) const {
        using std::acos; using std::sqrt; using std::abs; using std::sin;

        auto trace = R(0, 0) + R(1, 1) + R(2, 2);
        auto cos_angle = std::clamp((trace - T{1}) * T{0.5}, T{-1}, T{1});
        AlgebraType raw{R(2, 1) - R(1, 2), R(0, 2) - R(2, 0), R(1, 0) - R(0, 1)};

        if (cos_angle > T{1} - epsilon<T>()) {
            // Near R=I: acos'(x) = -1/sqrt(1-x²) itself diverges at x=1, so
            // computing angle=acos(cos_angle) here would corrupt the
            // derivative under Dual<T> before even reaching the (also
            // removable) θ/sin(θ) singularity below — same failure mode as
            // exp()'s old identity-at-origin shortcut, just one function
            // upstream. θ² ≈ 2(1-cos_angle) (small-angle Taylor, smooth in
            // R, no acos involved) sidesteps that; θ/(2 sin θ) ≈ 1/2 + θ²/12
            // is the matching Taylor coefficient. The old
            // `if (angle < eps) return AlgebraType{};` here had the exact
            // same bug as exp()'s: correct in value, zero derivative w.r.t.
            // R at the identity.
            T theta2 = (T{1} - cos_angle) * T{2};
            return raw * (T{0.5} + theta2 / T{12});
        }

        auto angle = acos(cos_angle);

        if (abs(angle - T{std::numbers::pi}) < epsilon<T>()) {
            // angle ≈ π: a genuine coordinate singularity of the axis-angle
            // chart itself (+πn̂ and -πn̂ represent the same rotation) — not
            // a removable Dual-derivative artifact like the two cases
            // above, so left as the existing value-only extraction. Not
            // part of this fix; a gradient-based optimizer essentially
            // never lands exactly here in practice.
            auto RpI = R + identity();
            int best = 0;
            T best_norm{0};
            for (int i = 0; i < 3; ++i) {
                auto col = RpI.col(i);
                auto n = col.norm();
                if (n > best_norm) { best_norm = n; best = i; }
            }
            auto axis = RpI.col(best).normalized();
            return axis * angle;
        }

        auto s = sin(angle);
        return raw * (angle / (T{2} * s));
    }

    // Action: rotate a point
    Vec3 act(const ElementType& R, const Vec3& p) const {
        return R * p;
    }

    // ── Convenience factories ──────────────────────────────────

    // Rotation around axis by angle (radians)
    ElementType from_axis_angle(const Vec3& axis, T angle) const {
        return exp(axis.normalized() * angle);
    }

    // Rotation around X/Y/Z axes
    ElementType rx(T angle) const { return from_axis_angle(Vec3{T{1}, T{0}, T{0}}, angle); }
    ElementType ry(T angle) const { return from_axis_angle(Vec3{T{0}, T{1}, T{0}}, angle); }
    ElementType rz(T angle) const { return from_axis_angle(Vec3{T{0}, T{0}, T{1}}, angle); }

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

// Concept checks — both the default double instantiation and, specifically,
// Dual<double> (this is the whole point of templating: a Scalar substitute
// that carries derivatives must still satisfy the same Group/LieGroup shape).
static_assert(Group<SO3<>>);
static_assert(LieGroup<SO3<>>);
static_assert(Group<SO3<Dual<double>>>);
static_assert(LieGroup<SO3<Dual<double>>>);

} // namespace spatium::algebra
