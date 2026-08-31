#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/concepts.hpp>
#  include <spatium/algebra/groups/so3.hpp>
#  include <spatium/algebra/matrix.hpp>
#  include <spatium/algebra/vector.hpp>
#endif

SPATIUM_EXPORT namespace spatium::inline algebra {

// SE(3) — Special Euclidean Group in 3D.
// Elements: 4x4 homogeneous matrices [R t; 0 1] (rotation + translation).
// Lie algebra se(3): 6-vectors (angular velocity, linear velocity).

struct SE3 {
    using ElementType = Matrix<double, 4, 4>;
    using AlgebraType = Vec<double, 6>; // [omega(3), velocity(3)]
    using ScalarType = double;

    SO3 so3;

    ElementType identity() const { return ElementType::identity(); }

    ElementType compose(const ElementType& a, const ElementType& b) const {
        return a * b;
    }

    ElementType inverse(const ElementType& a) const {
        // [R t; 0 1]^{-1} = [R^T -R^T t; 0 1]
        auto Rt = rotation_of(a).transpose();
        auto t = translation_of(a);
        auto neg_Rt_t = Rt * (t * -1.0);

        ElementType result;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                result(i, j) = Rt(i, j);
        result(0, 3) = neg_Rt_t[0];
        result(1, 3) = neg_Rt_t[1];
        result(2, 3) = neg_Rt_t[2];
        result(3, 3) = 1.0;
        return result;
    }

    // Exponential map: se(3) → SE(3)
    ElementType exp(const AlgebraType& xi) const {
        Vec3 omega{xi[0], xi[1], xi[2]};
        Vec3 vel{xi[3], xi[4], xi[5]};

        auto R = so3.exp(omega);
        auto angle = omega.norm();

        Vec3 t;
        if (angle < epsilon<double>()) {
            t = vel; // pure translation
        } else {
            // V = I + (1-cos θ)/θ² K + (θ - sin θ)/θ³ K²
            auto axis = omega / angle;
            auto K = skew(axis);
            using std::sin; using std::cos;
            auto a2 = angle * angle;
            auto a3 = a2 * angle;
            auto V = Mat3::identity()
                   + K * ((1.0 - cos(angle)) / a2)
                   + (K * K) * ((angle - sin(angle)) / a3);
            t = V * vel;
        }

        return from_Rt(R, t);
    }

    // Logarithmic map: SE(3) → se(3)
    AlgebraType log(const ElementType& T_mat) const {
        auto R = rotation_of(T_mat);
        auto t = translation_of(T_mat);
        auto omega = so3.log(R);
        auto angle = omega.norm();

        Vec3 vel;
        if (angle < epsilon<double>()) {
            vel = t; // pure translation
        } else {
            auto axis = omega / angle;
            auto K = skew(axis);
            using std::sin; using std::cos; using std::tan;
            // V^{-1} = I - θ/2 K + (1 - θ/(2 tan(θ/2))) K²
            auto half = angle * 0.5;
            auto V_inv = Mat3::identity()
                       - K * half
                       + (K * K) * (1.0 - half / tan(half));
            vel = V_inv * t;
        }

        return AlgebraType{omega[0], omega[1], omega[2], vel[0], vel[1], vel[2]};
    }

    // Action: transform a point
    Vec3 act(const ElementType& T_mat, const Vec3& p) const {
        Vec<double, 4> h{p[0], p[1], p[2], 1.0};
        auto r = T_mat * h;
        return Vec3{r[0], r[1], r[2]};
    }

    // ── Convenience ────────────────────────────────────────────

    static ElementType from_Rt(const Mat3& R, const Vec3& t) {
        ElementType m;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                m(i, j) = R(i, j);
        m(0, 3) = t[0]; m(1, 3) = t[1]; m(2, 3) = t[2];
        m(3, 3) = 1.0;
        return m;
    }

    static Mat3 rotation_of(const ElementType& T_mat) {
        Mat3 R;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                R(i, j) = T_mat(i, j);
        return R;
    }

    static Vec3 translation_of(const ElementType& T_mat) {
        return Vec3{T_mat(0, 3), T_mat(1, 3), T_mat(2, 3)};
    }

private:
    static Mat3 skew(const Vec3& v) {
        Mat3 K;
        K(0, 1) = -v[2]; K(0, 2) =  v[1];
        K(1, 0) =  v[2]; K(1, 2) = -v[0];
        K(2, 0) = -v[1]; K(2, 1) =  v[0];
        return K;
    }
};

// Concept checks
static_assert(Group<SE3>);
static_assert(LieGroup<SE3>);

} // namespace spatium::algebra
