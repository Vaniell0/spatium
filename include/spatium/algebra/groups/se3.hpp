#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/concepts.hpp>
#  include <spatium/algebra/groups/so3.hpp>
#  include <spatium/algebra/linear_solve.hpp>
#  include <spatium/algebra/matrix.hpp>
#  include <spatium/algebra/vector.hpp>
#endif

SPATIUM_EXPORT namespace spatium::inline algebra {

// SE(3) — Special Euclidean Group in 3D.
// Elements: 4x4 homogeneous matrices [R t; 0 1] (rotation + translation).
// Lie algebra se(3): 6-vectors (angular velocity, linear velocity).
//
// Templated on Scalar for the same reason SO3 is — see so3.hpp.
template<Scalar T = double>
struct SE3 {
    using ElementType = Matrix<T, 4, 4>;
    using AlgebraType = Vec<T, 6>; // [omega(3), velocity(3)]
    using ScalarType = T;
    using Vec3 = Vec<T, 3>;
    using Mat3 = Matrix<T, 3, 3>;

    SO3<T> so3;

    ElementType identity() const { return ElementType::identity(); }

    ElementType compose(const ElementType& a, const ElementType& b) const {
        return a * b;
    }

    ElementType inverse(const ElementType& a) const {
        // [R t; 0 1]^{-1} = [R^T -R^T t; 0 1]
        auto Rt = rotation_of(a).transpose();
        auto t = translation_of(a);
        auto neg_Rt_t = Rt * (t * T{-1});

        ElementType result;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                result(i, j) = Rt(i, j);
        result(0, 3) = neg_Rt_t[0];
        result(1, 3) = neg_Rt_t[1];
        result(2, 3) = neg_Rt_t[2];
        result(3, 3) = T{1};
        return result;
    }

    // Exponential map: se(3) → SE(3). t = V(ω)·vel — see
    // translation_jacobian() below for V and why it changed.
    ElementType exp(const AlgebraType& xi) const {
        Vec3 omega{xi[0], xi[1], xi[2]};
        Vec3 vel{xi[3], xi[4], xi[5]};

        auto R = so3.exp(omega);
        Vec3 t = translation_jacobian(omega) * vel;

        return from_Rt(R, t);
    }

    // Logarithmic map: SE(3) → se(3). vel = V(ω)⁻¹·t.
    AlgebraType log(const ElementType& T_mat) const {
        auto R = rotation_of(T_mat);
        auto t = translation_of(T_mat);
        auto omega = so3.log(R);

        auto V = translation_jacobian(omega);
        auto V_inv_result = invert(V);
        Mat3 V_inv = V_inv_result ? *V_inv_result : Mat3::identity();
        Vec3 vel = V_inv * t;

        return AlgebraType{omega[0], omega[1], omega[2], vel[0], vel[1], vel[2]};
    }

    // Action: transform a point
    Vec3 act(const ElementType& T_mat, const Vec3& p) const {
        Vec<T, 4> h{p[0], p[1], p[2], T{1}};
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
        m(3, 3) = T{1};
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

    // V(ω) = I + b(θ)K + c(θ)K², K = skew(ω) — UNNORMALIZED, unlike so3.hpp's
    // own exp()/log() which use skew of the axis-angle vector directly too.
    // b(θ) = (1-cos θ)/θ², c(θ) = (θ-sin θ)/θ³ — verified against a
    // brute-force 4×4 matrix exponential of the se(3) generator. Shared by
    // exp() (used as V itself) and log() (inverted via the library's own
    // invert(), rather than a second hand-derived V⁻¹ closed form — safer
    // given what follows).
    //
    // The previous code used skew(ω/θ) (the UNIT axis) with these same
    // coefficients in exp(), and a differently-shaped closed form for V⁻¹ in
    // log() — a real, pre-existing, silent correctness bug for any nonzero
    // rotation combined with a translation, not just the Dual-
    // differentiability gap this whole rewrite was actually chasing: no
    // prior test exercised exp()'s general branch against an independent
    // ground truth, only self-consistent exp/log roundtrips (pass under
    // either convention, since both sides used the same wrong one) and
    // act()-only checks that build T via from_Rt() directly, bypassing
    // exp() entirely.
    //
    // Taylor-expanded in θ² (not branching on the axis=ω/θ division) for
    // the same reason as so3.hpp's exp(): V is analytically smooth at ω=0
    // (V(0)=I exactly, with a well-defined first-order dependence on ω),
    // but the old `if (angle < eps) { t = vel; }` shortcut returned a
    // constant, ω-independent value there — correct at ω=0 itself, but
    // silently zero-derivative under Dual<T> at the identity, again the
    // single most common optimization starting point.
    static Mat3 translation_jacobian(const Vec3& omega) {
        using std::sin; using std::cos; using std::sqrt;
        T theta2 = omega.dot(omega);
        auto angle = sqrt(theta2);
        auto K = skew(omega);

        T b, c;
        if (angle < epsilon<T>()) {
            b = T{0.5} - theta2 / T{24};
            c = T{1} / T{6} - theta2 / T{120};
        } else {
            b = (T{1} - cos(angle)) / theta2;
            c = (angle - sin(angle)) / (theta2 * angle);
        }

        return Mat3::identity() + K * b + (K * K) * c;
    }
};

// Concept checks — default double instantiation and the Dual<double>
// instantiation that's the actual point of templating (see so3.hpp).
static_assert(Group<SE3<>>);
static_assert(LieGroup<SE3<>>);
static_assert(Group<SE3<Dual<double>>>);
static_assert(LieGroup<SE3<Dual<double>>>);

} // namespace spatium::algebra
