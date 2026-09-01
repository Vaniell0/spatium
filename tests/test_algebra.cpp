#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/algebra/concepts.hpp>
#include <spatium/algebra/functions.hpp>
#include <spatium/algebra/groups/so3.hpp>
#include <spatium/algebra/groups/se3.hpp>
#include <spatium/algebra/verify.hpp>
#include <numbers>

using namespace spatium;
using namespace spatium::algebra;
using Catch::Matchers::WithinAbs;

// ── Concept checks ─────────────────────────────────────────────

TEST_CASE("SO3 satisfies LieGroup", "[algebra]") {
    static_assert(Group<SO3<double>>);
    static_assert(LieGroup<SO3<double>>);
    SUCCEED();
}

TEST_CASE("SE3 satisfies LieGroup", "[algebra]") {
    static_assert(Group<SE3<double>>);
    static_assert(LieGroup<SE3<double>>);
    SUCCEED();
}

// ── SO(3) ──────────────────────────────────────────────────────

TEST_CASE("SO3 identity", "[algebra]") {
    SO3<double> so3;
    auto I = so3.identity();
    CHECK(I == Mat3::identity());
}

TEST_CASE("SO3 rotation around Z", "[algebra]") {
    SO3<double> so3;
    auto R = so3.rz(std::numbers::pi / 2);
    auto p = so3.act(R, Vec3{1, 0, 0});
    CHECK_THAT(p[0], WithinAbs(0.0, 1e-10));
    CHECK_THAT(p[1], WithinAbs(1.0, 1e-10));
    CHECK_THAT(p[2], WithinAbs(0.0, 1e-10));
}

TEST_CASE("SO3 inverse = transpose", "[algebra]") {
    SO3<double> so3;
    auto R = so3.rz(0.7);
    auto inv = so3.inverse(R);
    auto prod = so3.compose(R, inv);
    auto I = so3.identity();
    for (std::size_t i = 0; i < 9; ++i)
        CHECK_THAT(prod.data[i], WithinAbs(I.data[i], 1e-10));
}

TEST_CASE("SO3 exp/log roundtrip", "[algebra]") {
    SO3<double> so3;
    Vec3 omega{0.3, 0.5, -0.7}; // axis-angle
    auto R = so3.exp(omega);
    auto recovered = so3.log(R);
    CHECK_THAT(recovered[0], WithinAbs(omega[0], 1e-8));
    CHECK_THAT(recovered[1], WithinAbs(omega[1], 1e-8));
    CHECK_THAT(recovered[2], WithinAbs(omega[2], 1e-8));
}

TEST_CASE("SO3 exp zero = identity", "[algebra]") {
    SO3<double> so3;
    auto R = so3.exp(Vec3{0, 0, 0});
    CHECK(R == so3.identity());
}

TEST_CASE("SO3 composition = chained rotation", "[algebra]") {
    SO3<double> so3;
    auto Rx = so3.rx(std::numbers::pi / 4);
    auto Ry = so3.ry(std::numbers::pi / 4);
    auto Rxy = so3.compose(Ry, Rx);
    auto p = so3.act(Rxy, Vec3{1, 0, 0});
    // Should be rotated by both
    CHECK_THAT(p.norm(), WithinAbs(1.0, 1e-10)); // length preserved
}

TEST_CASE("SO3<Dual<double>> differentiates through rotation", "[algebra]") {
    // Rotating (1,0,0) by angle θ about Z gives (cos θ, sin θ, 0). Seeding θ
    // as a Dual variable recovers -sinθ/cosθ as the derivatives of the
    // rotated point's x/y coordinates with zero hand-derived Jacobian —
    // this is the actual point of templating SO3 on Scalar (pose-graph /
    // SLAM-style optimization needs exactly this).
    SO3<Dual<double>> so3;
    double theta0 = 0.6;
    auto theta = Dual<double>::variable(theta0);
    auto R = so3.rz(theta);
    Vec<Dual<double>, 3> p0{Dual<double>{1.0}, Dual<double>{0.0}, Dual<double>{0.0}};
    auto p = so3.act(R, p0);

    CHECK_THAT(p[0].value, WithinAbs(std::cos(theta0), 1e-10));
    CHECK_THAT(p[1].value, WithinAbs(std::sin(theta0), 1e-10));
    CHECK_THAT(p[0].deriv, WithinAbs(-std::sin(theta0), 1e-10));
    CHECK_THAT(p[1].deriv, WithinAbs(std::cos(theta0), 1e-10));
}

TEST_CASE("SO3<Dual<double>> log/exp roundtrip differentiates correctly", "[algebra]") {
    // log(exp(v)) == v exactly, so d(log(exp(v)))/dv0 must be exactly 1 by
    // the chain rule. Exercises log()'s acos()/sin() calls under Dual — real
    // motivation for adding Dual<T> acos()/tan() in this change, and a
    // non-tautological check: a wrong derivative formula anywhere in the
    // exp/log chain would break this composed identity, not just the value.
    SO3<Dual<double>> so3;
    double a1 = 0.3, a2 = -0.4, a3 = 0.5;
    auto v0 = Dual<double>::variable(a1);
    Vec<Dual<double>, 3> v{v0, Dual<double>{a2}, Dual<double>{a3}};
    auto R = so3.exp(v);
    auto recovered = so3.log(R);

    CHECK_THAT(recovered[0].value, WithinAbs(a1, 1e-8));
    CHECK_THAT(recovered[0].deriv, WithinAbs(1.0, 1e-6));
}

TEST_CASE("SO3<Dual<double>> exp differentiates AT v=0, not just away from it", "[algebra]") {
    // Regression test for a real bug: exp()'s old angle<eps special case
    // returned a v-INDEPENDENT identity() -- correct in VALUE exactly at
    // v=0, but silently zero-derivative there under Dual<T>, exactly at the
    // single most common optimization starting point (a standalone probe
    // confirmed gradient(residual, Vec3{0,0,0}) came back hard [0,0,0] pre-
    // fix on a rotation-averaging objective, vs. a real nonzero value one
    // step away). Checks against the known first-order Taylor expansion
    // exp(v) ~= I + skew(v): skew({1,0,0}) sets (1,2)=-1 and (2,1)=1, all
    // else 0 -- exactly the derivatives seeding v[0] should produce.
    SO3<Dual<double>> so3;
    Vec<Dual<double>, 3> v{Dual<double>::variable(0.0), Dual<double>{0.0}, Dual<double>{0.0}};
    auto R = so3.exp(v);

    CHECK_THAT(R(1, 2).deriv, WithinAbs(-1.0, 1e-9));
    CHECK_THAT(R(2, 1).deriv, WithinAbs(1.0, 1e-9));
    CHECK_THAT(R(0, 1).deriv, WithinAbs(0.0, 1e-9));
    CHECK_THAT(R(0, 0).value, WithinAbs(1.0, 1e-12)); // still exactly identity in value
}

TEST_CASE("SO3<Dual<double>> log/exp roundtrip differentiates AT v=0", "[algebra]") {
    // Same regression as the exp-at-origin test above, but through both
    // exp() AND log() at once (log()'s own near-identity branch had the
    // identical bug — acos'(1) diverges, and the old `angle<eps -> return
    // AlgebraType{}` shortcut zeroed the derivative there too). All three
    // components at exactly zero reproduces the original failure mode
    // exactly: a rotation-averaging optimizer starting from identity.
    SO3<Dual<double>> so3;
    auto v0 = Dual<double>::variable(0.0);
    Vec<Dual<double>, 3> v{v0, Dual<double>{0.0}, Dual<double>{0.0}};
    auto R = so3.exp(v);
    auto recovered = so3.log(R);

    CHECK_THAT(recovered[0].value, WithinAbs(0.0, 1e-12));
    CHECK_THAT(recovered[0].deriv, WithinAbs(1.0, 1e-6));
}

TEST_CASE("SE3 exp known-good: rotation + translation vs. brute-force matrix exponential", "[algebra]") {
    // Regression test for a real, pre-existing correctness bug found while
    // making exp()/log() Dual-differentiable: the old translation_jacobian
    // (V matrix) used skew(omega/angle) (the UNIT axis) with coefficients
    // that are only valid for skew(omega) (unnormalized) -- silently wrong
    // for any nonzero rotation combined with a translation. No prior test
    // caught it: "SE3 exp/log roundtrip" is self-consistent under either
    // convention, and the other "known-good" test builds T via from_Rt()
    // directly, bypassing exp()'s V matrix entirely. Expected values are an
    // independent ground truth: a 40-term brute-force Taylor sum of the
    // 4x4 se(3) generator's matrix exponential, unrelated to so3.hpp/
    // se3.hpp's own formulas.
    SE3<double> se3;
    Vec3 omega{0.4, -0.25, 0.15};
    Vec3 vel{1.0, 0.5, -0.3};
    auto T = se3.exp(Vec<double, 6>{omega[0], omega[1], omega[2], vel[0], vel[1], vel[2]});
    auto t = SE3<double>::translation_of(T);

    CHECK_THAT(t[0], WithinAbs(0.9748105589992107, 1e-10));
    CHECK_THAT(t[1], WithinAbs(0.6026315246925509, 1e-10));
    CHECK_THAT(t[2], WithinAbs(-0.06177561617697668, 1e-10));
}

TEST_CASE("SO3 verify group axioms", "[algebra]") {
    SO3<double> so3;
    std::array samples = {
        so3.rx(0.3),
        so3.ry(0.7),
        so3.rz(-0.5),
        so3.from_axis_angle(Vec3{1, 1, 1}, 0.8),
    };
    auto result = verify_matrix_group(so3, std::span{samples});
    if (!result) FAIL(result.failure);
    CHECK(result);
}

// ── SE(3) ──────────────────────────────────────────────────────

TEST_CASE("SE3 identity", "[algebra]") {
    SE3<double> se3;
    auto I = se3.identity();
    CHECK(I == Mat4::identity());
}

TEST_CASE("SE3 pure translation", "[algebra]") {
    SE3<double> se3;
    auto T = SE3<double>::from_Rt(Mat3::identity(), Vec3{10, 0, 0});
    auto p = se3.act(T, Vec3{1, 2, 3});
    CHECK(p == Vec3{11, 2, 3});
}

TEST_CASE("SE3 rotation + translation", "[algebra]") {
    SE3<double> se3;
    auto R = se3.so3.rz(std::numbers::pi / 2);
    auto T = SE3<double>::from_Rt(R, Vec3{1, 0, 0});
    auto p = se3.act(T, Vec3{1, 0, 0});
    CHECK_THAT(p[0], WithinAbs(1.0, 1e-10)); // rotation(1,0,0)→(0,1,0) + translate(1,0,0)
    CHECK_THAT(p[1], WithinAbs(1.0, 1e-10));
    CHECK_THAT(p[2], WithinAbs(0.0, 1e-10));
}

TEST_CASE("SE3 inverse", "[algebra]") {
    SE3<double> se3;
    auto R = se3.so3.from_axis_angle(Vec3{1, 0, 0}, 0.5);
    auto T = SE3<double>::from_Rt(R, Vec3{3, 4, 5});
    auto T_inv = se3.inverse(T);
    auto prod = se3.compose(T, T_inv);
    auto I = se3.identity();
    for (std::size_t i = 0; i < 16; ++i)
        CHECK_THAT(prod.data[i], WithinAbs(I.data[i], 1e-10));
}

TEST_CASE("SE3 exp/log roundtrip", "[algebra]") {
    SE3<double> se3;
    Vec<double, 6> xi{0.1, 0.2, 0.3, 1.0, 2.0, 3.0}; // omega + velocity
    auto T = se3.exp(xi);
    auto recovered = se3.log(T);
    for (int i = 0; i < 6; ++i)
        CHECK_THAT(recovered[i], WithinAbs(xi[i], 1e-6));
}

TEST_CASE("SE3 verify group axioms", "[algebra]") {
    SE3<double> se3;
    std::array samples = {
        SE3<double>::from_Rt(se3.so3.rx(0.3), Vec3{1, 0, 0}),
        SE3<double>::from_Rt(se3.so3.ry(0.5), Vec3{0, 2, 0}),
        SE3<double>::from_Rt(se3.so3.rz(-0.7), Vec3{0, 0, 3}),
    };
    auto result = verify_matrix_group(se3, std::span{samples});
    if (!result) FAIL(result.failure);
    CHECK(result);
}

// ── Generic algebra functions ─────────────────────────────────

TEST_CASE("power: SO3 rotation squared", "[algebra]") {
    SO3<double> so3;
    auto R = so3.rz(std::numbers::pi / 4);
    auto R2 = algebra::power(so3, R, 2);
    // π/4 * 2 = π/2 rotation
    auto p = so3.act(R2, Vec3{1, 0, 0});
    CHECK_THAT(p[0], WithinAbs(0.0, 1e-10));
    CHECK_THAT(p[1], WithinAbs(1.0, 1e-10));
}

TEST_CASE("power: n=0 gives identity", "[algebra]") {
    SO3<double> so3;
    auto R = so3.rz(0.7);
    auto I = algebra::power(so3, R, 0);
    CHECK(I == so3.identity());
}

TEST_CASE("power: negative exponent uses inverse", "[algebra]") {
    SO3<double> so3;
    auto R = so3.rz(std::numbers::pi / 3);
    auto R_inv = algebra::power(so3, R, -1);
    auto prod = so3.compose(R, R_inv);
    auto I = so3.identity();
    for (std::size_t i = 0; i < 9; ++i)
        CHECK_THAT(prod.data[i], WithinAbs(I.data[i], 1e-10));
}

TEST_CASE("commutator: SO3 commutator of Rx, Ry is non-trivial", "[algebra]") {
    SO3<double> so3;
    auto Rx = so3.rx(0.3);
    auto Ry = so3.ry(0.3);
    auto comm = algebra::commutator(so3, Rx, Ry);
    // Non-abelian: commutator should not be identity
    CHECK_FALSE(comm == so3.identity());
    // But it should still be a rotation (det=1, orthogonal)
    auto p = so3.act(comm, Vec3{1, 0, 0});
    CHECK_THAT(p.norm(), WithinAbs(1.0, 1e-10));
}

TEST_CASE("adjoint: SO3 adjoint action", "[algebra]") {
    SO3<double> so3;
    auto R = so3.rz(std::numbers::pi / 2);
    Vec3 v{1, 0, 0}; // x-axis angular velocity
    auto ad = algebra::adjoint(so3, R, v);
    // Rotating x-axis by 90° around z gives y-axis
    CHECK_THAT(ad[0], WithinAbs(0.0, 1e-4));
    CHECK_THAT(ad[1], WithinAbs(1.0, 1e-4));
    CHECK_THAT(ad[2], WithinAbs(0.0, 1e-4));
}
