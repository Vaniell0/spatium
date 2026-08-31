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
    static_assert(Group<SO3>);
    static_assert(LieGroup<SO3>);
    SUCCEED();
}

TEST_CASE("SE3 satisfies LieGroup", "[algebra]") {
    static_assert(Group<SE3>);
    static_assert(LieGroup<SE3>);
    SUCCEED();
}

// ── SO(3) ──────────────────────────────────────────────────────

TEST_CASE("SO3 identity", "[algebra]") {
    SO3 so3;
    auto I = so3.identity();
    CHECK(I == Mat3::identity());
}

TEST_CASE("SO3 rotation around Z", "[algebra]") {
    SO3 so3;
    auto R = so3.rz(std::numbers::pi / 2);
    auto p = so3.act(R, Vec3{1, 0, 0});
    CHECK_THAT(p[0], WithinAbs(0.0, 1e-10));
    CHECK_THAT(p[1], WithinAbs(1.0, 1e-10));
    CHECK_THAT(p[2], WithinAbs(0.0, 1e-10));
}

TEST_CASE("SO3 inverse = transpose", "[algebra]") {
    SO3 so3;
    auto R = so3.rz(0.7);
    auto inv = so3.inverse(R);
    auto prod = so3.compose(R, inv);
    auto I = so3.identity();
    for (std::size_t i = 0; i < 9; ++i)
        CHECK_THAT(prod.data[i], WithinAbs(I.data[i], 1e-10));
}

TEST_CASE("SO3 exp/log roundtrip", "[algebra]") {
    SO3 so3;
    Vec3 omega{0.3, 0.5, -0.7}; // axis-angle
    auto R = so3.exp(omega);
    auto recovered = so3.log(R);
    CHECK_THAT(recovered[0], WithinAbs(omega[0], 1e-8));
    CHECK_THAT(recovered[1], WithinAbs(omega[1], 1e-8));
    CHECK_THAT(recovered[2], WithinAbs(omega[2], 1e-8));
}

TEST_CASE("SO3 exp zero = identity", "[algebra]") {
    SO3 so3;
    auto R = so3.exp(Vec3{0, 0, 0});
    CHECK(R == so3.identity());
}

TEST_CASE("SO3 composition = chained rotation", "[algebra]") {
    SO3 so3;
    auto Rx = so3.rx(std::numbers::pi / 4);
    auto Ry = so3.ry(std::numbers::pi / 4);
    auto Rxy = so3.compose(Ry, Rx);
    auto p = so3.act(Rxy, Vec3{1, 0, 0});
    // Should be rotated by both
    CHECK_THAT(p.norm(), WithinAbs(1.0, 1e-10)); // length preserved
}

TEST_CASE("SO3 verify group axioms", "[algebra]") {
    SO3 so3;
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
    SE3 se3;
    auto I = se3.identity();
    CHECK(I == Mat4::identity());
}

TEST_CASE("SE3 pure translation", "[algebra]") {
    SE3 se3;
    auto T = SE3::from_Rt(Mat3::identity(), Vec3{10, 0, 0});
    auto p = se3.act(T, Vec3{1, 2, 3});
    CHECK(p == Vec3{11, 2, 3});
}

TEST_CASE("SE3 rotation + translation", "[algebra]") {
    SE3 se3;
    auto R = se3.so3.rz(std::numbers::pi / 2);
    auto T = SE3::from_Rt(R, Vec3{1, 0, 0});
    auto p = se3.act(T, Vec3{1, 0, 0});
    CHECK_THAT(p[0], WithinAbs(1.0, 1e-10)); // rotation(1,0,0)→(0,1,0) + translate(1,0,0)
    CHECK_THAT(p[1], WithinAbs(1.0, 1e-10));
    CHECK_THAT(p[2], WithinAbs(0.0, 1e-10));
}

TEST_CASE("SE3 inverse", "[algebra]") {
    SE3 se3;
    auto R = se3.so3.from_axis_angle(Vec3{1, 0, 0}, 0.5);
    auto T = SE3::from_Rt(R, Vec3{3, 4, 5});
    auto T_inv = se3.inverse(T);
    auto prod = se3.compose(T, T_inv);
    auto I = se3.identity();
    for (std::size_t i = 0; i < 16; ++i)
        CHECK_THAT(prod.data[i], WithinAbs(I.data[i], 1e-10));
}

TEST_CASE("SE3 exp/log roundtrip", "[algebra]") {
    SE3 se3;
    Vec<double, 6> xi{0.1, 0.2, 0.3, 1.0, 2.0, 3.0}; // omega + velocity
    auto T = se3.exp(xi);
    auto recovered = se3.log(T);
    for (int i = 0; i < 6; ++i)
        CHECK_THAT(recovered[i], WithinAbs(xi[i], 1e-6));
}

TEST_CASE("SE3 verify group axioms", "[algebra]") {
    SE3 se3;
    std::array samples = {
        SE3::from_Rt(se3.so3.rx(0.3), Vec3{1, 0, 0}),
        SE3::from_Rt(se3.so3.ry(0.5), Vec3{0, 2, 0}),
        SE3::from_Rt(se3.so3.rz(-0.7), Vec3{0, 0, 3}),
    };
    auto result = verify_matrix_group(se3, std::span{samples});
    if (!result) FAIL(result.failure);
    CHECK(result);
}

// ── Generic algebra functions ─────────────────────────────────

TEST_CASE("power: SO3 rotation squared", "[algebra]") {
    SO3 so3;
    auto R = so3.rz(std::numbers::pi / 4);
    auto R2 = algebra::power(so3, R, 2);
    // π/4 * 2 = π/2 rotation
    auto p = so3.act(R2, Vec3{1, 0, 0});
    CHECK_THAT(p[0], WithinAbs(0.0, 1e-10));
    CHECK_THAT(p[1], WithinAbs(1.0, 1e-10));
}

TEST_CASE("power: n=0 gives identity", "[algebra]") {
    SO3 so3;
    auto R = so3.rz(0.7);
    auto I = algebra::power(so3, R, 0);
    CHECK(I == so3.identity());
}

TEST_CASE("power: negative exponent uses inverse", "[algebra]") {
    SO3 so3;
    auto R = so3.rz(std::numbers::pi / 3);
    auto R_inv = algebra::power(so3, R, -1);
    auto prod = so3.compose(R, R_inv);
    auto I = so3.identity();
    for (std::size_t i = 0; i < 9; ++i)
        CHECK_THAT(prod.data[i], WithinAbs(I.data[i], 1e-10));
}

TEST_CASE("commutator: SO3 commutator of Rx, Ry is non-trivial", "[algebra]") {
    SO3 so3;
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
    SO3 so3;
    auto R = so3.rz(std::numbers::pi / 2);
    Vec3 v{1, 0, 0}; // x-axis angular velocity
    auto ad = algebra::adjoint(so3, R, v);
    // Rotating x-axis by 90° around z gives y-axis
    CHECK_THAT(ad[0], WithinAbs(0.0, 1e-4));
    CHECK_THAT(ad[1], WithinAbs(1.0, 1e-4));
    CHECK_THAT(ad[2], WithinAbs(0.0, 1e-4));
}
