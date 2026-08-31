#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/algebra/quaternion.hpp>
#include <cmath>
#include <numbers>

using namespace spatium;
using Catch::Matchers::WithinAbs;

TEST_CASE("Quaternion identity rotation", "[quaternion]") {
    Quat q;  // default = identity
    Vec3 v{1, 2, 3};
    auto r = q.rotate(v);
    CHECK_THAT(r[0], WithinAbs(1.0, 1e-12));
    CHECK_THAT(r[1], WithinAbs(2.0, 1e-12));
    CHECK_THAT(r[2], WithinAbs(3.0, 1e-12));
}

TEST_CASE("Quaternion 90deg rotation around Z", "[quaternion]") {
    auto q = Quat::from_axis_angle(Vec3{0, 0, 1}, std::numbers::pi / 2);
    auto r = q.rotate(Vec3{1, 0, 0});
    CHECK_THAT(r[0], WithinAbs(0.0, 1e-12));
    CHECK_THAT(r[1], WithinAbs(1.0, 1e-12));
    CHECK_THAT(r[2], WithinAbs(0.0, 1e-12));
}

TEST_CASE("Quaternion 180deg rotation around Y", "[quaternion]") {
    auto q = Quat::from_axis_angle(Vec3{0, 1, 0}, std::numbers::pi);
    auto r = q.rotate(Vec3{1, 0, 0});
    CHECK_THAT(r[0], WithinAbs(-1.0, 1e-10));
    CHECK_THAT(r[1], WithinAbs(0.0, 1e-10));
    CHECK_THAT(r[2], WithinAbs(0.0, 1e-10));
}

TEST_CASE("Quaternion composition", "[quaternion]") {
    auto q1 = Quat::from_axis_angle(Vec3{0, 0, 1}, std::numbers::pi / 2);
    auto q2 = Quat::from_axis_angle(Vec3{0, 0, 1}, std::numbers::pi / 2);
    auto q = q1 * q2;  // 180 deg around Z
    auto r = q.rotate(Vec3{1, 0, 0});
    CHECK_THAT(r[0], WithinAbs(-1.0, 1e-10));
    CHECK_THAT(r[1], WithinAbs(0.0, 1e-10));
}

TEST_CASE("Quaternion norm", "[quaternion]") {
    auto q = Quat::from_axis_angle(Vec3{1, 1, 1}, 1.0);
    CHECK_THAT(q.norm(), WithinAbs(1.0, 1e-12));
}

TEST_CASE("Quaternion inverse", "[quaternion]") {
    auto q = Quat::from_axis_angle(Vec3{1, 0, 0}, 0.7);
    auto qi = q.inverse();
    REQUIRE(qi.has_value());
    auto id = q * *qi;
    CHECK_THAT(id.w, WithinAbs(1.0, 1e-12));
    CHECK_THAT(id.x, WithinAbs(0.0, 1e-12));
}

TEST_CASE("Quaternion to/from matrix roundtrip", "[quaternion]") {
    auto q = Quat::from_axis_angle(Vec3{1, 2, 3}, 1.5);
    auto m = q.to_matrix();
    auto q2 = Quat::from_matrix(m);
    // May differ by sign (q and -q represent same rotation)
    auto diff = std::min(
        std::abs(q.w - q2.w) + std::abs(q.x - q2.x) + std::abs(q.y - q2.y) + std::abs(q.z - q2.z),
        std::abs(q.w + q2.w) + std::abs(q.x + q2.x) + std::abs(q.y + q2.y) + std::abs(q.z + q2.z)
    );
    CHECK_THAT(diff, WithinAbs(0.0, 1e-10));
}

TEST_CASE("Quaternion slerp endpoints", "[quaternion]") {
    auto a = Quat::from_axis_angle(Vec3{0, 0, 1}, 0.0);
    auto b = Quat::from_axis_angle(Vec3{0, 0, 1}, std::numbers::pi / 2);
    auto s0 = Quat::slerp(a, b, 0.0);
    auto s1 = Quat::slerp(a, b, 1.0);
    CHECK_THAT(s0.w, WithinAbs(a.w, 1e-10));
    CHECK_THAT(s1.w, WithinAbs(b.w, 1e-10));
}

TEST_CASE("Quaternion slerp midpoint", "[quaternion]") {
    auto a = Quat::from_axis_angle(Vec3{0, 0, 1}, 0.0);
    auto b = Quat::from_axis_angle(Vec3{0, 0, 1}, std::numbers::pi / 2);
    auto mid = Quat::slerp(a, b, 0.5);
    // Should be 45 degrees
    auto r = mid.rotate(Vec3{1, 0, 0});
    CHECK_THAT(r[0], WithinAbs(std::cos(std::numbers::pi / 4), 1e-10));
    CHECK_THAT(r[1], WithinAbs(std::sin(std::numbers::pi / 4), 1e-10));
}

TEST_CASE("Quaternion axis-angle roundtrip", "[quaternion]") {
    Vec3 axis{0.5, 0.7, 0.3};
    double angle = 1.2;
    auto q = Quat::from_axis_angle(axis, angle);
    auto [ax, an] = q.to_axis_angle();
    CHECK_THAT(an, WithinAbs(angle, 1e-10));
    // Axis should be parallel to original (normalized)
    auto n = axis.normalized();
    auto dot = std::abs(ax.dot(n));
    CHECK_THAT(dot, WithinAbs(1.0, 1e-10));
}

// ── Compile-time arithmetic guard ─────────────────────────────
// Pins down the constexpr surface of Quaternion: identity,
// Hamilton product, scalar mul, addition, conjugate,
// norm_squared, equality, and rotate(). Regressions that
// accidentally drop constexpr fail to compile here.

namespace {

constexpr Quat qid;            // identity (1, 0, 0, 0)
constexpr Quat qa{1.0, 2.0, 3.0, 4.0};
constexpr Quat qb{5.0, 6.0, 7.0, 8.0};

static_assert(qid.w == 1.0 && qid.x == 0.0 && qid.y == 0.0 && qid.z == 0.0);

// identity * anything = anything
static_assert(qid * qa == qa);
static_assert(qa * qid == qa);

// Hamilton product reference: (1,2,3,4)*(5,6,7,8)
//   w = 1*5 - 2*6 - 3*7 - 4*8 = 5 - 12 - 21 - 32 = -60
//   x = 1*6 + 2*5 + 3*8 - 4*7 = 6 + 10 + 24 - 28 = 12
//   y = 1*7 - 2*8 + 3*5 + 4*6 = 7 - 16 + 15 + 24 = 30
//   z = 1*8 + 2*7 - 3*6 + 4*5 = 8 + 14 - 18 + 20 = 24
static_assert((qa * qb) == Quat{-60.0, 12.0, 30.0, 24.0});

// scalar multiplication is symmetric
static_assert(qa * 2.0 == Quat{2.0, 4.0, 6.0, 8.0});
static_assert(2.0 * qa == Quat{2.0, 4.0, 6.0, 8.0});

// addition
static_assert(qa + qb == Quat{6.0, 8.0, 10.0, 12.0});

// conjugate flips the vector part
static_assert(qa.conjugate() == Quat{1.0, -2.0, -3.0, -4.0});

// norm_squared: 1 + 4 + 9 + 16 = 30
static_assert(qa.norm_squared() == 30.0);

// rotate is constexpr — identity quaternion leaves vectors alone
static_assert(qid.rotate(Vec3{1.0, 2.0, 3.0}) == Vec3{1.0, 2.0, 3.0});

} // namespace
