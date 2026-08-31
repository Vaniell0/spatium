#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/algebra/literals.hpp>

using namespace spatium;
using namespace spatium::literals;
using Catch::Matchers::WithinAbs;

TEST_CASE("_x literal creates X-axis vector", "[literals]") {
    Vec3 v = 3.0_x;
    CHECK_THAT(v[0], WithinAbs(3.0, 1e-15));
    CHECK_THAT(v[1], WithinAbs(0.0, 1e-15));
    CHECK_THAT(v[2], WithinAbs(0.0, 1e-15));
}

TEST_CASE("_y literal creates Y-axis vector", "[literals]") {
    Vec3 v = 2.0_y;
    CHECK_THAT(v[0], WithinAbs(0.0, 1e-15));
    CHECK_THAT(v[1], WithinAbs(2.0, 1e-15));
    CHECK_THAT(v[2], WithinAbs(0.0, 1e-15));
}

TEST_CASE("_z literal creates Z-axis vector", "[literals]") {
    Vec3 v = 1.0_z;
    CHECK_THAT(v[0], WithinAbs(0.0, 1e-15));
    CHECK_THAT(v[1], WithinAbs(0.0, 1e-15));
    CHECK_THAT(v[2], WithinAbs(1.0, 1e-15));
}

TEST_CASE("Combined axis literals produce correct Vec3", "[literals]") {
    Vec3 v = 3.0_x + 2.0_y + 1.0_z;
    CHECK_THAT(v[0], WithinAbs(3.0, 1e-15));
    CHECK_THAT(v[1], WithinAbs(2.0, 1e-15));
    CHECK_THAT(v[2], WithinAbs(1.0, 1e-15));
}

TEST_CASE("Integer axis literals work", "[literals]") {
    Vec3 v = 5_x + 3_y + 1_z;
    CHECK_THAT(v[0], WithinAbs(5.0, 1e-15));
    CHECK_THAT(v[1], WithinAbs(3.0, 1e-15));
    CHECK_THAT(v[2], WithinAbs(1.0, 1e-15));
}

TEST_CASE("Partial axis combination", "[literals]") {
    Vec3 v = 1.0_x + 2.0_y;
    CHECK_THAT(v[0], WithinAbs(1.0, 1e-15));
    CHECK_THAT(v[1], WithinAbs(2.0, 1e-15));
    CHECK_THAT(v[2], WithinAbs(0.0, 1e-15));
}

TEST_CASE("Axis literals are constexpr", "[literals]") {
    constexpr Vec3 v = 1.0_x;
    static_assert(v[0] == 1.0);
    static_assert(v[1] == 0.0);
    static_assert(v[2] == 0.0);
    SUCCEED();
}

TEST_CASE("_x2 / _y2 build Vec2", "[literals]") {
    Vec2 v = 1.0_x2 + 2.0_y2;
    CHECK_THAT(v[0], WithinAbs(1.0, 1e-15));
    CHECK_THAT(v[1], WithinAbs(2.0, 1e-15));
}

TEST_CASE("_w produces Vec4 with fourth component", "[literals]") {
    Vec4 w = 4.0_w;
    CHECK_THAT(w[0], WithinAbs(0.0, 1e-15));
    CHECK_THAT(w[1], WithinAbs(0.0, 1e-15));
    CHECK_THAT(w[2], WithinAbs(0.0, 1e-15));
    CHECK_THAT(w[3], WithinAbs(4.0, 1e-15));
}

TEST_CASE("Vec4 axis literals compose to full 4D vector", "[literals]") {
    Vec4 v = 1.0_x4 + 2.0_y4 + 3.0_z4 + 4.0_w;
    CHECK_THAT(v[0], WithinAbs(1.0, 1e-15));
    CHECK_THAT(v[1], WithinAbs(2.0, 1e-15));
    CHECK_THAT(v[2], WithinAbs(3.0, 1e-15));
    CHECK_THAT(v[3], WithinAbs(4.0, 1e-15));
}

TEST_CASE("Integer Vec2/Vec4 literals work", "[literals]") {
    Vec2 v2 = 1_x2 + 2_y2;
    CHECK(v2[0] == 1.0);
    CHECK(v2[1] == 2.0);

    Vec4 v4 = 1_x4 + 2_y4 + 3_z4 + 4_w;
    CHECK(v4[3] == 4.0);
}
