#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/geometry/hyperplane.hpp>

using namespace spatium;
using namespace spatium::geometry;
using Catch::Matchers::WithinAbs;

TEST_CASE("Hyperplane from normal and point", "[hyperplane]") {
    auto plane = Plane3::from_normal_and_point(Vec3{0.0, 0.0, 2.0}, Vec3{0.0, 0.0, 5.0});
    REQUIRE(plane.has_value());
    CHECK_THAT(plane->normal.norm(), WithinAbs(1.0, 1e-12));
    CHECK(plane->normal == Vec3{0.0, 0.0, 1.0});
    CHECK_THAT(plane->offset, WithinAbs(5.0, 1e-12));
}

TEST_CASE("Hyperplane from zero normal fails", "[hyperplane]") {
    auto plane = Plane3::from_normal_and_point(Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 2.0, 3.0});
    CHECK_FALSE(plane.has_value());
}

TEST_CASE("Plane from 3 points", "[hyperplane]") {
    auto plane = Plane3::from_points(
        Vec3{0.0, 0.0, 0.0},
        Vec3{1.0, 0.0, 0.0},
        Vec3{0.0, 1.0, 0.0}
    );
    REQUIRE(plane.has_value());
    // Normal should be (0,0,1) or (0,0,-1)
    CHECK_THAT(std::abs(plane->normal[2]), WithinAbs(1.0, 1e-12));
    CHECK_THAT(plane->offset, WithinAbs(0.0, 1e-12));
}

TEST_CASE("Plane from collinear points fails", "[hyperplane]") {
    auto plane = Plane3::from_points(
        Vec3{0.0, 0.0, 0.0},
        Vec3{1.0, 0.0, 0.0},
        Vec3{2.0, 0.0, 0.0}
    );
    CHECK_FALSE(plane.has_value());
}

TEST_CASE("Hyperplane signed distance", "[hyperplane]") {
    auto plane = *Plane3::from_normal_and_point(Vec3{0.0, 0.0, 1.0}, Vec3{0.0, 0.0, 3.0});

    CHECK_THAT(plane.signed_distance(Vec3{0.0, 0.0, 5.0}), WithinAbs(2.0, 1e-12));
    CHECK_THAT(plane.signed_distance(Vec3{0.0, 0.0, 1.0}), WithinAbs(-2.0, 1e-12));
    CHECK_THAT(plane.signed_distance(Vec3{0.0, 0.0, 3.0}), WithinAbs(0.0, 1e-12));
}

TEST_CASE("Hyperplane distance (unsigned)", "[hyperplane]") {
    auto plane = *Plane3::from_normal_and_point(Vec3{0.0, 0.0, 1.0}, Vec3{0.0, 0.0, 3.0});
    CHECK_THAT(plane.distance(Vec3{0.0, 0.0, 1.0}), WithinAbs(2.0, 1e-12));
}

TEST_CASE("Hyperplane project", "[hyperplane]") {
    auto plane = *Plane3::from_normal_and_point(Vec3{0.0, 0.0, 1.0}, Vec3{0.0, 0.0, 3.0});
    auto proj = plane.project(Vec3{5.0, 7.0, 10.0});
    CHECK_THAT(proj[0], WithinAbs(5.0, 1e-12));
    CHECK_THAT(proj[1], WithinAbs(7.0, 1e-12));
    CHECK_THAT(proj[2], WithinAbs(3.0, 1e-12));
}

TEST_CASE("Hyperplane contains", "[hyperplane]") {
    auto plane = *Plane3::from_normal_and_point(Vec3{0.0, 0.0, 1.0}, Vec3{0.0, 0.0, 3.0});
    CHECK(plane.contains(Vec3{5.0, 7.0, 3.0}));
    CHECK_FALSE(plane.contains(Vec3{5.0, 7.0, 4.0}));
}

TEST_CASE("Hyperplane side", "[hyperplane]") {
    auto plane = *Plane3::from_normal_and_point(Vec3{0.0, 0.0, 1.0}, Vec3{0.0, 0.0, 0.0});
    CHECK(plane.side(Vec3{0.0, 0.0, 1.0}) == 1);
    CHECK(plane.side(Vec3{0.0, 0.0, -1.0}) == -1);
    CHECK(plane.side(Vec3{0.0, 0.0, 0.0}) == 0);
}

TEST_CASE("2D Hyperplane (line)", "[hyperplane]") {
    auto line = Hyperplane<2>::from_points(Vec2{0.0, 0.0}, Vec2{1.0, 0.0});
    REQUIRE(line.has_value());
    // Normal should be perpendicular to x-axis
    CHECK_THAT(std::abs(line->normal[1]), WithinAbs(1.0, 1e-12));
    CHECK_THAT(line->distance(Vec2{0.0, 5.0}), WithinAbs(5.0, 1e-12));
}

TEST_CASE("Hyperplane concept satisfaction", "[hyperplane]") {
    static_assert(Shape<Plane3>);
    static_assert(DistanceQueryable<Plane3>);
    SUCCEED();
}
