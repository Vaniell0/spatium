#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/geometry/circle.hpp>
#include <numbers>

using namespace spatium;
using namespace spatium::geometry;
using Catch::Matchers::WithinAbs;

TEST_CASE("Circle circumference", "[circle]") {
    Circle3 c{Vec3{0.0, 0.0, 0.0}, 5.0, Vec3{0.0, 0.0, 1.0}};
    CHECK_THAT(c.circumference(), WithinAbs(2.0 * std::numbers::pi * 5.0, 1e-10));
}

TEST_CASE("Circle 2D distance", "[circle]") {
    Circle2 c{Vec2{0.0, 0.0}, 5.0, {}};

    // Point on circle
    CHECK_THAT(c.distance(Vec2{5.0, 0.0}), WithinAbs(0.0, 1e-12));

    // Point outside
    CHECK_THAT(c.distance(Vec2{8.0, 0.0}), WithinAbs(3.0, 1e-12));

    // Point inside
    CHECK_THAT(c.distance(Vec2{2.0, 0.0}), WithinAbs(3.0, 1e-12));
}

TEST_CASE("Circle 3D distance", "[circle]") {
    // Circle in XY plane, radius 1
    Circle3 c{Vec3{0.0, 0.0, 0.0}, 1.0, Vec3{0.0, 0.0, 1.0}};

    // Point on circle
    CHECK_THAT(c.distance(Vec3{1.0, 0.0, 0.0}), WithinAbs(0.0, 1e-12));

    // Point above center: distance = sqrt(1 + 0) = 1 (to nearest point on ring)
    // From (0,0,3) to any point on circle at distance sqrt(1+9)=sqrt(10)
    CHECK_THAT(c.distance(Vec3{0.0, 0.0, 3.0}), WithinAbs(std::sqrt(10.0), 1e-10));
}

TEST_CASE("Circle project 2D", "[circle]") {
    Circle2 c{Vec2{0.0, 0.0}, 5.0, {}};
    auto proj = c.project(Vec2{10.0, 0.0});
    CHECK_THAT(proj[0], WithinAbs(5.0, 1e-12));
    CHECK_THAT(proj[1], WithinAbs(0.0, 1e-12));
}

TEST_CASE("Circle bounding box 2D", "[circle]") {
    Circle2 c{Vec2{1.0, 2.0}, 3.0, {}};
    auto bb = c.bounding_box();
    CHECK_THAT(bb.min_corner[0], WithinAbs(-2.0, 1e-12));
    CHECK_THAT(bb.min_corner[1], WithinAbs(-1.0, 1e-12));
    CHECK_THAT(bb.max_corner[0], WithinAbs(4.0, 1e-12));
    CHECK_THAT(bb.max_corner[1], WithinAbs(5.0, 1e-12));
}

TEST_CASE("Disk area", "[circle]") {
    Disk3 d{Circle3{Vec3{0.0, 0.0, 0.0}, 3.0, Vec3{0.0, 0.0, 1.0}}};
    CHECK_THAT(d.area(), WithinAbs(std::numbers::pi * 9.0, 1e-10));
}

TEST_CASE("Disk contains 2D", "[circle]") {
    Disk2 d{Circle2{Vec2{0.0, 0.0}, 5.0, {}}};
    CHECK(d.contains(Vec2{0.0, 0.0}));
    CHECK(d.contains(Vec2{3.0, 4.0})); // on boundary
    CHECK_FALSE(d.contains(Vec2{4.0, 4.0}));
}

TEST_CASE("Disk contains 3D", "[circle]") {
    Disk3 d{Circle3{Vec3{0.0, 0.0, 0.0}, 5.0, Vec3{0.0, 0.0, 1.0}}};
    CHECK(d.contains(Vec3{3.0, 4.0, 0.0}));
    CHECK_FALSE(d.contains(Vec3{3.0, 4.0, 1.0})); // off plane
    CHECK_FALSE(d.contains(Vec3{6.0, 0.0, 0.0})); // outside radius
}

TEST_CASE("Disk project 2D", "[circle]") {
    Disk2 d{Circle2{Vec2{0.0, 0.0}, 5.0, {}}};

    // Inside: identity
    CHECK(d.project(Vec2{1.0, 1.0}) == Vec2{1.0, 1.0});

    // Outside: clamp to boundary
    auto proj = d.project(Vec2{10.0, 0.0});
    CHECK_THAT(proj[0], WithinAbs(5.0, 1e-12));
}

TEST_CASE("Circle/Disk concept satisfaction", "[circle]") {
    static_assert(Shape<Circle3>);
    static_assert(Measurable<Circle3>);
    static_assert(DistanceQueryable<Circle3>);
    static_assert(Bounded<Circle3>);

    static_assert(Shape<Disk3>);
    static_assert(ClosedShape<Disk3>);
    static_assert(Measurable<Disk3>);
    static_assert(DistanceQueryable<Disk3>);
    static_assert(Bounded<Disk3>);
    SUCCEED();
}
