#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/geometry/convex_hull.hpp>
#include <spatium/geometry/intersection.hpp>
#include <spatium/geometry/triangle.hpp>

using namespace spatium;
using namespace spatium::geometry;
using Catch::Matchers::WithinAbs;

TEST_CASE("Convex hull of square", "[convex_hull]") {
    std::vector<Vec2> points = {
        {0, 0}, {1, 0}, {1, 1}, {0, 1}, {0.5, 0.5}, // interior point
    };
    auto hull = convex_hull(points);
    REQUIRE(hull);
    CHECK(hull->size() == 4); // interior point excluded
    CHECK_THAT(hull->area(), WithinAbs(1.0, 1e-10));
}

TEST_CASE("Convex hull of collinear points", "[convex_hull]") {
    std::vector<Vec2> points = {{0, 0}, {1, 0}, {2, 0}};
    auto hull = convex_hull(points);
    REQUIRE(hull);
    // Three collinear points produce a degenerate hull (size ≤ 3) but
    // the call still succeeds — the input was non-empty.
    CHECK(hull->size() <= 3);
}

TEST_CASE("Convex hull of triangle", "[convex_hull]") {
    std::vector<Vec2> points = {{0, 0}, {4, 0}, {2, 3}};
    auto hull = convex_hull(points);
    REQUIRE(hull);
    CHECK(hull->size() == 3);
    CHECK_THAT(hull->area(), WithinAbs(6.0, 1e-10));
}

TEST_CASE("Convex hull of random cloud", "[convex_hull]") {
    std::vector<Vec2> points = {
        {0, 0}, {10, 0}, {10, 10}, {0, 10}, // outer
        {3, 3}, {7, 7}, {5, 2}, {2, 8},     // inner
    };
    auto hull = convex_hull(points);
    REQUIRE(hull);
    CHECK(hull->size() == 4);
    CHECK_THAT(hull->area(), WithinAbs(100.0, 1e-10));
}

TEST_CASE("Convex hull rejects fewer than 3 points", "[convex_hull][degenerate]") {
    CHECK_FALSE(convex_hull(std::vector<Vec2>{}).has_value());
    CHECK_FALSE(convex_hull(std::vector<Vec2>{{0, 0}}).has_value());
    CHECK_FALSE(convex_hull(std::vector<Vec2>{{0, 0}, {1, 1}}).has_value());

    auto err = convex_hull(std::vector<Vec2>{{0, 0}});
    REQUIRE_FALSE(err);
    CHECK(err.error().code == ErrorCode::DegenerateInput);
}

TEST_CASE("Triangle-Triangle intersection", "[intersection]") {
    Triangle3 t1(Vec3{0, 0, 0}, Vec3{2, 0, 0}, Vec3{1, 2, 0});
    Triangle3 t2(Vec3{1, 0, -1}, Vec3{1, 0, 1}, Vec3{1, 2, 0});

    auto result = intersect(t1, t2);
    REQUIRE(result.has_value());
    // Intersection should be a segment on the plane x=1
    CHECK_THAT(result->a[0], WithinAbs(1.0, 1e-8));
    CHECK_THAT(result->b[0], WithinAbs(1.0, 1e-8));
}

TEST_CASE("Triangle-Triangle no intersection", "[intersection]") {
    Triangle3 t1(Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0, 1, 0});
    Triangle3 t2(Vec3{0, 0, 5}, Vec3{1, 0, 5}, Vec3{0, 1, 5});
    CHECK_FALSE(intersect(t1, t2).has_value());
}

TEST_CASE("Segment is measurable", "[segment]") {
    Segment3 s{Vec3{0, 0, 0}, Vec3{3, 4, 0}};
    CHECK_THAT(s.measure(), WithinAbs(5.0, 1e-12));
    static_assert(Measurable<Segment<3>>);
}
