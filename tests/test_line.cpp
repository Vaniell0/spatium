#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/geometry/line.hpp>

using namespace spatium;
using namespace spatium::geometry;
using Catch::Matchers::WithinAbs;

TEST_CASE("Line construction", "[line]") {
    auto line = Line3::from(Vec3{0.0, 0.0, 0.0}, Vec3{3.0, 0.0, 0.0});
    REQUIRE(line.has_value());
    CHECK_THAT(line->direction.norm(), WithinAbs(1.0, 1e-12));
    CHECK(line->direction[0] == 1.0);
}

TEST_CASE("Line from zero direction fails", "[line]") {
    auto line = Line3::from(Vec3{0.0, 0.0, 0.0}, Vec3{0.0, 0.0, 0.0});
    CHECK_FALSE(line.has_value());
    CHECK(line.error().code == ErrorCode::DegenerateInput);
}

TEST_CASE("Line at", "[line]") {
    auto line = *Line3::from(Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0});
    auto p = line.at(5.0);
    CHECK(p == Vec3{1.0, 5.0, 0.0});
}

TEST_CASE("Line project and distance", "[line]") {
    auto line = *Line3::from(Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0});

    // Point directly on line
    CHECK_THAT(line.distance(Vec3{5.0, 0.0, 0.0}), WithinAbs(0.0, 1e-12));

    // Point off the line
    auto proj = line.project(Vec3{3.0, 4.0, 0.0});
    CHECK(proj == Vec3{3.0, 0.0, 0.0});
    CHECK_THAT(line.distance(Vec3{3.0, 4.0, 0.0}), WithinAbs(4.0, 1e-12));
}

TEST_CASE("Ray clamps at origin", "[ray]") {
    auto ray = *Ray3::from(Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0});

    // Behind origin: projects to origin
    auto proj = ray.project(Vec3{-5.0, 3.0, 0.0});
    CHECK(proj == Vec3{0.0, 0.0, 0.0});

    // In front: projects normally
    auto proj2 = ray.project(Vec3{5.0, 3.0, 0.0});
    CHECK(proj2 == Vec3{5.0, 0.0, 0.0});
}

TEST_CASE("Segment project clamps to endpoints", "[segment]") {
    Segment3 seg{Vec3{0.0, 0.0, 0.0}, Vec3{10.0, 0.0, 0.0}};

    // Before start
    CHECK(seg.project(Vec3{-5.0, 3.0, 0.0}) == Vec3{0.0, 0.0, 0.0});

    // After end
    CHECK(seg.project(Vec3{15.0, 3.0, 0.0}) == Vec3{10.0, 0.0, 0.0});

    // Middle
    auto proj = seg.project(Vec3{5.0, 3.0, 0.0});
    CHECK(proj == Vec3{5.0, 0.0, 0.0});
}

TEST_CASE("Segment length and midpoint", "[segment]") {
    Segment3 seg{Vec3{0.0, 0.0, 0.0}, Vec3{3.0, 4.0, 0.0}};
    CHECK_THAT(seg.length(), WithinAbs(5.0, 1e-12));
    CHECK(seg.midpoint() == Vec3{1.5, 2.0, 0.0});
}

TEST_CASE("Segment at (lerp)", "[segment]") {
    Segment3 seg{Vec3{0.0, 0.0, 0.0}, Vec3{10.0, 0.0, 0.0}};
    CHECK(seg.at(0.0) == Vec3{0.0, 0.0, 0.0});
    CHECK(seg.at(1.0) == Vec3{10.0, 0.0, 0.0});
    CHECK(seg.at(0.5) == Vec3{5.0, 0.0, 0.0});
}

TEST_CASE("Segment bounding box", "[segment]") {
    Segment3 seg{Vec3{1.0, 5.0, 3.0}, Vec3{4.0, 2.0, 7.0}};
    auto bb = seg.bounding_box();
    CHECK(bb.min_corner == Vec3{1.0, 2.0, 3.0});
    CHECK(bb.max_corner == Vec3{4.0, 5.0, 7.0});
}

TEST_CASE("Line/Ray/Segment concept satisfaction", "[line]") {
    static_assert(Shape<Line3>);
    static_assert(DistanceQueryable<Line3>);
    static_assert(Shape<Ray3>);
    static_assert(DistanceQueryable<Ray3>);
    static_assert(Shape<Segment3>);
    static_assert(DistanceQueryable<Segment3>);
    SUCCEED();
}
