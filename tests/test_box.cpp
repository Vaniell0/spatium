#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/geometry/box.hpp>

using namespace spatium;
using namespace spatium::geometry;
using Catch::Matchers::WithinAbs;

TEST_CASE("Box construction", "[box]") {
    Box3 box{Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 2.0, 3.0}};
    CHECK(box.min_corner == Vec3{0.0, 0.0, 0.0});
    CHECK(box.max_corner == Vec3{1.0, 2.0, 3.0});
}

TEST_CASE("Box from center and half-extents", "[box]") {
    auto box = Box3::from_center_half_extents(Vec3{5.0, 5.0, 5.0}, Vec3{1.0, 2.0, 3.0});
    CHECK(box.min_corner == Vec3{4.0, 3.0, 2.0});
    CHECK(box.max_corner == Vec3{6.0, 7.0, 8.0});
}

TEST_CASE("Box from points", "[box]") {
    std::array points = {
        Vec3{1.0, 5.0, 3.0},
        Vec3{-2.0, 0.0, 7.0},
        Vec3{4.0, 2.0, -1.0},
    };
    auto box = Box3::from_points(points);
    CHECK(box.min_corner == Vec3{-2.0, 0.0, -1.0});
    CHECK(box.max_corner == Vec3{4.0, 5.0, 7.0});
}

TEST_CASE("Box extents and centroid", "[box]") {
    Box3 box{Vec3{1.0, 2.0, 3.0}, Vec3{5.0, 6.0, 9.0}};
    CHECK(box.extents() == Vec3{4.0, 4.0, 6.0});
    CHECK(box.centroid() == Vec3{3.0, 4.0, 6.0});
}

TEST_CASE("Box measure", "[box]") {
    Box3 box{Vec3{0.0, 0.0, 0.0}, Vec3{2.0, 3.0, 4.0}};
    CHECK(box.measure() == 24.0); // volume

    Box2 box2{Vec2{0.0, 0.0}, Vec2{3.0, 4.0}};
    CHECK(box2.measure() == 12.0); // area
}

TEST_CASE("Box surface measure", "[box]") {
    Box2 box2{Vec2{0.0, 0.0}, Vec2{3.0, 4.0}};
    CHECK(box2.surface_measure() == 14.0); // perimeter

    Box3 box3{Vec3{0.0, 0.0, 0.0}, Vec3{2.0, 3.0, 4.0}};
    CHECK(box3.surface_measure() == 52.0); // 2*(6+12+8)
}

TEST_CASE("Box contains point", "[box]") {
    Box3 box{Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 1.0, 1.0}};
    CHECK(box.contains(Vec3{0.5, 0.5, 0.5}));
    CHECK(box.contains(Vec3{0.0, 0.0, 0.0})); // on boundary
    CHECK(box.contains(Vec3{1.0, 1.0, 1.0})); // on boundary
    CHECK_FALSE(box.contains(Vec3{1.1, 0.5, 0.5}));
    CHECK_FALSE(box.contains(Vec3{-0.1, 0.5, 0.5}));
}

TEST_CASE("Box contains box", "[box]") {
    Box3 outer{Vec3{0.0, 0.0, 0.0}, Vec3{10.0, 10.0, 10.0}};
    Box3 inner{Vec3{2.0, 2.0, 2.0}, Vec3{8.0, 8.0, 8.0}};
    Box3 partial{Vec3{5.0, 5.0, 5.0}, Vec3{15.0, 15.0, 15.0}};
    CHECK(outer.contains(inner));
    CHECK_FALSE(outer.contains(partial));
    CHECK_FALSE(inner.contains(outer));
}

TEST_CASE("Box intersection test", "[box]") {
    Box3 a{Vec3{0.0, 0.0, 0.0}, Vec3{2.0, 2.0, 2.0}};
    Box3 b{Vec3{1.0, 1.0, 1.0}, Vec3{3.0, 3.0, 3.0}};
    Box3 c{Vec3{5.0, 5.0, 5.0}, Vec3{6.0, 6.0, 6.0}};
    CHECK(a.intersects(b));
    CHECK_FALSE(a.intersects(c));
}

TEST_CASE("Box intersection result", "[box]") {
    Box3 a{Vec3{0.0, 0.0, 0.0}, Vec3{2.0, 2.0, 2.0}};
    Box3 b{Vec3{1.0, 1.0, 1.0}, Vec3{3.0, 3.0, 3.0}};

    auto result = a.intersection(b);
    REQUIRE(result.has_value());
    CHECK(result->min_corner == Vec3{1.0, 1.0, 1.0});
    CHECK(result->max_corner == Vec3{2.0, 2.0, 2.0});

    Box3 c{Vec3{5.0, 5.0, 5.0}, Vec3{6.0, 6.0, 6.0}};
    CHECK_FALSE(a.intersection(c).has_value());
}

TEST_CASE("Box union", "[box]") {
    Box3 a{Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 1.0, 1.0}};
    Box3 b{Vec3{2.0, 2.0, 2.0}, Vec3{3.0, 3.0, 3.0}};
    auto u = a.union_with(b);
    CHECK(u.min_corner == Vec3{0.0, 0.0, 0.0});
    CHECK(u.max_corner == Vec3{3.0, 3.0, 3.0});
}

TEST_CASE("Box distance", "[box]") {
    Box3 box{Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 1.0, 1.0}};

    // Inside: distance = 0
    CHECK_THAT(box.distance(Vec3{0.5, 0.5, 0.5}), WithinAbs(0.0, 1e-12));

    // Along axis: distance = gap
    CHECK_THAT(box.distance(Vec3{2.0, 0.5, 0.5}), WithinAbs(1.0, 1e-12));

    // Corner case: sqrt(1+1+1) from (2,2,2) to corner (1,1,1)
    CHECK_THAT(box.distance(Vec3{2.0, 2.0, 2.0}), WithinAbs(std::sqrt(3.0), 1e-12));
}

TEST_CASE("Box project", "[box]") {
    Box3 box{Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 1.0, 1.0}};

    // Inside: project is identity
    CHECK(box.project(Vec3{0.5, 0.5, 0.5}) == Vec3{0.5, 0.5, 0.5});

    // Outside: clamped
    CHECK(box.project(Vec3{2.0, -1.0, 0.5}) == Vec3{1.0, 0.0, 0.5});
}

TEST_CASE("Box concept satisfaction", "[box]") {
    static_assert(Shape<Box3>);
    static_assert(ClosedShape<Box3>);
    static_assert(Measurable<Box3>);
    static_assert(Bounded<Box3>);
    static_assert(DistanceQueryable<Box3>);
    SUCCEED();
}

TEST_CASE("Box 2D works", "[box]") {
    Box2 box{Vec2{-1.0, -1.0}, Vec2{1.0, 1.0}};
    CHECK(box.measure() == 4.0);
    CHECK(box.contains(Vec2{0.0, 0.0}));
    CHECK_FALSE(box.contains(Vec2{2.0, 0.0}));
}
