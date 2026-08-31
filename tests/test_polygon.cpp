#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/geometry/polygon.hpp>

using namespace spatium;
using namespace spatium::geometry;
using Catch::Matchers::WithinAbs;

static Polygon2 unit_square() {
    return {{{Vec2{0.0, 0.0}, Vec2{1.0, 0.0}, Vec2{1.0, 1.0}, Vec2{0.0, 1.0}}}};
}

TEST_CASE("Polygon area 2D (square)", "[polygon]") {
    CHECK_THAT(unit_square().area(), WithinAbs(1.0, 1e-12));
}

TEST_CASE("Polygon area 2D (irregular)", "[polygon]") {
    Polygon2 poly{{{Vec2{0.0, 0.0}, Vec2{4.0, 0.0}, Vec2{4.0, 3.0}}}};
    CHECK_THAT(poly.area(), WithinAbs(6.0, 1e-12)); // right triangle
}

TEST_CASE("Polygon perimeter", "[polygon]") {
    CHECK_THAT(unit_square().perimeter(), WithinAbs(4.0, 1e-12));
}

TEST_CASE("Polygon centroid", "[polygon]") {
    auto c = unit_square().centroid();
    CHECK_THAT(c[0], WithinAbs(0.5, 1e-12));
    CHECK_THAT(c[1], WithinAbs(0.5, 1e-12));
}

TEST_CASE("Polygon contains 2D", "[polygon]") {
    auto sq = unit_square();
    CHECK(sq.contains(Vec2{0.5, 0.5}));
    CHECK(sq.contains(Vec2{0.01, 0.01}));
    CHECK_FALSE(sq.contains(Vec2{1.5, 0.5}));
    CHECK_FALSE(sq.contains(Vec2{-0.1, 0.5}));
}

TEST_CASE("Polygon triangulate 2D", "[polygon]") {
    auto tris = unit_square().triangulate();
    CHECK(tris.size() == 2);

    // Total area should equal polygon area
    double total = 0.0;
    for (const auto& t : tris)
        total += t.area();
    CHECK_THAT(total, WithinAbs(1.0, 1e-10));
}

TEST_CASE("Polygon triangulate pentagon", "[polygon]") {
    // Regular-ish pentagon (convex)
    Polygon2 pent{{{
        Vec2{1.0, 0.0},
        Vec2{0.31, 0.95},
        Vec2{-0.81, 0.59},
        Vec2{-0.81, -0.59},
        Vec2{0.31, -0.95},
    }}};
    auto tris = pent.triangulate();
    CHECK(tris.size() == 3);
}

TEST_CASE("Polygon bounding box", "[polygon]") {
    auto bb = unit_square().bounding_box();
    CHECK(bb.min_corner == Vec2{0.0, 0.0});
    CHECK(bb.max_corner == Vec2{1.0, 1.0});
}

TEST_CASE("Polygon distance", "[polygon]") {
    auto sq = unit_square();
    // Point outside, nearest to bottom edge
    CHECK_THAT(sq.distance(Vec2{0.5, -1.0}), WithinAbs(1.0, 1e-12));
}

TEST_CASE("Polygon 3D area", "[polygon]") {
    Polygon3 poly{{{
        Vec3{0.0, 0.0, 0.0},
        Vec3{1.0, 0.0, 0.0},
        Vec3{1.0, 1.0, 0.0},
        Vec3{0.0, 1.0, 0.0},
    }}};
    CHECK_THAT(poly.area(), WithinAbs(1.0, 1e-10));
}

TEST_CASE("Polygon 3D normal", "[polygon]") {
    Polygon3 poly{{{
        Vec3{0.0, 0.0, 0.0},
        Vec3{1.0, 0.0, 0.0},
        Vec3{1.0, 1.0, 0.0},
        Vec3{0.0, 1.0, 0.0},
    }}};
    auto n = poly.normal();
    CHECK_THAT(std::abs(n[2]), WithinAbs(1.0, 1e-12));
}

TEST_CASE("Polygon concept satisfaction", "[polygon]") {
    static_assert(Shape<Polygon2>);
    static_assert(ClosedShape<Polygon2>);
    static_assert(Measurable<Polygon2>);
    static_assert(Bounded<Polygon2>);
    SUCCEED();
}
