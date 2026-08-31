#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/geometry/triangle.hpp>

using namespace spatium;
using namespace spatium::geometry;
using Catch::Matchers::WithinAbs;

static const Triangle3 unit_tri(
    Vec3{0.0, 0.0, 0.0},
    Vec3{1.0, 0.0, 0.0},
    Vec3{0.0, 1.0, 0.0}
);

TEST_CASE("Triangle area 3D", "[triangle]") {
    CHECK_THAT(unit_tri.area(), WithinAbs(0.5, 1e-12));
}

TEST_CASE("Triangle area 2D", "[triangle]") {
    Triangle2 tri(Vec2{0.0, 0.0}, Vec2{4.0, 0.0}, Vec2{0.0, 3.0});
    CHECK_THAT(tri.area(), WithinAbs(6.0, 1e-12));
}

TEST_CASE("Triangle perimeter", "[triangle]") {
    // Right triangle 3-4-5
    Triangle3 tri(Vec3{0.0, 0.0, 0.0}, Vec3{3.0, 0.0, 0.0}, Vec3{0.0, 4.0, 0.0});
    CHECK_THAT(tri.perimeter(), WithinAbs(12.0, 1e-12));
}

TEST_CASE("Triangle normal", "[triangle]") {
    auto n = unit_tri.normal();
    CHECK_THAT(n[2], WithinAbs(1.0, 1e-12));
    CHECK_THAT(n.norm(), WithinAbs(1.0, 1e-12));
}

TEST_CASE("Triangle centroid", "[triangle]") {
    auto c = unit_tri.centroid();
    CHECK_THAT(c[0], WithinAbs(1.0 / 3.0, 1e-12));
    CHECK_THAT(c[1], WithinAbs(1.0 / 3.0, 1e-12));
    CHECK_THAT(c[2], WithinAbs(0.0, 1e-12));
}

TEST_CASE("Triangle barycentric at vertices", "[triangle]") {
    auto b0 = unit_tri.barycentric(unit_tri[0]);
    CHECK_THAT(b0[0], WithinAbs(1.0, 1e-10));
    CHECK_THAT(b0[1], WithinAbs(0.0, 1e-10));
    CHECK_THAT(b0[2], WithinAbs(0.0, 1e-10));

    auto b1 = unit_tri.barycentric(unit_tri[1]);
    CHECK_THAT(b1[1], WithinAbs(1.0, 1e-10));

    auto b2 = unit_tri.barycentric(unit_tri[2]);
    CHECK_THAT(b2[2], WithinAbs(1.0, 1e-10));
}

TEST_CASE("Triangle barycentric at centroid", "[triangle]") {
    auto bc = unit_tri.barycentric(unit_tri.centroid());
    CHECK_THAT(bc[0], WithinAbs(1.0 / 3.0, 1e-10));
    CHECK_THAT(bc[1], WithinAbs(1.0 / 3.0, 1e-10));
    CHECK_THAT(bc[2], WithinAbs(1.0 / 3.0, 1e-10));
}

TEST_CASE("Triangle contains", "[triangle]") {
    CHECK(unit_tri.contains(Vec3{0.1, 0.1, 0.0}));
    CHECK(unit_tri.contains(unit_tri.centroid()));
    CHECK(unit_tri.contains(unit_tri[0]));     // vertex
    CHECK(unit_tri.contains(Vec3{0.5, 0.0, 0.0})); // on edge
    CHECK_FALSE(unit_tri.contains(Vec3{1.0, 1.0, 0.0})); // outside
    CHECK_FALSE(unit_tri.contains(Vec3{-0.1, 0.0, 0.0}));
}

TEST_CASE("Triangle project inside", "[triangle]") {
    // Point on the plane, inside triangle
    auto proj = unit_tri.project(Vec3{0.2, 0.2, 0.0});
    CHECK_THAT(proj[0], WithinAbs(0.2, 1e-10));
    CHECK_THAT(proj[1], WithinAbs(0.2, 1e-10));
    CHECK_THAT(proj[2], WithinAbs(0.0, 1e-10));
}

TEST_CASE("Triangle project above", "[triangle]") {
    // Point above centroid
    auto c = unit_tri.centroid();
    auto proj = unit_tri.project(Vec3{c[0], c[1], 5.0});
    CHECK_THAT(proj[0], WithinAbs(c[0], 1e-10));
    CHECK_THAT(proj[1], WithinAbs(c[1], 1e-10));
    CHECK_THAT(proj[2], WithinAbs(0.0, 1e-10));
}

TEST_CASE("Triangle project outside (nearest edge)", "[triangle]") {
    auto proj = unit_tri.project(Vec3{-1.0, 0.5, 0.0});
    // Should project to the edge between v0(0,0,0) and v2(0,1,0) at y=0.5
    CHECK_THAT(proj[0], WithinAbs(0.0, 1e-10));
    CHECK_THAT(proj[1], WithinAbs(0.5, 1e-10));
}

TEST_CASE("Triangle distance", "[triangle]") {
    // Point directly above centroid at height 3
    auto c = unit_tri.centroid();
    CHECK_THAT(unit_tri.distance(Vec3{c[0], c[1], 3.0}), WithinAbs(3.0, 1e-10));
}

TEST_CASE("Triangle subdivide", "[triangle]") {
    auto subs = unit_tri.subdivide();
    CHECK(subs.size() == 4);

    // Total area should equal original
    double total = 0.0;
    for (const auto& t : subs)
        total += t.area();
    CHECK_THAT(total, WithinAbs(unit_tri.area(), 1e-12));
}

TEST_CASE("Triangle bounding box", "[triangle]") {
    auto bb = unit_tri.bounding_box();
    CHECK(bb.min_corner == Vec3{0.0, 0.0, 0.0});
    CHECK(bb.max_corner == Vec3{1.0, 1.0, 0.0});
}

TEST_CASE("Triangle edges", "[triangle]") {
    auto e0 = unit_tri.edge(0); // opposite v0, connects v1-v2
    CHECK(e0.a == unit_tri[1]);
    CHECK(e0.b == unit_tri[2]);
}

TEST_CASE("Triangle supporting plane", "[triangle]") {
    auto plane = unit_tri.supporting_plane();
    REQUIRE(plane.has_value());
    // All vertices should lie on the plane
    for (const auto& v : unit_tri.vertices)
        CHECK(plane->contains(v));
}

TEST_CASE("Triangle concept satisfaction", "[triangle]") {
    static_assert(Shape<Triangle3>);
    static_assert(ClosedShape<Triangle3>);
    static_assert(Measurable<Triangle3>);
    static_assert(Bounded<Triangle3>);
    static_assert(DistanceQueryable<Triangle3>);
    SUCCEED();
}
