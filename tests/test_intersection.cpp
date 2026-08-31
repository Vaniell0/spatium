#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/geometry/intersection.hpp>
#include <spatium/geometry/distance.hpp>
#include <spatium/geometry/circle.hpp>

using namespace spatium;
using namespace spatium::geometry;
using Catch::Matchers::WithinAbs;

// ── Line-Hyperplane ────────────────────────────────────────────

TEST_CASE("Line-Plane intersection", "[intersection]") {
    auto line = *Line3::from(Vec3{0.0, 0.0, -5.0}, Vec3{0.0, 0.0, 1.0});
    auto plane = *Plane3::from_normal_and_point(Vec3{0.0, 0.0, 1.0}, Vec3{0.0, 0.0, 3.0});

    auto result = intersect(line, plane);
    REQUIRE(result.has_value());
    CHECK_THAT((*result)[0], WithinAbs(0.0, 1e-12));
    CHECK_THAT((*result)[1], WithinAbs(0.0, 1e-12));
    CHECK_THAT((*result)[2], WithinAbs(3.0, 1e-12));
}

TEST_CASE("Line-Plane parallel (no intersection)", "[intersection]") {
    auto line = *Line3::from(Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0});
    auto plane = *Plane3::from_normal_and_point(Vec3{0.0, 0.0, 1.0}, Vec3{0.0, 0.0, 3.0});

    auto result = intersect(line, plane);
    CHECK_FALSE(result.has_value());
    CHECK(result.error().code == ErrorCode::NoIntersection);
}

// ── Ray-Plane ──────────────────────────────────────────────────

TEST_CASE("Ray-Plane intersection", "[intersection]") {
    auto ray = *Ray3::from(Vec3{0.0, 0.0, 0.0}, Vec3{0.0, 0.0, 1.0});
    auto plane = *Plane3::from_normal_and_point(Vec3{0.0, 0.0, 1.0}, Vec3{0.0, 0.0, 5.0});

    auto result = intersect(ray, plane);
    REQUIRE(result.has_value());
    CHECK_THAT((*result)[2], WithinAbs(5.0, 1e-12));
}

TEST_CASE("Ray-Plane behind ray (no intersection)", "[intersection]") {
    auto ray = *Ray3::from(Vec3{0.0, 0.0, 10.0}, Vec3{0.0, 0.0, 1.0});
    auto plane = *Plane3::from_normal_and_point(Vec3{0.0, 0.0, 1.0}, Vec3{0.0, 0.0, 5.0});

    auto result = intersect(ray, plane);
    CHECK_FALSE(result.has_value());
}

// ── Segment-Plane ──────────────────────────────────────────────

TEST_CASE("Segment-Plane intersection", "[intersection]") {
    Segment3 seg{Vec3{0.0, 0.0, 0.0}, Vec3{0.0, 0.0, 10.0}};
    auto plane = *Plane3::from_normal_and_point(Vec3{0.0, 0.0, 1.0}, Vec3{0.0, 0.0, 5.0});

    auto result = intersect(seg, plane);
    REQUIRE(result.has_value());
    CHECK_THAT((*result)[2], WithinAbs(5.0, 1e-12));
}

TEST_CASE("Segment-Plane misses", "[intersection]") {
    Segment3 seg{Vec3{0.0, 0.0, 0.0}, Vec3{0.0, 0.0, 3.0}};
    auto plane = *Plane3::from_normal_and_point(Vec3{0.0, 0.0, 1.0}, Vec3{0.0, 0.0, 5.0});

    auto result = intersect(seg, plane);
    CHECK_FALSE(result.has_value());
}

// ── Line-Line 2D ───────────────────────────────────────────────

TEST_CASE("Line-Line 2D intersection", "[intersection]") {
    auto l1 = *Line<2>::from(Vec2{0.0, 0.0}, Vec2{1.0, 0.0});
    auto l2 = *Line<2>::from(Vec2{0.0, 0.0}, Vec2{0.0, 1.0});

    auto result = intersect(l1, l2);
    REQUIRE(result.has_value());
    CHECK_THAT((*result)[0], WithinAbs(0.0, 1e-12));
    CHECK_THAT((*result)[1], WithinAbs(0.0, 1e-12));
}

TEST_CASE("Line-Line 2D at offset", "[intersection]") {
    auto l1 = *Line<2>::from(Vec2{0.0, 1.0}, Vec2{1.0, 0.0});
    auto l2 = *Line<2>::from(Vec2{3.0, 0.0}, Vec2{0.0, 1.0});

    auto result = intersect(l1, l2);
    REQUIRE(result.has_value());
    CHECK_THAT((*result)[0], WithinAbs(3.0, 1e-12));
    CHECK_THAT((*result)[1], WithinAbs(1.0, 1e-12));
}

// ── Segment-Segment 2D ────────────────────────────────────────

TEST_CASE("Segment-Segment 2D intersection", "[intersection]") {
    Segment<2> s1{Vec2{0.0, 0.0}, Vec2{2.0, 2.0}};
    Segment<2> s2{Vec2{0.0, 2.0}, Vec2{2.0, 0.0}};

    auto result = intersect(s1, s2);
    REQUIRE(result.has_value());
    CHECK_THAT((*result)[0], WithinAbs(1.0, 1e-12));
    CHECK_THAT((*result)[1], WithinAbs(1.0, 1e-12));
}

TEST_CASE("Segment-Segment 2D no intersection", "[intersection]") {
    Segment<2> s1{Vec2{0.0, 0.0}, Vec2{1.0, 0.0}};
    Segment<2> s2{Vec2{0.0, 1.0}, Vec2{1.0, 1.0}};

    CHECK_FALSE(intersect(s1, s2).has_value());
}

// ── Ray-Triangle (Moller-Trumbore) ─────────────────────────────

TEST_CASE("Ray-Triangle hit", "[intersection]") {
    Triangle3 tri(Vec3{0.0, 0.0, 0.0}, Vec3{2.0, 0.0, 0.0}, Vec3{0.0, 2.0, 0.0});
    auto ray = *Ray3::from(Vec3{0.5, 0.5, 5.0}, Vec3{0.0, 0.0, -1.0});

    auto result = intersect(ray, tri);
    REQUIRE(result.has_value());
    CHECK_THAT((*result)[0], WithinAbs(0.5, 1e-12));
    CHECK_THAT((*result)[1], WithinAbs(0.5, 1e-12));
    CHECK_THAT((*result)[2], WithinAbs(0.0, 1e-12));
}

TEST_CASE("Ray-Triangle miss", "[intersection]") {
    Triangle3 tri(Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0});
    auto ray = *Ray3::from(Vec3{5.0, 5.0, 5.0}, Vec3{0.0, 0.0, -1.0});

    CHECK_FALSE(intersect(ray, tri).has_value());
}

TEST_CASE("Ray-Triangle behind", "[intersection]") {
    Triangle3 tri(Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0});
    auto ray = *Ray3::from(Vec3{0.25, 0.25, -5.0}, Vec3{0.0, 0.0, -1.0}); // pointing away

    CHECK_FALSE(intersect(ray, tri).has_value());
}

// ── Ray-Box ────────────────────────────────────────────────────

TEST_CASE("Ray-Box hit", "[intersection]") {
    Box3 box{Vec3{-1.0, -1.0, -1.0}, Vec3{1.0, 1.0, 1.0}};
    auto ray = *Ray3::from(Vec3{-5.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0});

    auto result = intersect(ray, box);
    REQUIRE(result.has_value());
    CHECK_THAT((*result)[0], WithinAbs(-1.0, 1e-12));
}

TEST_CASE("Ray-Box miss", "[intersection]") {
    Box3 box{Vec3{-1.0, -1.0, -1.0}, Vec3{1.0, 1.0, 1.0}};
    auto ray = *Ray3::from(Vec3{-5.0, 5.0, 0.0}, Vec3{1.0, 0.0, 0.0});

    CHECK_FALSE(intersect(ray, box).has_value());
}

// ── Plane-Plane ────────────────────────────────────────────────

TEST_CASE("Plane-Plane intersection", "[intersection]") {
    auto p1 = *Plane3::from_normal_and_point(Vec3{0.0, 0.0, 1.0}, Vec3{0.0, 0.0, 0.0});
    auto p2 = *Plane3::from_normal_and_point(Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 0.0, 0.0});

    auto result = intersect(p1, p2);
    REQUIRE(result.has_value());
    // Direction should be along X axis
    CHECK_THAT(std::abs(result->direction[0]), WithinAbs(1.0, 1e-12));
}

TEST_CASE("Plane-Plane parallel", "[intersection]") {
    auto p1 = *Plane3::from_normal_and_point(Vec3{0.0, 0.0, 1.0}, Vec3{0.0, 0.0, 0.0});
    auto p2 = *Plane3::from_normal_and_point(Vec3{0.0, 0.0, 1.0}, Vec3{0.0, 0.0, 5.0});

    CHECK_FALSE(intersect(p1, p2).has_value());
}

// ── Free-function distance ─────────────────────────────────────

TEST_CASE("Distance point-line", "[distance]") {
    auto line = *Line3::from(Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0});
    CHECK_THAT(distance(Vec3{0.0, 5.0, 0.0}, line), WithinAbs(5.0, 1e-12));
}

TEST_CASE("Distance line-line 3D", "[distance]") {
    auto l1 = *Line<3>::from(Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0});
    auto l2 = *Line<3>::from(Vec3{0.0, 0.0, 3.0}, Vec3{0.0, 1.0, 0.0});
    CHECK_THAT(distance(l1, l2), WithinAbs(3.0, 1e-12));
}

TEST_CASE("Distance segment-segment", "[distance]") {
    Segment3 s1{Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0}};
    Segment3 s2{Vec3{0.0, 0.0, 5.0}, Vec3{1.0, 0.0, 5.0}};
    CHECK_THAT(distance(s1, s2), WithinAbs(5.0, 1e-12));
}

// ── Generic BoundedRegion intersect (via subspace) ────────────

TEST_CASE("Triangle-Disk intersect via subspace", "[intersection][subspace]") {
    // Triangle in XY plane, disk in XZ plane — planes cross along X axis
    auto tri = Triangle3({0, 0, 0}, {4, 0, 0}, {0, 4, 0});
    Disk<3> disk{Circle<3>{Vec3{2, 0, 0}, 1.5, Vec3{0, 1, 0}}};

    auto r = intersect_via_subspace(tri, disk);
    REQUIRE(r.has_value());
    // Intersection line is along X axis (y=0, z=0).
    // Triangle clips to x ∈ [0, 4], disk clips to x ∈ [0.5, 3.5]
    CHECK_THAT(r->a[1], WithinAbs(0.0, 1e-10));
    CHECK_THAT(r->b[1], WithinAbs(0.0, 1e-10));
    CHECK_THAT(r->a[2], WithinAbs(0.0, 1e-10));
    CHECK_THAT(r->b[2], WithinAbs(0.0, 1e-10));
    CHECK(r->length() > 0.0);
    CHECK(r->length() < 3.01);  // bounded by disk diameter
}

TEST_CASE("Disk-Disk intersect via subspace", "[intersection][subspace]") {
    // Two disks on perpendicular planes, both centered at origin
    Disk<3> d1{Circle<3>{Vec3{0, 0, 0}, 1.0, Vec3{0, 0, 1}}};  // XY plane
    Disk<3> d2{Circle<3>{Vec3{0, 0, 0}, 1.0, Vec3{0, 1, 0}}};  // XZ plane

    auto r = intersect_via_subspace(d1, d2);
    REQUIRE(r.has_value());
    // Intersection is along X axis, both disks have radius 1
    CHECK_THAT(r->length(), WithinAbs(2.0, 1e-10));
}

TEST_CASE("Segment-Triangle intersect via subspace", "[intersection][subspace]") {
    // Segment crossing through a triangle
    auto tri = Triangle3({0, 0, 0}, {2, 0, 0}, {0, 2, 0});  // XY plane
    Segment<3> seg{Vec3{0.5, 0.5, -1}, Vec3{0.5, 0.5, 1}};  // vertical, through interior

    auto r = intersect_via_subspace(seg, tri);
    REQUIRE(r.has_value());
    CHECK_THAT((*r)[0], WithinAbs(0.5, 1e-10));
    CHECK_THAT((*r)[1], WithinAbs(0.5, 1e-10));
    CHECK_THAT((*r)[2], WithinAbs(0.0, 1e-10));
}

TEST_CASE("Segment-Triangle miss via subspace", "[intersection][subspace]") {
    auto tri = Triangle3({0, 0, 0}, {2, 0, 0}, {0, 2, 0});
    Segment<3> seg{Vec3{5, 5, -1}, Vec3{5, 5, 1}};

    auto r = intersect_via_subspace(seg, tri);
    CHECK_FALSE(r.has_value());
}

TEST_CASE("Parallel planes no intersection via subspace", "[intersection][subspace]") {
    auto tri = Triangle3({0, 0, 0}, {1, 0, 0}, {0, 1, 0});
    Disk<3> disk{Circle<3>{Vec3{0, 0, 5}, 1.0, Vec3{0, 0, 1}}};  // parallel to tri

    auto r = intersect_via_subspace(tri, disk);
    CHECK_FALSE(r.has_value());
}
