#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/geometry/clip.hpp>

using namespace spatium;
using namespace spatium::geometry;
using Catch::Matchers::WithinAbs;

// ── Point clip ────────────────────────────────────────────────

TEST_CASE("clip point inside triangle", "[clip]") {
    auto tri = Triangle3({0, 0, 0}, {1, 0, 0}, {0, 1, 0});
    auto r = clip(Vec3{0.2, 0.2, 0.0}, tri);
    REQUIRE(r.has_value());
    CHECK_THAT((*r)[0], WithinAbs(0.2, 1e-10));
}

TEST_CASE("clip point outside triangle", "[clip]") {
    auto tri = Triangle3({0, 0, 0}, {1, 0, 0}, {0, 1, 0});
    auto r = clip(Vec3{2.0, 2.0, 0.0}, tri);
    CHECK_FALSE(r.has_value());
}

TEST_CASE("clip point inside segment", "[clip]") {
    Segment3 seg{{0, 0, 0}, {1, 0, 0}};
    auto r = clip(Vec3{0.5, 0.0, 0.0}, seg);
    REQUIRE(r.has_value());
}

TEST_CASE("clip point outside segment", "[clip]") {
    Segment3 seg{{0, 0, 0}, {1, 0, 0}};
    auto r = clip(Vec3{2.0, 0.0, 0.0}, seg);
    CHECK_FALSE(r.has_value());
}

TEST_CASE("clip point inside disk", "[clip]") {
    Disk<3> disk{Circle<3>{Vec3{0, 0, 0}, 1.0, Vec3{0, 0, 1}}};
    auto r = clip(Vec3{0.5, 0.0, 0.0}, disk);
    REQUIRE(r.has_value());
}

TEST_CASE("clip point outside disk", "[clip]") {
    Disk<3> disk{Circle<3>{Vec3{0, 0, 0}, 1.0, Vec3{0, 0, 1}}};
    auto r = clip(Vec3{2.0, 0.0, 0.0}, disk);
    CHECK_FALSE(r.has_value());
}

TEST_CASE("clip point inside box", "[clip]") {
    Box3 box{Vec3{0, 0, 0}, Vec3{1, 1, 1}};
    auto r = clip(Vec3{0.5, 0.5, 0.5}, box);
    REQUIRE(r.has_value());
}

TEST_CASE("clip point outside box", "[clip]") {
    Box3 box{Vec3{0, 0, 0}, Vec3{1, 1, 1}};
    auto r = clip(Vec3{2.0, 0.5, 0.5}, box);
    CHECK_FALSE(r.has_value());
}

// ── Line clip ─────────────────────────────────────────────────

TEST_CASE("clip line to disk — chord", "[clip]") {
    Disk<3> disk{Circle<3>{Vec3{0, 0, 0}, 1.0, Vec3{0, 0, 1}}};
    auto line = *Line<3>::from(Vec3{-5, 0, 0}, Vec3{1, 0, 0});
    auto r = clip(line, disk);
    REQUIRE(r.has_value());
    CHECK_THAT(r->length(), WithinAbs(2.0, 1e-10));
    CHECK_THAT(r->a[0], WithinAbs(-1.0, 1e-10));
    CHECK_THAT(r->b[0], WithinAbs(1.0, 1e-10));
}

TEST_CASE("clip line to disk — miss", "[clip]") {
    Disk<3> disk{Circle<3>{Vec3{0, 0, 0}, 1.0, Vec3{0, 0, 1}}};
    auto line = *Line<3>::from(Vec3{-5, 2, 0}, Vec3{1, 0, 0});
    auto r = clip(line, disk);
    CHECK_FALSE(r.has_value());
}

TEST_CASE("clip line to disk — tangent", "[clip]") {
    Disk<3> disk{Circle<3>{Vec3{0, 0, 0}, 1.0, Vec3{0, 0, 1}}};
    auto line = *Line<3>::from(Vec3{-5, 1, 0}, Vec3{1, 0, 0});
    auto r = clip(line, disk);
    REQUIRE(r.has_value());
    CHECK_THAT(r->length(), WithinAbs(0.0, 1e-5));
}

TEST_CASE("clip line to triangle — through center", "[clip]") {
    auto tri = Triangle3({0, 0, 0}, {2, 0, 0}, {0, 2, 0});
    auto line = *Line<3>::from(Vec3{-1, 0.5, 0}, Vec3{1, 0, 0});
    auto r = clip(line, tri);
    REQUIRE(r.has_value());
    CHECK_THAT(r->a[0], WithinAbs(0.0, 1e-10));
    CHECK_THAT(r->b[0], WithinAbs(1.5, 1e-10));
}

TEST_CASE("clip line to segment — collinear", "[clip]") {
    Segment3 seg{{0, 0, 0}, {1, 0, 0}};
    auto line = *Line<3>::from(Vec3{-5, 0, 0}, Vec3{1, 0, 0});
    auto r = clip(line, seg);
    REQUIRE(r.has_value());
    CHECK_THAT(r->a[0], WithinAbs(0.0, 1e-10));
    CHECK_THAT(r->b[0], WithinAbs(1.0, 1e-10));
}

TEST_CASE("clip line to segment — not collinear", "[clip]") {
    Segment3 seg{{0, 0, 0}, {1, 0, 0}};
    auto line = *Line<3>::from(Vec3{0, 1, 0}, Vec3{1, 0, 0});
    auto r = clip(line, seg);
    CHECK_FALSE(r.has_value());
}

// ── Segment reclip ────────────────────────────────────────────

TEST_CASE("clip segment to disk — partial", "[clip]") {
    Disk<3> disk{Circle<3>{Vec3{0, 0, 0}, 1.0, Vec3{0, 0, 1}}};
    Segment<3> seg{Vec3{-2, 0, 0}, Vec3{0.5, 0, 0}};
    auto r = clip(seg, disk);
    REQUIRE(r.has_value());
    CHECK_THAT(r->a[0], WithinAbs(-1.0, 1e-10));
    CHECK_THAT(r->b[0], WithinAbs(0.5, 1e-10));
}

TEST_CASE("clip segment to disk — fully outside", "[clip]") {
    Disk<3> disk{Circle<3>{Vec3{0, 0, 0}, 1.0, Vec3{0, 0, 1}}};
    Segment<3> seg{Vec3{2, 0, 0}, Vec3{3, 0, 0}};
    auto r = clip(seg, disk);
    CHECK_FALSE(r.has_value());
}
