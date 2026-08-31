#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/geometry/transform.hpp>
#include <spatium/geometry/make.hpp>
#include <spatium/geometry/intersection.hpp>
#include <numbers>

using namespace spatium;
using namespace spatium::geometry;
using Catch::Matchers::WithinAbs;

TEST_CASE("Translation", "[transform]") {
    auto t = Transform3::translation(Vec3{10, 0, 0});
    auto result = t(Vec3{1, 2, 3});
    CHECK(result == Vec3{11, 2, 3});
}

TEST_CASE("Uniform scaling", "[transform]") {
    auto s = Transform3::uniform_scaling(2.0);
    auto result = s(Vec3{1, 2, 3});
    CHECK(result == Vec3{2, 4, 6});
}

TEST_CASE("Non-uniform scaling", "[transform]") {
    auto s = Transform3::scaling(Vec3{2, 3, 4});
    auto result = s(Vec3{1, 1, 1});
    CHECK(result == Vec3{2, 3, 4});
}

TEST_CASE("3D rotation around Z axis", "[transform]") {
    auto r = Transform3::rotation(Vec3{0, 0, 1}, std::numbers::pi / 2);
    auto result = r(Vec3{1, 0, 0});
    CHECK_THAT(result[0], WithinAbs(0.0, 1e-10));
    CHECK_THAT(result[1], WithinAbs(1.0, 1e-10));
    CHECK_THAT(result[2], WithinAbs(0.0, 1e-10));
}

TEST_CASE("2D rotation", "[transform]") {
    auto r = Transform2::rotation(std::numbers::pi / 2);
    auto result = r(Vec2{1, 0});
    CHECK_THAT(result[0], WithinAbs(0.0, 1e-10));
    CHECK_THAT(result[1], WithinAbs(1.0, 1e-10));
}

TEST_CASE("Transform composition", "[transform]") {
    auto s = Transform3::uniform_scaling(2.0);
    auto t = Transform3::translation(Vec3{10, 0, 0});
    // scale then translate
    auto composed = t * s;
    auto result = composed(Vec3{1, 0, 0});
    CHECK(result == Vec3{12, 0, 0});
}

TEST_CASE("Transform as Morphism pipe", "[transform]") {
    auto t = Transform3::translation(Vec3{5, 0, 0});
    Morphism<E3, E3> m = t;
    auto p = pt<E3>(Vec3{1, 0, 0});
    auto q = p | m;
    CHECK(q.raw() == Vec3{6, 0, 0});
}

TEST_CASE("Transform pipe operator", "[transform]") {
    auto t = Transform3::translation(Vec3{5, 0, 0});
    auto p = pt<E3>(Vec3{1, 0, 0});
    auto q = p | t;
    CHECK(q.raw() == Vec3{6, 0, 0});
}

// ── Result pipe-unwrap ───────────────────────────────────────

TEST_CASE("Result<Vec> pipe through transform: success", "[transform]") {
    auto t = translate(Vec3{10, 0, 0});
    Result<Vec3> ok = Vec3{1, 2, 3};

    auto r = ok | t;
    REQUIRE(r.has_value());
    CHECK(*r == Vec3{11, 2, 3});
}

TEST_CASE("Result<Vec> pipe through transform: error propagates", "[transform]") {
    auto t = translate(Vec3{10, 0, 0});
    Result<Vec3> err = std::unexpected(Error{ErrorCode::NoIntersection, "miss"});

    auto r = err | t;
    CHECK_FALSE(r.has_value());
}

TEST_CASE("Result<Point> pipe through transform", "[transform]") {
    auto t = translate(Vec3{5, 0, 0});
    Result<Point<E3>> ok = pt<E3>(Vec3{1, 0, 0});

    auto r = ok | t;
    REQUIRE(r.has_value());
    CHECK(r->raw() == Vec3{6, 0, 0});
}

// ── Lazy transform chains ────────────────────────────────────

TEST_CASE("Lazy transform: single leaf", "[transform]") {
    auto t = lazy(translate(Vec3{10, 0, 0}));
    auto result = t.apply(Vec3{1, 2, 3});
    CHECK(result == Vec3{11, 2, 3});
}

TEST_CASE("Lazy transform: chain of two", "[transform]") {
    auto chain = lazy(translate(Vec3{10, 0, 0})) * lazy(scale(2.0));
    // Apply: scale first (2,4,6), then translate (12,4,6)
    auto result = chain.apply(Vec3{1, 2, 3});
    CHECK(result == Vec3{12, 4, 6});
}

TEST_CASE("Lazy transform: chain of three", "[transform]") {
    auto chain = lazy(translate(Vec3{1, 0, 0})) * lazy(scale(2.0)) * lazy(translate(Vec3{10, 0, 0}));
    // translate(10,0,0) → (11,2,3), scale(2) → (22,4,6), translate(1,0,0) → (23,4,6)
    auto result = chain.apply(Vec3{1, 2, 3});
    CHECK(result == Vec3{23, 4, 6});
}

TEST_CASE("Lazy collapse equals eager composition", "[transform]") {
    auto eager = translate(Vec3{1, 0, 0}) * scale(2.0);
    auto lazy_chain = lazy(translate(Vec3{1, 0, 0})) * lazy(scale(2.0));
    auto collapsed = lazy_chain.collapse();

    auto p = Vec3{3, 4, 5};
    auto r1 = eager(p);
    auto r2 = collapsed(p);
    CHECK_THAT(r1[0], WithinAbs(r2[0], 1e-12));
    CHECK_THAT(r1[1], WithinAbs(r2[1], 1e-12));
    CHECK_THAT(r1[2], WithinAbs(r2[2], 1e-12));
}

TEST_CASE("Lazy transform: point pipe", "[transform]") {
    auto chain = lazy(translate(Vec3{10, 0, 0})) * lazy(scale(2.0));
    auto p = pt<E3>(Vec3{1, 2, 3});
    auto q = p | chain;
    CHECK(q.raw() == Vec3{12, 4, 6});
}

TEST_CASE("Lazy transform: Result<Vec> pipe", "[transform]") {
    auto chain = lazy(translate(Vec3{10, 0, 0})) * lazy(scale(2.0));
    Result<Vec3> ok = Vec3{1, 2, 3};
    auto r = ok | chain;
    REQUIRE(r.has_value());
    CHECK(*r == Vec3{12, 4, 6});
}

TEST_CASE("Intersection chain: line | plane | transform", "[transform]") {
    // line along X, plane at x=3 → intersection at (3,0,0), then translate +10 on X
    auto l = line(Vec3{0, 0, 0}, Vec3{1, 0, 0});
    auto pl = plane(Vec3{1, 0, 0}, Vec3{3, 0, 0});
    auto t = translate(Vec3{10, 0, 0});

    REQUIRE(l.has_value());
    REQUIRE(pl.has_value());

    auto result = l | *pl | t;
    REQUIRE(result.has_value());
    CHECK_THAT((*result)[0], WithinAbs(13.0, 1e-10));
    CHECK_THAT((*result)[1], WithinAbs(0.0, 1e-10));
    CHECK_THAT((*result)[2], WithinAbs(0.0, 1e-10));
}
