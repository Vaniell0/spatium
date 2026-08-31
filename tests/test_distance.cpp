#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/geometry/distance.hpp>
#include <spatium/geometry/make.hpp>

using namespace spatium;
using namespace spatium::geometry;
using Catch::Matchers::WithinAbs;

// ── Point-to-shape wrappers ──────────────────────────────────

TEST_CASE("distance(point, circle) 2D", "[distance]") {
    auto c = circle(Vec2{0, 0}, 5.0);
    CHECK_THAT(distance(Vec2{8, 0}, c), WithinAbs(3.0, 1e-12));
    CHECK_THAT(distance(Vec2{5, 0}, c), WithinAbs(0.0, 1e-12));
    CHECK_THAT(distance(Vec2{3, 0}, c), WithinAbs(2.0, 1e-12));
}

TEST_CASE("distance(point, disk) 2D", "[distance]") {
    auto d = disk(Vec2{0, 0}, 5.0);
    CHECK_THAT(distance(Vec2{8, 0}, d), WithinAbs(3.0, 1e-12));
    CHECK_THAT(distance(Vec2{3, 0}, d), WithinAbs(0.0, 1e-12));
}

TEST_CASE("distance(point, polygon) 2D", "[distance]") {
    auto p = poly<2, double>({{0, 0}, {4, 0}, {4, 3}, {0, 3}});
    CHECK_THAT(distance(Vec2{2, 1}, p), WithinAbs(0.0, 1e-12));
    CHECK_THAT(distance(Vec2{6, 0}, p), WithinAbs(2.0, 1e-12));
}

// ── Box-Box ─────────────────────────────────────────────────

TEST_CASE("distance(box, box) separated", "[distance]") {
    auto a = box(Vec2{0, 0}, Vec2{1, 1});
    auto b = box(Vec2{3, 0}, Vec2{4, 1});
    CHECK_THAT(distance(a, b), WithinAbs(2.0, 1e-12));
}

TEST_CASE("distance(box, box) touching", "[distance]") {
    auto a = box(Vec2{0, 0}, Vec2{1, 1});
    auto b = box(Vec2{1, 0}, Vec2{2, 1});
    CHECK_THAT(distance(a, b), WithinAbs(0.0, 1e-12));
}

TEST_CASE("distance(box, box) overlapping", "[distance]") {
    auto a = box(Vec2{0, 0}, Vec2{2, 2});
    auto b = box(Vec2{1, 1}, Vec2{3, 3});
    CHECK_THAT(distance(a, b), WithinAbs(0.0, 1e-12));
}

TEST_CASE("distance(box, box) diagonal gap", "[distance]") {
    auto a = box(Vec2{0, 0}, Vec2{1, 1});
    auto b = box(Vec2{4, 4}, Vec2{5, 5});
    CHECK_THAT(distance(a, b), WithinAbs(std::sqrt(18.0), 1e-12));
}

// ── Circle-Circle 2D ────────────────────────────────────────

TEST_CASE("distance(circle, circle) separated", "[distance]") {
    auto a = circle(Vec2{0, 0}, 1.0);
    auto b = circle(Vec2{5, 0}, 1.0);
    CHECK_THAT(distance(a, b), WithinAbs(3.0, 1e-12));
}

TEST_CASE("distance(circle, circle) touching", "[distance]") {
    auto a = circle(Vec2{0, 0}, 2.0);
    auto b = circle(Vec2{4, 0}, 2.0);
    CHECK_THAT(distance(a, b), WithinAbs(0.0, 1e-12));
}

TEST_CASE("distance(circle, circle) overlapping", "[distance]") {
    auto a = circle(Vec2{0, 0}, 3.0);
    auto b = circle(Vec2{1, 0}, 3.0);
    CHECK_THAT(distance(a, b), WithinAbs(0.0, 1e-12));
}

// ── Triangle-Triangle ───────────────────────────────────────

TEST_CASE("distance(triangle, triangle) separated", "[distance]") {
    auto a = tri(Vec2{0, 0}, Vec2{1, 0}, Vec2{0, 1});
    auto b = tri(Vec2{3, 0}, Vec2{4, 0}, Vec2{3, 1});
    CHECK_THAT(distance(a, b), WithinAbs(2.0, 1e-12));
}

TEST_CASE("distance(triangle, triangle) overlapping", "[distance]") {
    auto a = tri(Vec2{0, 0}, Vec2{2, 0}, Vec2{1, 2});
    auto b = tri(Vec2{1, 0}, Vec2{3, 0}, Vec2{2, 2});
    CHECK_THAT(distance(a, b), WithinAbs(0.0, 1e-12));
}

// ── Polygon-Polygon ─────────────────────────────────────────

TEST_CASE("distance(polygon, polygon) separated", "[distance]") {
    auto a = poly<2, double>({{0, 0}, {1, 0}, {1, 1}, {0, 1}});
    auto b = poly<2, double>({{3, 0}, {4, 0}, {4, 1}, {3, 1}});
    CHECK_THAT(distance(a, b), WithinAbs(2.0, 1e-12));
}

TEST_CASE("distance(polygon, polygon) contained", "[distance]") {
    auto outer = poly<2, double>({{0, 0}, {10, 0}, {10, 10}, {0, 10}});
    auto inner = poly<2, double>({{2, 2}, {4, 2}, {4, 4}, {2, 4}});
    CHECK_THAT(distance(outer, inner), WithinAbs(0.0, 1e-12));
}

// ── Segment-Triangle ────────────────────────────────────────

TEST_CASE("distance(segment, triangle) intersecting", "[distance]") {
    auto s = seg(Vec2{0.5, 0.25}, Vec2{0.5, 0.5});
    auto t = tri(Vec2{0, 0}, Vec2{1, 0}, Vec2{0, 1});
    CHECK_THAT(distance(s, t), WithinAbs(0.0, 1e-12));
}

TEST_CASE("distance(segment, triangle) separated", "[distance]") {
    auto s = seg(Vec2{3, 0}, Vec2{4, 0});
    auto t = tri(Vec2{0, 0}, Vec2{1, 0}, Vec2{0, 1});
    CHECK_THAT(distance(s, t), WithinAbs(2.0, 1e-12));
}
