#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/geometry/boolean.hpp>
#include <spatium/geometry/make.hpp>

using namespace spatium;
using namespace spatium::geometry;
using Catch::Matchers::WithinAbs;

// ── 2D Triangle-Triangle ──────────────────────────────────────

TEST_CASE("intersection_region 2D: overlapping triangles", "[boolean]") {
    Triangle<2> t1({0, 0}, {2, 0}, {1, 2});
    Triangle<2> t2({1, 0}, {3, 0}, {2, 2});

    auto r = intersection_region(t1, t2);
    REQUIRE(r.has_value());
    CHECK(r->vertices.size() >= 3);
    CHECK(r->area() > 0.0);
    // Intersection area must be less than both
    CHECK(r->area() < t1.area());
    CHECK(r->area() < t2.area());
}

TEST_CASE("intersection_region 2D: non-overlapping triangles", "[boolean]") {
    Triangle<2> t1({0, 0}, {1, 0}, {0, 1});
    Triangle<2> t2({5, 5}, {6, 5}, {5, 6});

    auto r = intersection_region(t1, t2);
    CHECK_FALSE(r.has_value());
}

TEST_CASE("intersection_region 2D: one inside other", "[boolean]") {
    Triangle<2> big({0, 0}, {10, 0}, {5, 10});
    Triangle<2> small({3, 1}, {5, 1}, {4, 3});

    auto r = intersection_region(big, small);
    REQUIRE(r.has_value());
    CHECK_THAT(r->area(), WithinAbs(small.area(), 1e-8));
}

// ── 3D Triangle-Triangle (coplanar) ───────────────────────────

TEST_CASE("intersection_region 3D: coplanar overlapping triangles", "[boolean]") {
    // Both in XY plane
    Triangle<3> t1({0, 0, 0}, {4, 0, 0}, {0, 4, 0});
    Triangle<3> t2({1, 1, 0}, {5, 1, 0}, {1, 5, 0});

    auto r = intersection_region(t1, t2);
    REQUIRE(r.has_value());
    CHECK(r->vertices.size() >= 3);
    CHECK(r->area() > 0.0);
}

// ── 3D Triangle-Disk (coplanar) ───────────────────────────────

TEST_CASE("intersection_region 3D: triangle and coplanar disk", "[boolean]") {
    Triangle<3> tri({0, 0, 0}, {4, 0, 0}, {0, 4, 0});
    Disk<3> disk{Circle<3>{Vec3{1, 1, 0}, 0.5, Vec3{0, 0, 1}}};

    auto r = intersection_region(tri, disk);
    REQUIRE(r.has_value());
    // Disk is fully inside triangle
    auto disk_area = std::numbers::pi * 0.5 * 0.5;
    CHECK_THAT(r->area(), WithinAbs(disk_area, disk_area * 0.05));  // 5% tolerance (polygon approx)
}

TEST_CASE("intersection_region 3D: disk partially outside triangle", "[boolean]") {
    Triangle<3> tri({0, 0, 0}, {2, 0, 0}, {0, 2, 0});
    Disk<3> disk{Circle<3>{Vec3{0, 0, 0}, 1.0, Vec3{0, 0, 1}}};

    auto r = intersection_region(tri, disk);
    REQUIRE(r.has_value());
    auto disk_area = std::numbers::pi;
    CHECK(r->area() < disk_area);
    CHECK(r->area() > 0.0);
}

// ── Symmetric difference area ─────────────────────────────────

TEST_CASE("symmetric_difference_area 2D: identity", "[boolean]") {
    Triangle<2> t1({0, 0}, {2, 0}, {1, 2});
    Triangle<2> t2({1, 0}, {3, 0}, {2, 2});

    auto inter = intersection_region(t1, t2);
    REQUIRE(inter.has_value());

    auto sd = symmetric_difference_area<2, double>(t1, t2);
    REQUIRE(sd.has_value());
    // area(A△B) = area(A) + area(B) - 2*area(A∩B)
    auto expected = t1.area() + t2.area() - 2.0 * inter->area();
    CHECK_THAT(*sd, WithinAbs(expected, 1e-10));
}

// ── Difference area ───────────────────────────────────────────

TEST_CASE("difference_area 2D: non-overlapping", "[boolean]") {
    Triangle<2> t1({0, 0}, {1, 0}, {0, 1});
    Triangle<2> t2({5, 5}, {6, 5}, {5, 6});

    auto d = difference_area<2, double>(t1, t2);
    REQUIRE(d.has_value());
    CHECK_THAT(*d, WithinAbs(t1.area(), 1e-10));
}

TEST_CASE("difference_area 2D: fully contained", "[boolean]") {
    Triangle<2> big({0, 0}, {10, 0}, {5, 10});
    Triangle<2> small({3, 1}, {5, 1}, {4, 3});

    auto d = difference_area<2, double>(big, small);
    REQUIRE(d.has_value());
    CHECK_THAT(*d, WithinAbs(big.area() - small.area(), 1e-8));
}

// ── Polygon boolean operators ─────────────────────────────────

TEST_CASE("operator& on polygons: overlapping squares", "[boolean]") {
    auto a = poly<2, double>({{0, 0}, {2, 0}, {2, 2}, {0, 2}});
    auto b = poly<2, double>({{1, 1}, {3, 1}, {3, 3}, {1, 3}});

    auto r = a & b;
    REQUIRE(r.has_value());
    CHECK_THAT(r->area(), WithinAbs(1.0, 1e-10));
}

TEST_CASE("operator& on polygons: non-overlapping", "[boolean]") {
    auto a = poly<2, double>({{0, 0}, {1, 0}, {1, 1}, {0, 1}});
    auto b = poly<2, double>({{5, 5}, {6, 5}, {6, 6}, {5, 6}});

    auto r = a & b;
    CHECK_FALSE(r.has_value());
}

TEST_CASE("operator- on polygons: partial overlap", "[boolean]") {
    auto a = poly<2, double>({{0, 0}, {4, 0}, {4, 4}, {0, 4}});
    auto b = poly<2, double>({{2, 0}, {6, 0}, {6, 4}, {2, 4}});

    auto r = a - b;
    REQUIRE(r.has_value());
    // A is 4x4=16, intersection is 2x4=8, difference is 8
    CHECK_THAT(r->area(), WithinAbs(8.0, 1e-8));
}

TEST_CASE("operator- on polygons: no overlap returns full A", "[boolean]") {
    auto a = poly<2, double>({{0, 0}, {1, 0}, {1, 1}, {0, 1}});
    auto b = poly<2, double>({{5, 5}, {6, 5}, {6, 6}, {5, 6}});

    auto r = a - b;
    REQUIRE(r.has_value());
    CHECK_THAT(r->area(), WithinAbs(1.0, 1e-10));
}

TEST_CASE("operator- on polygons: A inside B returns error", "[boolean]") {
    auto a = poly<2, double>({{1, 1}, {2, 1}, {2, 2}, {1, 2}});
    auto b = poly<2, double>({{0, 0}, {4, 0}, {4, 4}, {0, 4}});

    auto r = a - b;
    CHECK_FALSE(r.has_value());
}

TEST_CASE("operator+ on polygons: union is convex hull", "[boolean]") {
    auto a = poly<2, double>({{0, 0}, {2, 0}, {2, 2}, {0, 2}});
    auto b = poly<2, double>({{1, 1}, {3, 1}, {3, 3}, {1, 3}});

    auto r = a + b;
    // Convex hull of both squares
    REQUIRE(r);
    CHECK(r->vertices.size() >= 4);
    CHECK(r->area() >= 7.0);  // min: area(A) + area(B) - area(A∩B) = 4+4-1=7
}

TEST_CASE("operator+ on polygons: non-overlapping", "[boolean]") {
    auto a = poly<2, double>({{0, 0}, {1, 0}, {1, 1}, {0, 1}});
    auto b = poly<2, double>({{5, 0}, {6, 0}, {6, 1}, {5, 1}});

    auto r = a + b;
    REQUIRE(r);
    CHECK(r->area() >= 2.0);  // at least both areas
}

TEST_CASE("union_of named function", "[boolean]") {
    auto a = poly<2, double>({{0, 0}, {2, 0}, {2, 2}, {0, 2}});
    auto b = poly<2, double>({{1, 1}, {3, 1}, {3, 3}, {1, 3}});

    auto r = union_of(a, b);
    REQUIRE(r);
    CHECK(r->vertices.size() >= 4);
    CHECK(r->area() >= 7.0);
}
