#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/geometry/ray_surface.hpp>
#include <spatium/geometry/make.hpp>

using namespace spatium;
using namespace spatium::geometry;
using Catch::Matchers::WithinAbs;

// ── Quadric evaluation ──────────────────────────────────────

TEST_CASE("Quadric sphere: on surface evaluates to zero", "[ray_surface]") {
    auto q = Quadric<>::sphere(5.0);
    CHECK_THAT(q(Vec3{5, 0, 0}), WithinAbs(0.0, 1e-10));
    CHECK_THAT(q(Vec3{0, 5, 0}), WithinAbs(0.0, 1e-10));
}

TEST_CASE("Quadric sphere: inside is negative", "[ray_surface]") {
    auto q = Quadric<>::sphere(5.0);
    CHECK(q(Vec3{0, 0, 0}) < 0.0);
}

TEST_CASE("Quadric sphere: outside is positive", "[ray_surface]") {
    auto q = Quadric<>::sphere(5.0);
    CHECK(q(Vec3{10, 0, 0}) > 0.0);
}

// ── Ray-sphere intersection ─────────────────────────────────

TEST_CASE("Ray-sphere: through center, two hits", "[ray_surface]") {
    auto r = unwrap(ray(Vec3{-10, 0, 0}, Vec3{1, 0, 0}));
    auto q = Quadric<>::sphere(3.0);
    auto hits = ray_quadric(r, q);
    REQUIRE(hits.size() == 2);
    CHECK_THAT(hits[0].t, WithinAbs(7.0, 1e-10));   // enters at x=-3
    CHECK_THAT(hits[1].t, WithinAbs(13.0, 1e-10));  // exits at x=+3
    CHECK_THAT(hits[0].point[0], WithinAbs(-3.0, 1e-10));
    CHECK_THAT(hits[1].point[0], WithinAbs(3.0, 1e-10));
}

TEST_CASE("Ray-sphere: tangent, one hit", "[ray_surface]") {
    auto r = unwrap(ray(Vec3{-10, 3, 0}, Vec3{1, 0, 0}));
    auto q = Quadric<>::sphere(3.0);
    auto hits = ray_quadric(r, q);
    // Tangent: discriminant ≈ 0, two roots collapse
    REQUIRE(hits.size() >= 1);
    CHECK_THAT(hits[0].point[1], WithinAbs(3.0, 1e-10));
}

TEST_CASE("Ray-sphere: miss, no hits", "[ray_surface]") {
    auto r = unwrap(ray(Vec3{-10, 5, 0}, Vec3{1, 0, 0}));
    auto q = Quadric<>::sphere(3.0);
    auto hits = ray_quadric(r, q);
    CHECK(hits.empty());
}

TEST_CASE("Ray-sphere: from inside, one hit (forward only)", "[ray_surface]") {
    auto r = unwrap(ray(Vec3{0, 0, 0}, Vec3{1, 0, 0}));
    auto q = Quadric<>::sphere(3.0);
    auto hits = ray_quadric(r, q);
    // One root t<0 (behind), one t>0 (ahead)
    REQUIRE(hits.size() == 1);
    CHECK_THAT(hits[0].t, WithinAbs(3.0, 1e-10));
}

// ── Ray-sphere normals ──────────────────────────────────────

TEST_CASE("Ray-sphere: normal points outward", "[ray_surface]") {
    auto r = unwrap(ray(Vec3{-10, 0, 0}, Vec3{1, 0, 0}));
    auto q = Quadric<>::sphere(3.0);
    auto hits = ray_quadric(r, q);
    REQUIRE(hits.size() == 2);
    // Entry normal should point toward ray origin (-x)
    CHECK_THAT(hits[0].normal[0], WithinAbs(-1.0, 1e-8));
    // Exit normal should point away from center (+x)
    CHECK_THAT(hits[1].normal[0], WithinAbs(1.0, 1e-8));
}

// ── Proximity (miss distance) ───────────────────────────────

TEST_CASE("Ray-sphere proximity: miss distance", "[ray_surface]") {
    auto r = unwrap(ray(Vec3{-10, 5, 0}, Vec3{1, 0, 0}));
    auto q = Quadric<>::sphere(3.0);
    auto prox = ray_quadric_proximity(r, q);
    REQUIRE(prox.has_value());
    // Miss distance should be > 0
    CHECK(prox->miss_distance > 0.0);
    // Closest point should be near x=0 (point of closest approach to center)
    CHECK_THAT(prox->closest_point[0], WithinAbs(0.0, 1e-8));
}

TEST_CASE("Ray-sphere proximity: closer miss → smaller distance", "[ray_surface]") {
    auto q = Quadric<>::sphere(3.0);
    // Ray at y=4 (misses by 1)
    auto r1 = unwrap(ray(Vec3{-10, 4, 0}, Vec3{1, 0, 0}));
    auto p1 = ray_quadric_proximity(r1, q);
    // Ray at y=5 (misses by 2)
    auto r2 = unwrap(ray(Vec3{-10, 5, 0}, Vec3{1, 0, 0}));
    auto p2 = ray_quadric_proximity(r2, q);
    REQUIRE(p1.has_value());
    REQUIRE(p2.has_value());
    CHECK(p1->miss_distance < p2->miss_distance);
}

// ── Ray-quadric full ────────────────────────────────────────

TEST_CASE("ray_quadric_full: hit returns vector", "[ray_surface]") {
    auto r = unwrap(ray(Vec3{-10, 0, 0}, Vec3{1, 0, 0}));
    auto q = Quadric<>::sphere(3.0);
    auto result = ray_quadric_full(r, q);
    CHECK(std::holds_alternative<std::vector<RayHit<double>>>(result));
    CHECK(std::get<std::vector<RayHit<double>>>(result).size() == 2);
}

TEST_CASE("ray_quadric_full: miss returns proximity", "[ray_surface]") {
    auto r = unwrap(ray(Vec3{-10, 5, 0}, Vec3{1, 0, 0}));
    auto q = Quadric<>::sphere(3.0);
    auto result = ray_quadric_full(r, q);
    CHECK(std::holds_alternative<RayProximity<double>>(result));
}

// ── Cylinder ────────────────────────────────────────────────

TEST_CASE("Ray-cylinder: two hits", "[ray_surface]") {
    auto r = unwrap(ray(Vec3{-10, 0, 5}, Vec3{1, 0, 0}));
    auto q = Quadric<>::cylinder_z(3.0);
    auto hits = ray_quadric(r, q);
    REQUIRE(hits.size() == 2);
    CHECK_THAT(hits[0].point[0], WithinAbs(-3.0, 1e-10));
    CHECK_THAT(hits[1].point[0], WithinAbs(3.0, 1e-10));
    // z unchanged
    CHECK_THAT(hits[0].point[2], WithinAbs(5.0, 1e-10));
}

// ── Ellipsoid ───────────────────────────────────────────────

TEST_CASE("Ray-ellipsoid: two hits", "[ray_surface]") {
    auto r = unwrap(ray(Vec3{-10, 0, 0}, Vec3{1, 0, 0}));
    auto q = Quadric<>::ellipsoid(5.0, 3.0, 2.0);
    auto hits = ray_quadric(r, q);
    REQUIRE(hits.size() == 2);
    CHECK_THAT(hits[0].point[0], WithinAbs(-5.0, 1e-10));
    CHECK_THAT(hits[1].point[0], WithinAbs(5.0, 1e-10));
}

// ── Cone ────────────────────────────────────────────────────

TEST_CASE("Ray-cone: two hits through apex", "[ray_surface]") {
    auto r = unwrap(ray(Vec3{-10, 0, 0}, Vec3{1, 0, 0}));
    auto q = Quadric<>::cone_z();
    auto hits = ray_quadric(r, q);
    // Ray along x-axis at z=0 passes through apex (origin)
    // x² + y² - z² = 0 at y=0,z=0 → x²=0 → double root at x=0
    REQUIRE(hits.size() >= 1);
}


// ── Torus ────────────────────────────────────────────────────

TEST_CASE("Ray-torus: axial ray two hits", "[ray_surface]") {
    // Torus R=2, r=0.5 at origin, axis +Z
    Torus<> torus{.major_radius = 2.0, .minor_radius = 0.5};
    // Ray along +X through (2,0,0) — hits tube twice
    auto r = unwrap(ray(Vec3{-10, 0, 0}, Vec3{1, 0, 0}));
    auto hits = ray_torus(r, torus);
    REQUIRE(hits.size() == 4);
    // Expected x: -2.5, -1.5, 1.5, 2.5  → t-offsets 7.5, 8.5, 11.5, 12.5
    CHECK_THAT(hits[0].point[0], WithinAbs(-2.5, 1e-6));
    CHECK_THAT(hits[1].point[0], WithinAbs(-1.5, 1e-6));
    CHECK_THAT(hits[2].point[0], WithinAbs(1.5, 1e-6));
    CHECK_THAT(hits[3].point[0], WithinAbs(2.5, 1e-6));
}

TEST_CASE("Ray-torus: miss returns proximity", "[ray_surface]") {
    Torus<> torus{.major_radius = 2.0, .minor_radius = 0.3};
    // Ray parallel to X far above — misses
    auto r = unwrap(ray(Vec3{-10, 0, 5}, Vec3{1, 0, 0}));
    auto hits = ray_torus(r, torus);
    CHECK(hits.empty());
    auto prox = ray_torus_proximity(r, torus);
    REQUIRE(prox);
    CHECK(prox->miss_distance > 0.0);
}

TEST_CASE("Ray-torus: through hole", "[ray_surface]") {
    // Axial ray along +Z through the hole — no hit
    Torus<> torus{.major_radius = 3.0, .minor_radius = 0.5};
    auto r = unwrap(ray(Vec3{0, 0, -10}, Vec3{0, 0, 1}));
    auto hits = ray_torus(r, torus);
    CHECK(hits.empty());
}

TEST_CASE("Ray-torus: tilted axis", "[ray_surface]") {
    // Torus with Y-axis. Tube centerline = circle of R=2 in XZ plane.
    // Cross section at y=0 → two circles: R+r=2.5 and R-r=1.5.
    // Ray along +Z at x=2: hits outer circle at z=±1.5 (x²=4>2.25 → no inner hit).
    Vec3 axis{0, 1, 0};
    Torus<> torus{.axis = axis, .major_radius = 2.0, .minor_radius = 0.5};
    auto r = unwrap(ray(Vec3{2, 0, -10}, Vec3{0, 0, 1}));
    auto hits = ray_torus(r, torus);
    REQUIRE(hits.size() == 2);
    CHECK_THAT(hits[0].point[2], WithinAbs(-1.5, 1e-6));
    CHECK_THAT(hits[1].point[2], WithinAbs(1.5, 1e-6));
}

TEST_CASE("Ray-torus: tangent ray", "[ray_surface]") {
    // Ray tangent to outer equator: x=R+r, direction +Y, z=0
    Torus<> torus{.major_radius = 2.0, .minor_radius = 0.5};
    auto r = unwrap(ray(Vec3{2.5, -10, 0}, Vec3{0, 1, 0}));
    auto hits = ray_torus(r, torus);
    // Tangent: 1 or 2 very close hits
    CHECK(hits.size() >= 1);
}

