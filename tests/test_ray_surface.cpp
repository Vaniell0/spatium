#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/geometry/ray_surface.hpp>
#include <spatium/geometry/make.hpp>
#include <cmath>
#include <numbers>
#include <random>

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

TEST_CASE("Ray-torus: camera-distance sweep along the symmetric axis (solve_quartic resolvent-cubic regression)",
          "[ray_surface]") {
    // Exact repro documented in io/scene.hpp's make_torus(): a torus at
    // the origin (major_radius=1.4, minor_radius=0.35), camera stepping
    // along the ring's own symmetric axis with o=(0,-dist,0), d=(0,1,0),
    // always aimed exactly through the tube's true hit points at
    // y = ±1.05/±1.75. Before the resolvent-cubic fix, solve_quartic()
    // returned a single (t,y)=(inf,inf) "hit" at dist=5/9/15 and two
    // spurious roots collapsed near y=0 at dist=7/12/20, instead of the
    // four genuine crossings every one of these distances actually has.
    Torus<> torus{.major_radius = 1.4, .minor_radius = 0.35};
    for (double dist : {3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 11.6, 12.0, 15.0, 20.0}) {
        auto r = unwrap(ray(Vec3{0, -dist, 0}, Vec3{0, 1, 0}));
        auto hits = ray_torus(r, torus);
        INFO("dist = " << dist);
        REQUIRE(hits.size() == 4);
        for (auto& h : hits) {
            CHECK(std::isfinite(h.t));
            CHECK(std::isfinite(h.point[1]));
        }
        CHECK_THAT(hits[0].point[1], WithinAbs(-1.75, 1e-6));
        CHECK_THAT(hits[1].point[1], WithinAbs(-1.05, 1e-6));
        CHECK_THAT(hits[2].point[1], WithinAbs(1.05, 1e-6));
        CHECK_THAT(hits[3].point[1], WithinAbs(1.75, 1e-6));
    }
}

TEST_CASE("Ray-torus: off-axis hits land on the actual surface (solve_quartic resolvent-cubic regression)",
          "[ray_surface]") {
    // The symmetric-axis sweep above is the exact documented repro, but
    // make_torus()'s comment also reports the bug getting *worse* off
    // that one axis ("a sweep of viewing elevations away from it found
    // zero mathematically valid hits at every angle except the
    // exactly-symmetric one"). Cross-check a handful of hand-picked
    // off-axis rays directly against the implicit torus equation,
    // independent of solve_quartic(), the way make_torus()'s defensive
    // filter does.
    Torus<> torus{.major_radius = 1.4, .minor_radius = 0.35};
    double R = torus.major_radius, r = torus.minor_radius;
    auto on_surface = [&](const Vec3& p) {
        double lp2 = p.dot(p);
        double s = lp2 + R * R - r * r;
        double lhs = s * s;
        double rhs = 4 * R * R * (p[0] * p[0] + p[1] * p[1]);
        double scale = std::max({std::abs(lhs), std::abs(rhs), 1.0});
        return std::abs(lhs - rhs) < 1e-6 * scale;
    };
    struct Case { Vec3 origin, target; };
    Case cases[] = {
        {Vec3{6, 2, 1}, Vec3{0.2, -0.1, 0.3}},
        {Vec3{-5, 3, -2}, Vec3{-0.1, 0.2, -0.2}},
        {Vec3{4, -6, 2.5}, Vec3{0.0, 0.3, 0.1}},
        {Vec3{0, 7, -1.3}, Vec3{0.3, 0.0, -0.2}},
    };
    for (auto& c : cases) {
        Vec3 dir = c.target - c.origin;
        dir = dir / dir.norm();
        auto ry = unwrap(ray(c.origin, dir));
        auto hits = ray_torus(ry, torus);
        INFO("origin = (" << c.origin[0] << "," << c.origin[1] << "," << c.origin[2] << ")");
        REQUIRE(hits.size() >= 1);
        for (auto& h : hits) CHECK(on_surface(h.point));
    }
}

TEST_CASE("Ray-torus: fuzzed off-axis rays always land on the actual surface",
          "[ray_surface]") {
    // Same implicit-equation cross-check as the hand-picked case above,
    // swept over many random rays with a fixed seed for reproducibility.
    // This is the basis for io/scene.hpp's make_torus() comment's claim
    // that its defensive implicit-equation filter is not currently
    // observed to reject any hit solve_quartic() actually produces.
    Torus<> torus{.major_radius = 1.4, .minor_radius = 0.35};
    double R = torus.major_radius, r = torus.minor_radius;
    auto on_surface = [&](const Vec3& p) {
        double lp2 = p.dot(p);
        double s = lp2 + R * R - r * r;
        double lhs = s * s;
        double rhs = 4 * R * R * (p[0] * p[0] + p[1] * p[1]);
        double scale = std::max({std::abs(lhs), std::abs(rhs), 1.0});
        return std::abs(lhs - rhs) < 1e-6 * scale;
    };

    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> ang(0.0, 2 * std::numbers::pi);
    std::uniform_real_distribution<double> elev(-1.4, 1.4);
    std::uniform_real_distribution<double> distd(2.0, 8.0);

    int total_hits = 0;
    for (int trial = 0; trial < 300; ++trial) {
        double theta = ang(rng), dist = distd(rng), z = elev(rng);
        Vec3 origin{dist * std::cos(theta), dist * std::sin(theta), z};
        Vec3 target{elev(rng) * 0.3, elev(rng) * 0.3, elev(rng) * 0.3};
        Vec3 dir = target - origin;
        dir = dir / dir.norm();

        auto ry = unwrap(ray(origin, dir));
        auto hits = ray_torus(ry, torus);
        INFO("trial = " << trial << " theta = " << theta << " dist = " << dist << " z = " << z);
        for (auto& h : hits) {
            CHECK(std::isfinite(h.t));
            CHECK(on_surface(h.point));
        }
        total_hits += static_cast<int>(hits.size());
    }
    // Sanity: rays aimed roughly at the torus should hit it most of the
    // time -- a near-total absence of hits would mean the fuzz aim itself
    // is bad (not exercising the surface at all) rather than confirming
    // anything about solve_quartic().
    CHECK(total_hits > 300);
}

// ── Quadric Shape concept (project/distance/centroid) ───────

TEST_CASE("Quadric satisfies Shape and DistanceQueryable", "[ray_surface]") {
    static_assert(Shape<Quadric<double>>);
    static_assert(DistanceQueryable<Quadric<double>>);
    SUCCEED();
}

TEST_CASE("Quadric sphere: project pulls a point radially onto the surface", "[ray_surface]") {
    auto q = Quadric<>::sphere(5.0);
    auto p = q.project(Vec3{10, 0, 0});
    CHECK_THAT(p[0], WithinAbs(5.0, 1e-10));
    CHECK_THAT(p[1], WithinAbs(0.0, 1e-10));

    // Off-axis point: projection stays exactly on the sphere.
    auto p2 = q.project(Vec3{3, 4, 0}); // already on the |.|=5 shell
    CHECK_THAT(p2[0], WithinAbs(3.0, 1e-9));
    CHECK_THAT(p2[1], WithinAbs(4.0, 1e-9));

    auto p3 = q.project(Vec3{1, 1, 1});
    CHECK_THAT(p3.norm(), WithinAbs(5.0, 1e-9));
}

TEST_CASE("Quadric sphere: distance is zero on the surface, positive off it", "[ray_surface]") {
    auto q = Quadric<>::sphere(2.0);
    CHECK_THAT(q.distance(Vec3{2, 0, 0}), WithinAbs(0.0, 1e-9));
    CHECK_THAT(q.distance(Vec3{4, 0, 0}), WithinAbs(2.0, 1e-9));
    CHECK_THAT(q.distance(Vec3{0, 0, 0}), WithinAbs(2.0, 1e-9));
}

TEST_CASE("Quadric cylinder: project keeps the axial coordinate fixed", "[ray_surface]") {
    auto q = Quadric<>::cylinder_z(2.0);
    auto p = q.project(Vec3{5, 0, 7});
    CHECK_THAT(p[0], WithinAbs(2.0, 1e-9));
    CHECK_THAT(p[2], WithinAbs(7.0, 1e-9)); // Z is untouched -- exact for the cylinder too
}

TEST_CASE("Quadric sphere: is_round_sphere true, exposes center and radius", "[ray_surface]") {
    auto q = Quadric<>::sphere(3.5);
    CHECK(q.is_round_sphere());
    auto c = q.sphere_center();
    CHECK_THAT(c.norm(), WithinAbs(0.0, 1e-10));
    CHECK_THAT(q.sphere_radius(), WithinAbs(3.5, 1e-10));
}

TEST_CASE("Quadric ellipsoid/cylinder/cone are not round spheres", "[ray_surface]") {
    CHECK_FALSE(Quadric<>::ellipsoid(2.0, 1.0, 1.0).is_round_sphere());
    CHECK_FALSE(Quadric<>::cylinder_z(1.0).is_round_sphere());
    CHECK_FALSE(Quadric<>::cone_z().is_round_sphere());
}

TEST_CASE("Quadric: a uniform ellipsoid IS a round sphere", "[ray_surface]") {
    // ellipsoid(r, r, r) is mathematically identical to sphere(r).
    CHECK(Quadric<>::ellipsoid(4.0, 4.0, 4.0).is_round_sphere());
    CHECK_THAT(Quadric<>::ellipsoid(4.0, 4.0, 4.0).sphere_radius(), WithinAbs(4.0, 1e-10));
}

