#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <spatium/geometry/ray_surface.hpp>
#include <spatium/geometry/ray_hit.hpp>
#include <spatium/geometry/make.hpp>
#include <spatium/spatial/bvh.hpp>
#include <vector>
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

// The third instance of the degenerate-leading-coefficient defect, and
// the worst-presenting of them. Closed 2026-09-17 with the other two.
//
// Why this one was worse than a NaN. ray_quadric_proximity indexes
// roots[0] directly rather than iterating, so it read a root that a
// degenerate quadratic does not have. The NaN it got back had an
// imaginary part of exactly zero, and `miss` is abs() of that -- so the
// function returned Result *success* reporting a clean grazing hit, for
// a ray running down the middle of the tube that never approaches the
// wall. closest_t and the point were NaN and conspicuous; `miss` was 0.0
// and looked like an answer, in the one field a caller branches on.
TEST_CASE("ray_quadric_proximity: a ray down the axis is not a grazing hit",
          "[ray_surface]") {
    Quadric<double> cylinder{};
    cylinder.Q = {};
    cylinder.Q(0, 0) = 1.0;
    cylinder.Q(1, 1) = 1.0;
    cylinder.Q(3, 3) = -1.0;  // x^2 + y^2 = 1, infinite in z

    auto axis = unwrap(ray(Vec3{0, 0, -5}, Vec3{0, 0, 1}));
    auto prox = ray_quadric_proximity(axis, cylinder);

    // The right answer is a refusal, not a better number. This function's
    // whole model -- "re is the closest approach, |im| is the miss" --
    // assumes a genuine near-miss, and a ray inside the tube parallel to
    // its axis is not one: it never approaches and never recedes.
    //
    // This comment used to predict the fix would be a degeneracy check on
    // the t² coefficient *before* solve_quadratic, on the grounds that
    // checking after would catch the symptom and leave the semantics
    // wrong. The prediction was made when the return type could not say
    // "fewer roots", and it is superseded rather than merely unmet. The
    // check that shipped asks the *answer* -- fewer than two roots means
    // no conjugate pair to read a miss from -- and that is better than
    // re-deriving the degeneracy condition here, because it leaves one
    // place that decides what "degenerate" means. A coefficient test at
    // the call site would be a second copy of the solver's own rule, free
    // to drift from it.
    CHECK_FALSE(prox.has_value());
    if (!prox) CHECK(prox.error().code == ErrorCode::DegenerateInput);
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


// ── BoundedQuadric ────────────────────────────────────────────

TEST_CASE("Quadric::sphere(center, radius) puts the center in the linear terms",
          "[ray_surface]") {
    Vec<double, 3> c{2.0, -1.0, 0.5};
    auto q = Quadric<>::sphere(c, 1.5);
    CHECK(q.is_round_sphere());
    CHECK_THAT((q.sphere_center() - c).norm(), WithinAbs(0.0, 1e-12));
    CHECK_THAT(q.sphere_radius(), WithinAbs(1.5, 1e-12));
    // On-surface points evaluate to zero; the center evaluates to -r².
    CHECK_THAT(q(c + Vec<double, 3>{1.5, 0, 0}), WithinAbs(0.0, 1e-12));
    CHECK_THAT(q(c), WithinAbs(-2.25, 1e-12));
}

TEST_CASE("BoundedQuadric sphere: box is the surface's own extent", "[ray_surface]") {
    Vec<double, 3> c{2.0, 0.0, 0.0};
    auto s = BoundedQuadric<>::sphere(c, 0.5);
    auto b = s.bounding_box();
    CHECK_THAT(b.min_corner[0], WithinAbs(1.5, 1e-12));
    CHECK_THAT(b.max_corner[0], WithinAbs(2.5, 1e-12));
    CHECK_THAT((s.centroid() - c).norm(), WithinAbs(0.0, 1e-12));
}

TEST_CASE("BoundedQuadric sphere: a tangent hit is not rejected by its own box",
          "[ray_surface]") {
    // The silhouette of a sphere lands exactly on its bounding box, so an
    // exact containment test would discard every grazing ray -- the reason
    // within_clip() carries a tolerance rather than comparing directly.
    auto s = BoundedQuadric<>::sphere(1.0);
    Ray<3, double> grazing{Vec<double, 3>{-5.0, 1.0, 0.0}, Vec<double, 3>{1.0, 0.0, 0.0}};
    auto h = ray_hit(grazing, s);
    REQUIRE(h.has_value());
    CHECK_THAT(h->point[1], WithinAbs(1.0, 1e-6));
}

TEST_CASE("BoundedQuadric cylinder: the clip actually truncates", "[ray_surface]") {
    auto c = BoundedQuadric<>::cylinder_z(1.0, 0.0, 2.0);

    // Through the infinite surface, but above the extent: both algebraic
    // solutions are real and both must be rejected.
    Ray<3, double> above{Vec<double, 3>{-5.0, 0.0, 5.0}, Vec<double, 3>{1.0, 0.0, 0.0}};
    CHECK_FALSE(ray_hit(above, c).has_value());

    Ray<3, double> through{Vec<double, 3>{-5.0, 0.0, 1.0}, Vec<double, 3>{1.0, 0.0, 0.0}};
    auto h = ray_hit(through, c);
    REQUIRE(h.has_value());
    CHECK_THAT(h->point[0], WithinAbs(-1.0, 1e-9)); // near wall, not the far one

    auto b = c.bounding_box();
    CHECK_THAT(b.min_corner[2], WithinAbs(0.0, 1e-12));
    CHECK_THAT(b.max_corner[2], WithinAbs(2.0, 1e-12));
}

TEST_CASE("BoundedQuadric cylinder: a ray down the axis misses", "[ray_surface]") {
    // The ray travels inside the tube and never touches the wall, so the
    // honest answer is "no hit". Algebraically it is the degenerate case:
    // the direction lies in the quadric's null space, the t² coefficient
    // vanishes, and solve_quadratic divides by it. The roots come back NaN,
    // and NaN survives every naive rejection test -- including the clip box,
    // whose two comparisons are both false for NaN. Without a filter phrased
    // to *keep* only real, forward roots, this reports a hit at t = NaN.
    auto c = BoundedQuadric<>::cylinder_z(1.0, 0.0, 2.0);
    Ray<3, double> axis{Vec<double, 3>{0.0, 0.0, -1.0}, Vec<double, 3>{0.0, 0.0, 1.0}};
    CHECK_FALSE(ray_hit(axis, c).has_value());

    // Off-axis but still parallel to it: same degeneracy, and this one runs
    // the length of the tube without ever crossing the wall either.
    Ray<3, double> parallel{Vec<double, 3>{0.5, 0.0, -1.0}, Vec<double, 3>{0.0, 0.0, 1.0}};
    CHECK_FALSE(ray_hit(parallel, c).has_value());
}

// The other half of the NaN story, closed 2026-09-17. Was tagged
// [!shouldfail] so the suite stayed green while the defect stayed
// visible; when the solvers stopped dividing by a vanishing leading
// coefficient it started passing, Catch2 reported that as a failure, and
// the tag came off. See the matching case in test_polynomial.cpp.
TEST_CASE("BoundedQuadric cone: a ray along a generator keeps its hit",
          "[ray_surface]") {
    // A ray parallel to one of the cone's own generators meets the
    // surface exactly once -- a genuine, single, forward hit. Because it
    // is parallel to a generator the t² coefficient is exactly zero: one
    // quantity minus itself, zero to the bit, which is why an exact
    // comparison against zero is the right guard rather than a tolerance.
    // solve_quadratic used to divide by that zero and hand back two NaNs.
    // The leaf rejected them instead of reporting a hit at t = NaN (that
    // was the false positive, fixed earlier), but rejecting is not
    // finding, and the real intersection was silently gone.
    auto cone = BoundedQuadric<>::cone_z(-2.0, 2.0);  // straddles the apex
    double s = 1.0 / std::sqrt(2.0);
    Ray<3, double> along{Vec<double, 3>{0.0, 0.0, -1.0}, Vec<double, 3>{s, 0.0, s}};

    // By hand, on x² + y² = z²: (ts)² = (ts - 1)² gives t = 1/(2s),
    // landing at (0.5, 0, -0.5), comfortably inside the clip box.
    auto h = ray_hit(along, cone);
    REQUIRE(h.has_value());
    CHECK_THAT(h->t, WithinAbs(1.0 / (2.0 * s), 1e-9));
    CHECK_THAT(h->point[0], WithinAbs(0.5, 1e-9));
    CHECK_THAT(h->point[2], WithinAbs(-0.5, 1e-9));
}

TEST_CASE("BoundedQuadric cone: the box widens with distance from the apex",
          "[ray_surface]") {
    // x² + y² = z², so the radius at height z is |z|.
    auto k = BoundedQuadric<>::cone_z(0.0, 3.0);
    auto b = k.bounding_box();
    CHECK_THAT(b.max_corner[0], WithinAbs(3.0, 1e-12));
    CHECK_THAT(b.min_corner[0], WithinAbs(-3.0, 1e-12));

    // Straddling the apex: the wider end sets the width.
    auto straddling = BoundedQuadric<>::cone_z(-1.0, 4.0);
    CHECK_THAT(straddling.bounding_box().max_corner[0], WithinAbs(4.0, 1e-12));
}

TEST_CASE("Torus is bounded by construction, so it needs no clip wrapper",
          "[ray_surface]") {
    // The point of the type gaining these two members: a quadric needs a
    // BoundedQuadric wrapper because most quadrics are infinite, but R
    // and r make a torus finite by definition. So it satisfies Bounded
    // directly, and with its existing ray_hit overload that is all a BVH
    // asks for.
    STATIC_REQUIRE(Bounded<Torus<double>>);
    STATIC_REQUIRE(RayHittable<Torus<double>, double>);

    Torus<> t{.major_radius = 2.0, .minor_radius = 0.5};
    auto b = t.bounding_box();

    // Exact, not a loose sphere of radius R+r: perpendicular to the axis
    // the extent is R + r, but *along* it the torus is only as thick as
    // the tube.
    CHECK_THAT(b.max_corner[0], WithinAbs(2.5, 1e-12));
    CHECK_THAT(b.min_corner[0], WithinAbs(-2.5, 1e-12));
    CHECK_THAT(b.max_corner[2], WithinAbs(0.5, 1e-12));   // not 2.5
    CHECK_THAT(b.min_corner[2], WithinAbs(-0.5, 1e-12));

    // A tilted axis interpolates between the two, rather than falling
    // back to the loose bound.
    Torus<> tilted{.axis = Vec3{0, 1, 0}, .major_radius = 2.0, .minor_radius = 0.5};
    auto tb = tilted.bounding_box();
    CHECK_THAT(tb.max_corner[1], WithinAbs(0.5, 1e-12));  // now y is the thin axis
    CHECK_THAT(tb.max_corner[2], WithinAbs(2.5, 1e-12));
}

TEST_CASE("A BVH of exact tori hits them without a single triangle",
          "[ray_surface]") {
    std::vector<Torus<double>> tori;
    for (int i = 0; i < 8; ++i)
        tori.push_back(Torus<double>{.center = Vec3{i * 6.0, 0, 0},
                                     .major_radius = 2.0,
                                     .minor_radius = 0.5});

    auto tree = spatial::BVH<Torus<double>>::build(std::move(tori));

    // Straight at the fifth torus's near wall, along +x through its tube.
    auto r = unwrap(ray(Vec3{4 * 6.0, -10.0, 0.0}, Vec3{0, 1, 0}));
    auto hit = tree.ray_cast(r);
    REQUIRE(hit.has_value());
    // The tube's outer wall on that side sits at y = -2.5 relative to the
    // torus centre, so the first surface the ray meets is 7.5 away.
    CHECK_THAT(hit->t, WithinAbs(7.5, 1e-6));

    // And a ray down the hole of the first one meets nothing at all.
    auto through = unwrap(ray(Vec3{0, 0, -10.0}, Vec3{0, 0, 1}));
    CHECK_FALSE(tree.ray_cast(through).has_value());
}

TEST_CASE("BoundedQuadric satisfies what a BVH demands of a shape", "[ray_surface]") {
    // The whole point of the type: Bounded lets it into the tree, and
    // RayHittable makes the tree call the exact surface at the leaf
    // rather than a tessellation of it.
    STATIC_REQUIRE(Bounded<BoundedQuadric<double>>);
    STATIC_REQUIRE(RayHittable<BoundedQuadric<double>, double>);
    // A bare Quadric must NOT qualify: most quadrics are unbounded.
    STATIC_REQUIRE_FALSE(Bounded<Quadric<double>>);
}

// ── The instance leaf, oriented ──────────────────────────────

TEST_CASE("Instanced: the ray goes into the rotated shape's own space",
          "[ray_surface]") {
    // A cylinder along z, turned a quarter turn about x so its axis lies
    // along y. Chosen because the *unrotated* version of this exact test
    // is the degenerate case pinned above: a ray straight down the axis,
    // which the solvers currently lose. After the turn the ray is
    // perpendicular to the axis, so a wrong rotation does not merely
    // shift the answer -- it produces no hit at all.
    auto c = BoundedQuadric<>::cylinder_z(0.5, -1.0, 1.0);
    // Written out rather than exponentiated, so this pins ray_hit alone:
    // a quarter turn about x, mapping z to -y.
    Matrix<double, 3, 3> R{};
    R(0, 0) = 1.0; R(1, 2) = -1.0; R(2, 1) = 1.0;

    Instanced<BoundedQuadric<double>> inst{&c, Vec<double, 3>{}, 1.0, R};
    Ray<3, double> probe{Vec<double, 3>{0.0, 0.0, -5.0}, Vec<double, 3>{0.0, 0.0, 1.0}};

    auto h = ray_hit(probe, inst);
    REQUIRE(h.has_value());
    CHECK_THAT(h->t, WithinAbs(4.5, 1e-9));          // wall at z = -0.5
    CHECK_THAT(h->point[2], WithinAbs(-0.5, 1e-9));
    // The normal comes back out through R, so it points at the ray.
    CHECK_THAT(h->normal[2], WithinAbs(-1.0, 1e-9));
    CHECK_THAT(h->normal[0], WithinAbs(0.0, 1e-9));
    CHECK_THAT(h->normal[1], WithinAbs(0.0, 1e-9));

    // Unrotated, the same probe runs down the axis and finds nothing --
    // which is what makes the case above a test of the rotation rather
    // than of the cylinder.
    Instanced<BoundedQuadric<double>> upright{&c, Vec<double, 3>{}, 1.0};
    CHECK_FALSE(ray_hit(probe, upright).has_value());
}

TEST_CASE("Instanced: the box bounds the rotated shape, not the rotated box",
          "[ray_surface]") {
    // Getting this wrong is not a visibly wrong picture -- it is a box
    // that fails to contain its own shape, so rays that should hit are
    // culled before the leaf test ever runs, and the object develops
    // holes that look like a shading bug.
    auto c = BoundedQuadric<>::cylinder_z(0.5, -1.0, 1.0);
    // Written out rather than exponentiated, so this pins ray_hit alone:
    // a quarter turn about x, mapping z to -y.
    Matrix<double, 3, 3> R{};
    R(0, 0) = 1.0; R(1, 2) = -1.0; R(2, 1) = 1.0;

    Instanced<BoundedQuadric<double>> inst{&c, Vec<double, 3>{}, 1.0, R};
    auto b = inst.bounding_box();

    // The upright clip is x,y in [-0.5, 0.5] and z in [-1, 1]; after a
    // quarter turn about x the long axis is y.
    CHECK_THAT(b.min_corner[0], WithinAbs(-0.5, 1e-9));
    CHECK_THAT(b.max_corner[0], WithinAbs(0.5, 1e-9));
    CHECK_THAT(b.min_corner[1], WithinAbs(-1.0, 1e-9));
    CHECK_THAT(b.max_corner[1], WithinAbs(1.0, 1e-9));
    CHECK_THAT(b.min_corner[2], WithinAbs(-0.5, 1e-9));
    CHECK_THAT(b.max_corner[2], WithinAbs(0.5, 1e-9));

    // And it really does contain the hit the ray test reports.
    Ray<3, double> probe{Vec<double, 3>{0.0, 0.0, -5.0}, Vec<double, 3>{0.0, 0.0, 1.0}};
    auto h = ray_hit(probe, inst);
    REQUIRE(h.has_value());
    for (std::size_t i = 0; i < 3; ++i) {
        CHECK(h->point[i] >= b.min_corner[i] - 1e-9);
        CHECK(h->point[i] <= b.max_corner[i] + 1e-9);
    }
}
