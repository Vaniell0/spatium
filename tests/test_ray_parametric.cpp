#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/geometry/ray_parametric.hpp>
#include <spatium/geometry/ray_surface.hpp>
#include <spatium/geometry/make.hpp>
#include <spatium/spaces/parametric.hpp>
#include <numbers>

using namespace spatium;
using namespace spatium::geometry;
using Catch::Matchers::WithinAbs;

// ── Cross-check against analytical quartic on torus ──────────

TEST_CASE("ray_parametric: torus matches analytical ray_torus (axial)", "[ray_parametric]") {
    // Analytical torus: major=2, minor=0.5, axis +Z
    auto surf = make_torus<double>(2.0, 0.5);
    Torus<double> analytical{.major_radius = 2.0, .minor_radius = 0.5};

    // Ray along +X through (2,0,0) → 4 hits at x = ±2.5, ±1.5
    auto r = unwrap(ray(Vec3{-10, 0, 0}, Vec3{1, 0, 0}));
    auto hits = ray_parametric(r, surf);
    auto quartic = ray_torus(r, analytical);

    REQUIRE(quartic.size() == 4);
    REQUIRE(hits.size() == 4);

    // Compare t values (ordered) within Newton tolerance
    for (std::size_t i = 0; i < 4; ++i) {
        CHECK_THAT(hits[i].t, WithinAbs(quartic[i].t, 1e-4));
    }
}

TEST_CASE("ray_parametric: torus normals agree with analytical", "[ray_parametric]") {
    auto surf = make_torus<double>(2.0, 0.5);
    Torus<double> analytical{.major_radius = 2.0, .minor_radius = 0.5};

    auto r = unwrap(ray(Vec3{-10, 0, 0}, Vec3{1, 0, 0}));
    auto hits = ray_parametric(r, surf);
    auto quartic = ray_torus(r, analytical);
    REQUIRE(hits.size() == 4);

    for (std::size_t i = 0; i < 4; ++i) {
        // Normal direction should match up to sign (parameterization orientation)
        auto dot = hits[i].normal.dot(quartic[i].normal);
        CHECK(std::abs(dot) > 0.98);
    }
}

TEST_CASE("ray_parametric: torus miss, empty hits", "[ray_parametric]") {
    auto surf = make_torus<double>(2.0, 0.3);
    // Ray parallel to X, far above tube
    auto r = unwrap(ray(Vec3{-10, 0, 5}, Vec3{1, 0, 0}));
    auto hits = ray_parametric(r, surf);
    CHECK(hits.empty());
}

TEST_CASE("ray_parametric: torus through hole, no hit", "[ray_parametric]") {
    auto surf = make_torus<double>(3.0, 0.5);
    // Axial ray through center → no hit
    auto r = unwrap(ray(Vec3{0, 0, -10}, Vec3{0, 0, 1}));
    auto hits = ray_parametric(r, surf);
    CHECK(hits.empty());
}

// ── Möbius strip ─────────────────────────────────────────────

TEST_CASE("ray_parametric: Möbius strip hit on centerline", "[ray_parametric]") {
    auto mobius = make_mobius<double>(2.0, 0.5);
    // Ray straight down through the +X side (u=0, v=0) → hits strip at (2,0,0)
    auto r = unwrap(ray(Vec3{2, 0, 5}, Vec3{0, 0, -1}));
    auto hits = ray_parametric(r, mobius);
    REQUIRE(!hits.empty());
    CHECK_THAT(hits.front().point[0], WithinAbs(2.0, 1e-3));
    CHECK_THAT(hits.front().point[1], WithinAbs(0.0, 1e-3));
    CHECK_THAT(hits.front().point[2], WithinAbs(0.0, 1e-3));
}

TEST_CASE("ray_parametric: Möbius miss far from strip", "[ray_parametric]") {
    auto mobius = make_mobius<double>(2.0, 0.5);
    // Ray passing 10 units above — miss
    auto r = unwrap(ray(Vec3{-20, 0, 10}, Vec3{1, 0, 0}));
    auto hits = ray_parametric(r, mobius);
    CHECK(hits.empty());
}

// ── Domain / forward-only ────────────────────────────────────

TEST_CASE("ray_parametric: backward ray (origin past surface) rejected", "[ray_parametric]") {
    auto surf = make_torus<double>(2.0, 0.5);
    // Origin past the torus, direction +X → no forward hits
    auto r = unwrap(ray(Vec3{10, 0, 0}, Vec3{1, 0, 0}));
    auto hits = ray_parametric(r, surf);
    CHECK(hits.empty());
}

TEST_CASE("ray_parametric: first-hit convenience", "[ray_parametric]") {
    auto surf = make_torus<double>(2.0, 0.5);
    auto r = unwrap(ray(Vec3{-10, 0, 0}, Vec3{1, 0, 0}));
    auto first = ray_parametric_first(r, surf);
    REQUIRE(first.has_value());
    // Nearest forward hit is x = -2.5 (t = 7.5)
    CHECK_THAT(first->t, WithinAbs(7.5, 1e-4));
}

TEST_CASE("ray_parametric_first: miss returns unexpected", "[ray_parametric]") {
    auto surf = make_torus<double>(2.0, 0.5);
    auto r = unwrap(ray(Vec3{-10, 0, 5}, Vec3{1, 0, 0}));
    auto first = ray_parametric_first(r, surf);
    CHECK_FALSE(first.has_value());
}

// ── Custom parametric (bumpy sphere) ─────────────────────────

TEST_CASE("ray_parametric: custom f(u,v) — displaced sphere", "[ray_parametric]") {
    // Sphere of radius 1 with radial bump: r(u,v) = 1 + 0.1·sin(6u)·sin(4v)
    auto surf = parametric<double>(
        [](double u, double v) -> Vec<double, 3> {
            double r = 1.0 + 0.1 * std::sin(6.0 * u) * std::sin(4.0 * v);
            return {r * std::sin(v) * std::cos(u),
                    r * std::sin(v) * std::sin(u),
                    r * std::cos(v)};
        },
        typename ParametricSurface<double>::Domain{0.0, 2.0 * std::numbers::pi,
                                                    0.01, std::numbers::pi - 0.01},
        true, false);

    // Ray along +X toward origin from outside → hits twice
    auto r = unwrap(ray(Vec3{-5, 0, 0}, Vec3{1, 0, 0}));
    auto hits = ray_parametric(r, surf);
    // Between 2 hits (front/back) — allow extra if bump introduces ripples
    REQUIRE(hits.size() >= 2);
    // First hit near x = -1 ± 0.15 (bump amplitude + margin)
    CHECK(std::abs(hits.front().point[0] + 1.0) < 0.2);
}

// ── Performance fast-paths: uv_hint / cell_radius_hint ────────
// (2026-08-25 — the real ~1000-1700x per-ray slowdown found on Klein
// bottle/Möbius/bumpy-sphere renders; these prove the fast paths return
// the same answer as the exhaustive default, not just "something".)

TEST_CASE("ray_parametric_first: uv_hint seeded with the true answer matches the unhinted result",
          "[ray_parametric]") {
    auto surf = make_torus<double>(2.0, 0.5);
    auto r = unwrap(ray(Vec3{-10, 0, 0}, Vec3{1, 0, 0}));

    auto plain = ray_parametric_first(r, surf);
    REQUIRE(plain.has_value());

    auto hinted = ray_parametric_first(r, surf, RayParametricConfig<double>{},
                                        std::optional{std::make_pair(plain->u, plain->v)});
    REQUIRE(hinted.has_value());
    // t's residual tolerance is looser than u/v's: the hint seeds Newton
    // already essentially at the root, so it stops after the very first
    // ||F|| < cfg.tol check -- u/v match to the last bit, but t's exact
    // value is only as tight as that residual (cfg.tol = 1e-8) allows.
    CHECK_THAT(hinted->t, WithinAbs(plain->t, 1e-6));
    CHECK_THAT(hinted->u, WithinAbs(plain->u, 1e-9));
    CHECK_THAT(hinted->v, WithinAbs(plain->v, 1e-9));
}

TEST_CASE("ray_parametric_first: uv_hint that fails to converge falls back, doesn't fabricate a hit",
          "[ray_parametric]") {
    // Same "through the hole" geometry as the no-hit test above (torus
    // major=2, minor=0.5, ray offset in z past the tube) -- genuinely no
    // intersection exists anywhere on the surface, so a hint (however
    // it's chosen) must not manufacture a fake result once it fails and
    // falls through to the exhaustive scan, which also correctly finds
    // nothing.
    auto surf = make_torus<double>(2.0, 0.5);
    auto r = unwrap(ray(Vec3{-10, 0, 5}, Vec3{1, 0, 0}));

    auto hinted = ray_parametric_first(r, surf, RayParametricConfig<double>{},
                                        std::optional{std::make_pair(0.0, 0.0)});
    CHECK_FALSE(hinted.has_value());
}

TEST_CASE("ray_parametric: a precomputed cell_radius_hint matches the internally-computed default",
          "[ray_parametric]") {
    auto surf = make_torus<double>(2.0, 0.5);
    auto r = unwrap(ray(Vec3{-10, 0, 0}, Vec3{1, 0, 0}));

    auto default_radius = ray_parametric(r, surf);
    REQUIRE(default_radius.size() == 4);

    auto radius = estimate_cell_radius(surf);
    auto hinted_radius = ray_parametric(r, surf, RayParametricConfig<double>{},
                                         std::optional{radius});
    REQUIRE(hinted_radius.size() == 4);
    for (std::size_t i = 0; i < default_radius.size(); ++i)
        CHECK_THAT(hinted_radius[i].t, WithinAbs(default_radius[i].t, 1e-12));
}
