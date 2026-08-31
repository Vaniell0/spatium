#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/spaces/parametric.hpp>
#include <spatium/mesh/topology.hpp>
#include <spatium/mesh/geodesic.hpp>
#include <cmath>

using namespace spatium;
using namespace spatium::mesh;
using Catch::Matchers::WithinAbs;

// ── Torus ─────────────────────────────────────────────────────

TEST_CASE("ParametricSurface torus tessellate", "[parametric]") {
    auto torus = make_torus();
    auto mesh = tessellate(torus, 16, 8);
    CHECK(mesh.vertices.size() == 128);
    CHECK(mesh.faces.size() == 256);
}

TEST_CASE("ParametricSurface torus project", "[parametric]") {
    auto torus = make_torus();
    Vec3 on_torus{3.0, 0.0, 0.0};
    auto proj = torus.project(on_torus);
    CHECK_THAT((proj - on_torus).norm(), WithinAbs(0.0, 0.01));
}

TEST_CASE("ParametricSurface torus normal", "[parametric]") {
    auto torus = make_torus();
    Vec3 p{3.0, 0.0, 0.0};
    auto n = torus.normal(p);
    CHECK(n[0] > 0.5);
    CHECK_THAT(n.norm(), WithinAbs(1.0, 0.01));
}

TEST_CASE("ParametricSurface torus geodesic works", "[parametric]") {
    auto torus = make_torus();
    auto mesh = tessellate(torus, 16, 8);
    auto topo = MeshTopology<ParametricSurface<>>::build(mesh);
    auto field = geodesic_distances(topo, torus, uint32_t{0});
    CHECK(field.distances[0] == 0.0);
    for (auto d : field.distances)
        CHECK(d < 100.0);
}

// ── Cylinder ──────────────────────────────────────────────────

TEST_CASE("ParametricSurface cylinder tessellate", "[parametric]") {
    auto cyl = make_cylinder();
    auto mesh = tessellate(cyl, 12, 4);
    CHECK(mesh.vertices.size() == 60);
    CHECK(mesh.faces.size() == 96);
}

// ── Cone ──────────────────────────────────────────────────────

TEST_CASE("ParametricSurface cone tessellate", "[parametric]") {
    auto cone = make_cone();
    auto mesh = tessellate(cone, 12, 4);
    CHECK(mesh.faces.size() == 96);
}

// ── Möbius strip ──────────────────────────────────────────────

TEST_CASE("ParametricSurface mobius tessellate", "[parametric]") {
    auto mobius = make_mobius();
    auto mesh = tessellate(mobius, 24, 4);
    CHECK(mesh.vertices.size() == 125);
    CHECK(mesh.faces.size() == 192);
}

// ── Custom parametric ─────────────────────────────────────────

TEST_CASE("ParametricSurface custom function", "[parametric]") {
    ParametricSurface<> paraboloid(
        [](double u, double v) -> Vec3 { return {u, v, u*u + v*v}; },
        {-1.0, 1.0, -1.0, 1.0}
    );
    auto mesh = tessellate(paraboloid, 8, 8);
    CHECK(mesh.vertices.size() == 81);
    CHECK(mesh.faces.size() == 128);

    auto p = paraboloid(0.0, 0.0);
    CHECK_THAT(p[2], WithinAbs(0.0, 1e-15));
}

// ── Contains / exp / log ──────────────────────────────────────

TEST_CASE("ParametricSurface contains point on surface", "[parametric]") {
    auto torus = make_torus();
    Vec3 on{3.0, 0.0, 0.0};
    CHECK(torus.contains(on));
}

TEST_CASE("ParametricSurface exp then project stays on surface", "[parametric]") {
    auto torus = make_torus();
    Vec3 p{3.0, 0.0, 0.0};
    Vec3 v{0.0, 0.1, 0.0};
    auto q = torus.exp_map(p, v, 1.0);
    auto proj = torus.project(q);
    CHECK_THAT((q - proj).norm(), WithinAbs(0.0, 0.05));
}

// ── normal_at / parametrization_anisotropy ─────────────────────

TEST_CASE("ParametricSurface normal_at matches normal(p) on a torus", "[parametric]") {
    auto torus = make_torus();
    double u = 0.7, v = 1.3;
    auto direct = torus.normal_at(u, v);
    auto via_point = torus.normal(torus(u, v));
    // normal(p) round-trips through find_params()'s Newton search, so
    // allow a looser tolerance than normal_at's own precision.
    CHECK_THAT((direct - via_point).norm(), WithinAbs(0.0, 1e-6));
}

TEST_CASE("ParametricSurface anisotropy is >= 1 and ~1 for a near-square torus patch",
          "[parametric]") {
    // major_r == minor_r keeps |fu| and |fv| close (not identical -- fu
    // still scales with (R + r*cos(v)) -- but far from the thin-ring
    // extreme below) at v=0 where fu is largest.
    auto torus = make_torus(1.0, 1.0);
    double a = torus.parametrization_anisotropy(0.0, 0.0);
    CHECK(a >= 1.0);
    CHECK(a < 3.0);
}

TEST_CASE("ParametricSurface anisotropy grows with a thin-ring torus", "[parametric]") {
    // R >> r: going around the big loop (fu) covers much more R^3
    // distance per unit u than going around the tube (fv) does per unit
    // v, at every v away from the inner equator's near-cancellation.
    auto fat = make_torus(1.0, 1.0);
    auto thin = make_torus(20.0, 1.0);
    double a_fat = fat.parametrization_anisotropy(0.0, 0.0);
    double a_thin = thin.parametrization_anisotropy(0.0, 0.0);
    CHECK(a_thin > a_fat);
    CHECK(a_thin > 15.0);
}

TEST_CASE("ParametricSurface anisotropy diverges near a cone's apex", "[parametric]") {
    // Near v=height (the apex), the rim circumference -> 0 so fu -> 0
    // while fv (climbing the slant) stays roughly constant -- a genuine
    // degenerate point, not just a global aspect-ratio effect.
    auto cone = make_cone(1.0, 2.0);
    double a_base = cone.parametrization_anisotropy(0.0, 0.1);
    double a_near_apex = cone.parametrization_anisotropy(0.0, 1.99);
    CHECK(a_near_apex > a_base);
}
