#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/catch_approx.hpp>
#include <memory>
#include <spatium/algebra/noise.hpp>
#include <spatium/io/build.hpp>
#include <spatium/spaces/offset.hpp>
#include <spatium/spaces/sample.hpp>
#include <cmath>

using namespace spatium;
using Catch::Matchers::WithinAbs;

// ── algebra::PerlinNoise ─────────────────────────────────────────

TEST_CASE("PerlinNoise stays within [-1, 1] and is deterministic", "[noise]") {
    algebra::PerlinNoise n(42);
    for (int i = 0; i < 2000; ++i) {
        double x = i * 0.037, y = i * 0.11, z = i * 0.07;
        double v = n(x, y, z);
        CHECK(v >= -1.001);
        CHECK(v <= 1.001);
    }
    algebra::PerlinNoise n2(42);
    CHECK(n(1.234, 5.678, 9.012) == n2(1.234, 5.678, 9.012));
}

TEST_CASE("PerlinNoise: different seeds give different fields", "[noise]") {
    algebra::PerlinNoise a(1), b(2);
    CHECK(a(1.234, 5.678, 9.012) != b(1.234, 5.678, 9.012));
}

TEST_CASE("PerlinNoise is continuous", "[noise]") {
    algebra::PerlinNoise n(3);
    double a = n(2.0, 2.0, 2.0);
    double b = n(2.001, 2.0, 2.0);
    CHECK_THAT(a, WithinAbs(b, 0.05));
}

// ── spaces/offset.hpp ─────────────────────────────────────────────

TEST_CASE("offset_surface matches base + thickness*normal", "[offset]") {
    auto dough = make_torus<double>(2.0, 1.0);
    auto icing = offset_surface(dough, 0.05);

    double u = 0.3, v = 1.2;
    auto p_dough = dough.evaluate(u, v);
    auto n = dough.normal_at(u, v);
    auto p_icing = icing.evaluate(u, v);

    CHECK_THAT(p_icing[0], WithinAbs(p_dough[0] + 0.05 * n[0], 1e-9));
    CHECK_THAT(p_icing[1], WithinAbs(p_dough[1] + 0.05 * n[1], 1e-9));
    CHECK_THAT(p_icing[2], WithinAbs(p_dough[2] + 0.05 * n[2], 1e-9));
}

TEST_CASE("offset_surface still satisfies Surface -- project round-trips", "[offset]") {
    auto dough = make_torus<double>(2.0, 1.0);
    auto icing = offset_surface(dough, 0.03);
    auto p = icing.evaluate(0.7, 2.1);
    auto projected = icing.project(p);
    CHECK((projected - p).norm() < 1e-6);
}

TEST_CASE("offset_surface with a thickness field varies per point", "[offset]") {
    auto dough = make_torus<double>(2.0, 1.0);
    auto varying = offset_surface<double>(dough, [](double u, double) { return 0.01 + 0.1 * u; });
    auto near_zero = varying.evaluate(0.0, 0.0);
    auto near_big = varying.evaluate(3.0, 0.0);
    double d0 = (near_zero - dough.evaluate(0.0, 0.0)).norm();
    double d1 = (near_big - dough.evaluate(3.0, 0.0)).norm();
    CHECK(d1 > d0);
}

// ── spaces/sample.hpp ─────────────────────────────────────────────

TEST_CASE("sample_surface_uniform returns the requested count, on-surface", "[sample]") {
    auto dough = make_torus<double>(2.0, 1.0);
    auto pts = sample_surface_uniform(dough, 60, 7);
    REQUIRE(pts.size() == 60);
    for (const auto& s : pts) {
        auto projected = dough.project(s.position);
        CHECK((projected - s.position).norm() < 1e-6);
        CHECK_THAT(s.normal.norm(), WithinAbs(1.0, 1e-6));
    }
}

TEST_CASE("sample_surface_uniform is deterministic per seed", "[sample]") {
    auto dough = make_torus<double>(2.0, 1.0);
    auto a = sample_surface_uniform(dough, 20, 5);
    auto b = sample_surface_uniform(dough, 20, 5);
    REQUIRE(a.size() == b.size());
    for (std::size_t i = 0; i < a.size(); ++i)
        CHECK((a[i].position - b[i].position).norm() < 1e-12);
}

// ── io::build's Trace DSL ─────────────────────────────────────────

namespace bd = spatium::io::build;

TEST_CASE("Trace records node kinds and materializes the right shapes", "[build_dsl]") {
    bd::Trace<double> scene;
    auto dough = scene.torus(2.0, 1.0);
    auto icing = scene.offset(dough, 0.03);
    auto sprinkle = scene.cylinder(0.02, 0.1, 6, 2);
    auto sprinkles = scene.scatter(sprinkle, icing, 30);
    auto root = scene.compose({dough, icing, sprinkles});

    REQUIRE(scene.size() == 5);
    CHECK(scene.node(dough.index).kind == bd::Kind::Space);
    CHECK(scene.node(icing.index).kind == bd::Kind::Offset);
    CHECK(scene.node(sprinkles.index).kind == bd::Kind::Scatter);
    CHECK(scene.node(root.index).kind == bd::Kind::Compose);

    auto placed = bd::materialize(scene, root.index);
    REQUIRE(placed.size() == 3);
    CHECK(placed[0].mesh().vertex_count() > 0);
    CHECK(placed[1].mesh().vertex_count() == placed[0].mesh().vertex_count()); // icing shares dough's UV grid

    auto one_sprinkle = bd::materialize_mesh(scene, sprinkle.index);
    CHECK(placed[2].mesh().vertex_count() == 30 * one_sprinkle.vertex_count()); // 30 instances of the same item mesh

    // The point of the type: a materialized dough/icing object still
    // knows it is a surface, so nothing downstream is forced through
    // triangles.
    CHECK(placed[0].is_analytic());
    CHECK(placed[1].is_analytic());
    CHECK_FALSE(placed[2].is_analytic()); // Scatter is many objects, not one surface
    REQUIRE(placed[0].surface().has_value());
    auto s = *placed[0].surface();
    auto p = s.evaluate(0.3, 0.7);
    auto q = scene.node(dough.index).surface->evaluate(0.3, 0.7);
    CHECK_THAT((p - q).norm(), WithinAbs(0.0, 1e-12)); // same map, not a re-derivation
}

TEST_CASE("Placed::surface() carries the node's motion inside the map", "[build_dsl]") {
    using V3 = spatium::Vec<double, 3>;
    bd::Trace<double> scene;
    auto ring = scene.torus(2.0, 1.0).moving(
        [](const V3& p, double time) { return V3{p + V3{0.0, 0.0, 10.0 * time}}; });

    auto at0 = bd::materialize(scene, ring.index, 0.0);
    auto at1 = bd::materialize(scene, ring.index, 1.0);
    REQUIRE(at0[0].surface().has_value());
    REQUIRE(at1[0].surface().has_value());

    auto p0 = at0[0].surface()->evaluate(0.5, 0.5);
    auto p1 = at1[0].surface()->evaluate(0.5, 0.5);
    // Composed into the map rather than applied to vertices afterwards,
    // so the moved object is still a surface and not a deformed mesh.
    CHECK_THAT(p1[2] - p0[2], WithinAbs(10.0, 1e-9));

    // ...and the surface agrees with what the triangle path produces.
    auto m1 = at1[0].mesh();
    bool any_close = false;
    for (const auto& v : m1.vertices)
        if ((v - p1).norm() < 0.15) any_close = true;
    CHECK(any_close);
}

TEST_CASE("Trace::literal + cube() gives an escape hatch out of the analytic path", "[build_dsl]") {
    bd::Trace<double> scene;
    auto cube = scene.cube({0.5, 0.5, 0.5});
    CHECK(scene.node(cube.index).kind == bd::Kind::Literal);
    auto mesh = bd::materialize_mesh(scene, cube.index);
    CHECK(mesh.vertex_count() == 8);
    CHECK(mesh.face_count() == 12);
}

TEST_CASE("Trace .moving() applies a function of (point, t), no simulation state", "[build_dsl]") {
    bd::Trace<double> scene;
    auto cube = scene.cube({0.5, 0.5, 0.5})
                    .moving([](const Vec<double, 3>& p, double t) { return Vec<double, 3>{p * (1.0 - t)}; });

    auto at0 = bd::materialize_mesh(scene, cube.index, 0.0);
    auto at1 = bd::materialize_mesh(scene, cube.index, 1.0);
    for (const auto& v : at0.vertices) CHECK(v.norm() > 0.1); // untouched cube, real size
    for (const auto& v : at1.vertices) CHECK(v.norm() < 1e-9); // collapsed to the origin

    // materializing at an earlier t again gives the same result -- no
    // hidden state carried between calls.
    auto at0_again = bd::materialize_mesh(scene, cube.index, 0.0);
    for (std::size_t i = 0; i < at0.vertices.size(); ++i)
        CHECK((at0.vertices[i] - at0_again.vertices[i]).norm() < 1e-12);
}

TEST_CASE("resolve_surface throws for a non-Space/Offset node", "[build_dsl]") {
    bd::Trace<double> scene;
    auto cube = scene.cube();
    CHECK_THROWS_AS(bd::resolve_surface(scene, cube.index), std::logic_error);
}

TEST_CASE("moving() composes instead of replacing the previous motion", "[build_dsl]") {
    // Regression: .moving() used to assign the slot, so a second call
    // silently discarded the first and the chain quietly meant something
    // other than it reads as. Chained motion now applies in call order.
    using V3 = spatium::Vec<double, 3>;
    bd::Trace<double> scene;

    auto shift = [](const V3& p, double) { return V3{p + V3{1.0, 0.0, 0.0}}; };
    auto double_it = [](const V3& p, double) { return V3{p * 2.0}; };

    auto a = scene.cube({1.0, 1.0, 1.0}).moving(shift).moving(double_it);
    auto only_shift = scene.cube({1.0, 1.0, 1.0}).moving(shift);
    auto plain = scene.cube({1.0, 1.0, 1.0});

    auto composed = bd::materialize_mesh(scene, a.index, 0.0);
    auto shifted = bd::materialize_mesh(scene, only_shift.index, 0.0);
    auto base = bd::materialize_mesh(scene, plain.index, 0.0);

    REQUIRE(composed.vertex_count() == base.vertex_count());
    for (std::size_t i = 0; i < base.vertex_count(); ++i) {
        // double_it(shift(p)) == (p + x) * 2, not just (p + x) and not p * 2
        V3 expected{(base.vertices[i] + V3{1.0, 0.0, 0.0}) * 2.0};
        CHECK(composed.vertices[i][0] == Catch::Approx(expected[0]));
        CHECK(composed.vertices[i][1] == Catch::Approx(expected[1]));
        CHECK(composed.vertices[i][2] == Catch::Approx(expected[2]));
    }

    // ...and the result is genuinely not the old "last call wins" one.
    // Checked across the mesh rather than per vertex: individual vertices
    // can agree by coincidence (a vertex at x = -1 shifts to 0, and
    // doubling 0 is still 0), so only the whole mesh distinguishes them.
    bool differs_somewhere = false;
    for (std::size_t i = 0; i < base.vertex_count() && !differs_somewhere; ++i)
        for (std::size_t k = 0; k < 3; ++k)
            if (shifted.vertices[i][k] != Catch::Approx(composed.vertices[i][k]))
                differs_somewhere = true;
    CHECK(differs_somewhere);
}

TEST_CASE("A motion hook may own move-only state", "[build_dsl]") {
    // std::function required a copy-constructible callable, so a hook
    // could not own anything move-only and sharing heavy state meant
    // hand-rolling it. PointField is a move_only_function.
    using V3 = spatium::Vec<double, 3>;
    bd::Trace<double> scene;

    auto offset = std::make_unique<V3>(V3{0.0, 3.0, 0.0});
    auto node = scene.cube({1.0, 1.0, 1.0})
                    .moving([off = std::move(offset)](const V3& p, double) { return V3{p + *off}; });

    auto moved = bd::materialize_mesh(scene, node.index, 0.0);
    auto plain = bd::materialize_mesh(scene, scene.cube({1.0, 1.0, 1.0}).index, 0.0);
    REQUIRE(moved.vertex_count() == plain.vertex_count());
    for (std::size_t i = 0; i < plain.vertex_count(); ++i)
        CHECK(moved.vertices[i][1] == Catch::Approx(plain.vertices[i][1] + 3.0));
}
