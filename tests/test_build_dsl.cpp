#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/catch_approx.hpp>
#include <any>
#include <concepts>
#include <memory>
#include <spatium/algebra/noise.hpp>
#include <spatium/geometry/concepts.hpp>
#include <spatium/algebra/quaternion.hpp>
#include <spatium/io/build.hpp>
#include <spatium/spaces/offset.hpp>
#include <spatium/spaces/sample.hpp>
#include <cmath>
#include <numbers>
#include <stdexcept>

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

// ── Closure, and the two operations it separates ─────────────────

TEST_CASE("is_closed separates a surface that ends from one that comes back",
          "[parametric]") {
    constexpr double pi = std::numbers::pi;

    CHECK(is_closed(make_torus<double>(2.0, 1.0)));           // periodic both ways
    CHECK_FALSE(is_closed(make_cylinder<double>(1.0, 2.0)));  // two open rims
    CHECK_FALSE(is_closed(make_cone<double>(1.0, 2.0)));      // apex closes one end, the base rim does not

    // The sphere is the case the periodicity flags alone get wrong: v
    // runs over a plain interval and is not periodic, yet both of its
    // edges are poles -- single points the surface passes through rather
    // than rims where it stops.
    auto sphere_fn = [](double u, double v) -> Vec<double, 3> {
        return {std::sin(v) * std::cos(u), std::sin(v) * std::sin(u), std::cos(v)};
    };
    auto sphere = ParametricSurface<double>(sphere_fn, {0.0, 2 * pi, 0.0, pi}, true, false);
    CHECK_FALSE(sphere.periodic_v());  // the flags would have said "open"
    CHECK(is_closed(sphere));          // the geometry says otherwise, and is right

    // Half of it: the poles are still poles, but the u edges are now a
    // genuine rim, so the surface as a whole is not closed.
    auto half = ParametricSurface<double>(sphere_fn, {0.0, pi, 0.0, pi}, false, false);
    CHECK_FALSE(is_closed(half));
}

TEST_CASE("offset() refuses a base that is not a surface at all", "[build_dsl]") {
    bd::Trace<double> scene;
    auto cube = scene.cube();
    // Used to be a std::logic_error out of resolve_surface, raised from
    // inside materialize() and only if the scene was ever materialized.
    CHECK_THROWS_AS(scene.offset(cube, 0.1), std::invalid_argument);
    CHECK(scene.size() == 1);  // and no half-built node was left behind
}

TEST_CASE("offset() refuses an open base; offset_shell() takes it", "[build_dsl]") {
    bd::Trace<double> scene;
    auto tube = scene.cylinder(1.0, 2.0);  // a tube with two open rims
    CHECK_THROWS_AS(scene.offset(tube, 0.1), std::invalid_argument);
    CHECK(scene.size() == 1);

    // The same base is fine once the caller has said what the rim means.
    auto shell = scene.offset_shell(tube, 0.1, bd::EdgeRule::ZeroThickness);
    CHECK(scene.node(shell.index).kind == bd::Kind::Offset);
    CHECK(scene.node(shell.index).edge == bd::EdgeRule::ZeroThickness);
}

TEST_CASE("offset() and offset_shell() agree wherever both are legal", "[build_dsl]") {
    // The split names two operations apart; it does not change either.
    // Given a base with no rim, there is nothing for the rule to do and
    // the two must produce identical geometry.
    bd::Trace<double> scene;
    auto dough = scene.torus(2.0, 1.0);
    auto plain = scene.offset(dough, 0.05);
    auto shell = scene.offset_shell(dough, 0.05, bd::EdgeRule::ZeroThickness);

    auto a = bd::materialize_mesh(scene, plain.index);
    auto b = bd::materialize_mesh(scene, shell.index);
    REQUIRE(a.vertex_count() == b.vertex_count());
    REQUIRE(a.vertex_count() > 0);
    for (std::size_t i = 0; i < a.vertex_count(); ++i)
        CHECK_THAT((a.vertices[i] - b.vertices[i]).norm(), WithinAbs(0.0, 1e-12));
}

TEST_CASE("scatter() refuses a target that is not a surface", "[build_dsl]") {
    // Same hole, same fix: placement happens *on* a space, and a
    // precomputed mesh is not one. An open target stays legal, though --
    // sampling a band by its own area element is well posed.
    bd::Trace<double> scene;
    auto item = scene.cylinder(0.02, 0.1);
    auto slab = scene.cube();
    CHECK_THROWS_AS(scene.scatter(item, slab, 10), std::invalid_argument);

    auto band = scene.cylinder(1.0, 2.0);
    CHECK_NOTHROW(scene.scatter(item, band, 10));
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

// ── The exact form a node keeps, and when it loses it ────────────

TEST_CASE("A torus node keeps the shape it is, not just its (u,v) map",
          "[build_dsl]") {
    bd::Trace<double> scene;
    auto dough = scene.torus(2.0, 1.0);

    REQUIRE(bd::materialize(scene, dough.index).size() == 1);
    auto obj = bd::materialize(scene, dough.index)[0];
    REQUIRE(obj.is_exact());
    CHECK(obj.exact_type() == typeid(geometry::Torus<double>));

    const auto* t = obj.exact_as<geometry::Torus<double>>();
    REQUIRE(t != nullptr);
    CHECK_THAT(t->major_radius, WithinAbs(2.0, 1e-12));
    CHECK_THAT(t->minor_radius, WithinAbs(1.0, 1e-12));

    // Asking for the wrong shape gives nothing rather than garbage, so a
    // renderer's bucketing loop can be a cast attempt.
    CHECK(obj.exact_as<geometry::BoundedQuadric<double>>() == nullptr);
}

TEST_CASE("The exact form agrees with the map it sits beside", "[build_dsl]") {
    // The invariant that makes the exact form usable at all. If these
    // two ever disagree, a render is correct for the shape and wrong for
    // the scene, with nothing to notice -- so it gets checked against
    // the torus's own implicit equation rather than assumed from the
    // fact that both were built from the same two numbers.
    bd::Trace<double> scene;
    auto dough = scene.torus(2.0, 1.0);
    auto obj = bd::materialize(scene, dough.index)[0];
    const auto* t = obj.exact_as<geometry::Torus<double>>();
    REQUIRE(t != nullptr);
    auto surf = obj.surface();
    REQUIRE(surf.has_value());

    const double R = t->major_radius, r = t->minor_radius;
    for (double u : {0.0, 0.7, 2.1, 4.5}) {
        for (double v : {0.0, 1.3, 3.0, 5.5}) {
            auto p = surf->evaluate(u, v);
            // (|p|^2 + R^2 - r^2)^2 - 4 R^2 (x^2 + y^2) = 0
            double n2 = p.dot(p);
            double lhs = (n2 + R * R - r * r) * (n2 + R * R - r * r)
                       - 4.0 * R * R * (p[0] * p[0] + p[1] * p[1]);
            CHECK_THAT(lhs, WithinAbs(0.0, 1e-9));
        }
    }
}

TEST_CASE("space() records no exact form, because it has none", "[build_dsl]") {
    // An arbitrary (u,v) -> R^3 formula has no closed form to record,
    // even when it happens to be a torus underneath: the DSL only knows
    // what the caller told it.
    bd::Trace<double> scene;
    auto same_shape = scene.space(make_torus<double>(2.0, 1.0));
    auto obj = bd::materialize(scene, same_shape.index)[0];
    CHECK_FALSE(obj.is_exact());
    CHECK(obj.is_analytic());   // still a surface, just not a named one
}

TEST_CASE("moving() drops the exact form, isometry or not", "[build_dsl]") {
    using V3 = spatium::Vec<double, 3>;

    // A pure translation *is* an isometry and does send a torus to a
    // torus, so this case could in principle be preserved. It is not,
    // deliberately: deciding that means asking an opaque callable what
    // it does, which is the one question a callable cannot answer. The
    // conservative rule costs the exact path on a moving node and keeps
    // the invariant true, which is the cheaper mistake.
    bd::Trace<double> scene;
    auto moved = scene.torus(2.0, 1.0)
                     .moving([](const V3& p, double) { return V3{p + V3{0, 0, 5}}; });
    CHECK_FALSE(bd::materialize(scene, moved.index)[0].is_exact());

    // And a genuinely non-rigid motion, where keeping it would be wrong.
    bd::Trace<double> squashed_scene;
    auto squashed = squashed_scene.torus(2.0, 1.0)
                        .moving([](const V3& p, double) { return V3{p * 0.5}; });
    CHECK_FALSE(bd::materialize(squashed_scene, squashed.index)[0].is_exact());
}

TEST_CASE("A cylinder node's exact form is clipped to the map's extent",
          "[build_dsl]") {
    // make_cylinder puts v in [0, height], so an unclipped quadric would
    // describe an infinite tube the map never covers.
    bd::Trace<double> scene;
    auto tube = scene.cylinder(0.5, 3.0);
    auto obj = bd::materialize(scene, tube.index)[0];
    REQUIRE(obj.is_exact());

    const auto* q = obj.exact_as<geometry::BoundedQuadric<double>>();
    REQUIRE(q != nullptr);
    auto b = q->bounding_box();
    CHECK_THAT(b.min_corner[2], WithinAbs(0.0, 1e-12));
    CHECK_THAT(b.max_corner[2], WithinAbs(3.0, 1e-12));
    CHECK_THAT(b.max_corner[0], WithinAbs(0.5, 1e-12));
}

// Where the copy-constructibility wall actually stands. Not a bug and
// not a [!shouldfail] -- it is a real, currently-unhit limitation, and
// the point of testing it is that it stops being discovered as a
// template error out of <any> that says nothing about this design. See
// ROADMAP's "Copy-constructibility, leaking upward from std::function".
//
// The trigger for fixing it is structural fields: a Field leaf owning
// state meets this wall by construction, which is a named event rather
// than a hope that someone remembers.
namespace {
struct MoveOnlyShape {
    using ScalarType = double;
    using PointType = spatium::Vec<double, 3>;
    static constexpr std::size_t ambient_dimension = 3;
    std::unique_ptr<int> owned;           // e.g. a device handle, a voxel grid
    PointType centroid() const { return {}; }
    spatium::geometry::Box<3, double> bounding_box() const { return {}; }
};
} // namespace

TEST_CASE("The exact slot takes copyable shapes, and that is the boundary",
          "[build_dsl]") {
    // Why every shape in the tree fits today.
    STATIC_REQUIRE(std::copy_constructible<geometry::Torus<double>>);
    STATIC_REQUIRE(std::copy_constructible<geometry::BoundedQuadric<double>>);

    // And where it stops: std::any requires copy-constructibility of what
    // it stores, so a shape owning move-only state cannot be recorded --
    // even though it is otherwise a perfectly good Bounded shape.
    STATIC_REQUIRE(geometry::Bounded<MoveOnlyShape>);
    STATIC_REQUIRE_FALSE(std::constructible_from<std::any, MoveOnlyShape>);
}

// ── Which level a renderer should use ────────────────────────────

TEST_CASE("The render level is inferred from what the node is", "[build_dsl]") {
    bd::Trace<double> scene;
    auto dough    = scene.torus(2.0, 1.0);
    auto freeform = scene.space(make_torus<double>(2.0, 1.0));
    auto block    = scene.cube();

    CHECK(bd::materialize(scene, dough.index)[0].render_level()    == bd::RenderLevel::Exact);
    CHECK(bd::materialize(scene, freeform.index)[0].render_level() == bd::RenderLevel::Tessellated);
    CHECK(bd::materialize(scene, block.index)[0].render_level()    == bd::RenderLevel::Tessellated);
}

TEST_CASE("Newton is never inferred, only asked for", "[build_dsl]") {
    // At 23 612 ns/ray against a torus's 20.9, a default that reached for
    // Newton would hand a caller three orders of magnitude they never
    // requested. Every inference path lands somewhere else.
    bd::Trace<double> scene;
    auto dough    = scene.torus(2.0, 1.0);
    auto freeform = scene.space(make_torus<double>(2.0, 1.0));
    auto icing    = scene.offset(dough, 0.05);
    for (auto h : {dough, freeform, icing})
        CHECK(bd::materialize(scene, h.index)[0].render_level() != bd::RenderLevel::Newton);

    // ...and asking is enough to get it, on a node that has a map.
    freeform.rendered_as(bd::RenderLevel::Newton);
    CHECK(bd::materialize(scene, freeform.index)[0].render_level() == bd::RenderLevel::Newton);
}

TEST_CASE("An explicit level overrides the inference", "[build_dsl]") {
    bd::Trace<double> scene;
    auto dough = scene.torus(2.0, 1.0);
    REQUIRE(bd::materialize(scene, dough.index)[0].render_level() == bd::RenderLevel::Exact);

    // A shape that *can* be exact may still be asked to tessellate --
    // matching a mixed scene's other objects, say.
    dough.rendered_as(bd::RenderLevel::Tessellated);
    CHECK(bd::materialize(scene, dough.index)[0].render_level() == bd::RenderLevel::Tessellated);
    CHECK(bd::materialize(scene, dough.index)[0].is_exact());  // the form is still there
}

TEST_CASE("rendered_as refuses a level the node cannot serve", "[build_dsl]") {
    bd::Trace<double> scene;

    // No exact form to be exact with.
    auto freeform = scene.space(make_torus<double>(2.0, 1.0));
    CHECK_THROWS_AS(freeform.rendered_as(bd::RenderLevel::Exact), std::invalid_argument);

    // A mesh has no (u,v) map for Newton to iterate on.
    auto block = scene.cube();
    CHECK_THROWS_AS(block.rendered_as(bd::RenderLevel::Newton), std::invalid_argument);

    // A group is many objects, so it has no single mesh.
    auto group = scene.compose({freeform, block});
    CHECK_THROWS_AS(group.rendered_as(bd::RenderLevel::Tessellated), std::invalid_argument);

    // Refusal leaves the node as it was, rather than half-set.
    CHECK(bd::materialize(scene, freeform.index)[0].render_level() == bd::RenderLevel::Tessellated);
}

TEST_CASE("Losing the exact form drops the level back with it", "[build_dsl]") {
    using V3 = spatium::Vec<double, 3>;
    bd::Trace<double> scene;
    auto moving_dough = scene.torus(2.0, 1.0)
                            .moving([](const V3& p, double t) { return V3{p + V3{0, 0, t}}; });
    // .moving() cleared the exact form, so the inference has nothing to
    // reach for and says so rather than promising a closed-form hit that
    // no longer describes this node.
    CHECK_FALSE(bd::materialize(scene, moving_dough.index)[0].is_exact());
    CHECK(bd::materialize(scene, moving_dough.index)[0].render_level()
          == bd::RenderLevel::Tessellated);
}

// ── Rotation in a placement ───────────────────────────────────────

TEST_CASE("rotated() stays a placement and placement_at recovers the turn",
          "[build_dsl]") {
    using V3 = spatium::Vec<double, 3>;
    const double quarter = std::numbers::pi / 2;

    // Axis-angle overload: pi/2 about z. R maps x -> y.
    auto f = rotated(bd::PointField<double>::point(),
                     [quarter](double) { return V3{0.0, 0.0, quarter}; });

    // A rotation is affine in the point, exactly like a scale, so it does
    // not cost the node its instanceability.
    CHECK(f.is_placement());

    auto pl = f.placement_at(bd::MotionEnv<double>{V3{}, 0.0});
    V3 turned{pl.rotation * V3{1.0, 0.0, 0.0}};
    CHECK_THAT(turned[0], WithinAbs(0.0, 1e-12));
    CHECK_THAT(turned[1], WithinAbs(1.0, 1e-12));
    CHECK_THAT(turned[2], WithinAbs(0.0, 1e-12));
    CHECK_THAT(pl.scale, WithinAbs(1.0, 1e-12));
}

TEST_CASE("A translation below a rotation is itself rotated", "[build_dsl]") {
    // The claim Placement's comment makes: translation needed no change
    // when rotation arrived, because "evaluate the motion at p = 0"
    // already reports where the origin ends up -- through every rotation
    // on the way out. Worth pinning, because the plausible-looking
    // alternative (translation as a separate additive term) is wrong here
    // and wrong silently.
    using V3 = spatium::Vec<double, 3>;
    const double quarter = std::numbers::pi / 2;

    auto inner = bd::PointField<double>::point() +
                 bd::PointField<double>::constant(V3{1.0, 0.0, 0.0});
    auto f = rotated(std::move(inner),
                     [quarter](double) { return V3{0.0, 0.0, quarter}; });

    REQUIRE(f.is_placement());
    auto pl = f.placement_at(bd::MotionEnv<double>{V3{}, 0.0});

    // R * (0 + x_hat) = y_hat, not x_hat.
    CHECK_THAT(pl.translation[0], WithinAbs(0.0, 1e-12));
    CHECK_THAT(pl.translation[1], WithinAbs(1.0, 1e-12));
    CHECK_THAT(pl.translation[2], WithinAbs(0.0, 1e-12));

    // And the placement as a whole still agrees with evaluating the field
    // directly, which is the property everything downstream relies on.
    V3 p{0.3, -0.7, 0.2};
    V3 direct = f(p, 0.0);
    V3 via_placement{pl.rotation * V3{p * pl.scale} + pl.translation};
    for (std::size_t i = 0; i < 3; ++i)
        CHECK_THAT(via_placement[i], WithinAbs(direct[i], 1e-12));
}

TEST_CASE("cook() puts a scattered item exactly where materialize_mesh draws it",
          "[build_dsl]") {
    // The agreement that had no test, and the bug that cost. cook() used
    // to keep only site.position and drop the site's {t1, t2, normal}
    // frame, so every scattered object was placed unrotated while
    // materialize_mesh baked the frame into its vertices. Two renderings
    // of one scene, disagreeing -- invisible only because nothing
    // consumed the cooked one yet.
    //
    // The node also carries a placement with all three parts, so the
    // composition order is pinned too: the site position has to pass
    // through the node's scale and rotation, not be added outside them.
    using V3 = spatium::Vec<double, 3>;
    const double quarter = std::numbers::pi / 2;
    const double t = 0.0;

    bd::Trace<double> scene;
    auto dough = scene.torus(2.0, 1.0);
    auto sprinkle = scene.cylinder(0.02, 0.1, 6, 2);
    auto motion = rotated(bd::PointField<double>::point(),
                          [quarter](double) { return V3{quarter, 0.0, 0.0}; });
    motion = scaled(std::move(motion), [](double) { return 0.5; });
    motion += bd::PointField<double>::constant(V3{0.0, 3.0, 0.0});
    auto sprinkles = scene.scatter(sprinkle, dough, 8).moving(std::move(motion));

    auto baked  = bd::materialize_mesh(scene, sprinkles.index, t);
    auto cooked = bd::cook(scene, sprinkles.index, t);

    REQUIRE(cooked.object_count() == 8);
    REQUIRE(cooked.shape_count() == 1);          // one geometry, eight placements
    const auto& rest = cooked.shapes()[cooked.objects()[0].shape].geometry;
    REQUIRE(baked.vertex_count() == 8 * rest.vertex_count());

    for (std::size_t o = 0; o < cooked.objects().size(); ++o) {
        const auto& obj = cooked.objects()[o];
        REQUIRE(obj.instanceable);
        for (std::size_t v = 0; v < rest.vertex_count(); ++v) {
            V3 world{obj.rotation_q.to_matrix() * V3{rest.vertices[v] * obj.scale} + obj.translation};
            const auto& want = baked.vertices[o * rest.vertex_count() + v];
            for (std::size_t k = 0; k < 3; ++k)
                CHECK_THAT(world[k], WithinAbs(want[k], 1e-9));
        }
    }
}

TEST_CASE("cook() hands back geometry a renderer can actually use for a deformation",
          "[build_dsl]") {
    // The regression this fix exists for. cook() used to store the *rest*
    // geometry for a refused (deforming) object and no transform capable
    // of expressing its motion -- because no such transform exists, which
    // is what makes it a deformation. A renderer driven by Cooked alone
    // would therefore have drawn the donut demo's exploding cube
    // unexploded, and nothing noticed for as long as nothing consumed a
    // cooked scene.
    //
    // So for a refused object the stored geometry is the *placed* one,
    // and the test is that it actually moved.
    using V3 = spatium::Vec<double, 3>;
    bd::Trace<double> scene;

    // Reads the point in a way no placement can express: each vertex is
    // pushed along its own direction, so the shape genuinely changes.
    auto burst = scene.cube({1.0, 1.0, 1.0}).moving([](const V3& p, double time) {
        return V3{p + V3{p * (time * 0.5)}};
    });

    auto cooked = bd::cook(scene, burst.index, 2.0);
    REQUIRE(cooked.object_count() == 1);
    REQUIRE_FALSE(cooked.objects()[0].instanceable);
    CHECK(cooked.opaque_refused() == 1);

    const auto& geometry = cooked.shapes()[cooked.objects()[0].shape].geometry;
    const auto& obj = cooked.objects()[0];

    // At t = 2 every vertex is twice as far out as at rest. A corner of
    // the unit cube sits at |(1,1,1)| = sqrt(3); doubled, 2*sqrt(3).
    double farthest = 0.0;
    for (const auto& v : geometry.vertices) {
        V3 world{obj.rotation_q.to_matrix() * V3{v * obj.scale} + obj.translation};
        farthest = std::max(farthest, world.norm());
    }
    CHECK_THAT(farthest, WithinAbs(2.0 * std::sqrt(3.0), 1e-9));

    // And the rest shape really is smaller, so the check above could fail.
    auto rest = bd::cook(scene, burst.index, 0.0);
    double rest_farthest = 0.0;
    for (const auto& v : rest.shapes()[0].geometry.vertices)
        rest_farthest = std::max(rest_farthest, V3{v}.norm());
    CHECK_THAT(rest_farthest, WithinAbs(std::sqrt(3.0), 1e-9));
}

TEST_CASE("cook() resolves colour and emission per object, not per node",
          "[build_dsl]") {
    // cook() used to copy the node's raw Material and never evaluate
    // color_fn or emissive_fn at all, so a renderer reading a cooked scene
    // would have lost every per-particle colour and the whole glow on the
    // donut's letterforms -- silently, since a default Material is a
    // perfectly plausible grey.
    //
    // Scattered objects are the case that matters: they share one node and
    // must still differ, because the field is evaluated where each one
    // ended up.
    using V3 = spatium::Vec<double, 3>;
    bd::Trace<double> scene;
    auto dough = scene.torus(2.0, 1.0);
    auto speck = scene.cube({0.02, 0.02, 0.02});

    // Colour from height, so two objects at different heights differ.
    auto by_height = bd::PointField<double>{
        [](const V3& p, double) { return V3{0.5 + 0.5 * p[2], 0.0, 0.0}; }};
    auto glow = bd::PointField<double>{
        [](const V3& p, double) { return V3{0.0, p[2] > 0.0 ? 1.0 : 0.0, 0.0}; }};

    auto sprinkled = scene.scatter(speck, dough, 24, 7)
                         .colored(std::move(by_height))
                         .glowing(std::move(glow));

    auto cooked = bd::cook(scene, sprinkled.index, 0.0);
    REQUIRE(cooked.object_count() == 24);

    bool any_different = false;
    bool any_glowing = false, any_dark = false;
    for (const auto& o : cooked.objects()) {
        if (std::abs(o.material.base_color[0] - cooked.objects()[0].material.base_color[0]) > 1e-9)
            any_different = true;
        if (o.material.emissive[1] > 0.5) any_glowing = true;
        else any_dark = true;
    }
    CHECK(any_different);   // the colour field was evaluated, and per object
    CHECK(any_glowing);     // so was the emission field
    CHECK(any_dark);        // and it distinguishes objects, rather than being constant
}

TEST_CASE("Every factory's exact form agrees with the chart recorded beside it",
          "[build_dsl]") {
    // The invariant `docs/api-reference.md` already states -- "the exact
    // form must agree with the map" -- checked where the two are first
    // written down together, which is where it is cheapest to check and
    // was the one place nothing looked.
    //
    // `.moving()` has always honoured it by clearing the exact form when
    // a motion could move them apart. A *factory* can violate it at the
    // moment of creation, and `flake()` did: its chart was a spherical
    // cap ending where its radius reached the slab's half-width, while
    // its clip box ran on to the slab's floor, so the exact form carried
    // a skirt the tessellation had never heard of.
    //
    // Comparing bounding boxes rather than surfaces is deliberate. It is
    // cheap, it needs nothing that can itself be wrong, and it catches
    // the failure that matters -- one of the pair covering ground the
    // other does not. A divergence here is a bug until someone adds a
    // documented exception, and there are none.
    using V3 = spatium::Vec<double, 3>;

    auto chart_box = [](const bd::Trace<double>& tr, std::size_t idx) {
        auto s = bd::resolve_surface(tr, idx);
        auto [u0, u1, v0, v1] = s.domain();
        V3 lo{1e30, 1e30, 1e30}, hi{-1e30, -1e30, -1e30};
        constexpr int N = 96;
        for (int i = 0; i <= N; ++i)
            for (int j = 0; j <= N; ++j) {
                auto p = s.evaluate(u0 + (u1 - u0) * i / N, v0 + (v1 - v0) * j / N);
                for (int k = 0; k < 3; ++k) {
                    lo[k] = std::min(lo[k], p[k]);
                    hi[k] = std::max(hi[k], p[k]);
                }
            }
        return std::pair{lo, hi};
    };

    auto agrees = [&](const char* what, const bd::Trace<double>& tr, std::size_t idx,
                      const V3& elo, const V3& ehi, double tol) {
        INFO(what);
        auto [lo, hi] = chart_box(tr, idx);
        for (int k = 0; k < 3; ++k) {
            CHECK_THAT(lo[k], WithinAbs(elo[k], tol));
            CHECK_THAT(hi[k], WithinAbs(ehi[k], tol));
        }
    };

    bd::Trace<double> tr;

    auto t = tr.torus(2.0, 1.0);
    {
        auto* e = std::any_cast<spatium::geometry::Torus<double>>(&tr.node(t.index).exact);
        REQUIRE(e != nullptr);
        auto b = e->bounding_box();
        agrees("torus", tr, t.index, V3{b.min_corner}, V3{b.max_corner}, 1e-9);
    }

    auto c = tr.cylinder(0.011, 0.062, 24, 4);
    {
        auto* e = std::any_cast<spatium::geometry::BoundedQuadric<double>>(&tr.node(c.index).exact);
        REQUIRE(e != nullptr);
        auto b = e->bounding_box();
        agrees("cylinder", tr, c.index, V3{b.min_corner}, V3{b.max_corner}, 1e-9);
    }

    auto sp = tr.sphere(0.5);
    {
        auto* e = std::any_cast<spatium::geometry::BoundedQuadric<double>>(&tr.node(sp.index).exact);
        REQUIRE(e != nullptr);
        auto b = e->bounding_box();
        agrees("sphere", tr, sp.index, V3{b.min_corner}, V3{b.max_corner}, 1e-9);
    }

    // The one that was wrong, and the asymmetric one: a square slab, so
    // `rim` is unambiguous, and a rectangular one, where `rim` is the
    // narrower half-width and the clip must follow it rather than the
    // caller's larger number.
    for (auto half : {V3{0.010, 0.010, 0.003}, V3{0.020, 0.008, 0.004}}) {
        auto f = tr.flake(half);
        auto* e = std::any_cast<spatium::geometry::BoundedQuadric<double>>(&tr.node(f.index).exact);
        REQUIRE(e != nullptr);
        auto b = e->bounding_box();
        agrees("flake", tr, f.index, V3{b.min_corner}, V3{b.max_corner}, 1e-6);
    }
}

TEST_CASE("Same t, different origin, different position -- and still a placement",
          "[build_dsl]") {
    // The test without which `origin` is stored and not read. A field can
    // take the new member, the plumbing can carry it, every existing case
    // can keep passing, and nothing anywhere would notice that the value
    // never reaches the expression -- because every test until now used
    // one instance, where an origin of zero and an origin that is ignored
    // look identical.
    //
    // So: two instances of one node, one time, one motion field, and the
    // only thing that differs between them is where they started.
    using V3 = spatium::Vec<double, 3>;

    // Fly away from the origin along your own direction, at a speed that
    // does not depend on the point being moved.
    auto burst = bd::VecField<double>::opaque_per_instance(
                     [](const V3& origin, double time) { return V3{origin * (time * 3.0)}; }) +
                 bd::VecField<double>::point();

    // Reading the origin must NOT make this a deformation: an origin is
    // one value per object, so the motion is still affine in the point.
    // This is the load-bearing half -- a per-instance motion that refused
    // instancing would be worth nothing.
    REQUIRE(burst.is_placement());

    const double t = 2.0;
    auto at = [&](const V3& origin) {
        return burst.placement_at(bd::MotionEnv<double>{V3{}, t, origin}).translation;
    };

    V3 a = at(V3{1.0, 0.0, 0.0});
    V3 b = at(V3{0.0, -2.0, 0.0});

    CHECK_THAT(a[0], WithinAbs(6.0, 1e-12));   // 1 * 2 * 3
    CHECK_THAT(a[1], WithinAbs(0.0, 1e-12));
    CHECK_THAT(b[1], WithinAbs(-12.0, 1e-12)); // -2 * 2 * 3
    CHECK_THAT(b[0], WithinAbs(0.0, 1e-12));

    // Stated as the property rather than as two numbers, because the
    // numbers could both be right while the mechanism is wrong.
    CHECK((a - b).norm() > 1.0);

    // Same origin, same answer: the field is a function of (origin, t) and
    // of nothing else.
    V3 again = at(V3{1.0, 0.0, 0.0});
    for (int k = 0; k < 3; ++k) CHECK_THAT(again[k], WithinAbs(a[k], 1e-15));
}

TEST_CASE("A scatter's instances fly apart, and materialize_mesh agrees with cook()",
          "[build_dsl]") {
    // The same property one level up, where it has to survive the plumbing
    // rather than just the field: one Scatter node, one motion, and every
    // instance going somewhere different because of where it sat.
    //
    // And the agreement between the two renderings re-checked *with* a
    // per-instance motion in play -- that test has guarded this pair since
    // the site frame was dropped, and the per-site motion path is a new
    // way for them to disagree.
    using V3 = spatium::Vec<double, 3>;
    bd::Trace<double> scene;
    auto ball = scene.sphere(1.0);
    auto speck = scene.cube({0.02, 0.02, 0.02});

    auto fly = bd::VecField<double>::opaque_per_instance(
                   [](const V3& origin, double time) { return V3{origin * (time * 2.0)}; }) +
               bd::VecField<double>::point();

    const double t = 1.0;
    auto burst = scene.scatter(speck, ball, 16, 3).moving(std::move(fly));

    auto cooked = bd::cook(scene, burst.index, t);
    REQUIRE(cooked.object_count() == 16);
    REQUIRE(cooked.opaque_refused() == 0);          // per-instance, still instanceable
    REQUIRE(cooked.shapes().size() == 1);           // one geometry for all of them

    // Every instance ended up somewhere different, and further out than
    // it started: the sphere has radius 1, and t * 2 sends each seat to
    // roughly three times its own distance.
    for (const auto& o : cooked.objects()) CHECK(o.translation.norm() > 2.0);
    for (std::size_t i = 1; i < cooked.objects().size(); ++i)
        CHECK((cooked.objects()[i].translation - cooked.objects()[0].translation).norm() > 1e-6);

    // Vertex for vertex against the other path.
    auto baked = bd::materialize_mesh(scene, burst.index, t);
    const auto& rest = cooked.shapes()[0].geometry;
    REQUIRE(baked.vertex_count() == 16 * rest.vertex_count());
    for (std::size_t o = 0; o < cooked.objects().size(); ++o) {
        const auto& obj = cooked.objects()[o];
        for (std::size_t v = 0; v < rest.vertex_count(); ++v) {
            V3 world{obj.rotation_q.to_matrix() * V3{rest.vertices[v] * obj.scale} + obj.translation};
            const auto& want = baked.vertices[o * rest.vertex_count() + v];
            for (std::size_t k = 0; k < 3; ++k) CHECK_THAT(world[k], WithinAbs(want[k], 1e-9));
        }
    }
}

TEST_CASE("A rotation survives a quaternion round trip to within a few ulp",
          "[build_dsl]") {
    // Pinning a magnitude, not a property. Storing a rotation as a
    // quaternion and rebuilding the matrix is not bit-exact -- measured
    // over the rotations this DSL actually produces (a scatter frame
    // composed with an SO3::exp placement), *zero* of four thousand round
    // trips come back identical, and the worst element moves by about
    // 8 ulp.
    //
    // That is small enough to be invisible and large enough that a frame
    // hash cannot survive it, which is worth knowing before a compaction
    // rather than after a "the picture changed" panic. What this test
    // guards is the *size* of that error: if a future change to
    // from_matrix or to_matrix makes it 1e-8 instead of 1e-15, the
    // rotations are still rotations and nothing else would notice.
    using V3 = spatium::Vec<double, 3>;
    bd::Trace<double> tr;
    auto ball = tr.sphere(1.0);
    auto speck = tr.cube({0.02, 0.02, 0.02});
    auto spun = tr.scatter(speck, ball, 512, 11, 0.5, bd::SeatAxis::X)
                    .moving(rotated(bd::VecField<double>::point(), [](double t) {
                        return V3{0.3 + t, 0.7, -0.4 * t};
                    }));
    auto cooked = bd::cook(tr, spun.index, 1.0);
    REQUIRE(cooked.object_count() == 512);

    // Every object's carried quaternion really is its own rotation, and
    // not a default that happens to sit next to one. At N = 1 those look
    // alike -- an identity quaternion beside an identity matrix agrees
    // perfectly -- which is why the scatter below has hundreds of
    // instances at genuinely different orientations.
    for (const auto& o : cooked.objects()) {
        auto from_stored = o.rotation_q.to_matrix();
        // Still a rotation after the trip: orthonormal columns. The
        // matrix it used to be compared against is gone, so what is left
        // to check is the property rather than the equality.
        V3 c0{from_stored(0, 0), from_stored(1, 0), from_stored(2, 0)};
        V3 c1{from_stored(0, 1), from_stored(1, 1), from_stored(2, 1)};
        CHECK_THAT(c0.norm(), WithinAbs(1.0, 1e-12));
        CHECK_THAT(c0.dot(c1), WithinAbs(0.0, 1e-12));
    }

    // Two quantities, and the second is the one that matters. A matrix
    // element error is dimensionless; what actually moves is a vertex,
    // and it moves by that error times the object's own world radius. A
    // scene built at radius 10^3 has a thousand times the displacement
    // for the same matrix error, and a test that only watched the matrix
    // would call that unchanged.
    double worst = 0.0, worst_shift = 0.0;
    for (const auto& o : cooked.objects()) {
        // Round-tripped once more, from the stored quaternion, so the
        // error being measured is the one a second conversion would add
        // on top of what is already stored.
        auto stored = o.rotation_q.to_matrix();
        auto back = spatium::Quaternion<double>::from_matrix(stored).to_matrix();
        double err = 0.0;
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = 0; j < 3; ++j)
                err = std::max(err, std::abs(back(i, j) - stored(i, j)));
        worst = std::max(worst, err);

        double rest = 0.0;
        for (const auto& v : cooked.shapes()[o.shape].geometry.vertices)
            rest = std::max(rest, V3{v}.norm());
        worst_shift = std::max(worst_shift, err * rest * std::abs(o.scale));

        // Still a rotation afterwards: orthonormal columns, determinant
        // +1. A round trip that drifted off SO(3) would mirror or shear
        // an instance, which is a different failure from losing a few ulp.
        V3 c0{back(0, 0), back(1, 0), back(2, 0)};
        V3 c1{back(0, 1), back(1, 1), back(2, 1)};
        CHECK_THAT(c0.norm(), WithinAbs(1.0, 1e-12));
        CHECK_THAT(c0.dot(c1), WithinAbs(0.0, 1e-12));
    }
    // Both thresholds are alarms rather than values: each sits about two
    // orders above what is measured, so a real degradation trips them
    // while ordinary noise does not. A threshold set far above the
    // measurement -- 1e-8, say -- would let the error grow by a factor of
    // ten million and still pass, which is a test about catastrophes and
    // not about drift.
    INFO("worst element error " << worst << ", worst vertex shift " << worst_shift);
    CHECK(worst < 1e-13);          // measured ~1.7e-15
    CHECK(worst_shift < 1e-12);    // measured ~1e-16 here; the donut's worst object is 10.3 units
}

TEST_CASE("Cooking preserves the object count, whatever the rotations do",
          "[build_dsl]") {
    // A count, checked separately from any picture, because a picture
    // cannot check it. If a rotation were lost -- stored as identity, or
    // dropped somewhere in the compaction -- instances would coalesce
    // toward each other and the render would look *nearly* the same: a
    // slightly denser cluster is not a visible defect, and neither a
    // frame hash nor a pixel diff distinguishes "forty thousand
    // instances" from "forty thousand instances, two of them on top of
    // each other".
    //
    // So: one object per site, every site distinct, and the distinctness
    // asserted on positions rather than inferred from the count.
    using V3 = spatium::Vec<double, 3>;
    bd::Trace<double> scene;
    auto ball = scene.sphere(1.0);
    auto speck = scene.cube({0.03, 0.03, 0.03});

    constexpr std::size_t N = 200;
    auto sown = scene.scatter(speck, ball, N, 17, 0.6, bd::SeatAxis::X)
                    .moving(rotated(bd::VecField<double>::point(),
                                    [](double t) { return V3{0.2, 0.5 + t, -0.3}; }));

    auto cooked = bd::cook(scene, sown.index, 0.5);
    REQUIRE(cooked.object_count() == N);
    REQUIRE(cooked.shape_count() == 1);
    CHECK(cooked.shapes()[0].instances == N);

    // No two objects share a position, and no two share an orientation.
    // Coalescing shows up here and nowhere else.
    std::size_t coincident = 0, same_rotation = 0;
    for (std::size_t i = 0; i < cooked.objects().size(); ++i)
        for (std::size_t j = i + 1; j < cooked.objects().size(); ++j) {
            if ((cooked.objects()[i].translation - cooked.objects()[j].translation).norm() < 1e-9)
                ++coincident;
            auto a = cooked.objects()[i].rotation_q, b = cooked.objects()[j].rotation_q;
            if (std::abs(a.w - b.w) < 1e-12 && std::abs(a.x - b.x) < 1e-12 &&
                std::abs(a.y - b.y) < 1e-12 && std::abs(a.z - b.z) < 1e-12)
                ++same_rotation;
        }
    CHECK(coincident == 0);
    CHECK(same_rotation == 0);

    // And the rotations are not all the identity, which is the specific
    // way a lost rotation would present.
    std::size_t identities = 0;
    for (const auto& o : cooked.objects())
        if (std::abs(o.rotation_q.w - 1.0) < 1e-12) ++identities;
    CHECK(identities == 0);
}
