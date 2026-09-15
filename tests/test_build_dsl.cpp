#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/catch_approx.hpp>
#include <any>
#include <concepts>
#include <memory>
#include <spatium/algebra/noise.hpp>
#include <spatium/geometry/concepts.hpp>
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
