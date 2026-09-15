// The Chart concept and the chart_of() extension point.
//
// These tests exist to make "the DSL is open to new spaces" a property
// rather than a promise. The claim has two halves and each gets a test:
// a space that ships with the library but had no parametrization until
// now (Sphere<2, T>) goes through the DSL end to end, and a space
// defined *outside* spatium:: -- in a namespace this header knows
// nothing about -- does the same, reached only by ADL.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/io/build.hpp>
#include <spatium/spaces/chart.hpp>
#include <cmath>
#include <concepts>

using namespace spatium;
using Catch::Matchers::WithinAbs;

namespace bd = spatium::io::build;

// ── The concept names what the DSL operations already required ───

TEST_CASE("ParametricSurface is the erased Chart", "[chart]") {
    STATIC_REQUIRE(Chart<ParametricSurface<double>>);
    STATIC_REQUIRE(Chartable<ParametricSurface<double>, double>);
}

TEST_CASE("a space with no (u,v) map is not a Chart, and says so", "[chart]") {
    // Sphere<2> is a perfectly good RiemannianManifold. It is not a
    // chart, and before chart_of() that is precisely why it could not
    // enter the DSL -- not because the node type was closed.
    STATIC_REQUIRE_FALSE(Chart<Sphere<2, double>>);
    STATIC_REQUIRE(Chartable<Sphere<2, double>, double>);

    // Sphere<3> lives in R^4. There is no (u,v) -> R^3 chart to write
    // for it, and none is written.
    STATIC_REQUIRE_FALSE(Chartable<Sphere<3, double>, double>);
}

// ── Sphere<2>'s first chart ──────────────────────────────────────

TEST_CASE("chart_of(Sphere<2>) puts every point on the sphere", "[chart]") {
    const double r = 2.5;
    auto chart = chart_of(Sphere<2, double>{r});

    auto [u_min, u_max, v_min, v_max] = chart.domain();
    const double pi = std::acos(-1.0);
    REQUIRE_THAT(u_min, WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(u_max, WithinAbs(2.0 * pi, 1e-12));
    REQUIRE_THAT(v_min, WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(v_max, WithinAbs(pi, 1e-12));

    // The seam in u is a seam; the v edges are poles, not a boundary,
    // and the flag cannot say so -- which is the whole reason
    // is_closed() reads the geometry instead of these two flags.
    REQUIRE(chart.periodic_u());
    REQUIRE_FALSE(chart.periodic_v());
    REQUIRE(is_closed(chart));

    for (int i = 0; i <= 16; ++i) {
        for (int j = 0; j <= 16; ++j) {
            double u = u_min + (u_max - u_min) * i / 16.0;
            double v = v_min + (v_max - v_min) * j / 16.0;
            REQUIRE_THAT(chart.evaluate(u, v).norm(), WithinAbs(r, 1e-12));
        }
    }

    // The poles: both v edges collapse to a single point.
    auto north = chart.evaluate(0.0, v_min);
    auto south = chart.evaluate(0.0, v_max);
    REQUIRE_THAT(north[2], WithinAbs(r, 1e-12));
    REQUIRE_THAT(south[2], WithinAbs(-r, 1e-12));
    REQUIRE_THAT((chart.evaluate(1.3, v_min) - north).norm(), WithinAbs(0.0, 1e-12));
    REQUIRE_THAT((chart.evaluate(1.3, v_max) - south).norm(), WithinAbs(0.0, 1e-12));

    // area_element is r^2 sin v: maximal at the equator, zero at a pole.
    REQUIRE_THAT(chart.area_element(0.0, pi / 2), WithinAbs(r * r, 1e-6));
    REQUIRE_THAT(chart.area_element(0.0, v_min), WithinAbs(0.0, 1e-6));
}

TEST_CASE("Sphere<2> goes through the DSL and materializes", "[chart][dsl]") {
    const double r = 1.75;
    bd::Trace<double> trace;
    auto s = trace.space(Sphere<2, double>{r}, 32, 16);

    REQUIRE(trace.node(s.index).kind == bd::Kind::Space);
    REQUIRE(trace.node(s.index).surface.has_value());

    auto placed = bd::materialize(trace, s.index, 0.0);
    REQUIRE(placed.size() == 1);

    auto m = placed[0].mesh();
    REQUIRE(m.vertex_count() > 0);
    REQUIRE(m.face_count() > 0);

    for (const auto& p : m.vertices)
        REQUIRE_THAT(p.norm(), WithinAbs(r, 1e-9));
}

TEST_CASE("a sphere node tessellates, because no renderer hits it exactly", "[chart][dsl]") {
    bd::Trace<double> trace;
    auto s = trace.space(Sphere<2, double>{1.0}, 16, 8);

    // chart_of gives the node a map, not an exact analytic form. The
    // `exact` slot means "a renderer can hit this without triangles",
    // and nothing in the library can do that for a Sphere today. Leaving
    // it empty is what keeps render_level() honest -- claiming Exact for
    // a shape no renderer can hit exactly would be a picture that is
    // right about the shape and wrong about the scene.
    REQUIRE_FALSE(trace.node(s.index).exact.has_value());

    auto placed = bd::materialize(trace, s.index, 0.0);
    REQUIRE(placed[0].render_level() == bd::RenderLevel::Tessellated);
}

// ── A space from outside spatium::, reached only by ADL ──────────

namespace someone_elses_library {

// A vertical ruled band around the z axis. Deliberately not derived from
// anything of ours and not mentioned in any spatium header -- the only
// thing connecting it to the DSL is the free function below.
struct Ribbon {
    using ScalarType = double;
    double radius = 1.0;
    double height = 2.0;
};

spatium::ParametricSurface<double> chart_of(const Ribbon& rb) {
    const double pi = std::acos(-1.0);
    return spatium::ParametricSurface<double>(
        [rb](double u, double v) -> spatium::Vec<double, 3> {
            return {rb.radius * std::cos(u), rb.radius * std::sin(u), v};
        },
        {0.0, 2.0 * pi, 0.0, rb.height},
        /*periodic_u=*/true, /*periodic_v=*/false);
}

} // namespace someone_elses_library

TEST_CASE("a foreign space enters the DSL through ADL alone", "[chart][dsl]") {
    using someone_elses_library::Ribbon;

    // No spatium header names Ribbon, and no registry was written to.
    STATIC_REQUIRE(Chartable<Ribbon, double>);

    bd::Trace<double> trace;
    auto h = trace.space(Ribbon{.radius = 3.0, .height = 4.0}, 24, 6);

    auto placed = bd::materialize(trace, h.index, 0.0);
    REQUIRE(placed.size() == 1);

    auto m = placed[0].mesh();
    REQUIRE(m.vertex_count() > 0);
    for (const auto& p : m.vertices) {
        REQUIRE_THAT(std::hypot(p[0], p[1]), WithinAbs(3.0, 1e-9));
        REQUIRE(p[2] >= -1e-9);
        REQUIRE(p[2] <= 4.0 + 1e-9);
    }
}

TEST_CASE("the DSL operations work on a foreign space too", "[chart][dsl]") {
    // The point of entering as a chart rather than as a special case:
    // every operation that consumed ParametricSurface consumes this.
    using someone_elses_library::Ribbon;

    bd::Trace<double> trace;
    auto base  = trace.space(Ribbon{.radius = 2.0, .height = 1.0}, 24, 6);
    auto shell = trace.offset_shell(base, 0.25, bd::EdgeRule::ZeroThickness);
    auto dots  = trace.scatter(trace.cube(Vec<double, 3>{0.05, 0.05, 0.05}), base, 32);

    auto shell_mesh = bd::materialize(trace, shell.index, 0.0)[0].mesh();
    REQUIRE(shell_mesh.vertex_count() > 0);
    for (const auto& p : shell_mesh.vertices)
        REQUIRE_THAT(std::hypot(p[0], p[1]), WithinAbs(2.25, 1e-9));

    auto dots_mesh = bd::materialize(trace, dots.index, 0.0)[0].mesh();
    REQUIRE(dots_mesh.vertex_count() > 0);
}
