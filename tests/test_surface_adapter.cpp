#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/geometry/surface_adapter.hpp>
#include <spatium/geometry/triangle.hpp>
#include <spatium/geometry/circle.hpp>
#include <spatium/geometry/hyperplane.hpp>
#include <spatium/geometry/ray_surface.hpp>
#include <spatium/mesh/mesh.hpp>
#include <spatium/mesh/subdivision.hpp>
#include <spatium/point.hpp>
#include <spatium/spaces/sphere.hpp>

using namespace spatium;
using namespace spatium::geometry;
using Catch::Matchers::WithinAbs;

// ── Triangle as Surface ────────────────────────────────────────

TEST_CASE("Triangle as Surface satisfies concepts", "[surface_adapter]") {
    auto tri = Triangle3(Vec3{0,0,0}, Vec3{2,0,0}, Vec3{0,2,0});
    auto surface = as_surface(tri);

    static_assert(TopologicalSpace<decltype(surface)>);
    static_assert(MetricSpace<decltype(surface)>);
    static_assert(Manifold<decltype(surface)>);
    static_assert(RiemannianManifold<decltype(surface)>);
    static_assert(Surface<decltype(surface)>);
    SUCCEED();
}

TEST_CASE("Triangle surface contains", "[surface_adapter]") {
    auto surface = as_surface(Triangle3(Vec3{0,0,0}, Vec3{2,0,0}, Vec3{0,2,0}));
    CHECK(surface.contains(Vec3{0.5, 0.5, 0.0}));
    CHECK_FALSE(surface.contains(Vec3{0.5, 0.5, 1.0})); // off-plane
}

TEST_CASE("Triangle surface distance", "[surface_adapter]") {
    auto surface = as_surface(Triangle3(Vec3{0,0,0}, Vec3{2,0,0}, Vec3{0,2,0}));
    auto d = surface.distance(Vec3{0, 0, 0}, Vec3{1, 0, 0});
    CHECK_THAT(d, WithinAbs(1.0, 1e-10));
}

TEST_CASE("Triangle surface exp/log roundtrip", "[surface_adapter]") {
    auto surface = as_surface(Triangle3(Vec3{0,0,0}, Vec3{4,0,0}, Vec3{0,4,0}));
    Vec3 p{1, 1, 0};
    Vec3 q{2, 1, 0};

    auto v = surface.log_map(p, q);
    auto recovered = surface.exp_map(p, v, 1.0);
    CHECK_THAT(surface.distance(q, recovered), WithinAbs(0.0, 1e-8));
}

TEST_CASE("Triangle surface project stays on surface", "[surface_adapter]") {
    auto surface = as_surface(Triangle3(Vec3{0,0,0}, Vec3{2,0,0}, Vec3{0,2,0}));
    auto proj = surface.project(Vec3{0.5, 0.5, 5.0}); // above triangle
    CHECK_THAT(proj[2], WithinAbs(0.0, 1e-10));
}

// ── Mesh on Triangle Surface ───────────────────────────────────

TEST_CASE("Mesh on triangle surface", "[surface_adapter]") {
    auto surface = as_surface(Triangle3(Vec3{0,0,0}, Vec3{4,0,0}, Vec3{0,4,0}));

    // Create a simple mesh on this triangle
    mesh::Mesh<decltype(surface)> m;
    m.vertices = {Vec3{0,0,0}, Vec3{2,0,0}, Vec3{0,2,0}, Vec3{2,2,0}};
    m.faces = {{0, 1, 2}, {1, 3, 2}};

    auto sub = mesh::subdivide_once(m, surface);
    CHECK(sub.face_count() == 8);

    // All vertices should be on the triangle plane (z=0)
    for (const auto& v : sub.vertices)
        CHECK_THAT(v[2], WithinAbs(0.0, 1e-10));
}

// Hyperplane doesn't have normal() method (normal is a field),
// so as_surface(plane) is not directly supported.
// Use Euclidean<3> as Surface instead for infinite flat surfaces.

// ── Direct navigation on surface ───────────────────────────────

TEST_CASE("Navigate on triangle surface", "[surface_adapter]") {
    auto surface = as_surface(Triangle3(Vec3{0,0,0}, Vec3{4,0,0}, Vec3{0,4,0}));
    Vec3 p{1, 1, 0};
    Vec3 q{2, 1, 0};

    auto d = surface.distance(p, q);
    CHECK_THAT(d, WithinAbs(1.0, 1e-10));

    auto tangent = surface.log_map(p, q);
    auto mid = surface.exp_map(p, tangent, 0.5);
    CHECK_THAT(mid[0], WithinAbs(1.5, 1e-8));
}

// ── Quadric sphere as Surface: exact geodesics ──────────────────
//
// ShapeSurface wraps a Quadric::sphere() with HasExactGeodesic, so its
// exp_map/log_map delegate to spaces::Sphere<2,T>'s closed-form
// great-circle formula instead of the generic tangent-plane-projection
// approximation used for every other curved shape (ellipsoid, cylinder,
// cone, Torus). These tests prove that delegation actually happens and
// actually matters numerically.

TEST_CASE("Quadric sphere surface satisfies Surface concept", "[surface_adapter]") {
    auto surface = as_surface(Quadric<double>::sphere(2.0));

    static_assert(HasExactGeodesic<Quadric<double>>);
    static_assert(TopologicalSpace<decltype(surface)>);
    static_assert(MetricSpace<decltype(surface)>);
    static_assert(Manifold<decltype(surface)>);
    static_assert(RiemannianManifold<decltype(surface)>);
    static_assert(Surface<decltype(surface)>);
    SUCCEED();
}

TEST_CASE("Quadric sphere exp_map matches spaces::Sphere<2,T> to near machine precision",
          "[surface_adapter]") {
    const double radius = 2.0;
    auto surface = as_surface(Quadric<double>::sphere(radius));

    Vec3 p{radius, 0, 0};
    Vec3 v{0, radius * 0.8, 0}; // sizable tangent step, not infinitesimal
    double t = 1.0;

    // Reference answer: spaces::Sphere<2,double> directly, via the same
    // translate-to-origin / scale-to-unit-radius change of coordinates
    // that ShapeSurface now performs internally for a round sphere.
    Sphere<2, double> unit_sphere{};
    Vec3 expected = unit_sphere.exp_map(p / radius, v / radius, t) * radius;

    Vec3 actual = surface.exp_map(p, v, t);
    CHECK_THAT((actual - expected).norm(), WithinAbs(0.0, 1e-10));

    // The result must also stay exactly on the sphere.
    CHECK_THAT(actual.norm(), WithinAbs(radius, 1e-9));
}

TEST_CASE("Quadric sphere exp_map is measurably better than the old tangent-plane "
          "approximation at a comparable step size", "[surface_adapter]") {
    const double radius = 2.0;
    auto surface = as_surface(Quadric<double>::sphere(radius));

    Vec3 p{radius, 0, 0};
    Vec3 v{0, radius * 0.8, 0};
    double t = 1.0;

    Sphere<2, double> unit_sphere{};
    Vec3 expected = unit_sphere.exp_map(p / radius, v / radius, t) * radius;

    // New exact path (what ShapeSurface::exp_map now returns for a
    // round-sphere Quadric).
    Vec3 exact_result = surface.exp_map(p, v, t);
    double exact_error = (exact_result - expected).norm();

    // Old approximate path: every curved shape (sphere included) used
    // to just walk the straight line p + v*t and project the endpoint
    // back onto the surface. Quadric::project() reproduces exactly that
    // step for a sphere, so calling it directly reconstructs what
    // ShapeSurface::exp_map returned before this change.
    Vec3 approx_result = Quadric<double>::sphere(radius).project(p + v * t);
    double approx_error = (approx_result - expected).norm();

    CHECK(exact_error < 1e-9);
    CHECK(approx_error > 1e-3); // real, nonzero curvature error
    CHECK(approx_error > 1000 * exact_error); // improvement is orders of magnitude
}

TEST_CASE("Quadric sphere log_map matches spaces::Sphere<2,T> and roundtrips via exp_map",
          "[surface_adapter]") {
    const double radius = 1.5;
    auto surface = as_surface(Quadric<double>::sphere(radius));

    Vec3 p{radius, 0, 0};
    Vec3 q{0, radius, 0}; // a quarter-turn away

    Sphere<2, double> unit_sphere{};
    Vec3 expected_v = unit_sphere.log_map(p / radius, q / radius) * radius;

    Vec3 v = surface.log_map(p, q);
    CHECK_THAT((v - expected_v).norm(), WithinAbs(0.0, 1e-10));

    Vec3 recovered = surface.exp_map(p, v, 1.0);
    CHECK_THAT((recovered - q).norm(), WithinAbs(0.0, 1e-9));
}

TEST_CASE("Quadric ellipsoid surface stays on the approximate path", "[surface_adapter]") {
    // ellipsoid(2,1,1) is not a round sphere, so ShapeSurface must fall
    // back to the tangent-plane-projection approximation -- no
    // HasExactGeodesic special case fires for it.
    auto surface = as_surface(Quadric<double>::ellipsoid(2.0, 1.0, 1.0));

    Vec3 p{2, 0, 0};
    Vec3 v{0, 0, 1};
    Vec3 result = surface.exp_map(p, v, 1.0);

    // Matches the plain project(p + v*t) formula -- the approximate path.
    Vec3 approx = Quadric<double>::ellipsoid(2.0, 1.0, 1.0).project(p + v);
    CHECK_THAT((result - approx).norm(), WithinAbs(0.0, 1e-12));
}
