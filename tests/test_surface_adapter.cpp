#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/geometry/surface_adapter.hpp>
#include <spatium/geometry/triangle.hpp>
#include <spatium/geometry/circle.hpp>
#include <spatium/geometry/hyperplane.hpp>
#include <spatium/mesh/mesh.hpp>
#include <spatium/mesh/subdivision.hpp>
#include <spatium/point.hpp>

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
