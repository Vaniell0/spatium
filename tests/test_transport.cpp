#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/mesh/transport.hpp>
#include <spatium/mesh/primitives.hpp>
#include <spatium/mesh/subdivision.hpp>
#include <spatium/mesh/topology.hpp>
#include <spatium/spaces/sphere.hpp>
#include <spatium/spaces/euclidean.hpp>
#include <cmath>

using namespace spatium;
using namespace spatium::mesh;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

using S2 = Sphere<2>;
using E3 = Euclidean<3>;

// ── Norm preservation ─────────────────────────────────────────

TEST_CASE("Transport preserves norm on S2", "[transport]") {
    S2 sphere;
    auto mesh = subdivide(icosahedron(sphere), sphere, 2);
    auto topo = MeshTopology<S2>::build(mesh);

    auto path = shortest_path(topo, sphere, uint32_t{0}, uint32_t{topo.vertex_count() / 2});
    REQUIRE(path.vertices.size() >= 2);

    // Tangent vector at vertex 0: use log_map to a neighbor
    auto v0 = mesh.vertices[path.vertices[0]];
    auto v1 = mesh.vertices[path.vertices[1]];
    auto tangent = sphere.log_map(v0, v1);

    // Normalize to unit tangent
    auto norm_before = tangent.norm();
    auto transported = parallel_transport(sphere, mesh, path, tangent);
    auto norm_after = transported.norm();

    // Schild's ladder preserves norm to second order
    CHECK_THAT(norm_after, WithinRel(norm_before, 0.1));  // 10% tolerance for coarse mesh
}

// ── Euclidean: vector unchanged ───────────────────────────────

TEST_CASE("Transport in Euclidean space preserves vector", "[transport]") {
    E3 space;
    Mesh<E3> mesh;
    mesh.vertices = {
        Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{2, 0, 0},
        Vec3{0.5, 1, 0}, Vec3{1.5, 1, 0}
    };
    mesh.faces = {{0, 1, 3}, {1, 2, 4}, {1, 4, 3}};

    auto topo = MeshTopology<E3>::build(mesh);
    auto path = shortest_path(topo, space, uint32_t{0}, uint32_t{2});
    REQUIRE(path.vertices.size() >= 2);

    Vec3 tangent{0, 1, 0};  // perpendicular to path
    auto transported = parallel_transport(space, mesh, path, tangent);

    // In flat space, parallel transport preserves direction exactly
    CHECK_THAT(transported[0], WithinAbs(0.0, 1e-10));
    CHECK_THAT(transported[1], WithinAbs(1.0, 1e-10));
    CHECK_THAT(transported[2], WithinAbs(0.0, 1e-10));
}

// ── Single edge transport ─────────────────────────────────────

TEST_CASE("Transport along single edge preserves norm", "[transport]") {
    S2 sphere;
    auto mesh = icosahedron(sphere);
    auto topo = MeshTopology<S2>::build(mesh);

    auto nb = topo.neighbors(0);
    GeodesicPath<S2> path;
    path.vertices = {0, nb[0]};
    path.total_length = sphere.distance(mesh.vertices[0], mesh.vertices[nb[0]]);

    auto v = sphere.log_map(mesh.vertices[0], mesh.vertices[nb[1]]);
    auto norm_before = v.norm();
    auto transported = parallel_transport(sphere, mesh, path, v);
    auto norm_after = transported.norm();

    CHECK_THAT(norm_after, WithinRel(norm_before, 0.05));
}

// ── Round trip ────────────────────────────────────────────────

TEST_CASE("Transport round trip approximate recovery", "[transport]") {
    S2 sphere;
    auto mesh = subdivide(icosahedron(sphere), sphere, 2);
    auto topo = MeshTopology<S2>::build(mesh);

    uint32_t a = 0, b = 5;
    auto path_ab = shortest_path(topo, sphere, a, b);
    auto path_ba = shortest_path(topo, sphere, b, a);
    REQUIRE(path_ab.vertices.size() >= 2);

    auto v0 = mesh.vertices[a];
    auto tangent = sphere.log_map(v0, mesh.vertices[topo.neighbors(a)[0]]);

    auto at_b = parallel_transport(sphere, mesh, path_ab, tangent);
    auto back = parallel_transport(sphere, mesh, path_ba, at_b);

    // On a curved surface, round-trip picks up holonomy.
    // But for a short path, should be close to original.
    auto diff = (back - tangent).norm();
    auto orig_norm = tangent.norm();
    // Allow generous tolerance — Schild's ladder on coarse mesh is approximate
    CHECK(diff / orig_norm < 0.5);
}

// ── Zero vector ───────────────────────────────────────────────

TEST_CASE("Transport zero vector stays zero", "[transport]") {
    S2 sphere;
    auto mesh = icosahedron(sphere);
    auto topo = MeshTopology<S2>::build(mesh);

    auto path = shortest_path(topo, sphere, uint32_t{0}, uint32_t{5});
    Vec<double, 3> zero{0, 0, 0};
    auto result = parallel_transport(sphere, mesh, path, zero);

    CHECK_THAT(result.norm(), WithinAbs(0.0, 1e-10));
}

// ── Trivial path ──────────────────────────────────────────────

TEST_CASE("Transport on trivial path returns original", "[transport]") {
    S2 sphere;
    auto mesh = icosahedron(sphere);

    GeodesicPath<S2> path;
    path.vertices = {0};
    path.total_length = 0.0;

    Vec<double, 3> tangent{0.1, 0.2, 0.0};
    auto result = parallel_transport(sphere, mesh, path, tangent);

    CHECK_THAT(result[0], WithinAbs(0.1, 1e-15));
    CHECK_THAT(result[1], WithinAbs(0.2, 1e-15));
    CHECK_THAT(result[2], WithinAbs(0.0, 1e-15));
}
