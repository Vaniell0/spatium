#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/mesh/geodesic.hpp>
#include <spatium/mesh/primitives.hpp>
#include <spatium/mesh/subdivision.hpp>
#include <spatium/spaces/sphere.hpp>
#include <spatium/spaces/euclidean.hpp>
#include <cmath>
#include <set>

using namespace spatium;
using namespace spatium::mesh;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

using S2 = Sphere<2>;
using E3 = Euclidean<3>;

// ── Distance field basics ─────────────────────────────────────

TEST_CASE("Geodesic self-distance is zero", "[geodesic]") {
    S2 sphere;
    auto mesh = icosahedron(sphere);
    auto topo = MeshTopology<S2>::build(mesh);
    auto field = geodesic_distances(topo, sphere, uint32_t{0});
    CHECK(field.distances[0] == 0.0);
}

TEST_CASE("Geodesic adjacent distance equals edge weight", "[geodesic]") {
    S2 sphere;
    auto mesh = icosahedron(sphere);
    auto topo = MeshTopology<S2>::build(mesh);
    auto field = geodesic_distances(topo, sphere, uint32_t{0});

    for (auto v : topo.neighbors(0)) {
        double expected = sphere.distance(mesh.vertices[0], mesh.vertices[v]);
        CHECK_THAT(field.distances[v], WithinAbs(expected, 1e-12));
    }
}

TEST_CASE("Geodesic all vertices reachable on closed mesh", "[geodesic]") {
    S2 sphere;
    auto mesh = icosahedron(sphere);
    auto topo = MeshTopology<S2>::build(mesh);
    auto field = geodesic_distances(topo, sphere, uint32_t{0});

    for (uint32_t v = 0; v < topo.vertex_count(); ++v)
        CHECK(field.distances[v] < 1e10);
}

TEST_CASE("Geodesic triangle inequality", "[geodesic]") {
    S2 sphere;
    auto mesh = icosahedron(sphere);
    auto topo = MeshTopology<S2>::build(mesh);
    auto field = geodesic_distances(topo, sphere, uint32_t{0});

    for (uint32_t e = 0; e < topo.edge_count(); ++e) {
        auto& edge = topo.edge(e);
        double w = sphere.distance(mesh.vertices[edge.v0], mesh.vertices[edge.v1]);
        CHECK(field.distances[edge.v0] <= field.distances[edge.v1] + w + 1e-10);
        CHECK(field.distances[edge.v1] <= field.distances[edge.v0] + w + 1e-10);
    }
}

// ── Sphere-specific ───────────────────────────────────────────

TEST_CASE("Geodesic S2 max distance approaches pi", "[geodesic]") {
    S2 sphere;
    // Subdivide for better approximation
    auto mesh = subdivide(icosahedron(sphere), sphere, 2);
    auto topo = MeshTopology<S2>::build(mesh);
    auto field = geodesic_distances(topo, sphere, uint32_t{0});

    double max_dist = *std::max_element(field.distances.begin(), field.distances.end());
    // On unit sphere, max geodesic distance = pi
    CHECK(max_dist > 2.5);
    CHECK(max_dist < M_PI + 0.5);  // within mesh discretization
}

TEST_CASE("Geodesic S2 improves with subdivision", "[geodesic]") {
    S2 sphere;
    auto mesh_coarse = icosahedron(sphere);
    auto mesh_fine = subdivide(icosahedron(sphere), sphere, 2);

    auto topo_c = MeshTopology<S2>::build(mesh_coarse);
    auto topo_f = MeshTopology<S2>::build(mesh_fine);

    // Pick vertex 0 on coarse mesh (exists in both)
    auto field_c = geodesic_distances(topo_c, sphere, uint32_t{0});
    auto field_f = geodesic_distances(topo_f, sphere, uint32_t{0});

    // Fine mesh max distance should be closer to pi
    double max_c = *std::max_element(field_c.distances.begin(), field_c.distances.end());
    double max_f = *std::max_element(field_f.distances.begin(), field_f.distances.end());
    CHECK(std::abs(max_f - M_PI) < std::abs(max_c - M_PI));
}

// ── Shortest path ─────────────────────────────────────────────

TEST_CASE("Geodesic shortest path source to self", "[geodesic]") {
    S2 sphere;
    auto mesh = icosahedron(sphere);
    auto topo = MeshTopology<S2>::build(mesh);
    auto path = shortest_path(topo, sphere, uint32_t{0}, uint32_t{0});
    CHECK(path.vertices.size() == 1);
    CHECK(path.vertices[0] == 0);
    CHECK(path.total_length == 0.0);
}

TEST_CASE("Geodesic shortest path adjacent vertices", "[geodesic]") {
    S2 sphere;
    auto mesh = icosahedron(sphere);
    auto topo = MeshTopology<S2>::build(mesh);
    auto nb = topo.neighbors(0);
    uint32_t target = nb[0];

    auto path = shortest_path(topo, sphere, uint32_t{0}, target);
    CHECK(path.vertices.size() == 2);
    CHECK(path.vertices.front() == 0);
    CHECK(path.vertices.back() == target);

    double expected = sphere.distance(mesh.vertices[0], mesh.vertices[target]);
    CHECK_THAT(path.total_length, WithinAbs(expected, 1e-12));
}

TEST_CASE("Geodesic shortest path uses valid edges", "[geodesic]") {
    S2 sphere;
    auto mesh = subdivide(icosahedron(sphere), sphere, 1);
    auto topo = MeshTopology<S2>::build(mesh);

    auto path = shortest_path(topo, sphere, uint32_t{0}, uint32_t{topo.vertex_count() - 1});
    REQUIRE(path.vertices.size() >= 2);

    for (std::size_t i = 0; i + 1 < path.vertices.size(); ++i) {
        auto u = path.vertices[i], v = path.vertices[i + 1];
        auto nb = topo.neighbors(u);
        bool found = false;
        for (auto n : nb)
            if (n == v) { found = true; break; }
        CHECK(found);
    }
}

TEST_CASE("Geodesic shortest path length matches distance field", "[geodesic]") {
    S2 sphere;
    auto mesh = subdivide(icosahedron(sphere), sphere, 1);
    auto topo = MeshTopology<S2>::build(mesh);

    uint32_t target = topo.vertex_count() / 2;
    auto path = shortest_path(topo, sphere, uint32_t{0}, target);
    auto field = geodesic_distances(topo, sphere, uint32_t{0});
    CHECK_THAT(path.total_length, WithinAbs(field.distances[target], 1e-10));
}

// ── Multi-source ──────────────────────────────────────────────

TEST_CASE("Geodesic multi-source closer to nearest", "[geodesic]") {
    S2 sphere;
    auto mesh = subdivide(icosahedron(sphere), sphere, 1);
    auto topo = MeshTopology<S2>::build(mesh);

    std::array<uint32_t, 2> sources = {0, topo.vertex_count() / 2};
    auto field = geodesic_distances(topo, sphere, std::span<const uint32_t>(sources));

    // Each source has distance 0
    CHECK(field.distances[sources[0]] == 0.0);
    CHECK(field.distances[sources[1]] == 0.0);

    // Multi-source distance <= single-source distance for each source
    auto field0 = geodesic_distances(topo, sphere, sources[0]);
    auto field1 = geodesic_distances(topo, sphere, sources[1]);

    for (uint32_t v = 0; v < topo.vertex_count(); ++v) {
        double min_single = std::min(field0.distances[v], field1.distances[v]);
        CHECK(field.distances[v] <= min_single + 1e-10);
    }
}

// ── Euclidean ─────────────────────────────────────────────────

TEST_CASE("Geodesic Euclidean flat mesh", "[geodesic]") {
    // Build a small flat Euclidean mesh (triangle strip)
    E3 space;
    Mesh<E3> mesh;
    mesh.vertices = {
        Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0.5, 1, 0},
        Vec3{1.5, 1, 0}, Vec3{2, 0, 0}
    };
    mesh.faces = {{0, 1, 2}, {1, 3, 2}, {1, 4, 3}};

    auto topo = MeshTopology<E3>::build(mesh);
    auto field = geodesic_distances(topo, space, uint32_t{0});

    // Distance to vertex 1 should be 1.0 (direct edge)
    CHECK_THAT(field.distances[1], WithinAbs(1.0, 1e-12));
}
