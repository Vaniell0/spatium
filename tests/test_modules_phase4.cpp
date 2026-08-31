// Phase 4 verification — mesh algorithms via `import spatium.mesh`.
// Note: primitives.hpp, operations.hpp, differential.hpp, heat_geodesic.hpp
// stay header-only (cycle / SSE-BMI conflict). This test does not touch them.
#include <catch2/catch_test_macros.hpp>
#include <cmath>

import spatium.core;
import spatium.algebra;
import spatium.spaces;
import spatium.mesh;

using spatium::Vec3;
using spatium::E3;
using spatium::mesh::Mesh;
using spatium::mesh::MeshTopology;
using spatium::mesh::subdivide_once;

TEST_CASE("module spatium.mesh: subdivide_once + topology", "[modules][phase4][mesh]") {
    // Unit tetrahedron
    Mesh<E3> m;
    m.vertices = {Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0},
                  Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 0.0, 1.0}};
    m.faces = {{0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {0, 2, 3}};

    auto sub = subdivide_once(m, E3{});
    REQUIRE(sub.faces.size() == m.faces.size() * 4);
    REQUIRE(sub.vertices.size() > m.vertices.size());

    auto topo = MeshTopology<E3>::build(m);
    REQUIRE(topo.edge_count() == 6);  // tetrahedron has 6 edges
}

TEST_CASE("module spatium.mesh: Dijkstra geodesic on tetrahedron", "[modules][phase4][mesh]") {
    Mesh<E3> m;
    m.vertices = {Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0},
                  Vec3{0.0, 1.0, 0.0}, Vec3{0.0, 0.0, 1.0}};
    m.faces = {{0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {0, 2, 3}};
    auto topo = MeshTopology<E3>::build(m);

    auto df = spatium::mesh::geodesic_distances(topo, E3{}, 0u);
    REQUIRE(df.distances.size() == 4);
    REQUIRE(df.distances[0] == 0.0);
    REQUIRE(std::abs(df.distances[1] - 1.0) < 1e-12);
    REQUIRE(std::abs(df.distances[2] - 1.0) < 1e-12);
    REQUIRE(std::abs(df.distances[3] - 1.0) < 1e-12);
}
