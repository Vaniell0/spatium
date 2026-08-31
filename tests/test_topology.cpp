#include <catch2/catch_test_macros.hpp>
#include <spatium/mesh/topology.hpp>
#include <spatium/mesh/primitives.hpp>
#include <spatium/mesh/subdivision.hpp>
#include <spatium/spaces/sphere.hpp>
#include <spatium/spaces/euclidean.hpp>
#include <memory>
#include <set>

using namespace spatium;
using namespace spatium::mesh;

using S2 = Sphere<2>;

// ── Icosahedron ───────────────────────────────────────────────

TEST_CASE("Topology icosahedron edge count", "[topology]") {
    S2 sphere;
    auto mesh = icosahedron(sphere);
    auto topo = MeshTopology<S2>::build(mesh);
    // Icosahedron: V=12, F=20, E=30 (Euler: V-E+F=2)
    CHECK(topo.edge_count() == 30);
}

TEST_CASE("Topology icosahedron vertex count", "[topology]") {
    S2 sphere;
    auto mesh = icosahedron(sphere);
    auto topo = MeshTopology<S2>::build(mesh);
    CHECK(topo.vertex_count() == 12);
}

TEST_CASE("Topology icosahedron vertex degree", "[topology]") {
    S2 sphere;
    auto mesh = icosahedron(sphere);
    auto topo = MeshTopology<S2>::build(mesh);
    // Every vertex of an icosahedron has degree 5
    for (uint32_t v = 0; v < topo.vertex_count(); ++v)
        CHECK(topo.neighbors(v).size() == 5);
}

TEST_CASE("Topology icosahedron no boundary", "[topology]") {
    S2 sphere;
    auto mesh = icosahedron(sphere);
    auto topo = MeshTopology<S2>::build(mesh);
    CHECK_FALSE(topo.has_boundary());
}

TEST_CASE("Topology icosahedron all edges have 2 faces", "[topology]") {
    S2 sphere;
    auto mesh = icosahedron(sphere);
    auto topo = MeshTopology<S2>::build(mesh);
    for (uint32_t e = 0; e < topo.edge_count(); ++e) {
        auto& edge = topo.edge(e);
        CHECK(edge.adj_faces[0] != no_face);
        CHECK(edge.adj_faces[1] != no_face);
    }
}

TEST_CASE("Topology icosahedron edge ordering", "[topology]") {
    S2 sphere;
    auto mesh = icosahedron(sphere);
    auto topo = MeshTopology<S2>::build(mesh);
    for (uint32_t e = 0; e < topo.edge_count(); ++e)
        CHECK(topo.edge(e).v0 < topo.edge(e).v1);
}

// ── Tetrahedron ───────────────────────────────────────────────

TEST_CASE("Topology tetrahedron edge count", "[topology]") {
    S2 sphere;
    auto mesh = tetrahedron(sphere);
    auto topo = MeshTopology<S2>::build(mesh);
    // Tetrahedron: V=4, F=4, E=6
    CHECK(topo.edge_count() == 6);
}

TEST_CASE("Topology tetrahedron vertex degree", "[topology]") {
    S2 sphere;
    auto mesh = tetrahedron(sphere);
    auto topo = MeshTopology<S2>::build(mesh);
    for (uint32_t v = 0; v < topo.vertex_count(); ++v)
        CHECK(topo.neighbors(v).size() == 3);
}

// ── Euler formula ─────────────────────────────────────────────

TEST_CASE("Topology Euler formula icosahedron", "[topology]") {
    S2 sphere;
    auto mesh = icosahedron(sphere);
    auto topo = MeshTopology<S2>::build(mesh);
    int V = topo.vertex_count();
    int E = topo.edge_count();
    int F = static_cast<int>(mesh.faces.size());
    CHECK(V - E + F == 2);
}

TEST_CASE("Topology Euler formula subdivided", "[topology]") {
    S2 sphere;
    auto mesh = subdivide(icosahedron(sphere), sphere, 2);
    auto topo = MeshTopology<S2>::build(mesh);
    int V = topo.vertex_count();
    int E = topo.edge_count();
    int F = static_cast<int>(mesh.faces.size());
    CHECK(V - E + F == 2);
}

// ── Face edges ────────────────────────────────────────────────

TEST_CASE("Topology face edges valid", "[topology]") {
    S2 sphere;
    auto mesh = icosahedron(sphere);
    auto topo = MeshTopology<S2>::build(mesh);
    for (uint32_t f = 0; f < mesh.faces.size(); ++f) {
        auto& fe = topo.face_edges(f);
        for (int k = 0; k < 3; ++k)
            CHECK(fe[k] < topo.edge_count());
    }
}

// ── Neighbor symmetry ─────────────────────────────────────────

TEST_CASE("Topology neighbor symmetry", "[topology]") {
    S2 sphere;
    auto mesh = icosahedron(sphere);
    auto topo = MeshTopology<S2>::build(mesh);
    for (uint32_t v = 0; v < topo.vertex_count(); ++v) {
        for (auto n : topo.neighbors(v)) {
            auto nb = topo.neighbors(n);
            bool found = false;
            for (auto x : nb)
                if (x == v) { found = true; break; }
            CHECK(found);
        }
    }
}

// ── Mesh ownership ───────────────────────────────────────────

TEST_CASE("Topology holds mesh copy", "[topology]") {
    S2 sphere;
    auto mesh = icosahedron(sphere);
    auto topo = MeshTopology<S2>::build(mesh);
    CHECK(topo.mesh().vertices.size() == mesh.vertices.size());
    CHECK(topo.mesh().faces.size() == mesh.faces.size());
}

TEST_CASE("Topology shared_ptr avoids copy", "[topology]") {
    S2 sphere;
    auto mesh_ptr = std::make_shared<const Mesh<S2>>(icosahedron(sphere));
    auto topo = MeshTopology<S2>::build(mesh_ptr);
    CHECK(&topo.mesh() == mesh_ptr.get());
}
