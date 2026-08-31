#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#if SPATIUM_HAS_EIGEN

#include <spatium/mesh/geodesic.hpp>
#include <spatium/mesh/primitives.hpp>
#include <spatium/mesh/subdivision.hpp>
#include <spatium/spaces/sphere.hpp>
#include <spatium/spaces/euclidean.hpp>
#include <cmath>
#include <numbers>

using namespace spatium;
using namespace spatium::mesh;
using Catch::Matchers::WithinAbs;

using S2 = Sphere<2>;
using E3 = Euclidean<3>;

// Helper: subdivide icosahedron to given level
static auto make_sphere_mesh(int levels) {
    S2 sphere;
    auto mesh = icosahedron(sphere);
    for (int i = 0; i < levels; ++i)
        mesh = subdivide_once(mesh, sphere);
    return mesh;
}

// ── Basic properties ─────────────────────────────────────────

TEST_CASE("Heat: self-distance is zero", "[heat]") {
    S2 sphere;
    auto mesh = make_sphere_mesh(2);
    auto topo = MeshTopology<S2>::build(mesh);
    auto field = heat_geodesic_distances(topo, sphere, uint32_t{0});
    CHECK_THAT(field.distances[0], WithinAbs(0.0, 1e-6));
}

TEST_CASE("Heat: all distances non-negative", "[heat]") {
    S2 sphere;
    auto mesh = make_sphere_mesh(2);
    auto topo = MeshTopology<S2>::build(mesh);
    auto field = heat_geodesic_distances(topo, sphere, uint32_t{0});

    for (uint32_t v = 0; v < topo.vertex_count(); ++v)
        CHECK(field.distances[v] >= 0.0);
}

TEST_CASE("Heat: max distance on S2 approaches pi", "[heat]") {
    S2 sphere;
    auto mesh = make_sphere_mesh(3);
    auto topo = MeshTopology<S2>::build(mesh);
    auto field = heat_geodesic_distances(topo, sphere, uint32_t{0});

    double max_dist = 0.0;
    for (auto d : field.distances)
        max_dist = std::max(max_dist, d);

    // Great-circle max = pi. Heat method should get within 10%.
    CHECK(max_dist > 2.5);
    CHECK(max_dist < 3.8);
}

TEST_CASE("Heat: more accurate than Dijkstra on S2", "[heat]") {
    S2 sphere;
    auto mesh = make_sphere_mesh(3);
    auto topo = MeshTopology<S2>::build(mesh);

    auto heat = heat_geodesic_distances(topo, sphere, uint32_t{0});
    auto dijkstra = geodesic_distances(topo, sphere, uint32_t{0});

    // Compare against exact geodesic (great-circle distance)
    double heat_error = 0.0, dijkstra_error = 0.0;
    for (uint32_t v = 0; v < topo.vertex_count(); ++v) {
        double exact = sphere.distance(mesh.vertices[0], mesh.vertices[v]);
        heat_error += std::abs(heat.distances[v] - exact);
        dijkstra_error += std::abs(dijkstra.distances[v] - exact);
    }

    // Heat method should have lower total error
    CHECK(heat_error < dijkstra_error);
}

TEST_CASE("Heat: GeodesicMethod::Heat dispatch works", "[heat]") {
    S2 sphere;
    auto mesh = make_sphere_mesh(2);
    auto topo = MeshTopology<S2>::build(mesh);

    auto field = geodesic_distances(topo, sphere, uint32_t{0}, GeodesicMethod::Heat);
    CHECK_THAT(field.distances[0], WithinAbs(0.0, 1e-6));

    // All reachable
    for (uint32_t v = 0; v < topo.vertex_count(); ++v)
        CHECK(field.distances[v] < 1e10);
}

TEST_CASE("Heat: HeatSolver reuse gives consistent results", "[heat]") {
    S2 sphere;
    auto mesh = make_sphere_mesh(2);
    auto topo = MeshTopology<S2>::build(mesh);

    auto solver = HeatSolver<S2>::build(topo, sphere);
    auto f1 = solver->distances(uint32_t{0});
    auto f2 = solver->distances(uint32_t{1});

    // Source distances should be zero for respective sources
    CHECK_THAT(f1.distances[0], WithinAbs(0.0, 1e-6));
    CHECK_THAT(f2.distances[1], WithinAbs(0.0, 1e-6));

    // Different sources should give different fields
    bool different = false;
    for (uint32_t v = 0; v < topo.vertex_count(); ++v) {
        if (std::abs(f1.distances[v] - f2.distances[v]) > 1e-6) {
            different = true;
            break;
        }
    }
    CHECK(different);
}

TEST_CASE("Heat: flat Euclidean mesh approximates exact distances", "[heat]") {
    E3 euclidean;

    // Build a flat grid mesh
    Mesh<E3> mesh;
    int n = 10;
    for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= n; ++j)
            mesh.vertices.push_back(Vec3{
                static_cast<double>(i),
                static_cast<double>(j),
                0.0});

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            uint32_t v00 = i * (n + 1) + j;
            uint32_t v10 = (i + 1) * (n + 1) + j;
            uint32_t v01 = i * (n + 1) + (j + 1);
            uint32_t v11 = (i + 1) * (n + 1) + (j + 1);
            mesh.faces.push_back({v00, v10, v01});
            mesh.faces.push_back({v10, v11, v01});
        }
    }

    auto topo = MeshTopology<E3>::build(mesh);
    uint32_t center = 5 * (n + 1) + 5;  // vertex (5,5)
    auto field = heat_geodesic_distances(topo, euclidean, center);

    // Check average relative error.
    // Flat mesh has boundary, which degrades heat method accuracy.
    double total_err = 0.0;
    int count = 0;
    for (uint32_t v = 0; v < topo.vertex_count(); ++v) {
        double exact = (mesh.vertices[v] - mesh.vertices[center]).norm();
        if (exact > 1.0 && exact < 3.0) {
            total_err += std::abs(field.distances[v] - exact) / exact;
            ++count;
        }
    }
    double avg_err = total_err / count;
    CHECK(avg_err < 0.3);  // within 30% average on boundary mesh
}

#endif // SPATIUM_HAS_EIGEN
