#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/mesh/voronoi.hpp>
#include <spatium/mesh/geodesic.hpp>
#include <spatium/mesh/primitives.hpp>
#include <spatium/mesh/subdivision.hpp>
#include <spatium/spaces/sphere.hpp>
#include <algorithm>
#include <set>

using namespace spatium;
using namespace spatium::mesh;
using Catch::Matchers::WithinAbs;

using S2 = Sphere<2>;

// ── Basic properties ──────────────────────────────────────────

TEST_CASE("Voronoi site vertices have distance zero", "[voronoi]") {
    S2 sphere;
    auto mesh = subdivide(icosahedron(sphere), sphere, 1);
    auto topo = MeshTopology<S2>::build(mesh);

    std::array<uint32_t, 3> sites = {0, 10, 20};
    auto vd = geodesic_voronoi(topo, sphere, std::span<const uint32_t>(sites));

    for (uint32_t i = 0; i < sites.size(); ++i) {
        CHECK(vd.distances[sites[i]] == 0.0);
        CHECK(vd.labels[sites[i]] == i);
    }
}

TEST_CASE("Voronoi all vertices assigned", "[voronoi]") {
    S2 sphere;
    auto mesh = subdivide(icosahedron(sphere), sphere, 1);
    auto topo = MeshTopology<S2>::build(mesh);

    std::array<uint32_t, 2> sites = {0, topo.vertex_count() / 2};
    auto vd = geodesic_voronoi(topo, sphere, std::span<const uint32_t>(sites));

    for (uint32_t v = 0; v < topo.vertex_count(); ++v) {
        CHECK(vd.labels[v] != no_vertex);
        CHECK(vd.distances[v] < 1e10);
    }
}

TEST_CASE("Voronoi single site all same label", "[voronoi]") {
    S2 sphere;
    auto mesh = icosahedron(sphere);
    auto topo = MeshTopology<S2>::build(mesh);

    std::array<uint32_t, 1> sites = {0};
    auto vd = geodesic_voronoi(topo, sphere, std::span<const uint32_t>(sites));

    for (uint32_t v = 0; v < topo.vertex_count(); ++v)
        CHECK(vd.labels[v] == 0);
}

TEST_CASE("Voronoi two sites non-empty cells", "[voronoi]") {
    S2 sphere;
    auto mesh = subdivide(icosahedron(sphere), sphere, 2);
    auto topo = MeshTopology<S2>::build(mesh);

    std::array<uint32_t, 2> sites = {0, topo.vertex_count() - 1};
    auto vd = geodesic_voronoi(topo, sphere, std::span<const uint32_t>(sites));

    uint32_t count0 = 0, count1 = 0;
    for (auto l : vd.labels) {
        if (l == 0) count0++;
        else if (l == 1) count1++;
    }
    CHECK(count0 > 0);
    CHECK(count1 > 0);
}

TEST_CASE("Voronoi distance matches single-source Dijkstra", "[voronoi]") {
    S2 sphere;
    auto mesh = subdivide(icosahedron(sphere), sphere, 1);
    auto topo = MeshTopology<S2>::build(mesh);

    std::array<uint32_t, 2> sites = {0, 20};
    auto vd = geodesic_voronoi(topo, sphere, std::span<const uint32_t>(sites));

    auto field0 = geodesic_distances(topo, sphere, sites[0]);
    auto field1 = geodesic_distances(topo, sphere, sites[1]);

    for (uint32_t v = 0; v < topo.vertex_count(); ++v) {
        double min_dist = std::min(field0.distances[v], field1.distances[v]);
        CHECK_THAT(vd.distances[v], WithinAbs(min_dist, 1e-10));
    }
}

// ── Face labels ───────────────────────────────────────────────

TEST_CASE("Voronoi face labels valid", "[voronoi]") {
    S2 sphere;
    auto mesh = subdivide(icosahedron(sphere), sphere, 1);
    auto topo = MeshTopology<S2>::build(mesh);

    std::array<uint32_t, 3> sites = {0, 10, 20};
    auto vd = geodesic_voronoi(topo, sphere, std::span<const uint32_t>(sites));
    auto fl = face_labels(vd, mesh);

    CHECK(fl.size() == mesh.faces.size());
    for (auto l : fl)
        CHECK((l < sites.size() || l == no_vertex));
}

TEST_CASE("Voronoi single site all faces same label", "[voronoi]") {
    S2 sphere;
    auto mesh = icosahedron(sphere);
    auto topo = MeshTopology<S2>::build(mesh);

    std::array<uint32_t, 1> sites = {0};
    auto vd = geodesic_voronoi(topo, sphere, std::span<const uint32_t>(sites));
    auto fl = face_labels(vd, mesh);

    for (auto l : fl)
        CHECK(l == 0);
}

// ── N sites all non-empty ─────────────────────────────────────

TEST_CASE("Voronoi 5 sites all cells non-empty", "[voronoi]") {
    S2 sphere;
    auto mesh = subdivide(icosahedron(sphere), sphere, 2);
    auto topo = MeshTopology<S2>::build(mesh);

    auto nv = topo.vertex_count();
    std::array<uint32_t, 5> sites = {
        0, nv / 5, nv * 2 / 5, nv * 3 / 5, nv * 4 / 5
    };
    auto vd = geodesic_voronoi(topo, sphere, std::span<const uint32_t>(sites));

    std::array<uint32_t, 5> counts{};
    for (auto l : vd.labels) {
        REQUIRE(l < 5);
        counts[l]++;
    }
    for (auto c : counts)
        CHECK(c > 0);
}
