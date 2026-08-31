#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/spatial/bvh.hpp>
#include <spatium/geometry/triangle.hpp>
#include <spatium/geometry/line.hpp>

using namespace spatium;
using namespace spatium::geometry;
using namespace spatium::spatial;
using Catch::Matchers::WithinAbs;

// ── Helpers ───────────────────────────────────────────────────

static std::vector<Triangle3> make_grid_triangles(int nx, int ny) {
    std::vector<Triangle3> tris;
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            double fx = x, fy = y;
            tris.push_back(Triangle3({fx, fy, 0}, {fx + 1, fy, 0}, {fx, fy + 1, 0}));
            tris.push_back(Triangle3({fx + 1, fy, 0}, {fx + 1, fy + 1, 0}, {fx, fy + 1, 0}));
        }
    }
    return tris;
}

// ── Build ─────────────────────────────────────────────────────

TEST_CASE("BVH build from triangles", "[bvh]") {
    auto tris = make_grid_triangles(5, 5);
    REQUIRE(tris.size() == 50);

    auto bvh = BVH<Triangle3>::build(tris);
    CHECK(bvh.node_count() > 0);
    CHECK(bvh.shapes().size() == 50);
}

TEST_CASE("BVH build empty", "[bvh]") {
    auto bvh = BVH<Triangle3>::build({});
    CHECK(bvh.node_count() == 0);
    CHECK(bvh.shapes().empty());
}

TEST_CASE("BVH build single element", "[bvh]") {
    std::vector<Triangle3> tris = {Triangle3({0, 0, 0}, {1, 0, 0}, {0, 1, 0})};
    auto bvh = BVH<Triangle3>::build(tris);
    CHECK(bvh.node_count() == 1);
    CHECK(bvh.shapes().size() == 1);
}

// ── ray_cast ──────────────────────────────────────────────────

TEST_CASE("BVH ray_cast hit", "[bvh]") {
    auto tris = make_grid_triangles(5, 5);
    auto bvh = BVH<Triangle3>::build(tris);

    // Shoot ray down at (2.5, 2.5)
    auto ray = *Ray3::from(Vec3{2.5, 2.5, 10.0}, Vec3{0.0, 0.0, -1.0});
    auto hit = bvh.ray_cast(ray);
    REQUIRE(hit.has_value());
    CHECK_THAT(hit->point[2], WithinAbs(0.0, 1e-10));
    CHECK_THAT(hit->t, WithinAbs(10.0, 1e-10));
}

TEST_CASE("BVH ray_cast miss", "[bvh]") {
    auto tris = make_grid_triangles(5, 5);
    auto bvh = BVH<Triangle3>::build(tris);

    // Shoot ray above and parallel
    auto ray = *Ray3::from(Vec3{2.5, 2.5, 10.0}, Vec3{1.0, 0.0, 0.0});
    auto hit = bvh.ray_cast(ray);
    CHECK_FALSE(hit.has_value());
}

TEST_CASE("BVH ray_cast empty", "[bvh]") {
    auto bvh = BVH<Triangle3>::build({});
    auto ray = *Ray3::from(Vec3{0, 0, 1}, Vec3{0, 0, -1});
    CHECK_FALSE(bvh.ray_cast(ray).has_value());
}

TEST_CASE("BVH ray_cast finds closest hit", "[bvh]") {
    // Two triangles at different Z: z=0 and z=-5
    std::vector<Triangle3> tris = {
        Triangle3({-1, -1, 0}, {1, -1, 0}, {0, 1, 0}),
        Triangle3({-1, -1, -5}, {1, -1, -5}, {0, 1, -5}),
    };
    auto bvh = BVH<Triangle3>::build(tris);
    auto ray = *Ray3::from(Vec3{0, 0, 10}, Vec3{0, 0, -1});
    auto hit = bvh.ray_cast(ray);
    REQUIRE(hit.has_value());
    CHECK_THAT(hit->point[2], WithinAbs(0.0, 1e-10));  // closest
}

// ── ray_test ──────────────────────────────────────────────────

TEST_CASE("BVH ray_test hit", "[bvh]") {
    auto tris = make_grid_triangles(5, 5);
    auto bvh = BVH<Triangle3>::build(tris);
    auto ray = *Ray3::from(Vec3{2.5, 2.5, 10.0}, Vec3{0.0, 0.0, -1.0});
    CHECK(bvh.ray_test(ray));
}

TEST_CASE("BVH ray_test miss", "[bvh]") {
    auto tris = make_grid_triangles(5, 5);
    auto bvh = BVH<Triangle3>::build(tris);
    auto ray = *Ray3::from(Vec3{100.0, 100.0, 10.0}, Vec3{0.0, 0.0, -1.0});
    CHECK_FALSE(bvh.ray_test(ray));
}

// ── nearest ───────────────────────────────────────────────────

TEST_CASE("BVH nearest point", "[bvh]") {
    auto tris = make_grid_triangles(5, 5);
    auto bvh = BVH<Triangle3>::build(tris);

    // Point directly above (2.5, 2.5, 3)
    auto result = bvh.nearest(Vec3{2.5, 2.5, 3.0});
    REQUIRE(result.has_value());
    CHECK_THAT(result->distance, WithinAbs(3.0, 1e-10));
    CHECK_THAT(result->point[2], WithinAbs(0.0, 1e-10));
}

TEST_CASE("BVH nearest point on surface", "[bvh]") {
    auto tris = make_grid_triangles(5, 5);
    auto bvh = BVH<Triangle3>::build(tris);

    auto result = bvh.nearest(Vec3{2.5, 2.5, 0.0});
    REQUIRE(result.has_value());
    CHECK_THAT(result->distance, WithinAbs(0.0, 1e-10));
}

TEST_CASE("BVH nearest empty", "[bvh]") {
    auto bvh = BVH<Triangle3>::build({});
    CHECK_FALSE(bvh.nearest(Vec3{0, 0, 0}).has_value());
}

// ── query_box ─────────────────────────────────────────────────

TEST_CASE("BVH query_box all", "[bvh]") {
    auto tris = make_grid_triangles(5, 5);
    auto bvh = BVH<Triangle3>::build(tris);

    Box3 query{Vec3{-1, -1, -1}, Vec3{6, 6, 1}};
    auto result = bvh.query_box(query);
    CHECK(result.size() == 50);
}

TEST_CASE("BVH query_box partial", "[bvh]") {
    auto tris = make_grid_triangles(5, 5);
    auto bvh = BVH<Triangle3>::build(tris);

    // Query a small region — should return subset
    Box3 query{Vec3{0.5, 0.5, -0.5}, Vec3{1.5, 1.5, 0.5}};
    auto result = bvh.query_box(query);
    CHECK(result.size() > 0);
    CHECK(result.size() < 50);
}

TEST_CASE("BVH query_box miss", "[bvh]") {
    auto tris = make_grid_triangles(5, 5);
    auto bvh = BVH<Triangle3>::build(tris);

    Box3 query{Vec3{100, 100, 100}, Vec3{200, 200, 200}};
    auto result = bvh.query_box(query);
    CHECK(result.empty());
}

TEST_CASE("BVH query_box empty tree", "[bvh]") {
    auto bvh = BVH<Triangle3>::build({});
    Box3 query{Vec3{0, 0, 0}, Vec3{1, 1, 1}};
    CHECK(bvh.query_box(query).empty());
}

// ── Hit barycentric + normal ─────────────────────────────────

TEST_CASE("BVH ray_cast Triangle3 Hit barycentric + normal", "[bvh]") {
    std::vector<Triangle3> tris{Triangle3({0,0,0}, {2,0,0}, {0,2,0})};
    auto bvh = BVH<Triangle3>::build(tris);

    // Centroid ray: u ≈ 1/3, v ≈ 1/3, normal = +Z
    Ray3 ray{Vec3{2.0/3, 2.0/3, 1}, Vec3{0, 0, -1}};
    auto hit = bvh.ray_cast(ray);
    REQUIRE(hit);
    CHECK_THAT(hit->u, WithinAbs(1.0/3, 1e-9));
    CHECK_THAT(hit->v, WithinAbs(1.0/3, 1e-9));
    CHECK_THAT(hit->normal[0], WithinAbs(0.0, 1e-9));
    CHECK_THAT(hit->normal[1], WithinAbs(0.0, 1e-9));
    CHECK_THAT(hit->normal[2], WithinAbs(1.0, 1e-9));

    // Corner at vertex 1 → u=1, v=0
    Ray3 r1{Vec3{2, 0, 1}, Vec3{0, 0, -1}};
    auto h1 = bvh.ray_cast(r1);
    REQUIRE(h1);
    CHECK_THAT(h1->u, WithinAbs(1.0, 1e-9));
    CHECK_THAT(h1->v, WithinAbs(0.0, 1e-9));

    // Corner at vertex 2 → u=0, v=1
    Ray3 r2{Vec3{0, 2, 1}, Vec3{0, 0, -1}};
    auto h2 = bvh.ray_cast(r2);
    REQUIRE(h2);
    CHECK_THAT(h2->u, WithinAbs(0.0, 1e-9));
    CHECK_THAT(h2->v, WithinAbs(1.0, 1e-9));
}
