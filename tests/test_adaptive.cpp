#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/mesh/subdivision.hpp>
#include <spatium/mesh/quality.hpp>
#include <spatium/mesh/primitives.hpp>
#include <spatium/mesh/topology.hpp>
#include <spatium/spaces/sphere.hpp>
#include <spatium/spaces/euclidean.hpp>
#include <cmath>
#include <numbers>

using namespace spatium;
using namespace spatium::mesh;
using Catch::Matchers::WithinAbs;

using S2 = Sphere<2>;
using E3 = Euclidean<3>;

// ── Adaptive subdivision ──────────────────────────────────────

TEST_CASE("Adaptive subdivide flat grid → no refinement", "[adaptive]") {
    auto mesh = grid_mesh(4, 4);
    E3 space;
    // Very small threshold but flat surface → normals all parallel → no splits
    auto refined = subdivide_adaptive(mesh, space, 0.01, 3);
    CHECK(refined.faces.size() == mesh.faces.size());
}

TEST_CASE("Adaptive subdivide icosphere → refines all", "[adaptive]") {
    S2 sphere;
    auto mesh = icosahedron(sphere);
    // Low threshold → all faces refined
    auto refined = subdivide_adaptive(mesh, sphere, 0.1, 1);
    CHECK(refined.faces.size() > mesh.faces.size());
}

TEST_CASE("Adaptive subdivide high threshold → no refinement", "[adaptive]") {
    S2 sphere;
    auto mesh = icosahedron(sphere);
    // Very high threshold → no faces refined
    auto refined = subdivide_adaptive(mesh, sphere, 10.0, 1);
    CHECK(refined.faces.size() == mesh.faces.size());
}

TEST_CASE("Adaptive subdivide fewer faces than uniform", "[adaptive]") {
    S2 sphere;
    auto mesh = subdivide(icosahedron(sphere), sphere, 1);  // 80 faces, mixed curvature
    // Threshold 0.5 rad (~29 deg) — only highly curved faces subdivide
    auto adaptive = subdivide_adaptive(mesh, sphere, 0.5, 1);
    auto uniform = subdivide_once(mesh, sphere);  // always 320 faces
    // Adaptive should produce <= uniform (some faces skip)
    CHECK(adaptive.faces.size() <= uniform.faces.size());
}

TEST_CASE("Adaptive subdivide valid mesh", "[adaptive]") {
    S2 sphere;
    auto mesh = icosahedron(sphere);
    auto refined = subdivide_adaptive(mesh, sphere, 0.2, 2);
    // All face indices valid
    for (auto& [a, b, c] : refined.faces) {
        CHECK(a < refined.vertices.size());
        CHECK(b < refined.vertices.size());
        CHECK(c < refined.vertices.size());
    }
}

TEST_CASE("Adaptive subdivide vertices on surface", "[adaptive]") {
    S2 sphere;
    auto mesh = icosahedron(sphere);
    auto refined = subdivide_adaptive(mesh, sphere, 0.2, 2);
    for (auto& v : refined.vertices)
        CHECK_THAT(v.norm(), WithinAbs(1.0, 1e-10));
}

// ── Mesh quality ──────────────────────────────────────────────

TEST_CASE("Mesh quality icosahedron", "[quality]") {
    S2 sphere;
    auto mesh = icosahedron(sphere);
    auto q = mesh_quality(mesh);
    CHECK(q.min_angle > 0.1);
    CHECK(q.max_angle < std::numbers::pi);
    CHECK(q.avg_aspect_ratio > 0.9);
    CHECK(q.degenerate_count == 0);
}

TEST_CASE("Mesh quality grid", "[quality]") {
    auto mesh = grid_mesh(4, 4);
    auto q = mesh_quality(mesh);
    CHECK(q.degenerate_count == 0);
    CHECK(q.avg_aspect_ratio < 2.0);  // right triangles from grid
}

TEST_CASE("Face aspect ratio equilateral", "[quality]") {
    S2 sphere;
    auto mesh = icosahedron(sphere);
    for (std::size_t fi = 0; fi < mesh.faces.size(); ++fi) {
        auto ar = face_aspect_ratio<S2>(mesh, fi);
        CHECK(ar > 0.9);
        CHECK(ar < 1.5);  // icosahedron faces are nearly equilateral
    }
}
