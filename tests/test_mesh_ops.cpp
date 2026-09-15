#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/mesh/conform.hpp>
#include <spatium/mesh/operations.hpp>
#include <spatium/mesh/primitives.hpp>
#include <spatium/mesh/scatter.hpp>
#include <spatium/mesh/subdivision.hpp>
#include <spatium/spaces/parametric.hpp>
#include <spatium/spaces/sphere.hpp>
#include <spatium/spaces/euclidean.hpp>
#include <spatium/io/obj.hpp>
#include <filesystem>

using namespace spatium;
using namespace spatium::mesh;
using Catch::Matchers::WithinAbs;

using S2 = Sphere<2>;
using E3 = Euclidean<3>;

// ── Merge ─────────────────────────────────────────────────────

TEST_CASE("Mesh merge combines vertices and faces", "[mesh_ops]") {
    S2 sphere;
    auto a = icosahedron(sphere);
    auto b = tetrahedron(sphere);
    auto merged = merge(a, b);
    CHECK(merged.vertices.size() == a.vertices.size() + b.vertices.size());
    CHECK(merged.faces.size() == a.faces.size() + b.faces.size());
}

TEST_CASE("Mesh merge face indices valid", "[mesh_ops]") {
    S2 sphere;
    auto merged = merge(icosahedron(sphere), tetrahedron(sphere));
    for (auto& [a, b, c] : merged.faces) {
        CHECK(a < merged.vertices.size());
        CHECK(b < merged.vertices.size());
        CHECK(c < merged.vertices.size());
    }
}

// ── Flip normals ──────────────────────────────────────────────

TEST_CASE("Mesh flip normals reverses winding", "[mesh_ops]") {
    S2 sphere;
    auto mesh = icosahedron(sphere);
    auto flipped = flip_normals(mesh);
    CHECK(flipped.faces.size() == mesh.faces.size());
    for (std::size_t i = 0; i < mesh.faces.size(); ++i) {
        CHECK(flipped.faces[i][0] == mesh.faces[i][0]);
        CHECK(flipped.faces[i][1] == mesh.faces[i][2]);
        CHECK(flipped.faces[i][2] == mesh.faces[i][1]);
    }
}

// ── Transform ─────────────────────────────────────────────────

TEST_CASE("Mesh transform with function", "[mesh_ops]") {
    auto mesh = grid_mesh(2, 2);
    auto scaled = transform<E3>(mesh, [](const Vec3& v) -> Vec3 {
        return Vec3{v * 2.0};
    });
    CHECK_THAT(scaled.vertices[0][0], WithinAbs(0.0, 1e-15));
    // Corner (1,1,0) → (2,2,0)
    CHECK_THAT(scaled.vertices.back()[0], WithinAbs(2.0, 1e-10));
    CHECK_THAT(scaled.vertices.back()[1], WithinAbs(2.0, 1e-10));
}

// ── Normals ───────────────────────────────────────────────────

TEST_CASE("Mesh face normals unit length", "[mesh_ops]") {
    S2 sphere;
    auto mesh = icosahedron(sphere);
    auto normals = compute_face_normals(mesh);
    CHECK(normals.size() == mesh.faces.size());
    for (auto& n : normals)
        CHECK_THAT(n.norm(), WithinAbs(1.0, 1e-10));
}

TEST_CASE("Mesh vertex normals unit length", "[mesh_ops]") {
    S2 sphere;
    auto mesh = icosahedron(sphere);
    auto normals = compute_vertex_normals(mesh);
    CHECK(normals.size() == mesh.vertices.size());
    for (auto& n : normals)
        CHECK_THAT(n.norm(), WithinAbs(1.0, 1e-10));
}

// ── Centered ──────────────────────────────────────────────────

TEST_CASE("Mesh centered at origin", "[mesh_ops]") {
    auto mesh = grid_mesh(4, 4, 2.0, 2.0);
    auto c = centered(mesh);
    Vec3 sum{};
    for (auto& v : c.vertices) sum = Vec3{sum + v};
    auto centroid = Vec3{sum / static_cast<double>(c.vertices.size())};
    CHECK_THAT(centroid.norm(), WithinAbs(0.0, 1e-10));
}

// ── Grid mesh ─────────────────────────────────────────────────

TEST_CASE("Grid mesh vertex and face count", "[mesh_ops]") {
    auto mesh = grid_mesh(4, 3);
    CHECK(mesh.vertices.size() == 20);  // 5*4
    CHECK(mesh.faces.size() == 24);     // 4*3*2
}

// ── UV sphere ─────────────────────────────────────────────────

TEST_CASE("UV sphere mesh vertex count", "[mesh_ops]") {
    auto mesh = uv_sphere_mesh(8, 4);
    // 2 poles + (4-1)*8 = 2 + 24 = 26
    CHECK(mesh.vertices.size() == 26);
}

TEST_CASE("UV sphere mesh closed", "[mesh_ops]") {
    auto mesh = uv_sphere_mesh(12, 6);
    // Euler: V-E+F=2 for closed mesh
    // V = 2 + (6-1)*12 = 62
    // F = 12 + (6-2)*12*2 + 12 = 12+96+12 = 120
    // E = V+F-2 = 180
    CHECK(mesh.vertices.size() == 62);
    CHECK(mesh.faces.size() == 120);
}

// ── Box mesh ──────────────────────────────────────────────────

TEST_CASE("Box mesh 8 vertices 12 faces", "[mesh_ops]") {
    auto mesh = box_mesh();
    CHECK(mesh.vertices.size() == 8);
    CHECK(mesh.faces.size() == 12);
}

// ── OBJ round-trip ────────────────────────────────────────────

TEST_CASE("OBJ save and load round-trip", "[obj]") {
    auto mesh = box_mesh();
    auto path = std::filesystem::temp_directory_path() / "spatium_test.obj";

    auto save_result = io::save_obj(mesh, path);
    REQUIRE(save_result.has_value());

    auto load_result = io::load_obj(path);
    REQUIRE(load_result.has_value());

    auto& loaded = *load_result;
    CHECK(loaded.vertices.size() == mesh.vertices.size());
    CHECK(loaded.faces.size() == mesh.faces.size());

    std::filesystem::remove(path);
}

TEST_CASE("OBJ load nonexistent file fails", "[obj]") {
    auto result = io::load_obj("/nonexistent/path.obj");
    CHECK_FALSE(result.has_value());
}

// ── mesh/conform.hpp and mesh/scatter.hpp ────────────────────────
//
// These two had no caller anywhere in the tree -- not a demo, not a
// test, not even the umbrella header, which does not include them. So
// their bodies are templates that were never instantiated: they parsed,
// and nothing more was ever checked. The roadmap claimed they "pass
// their own tests", which was simply untrue; there were none.
//
// Kept rather than deleted (they are a real alternative to the analytic
// placement path in spaces/sample.hpp, for targets that are a mesh
// rather than a ParametricSurface), so the cheap half of "keep" is owed:
// instantiate them once and check the result is the shape it claims.
// This is not a demo and is not meant to grow into one.

TEST_CASE("conform_to_surface drapes a guide mesh onto a target", "[conform]") {
    auto torus = make_torus<double>(2.0, 1.0);
    auto guide = box_mesh<double>(Vec<double, 3>{0.2, 0.2, 0.2});

    auto draped = conform_to_surface(guide, torus, 0.05);

    // Topology carried over unchanged -- it drapes, it does not remesh.
    CHECK(draped.vertex_count() == guide.vertex_count());
    CHECK(draped.face_count() == guide.face_count());

    // Every vertex ends up one thickness off the target, measured
    // against the target's own projection rather than against the
    // guide's starting position.
    for (const auto& v : draped.vertices) {
        auto on_surface = torus.project(v);
        CHECK_THAT((v - on_surface).norm(), WithinAbs(0.05, 1e-6));
    }
}

TEST_CASE("scatter_on_surface spreads points by geodesic distance", "[scatter]") {
    using S = Sphere<2, double>;
    S sphere{1.0};

    auto euclidean = uv_sphere_mesh<double>(24, 12, 1.0);
    Mesh<S> m;
    m.vertices.assign(euclidean.vertices.begin(), euclidean.vertices.end());
    m.faces = euclidean.faces;

    auto pts = scatter_on_surface(m, sphere, std::size_t{20});
    REQUIRE(pts.size() == 20);

    for (const auto& p : pts) {
        CHECK_THAT(p.position.norm(), WithinAbs(1.0, 1e-9));  // on the sphere
        CHECK_THAT(p.normal.norm(), WithinAbs(1.0, 1e-9));    // oriented, ready to place an object
    }

    // Farthest-point sampling, so no two sites coincide -- the property
    // that distinguishes this from independent random placement, which
    // clusters and gaps.
    for (std::size_t i = 0; i < pts.size(); ++i)
        for (std::size_t j = i + 1; j < pts.size(); ++j)
            CHECK((pts[i].position - pts[j].position).norm() > 1e-9);
}
