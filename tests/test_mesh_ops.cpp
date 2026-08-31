#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/mesh/operations.hpp>
#include <spatium/mesh/primitives.hpp>
#include <spatium/mesh/subdivision.hpp>
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
