#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/spaces/implicit.hpp>
#include <spatium/mesh/primitives.hpp>
#include <spatium/io/stl.hpp>
#include <cmath>
#include <filesystem>

using namespace spatium;
using Catch::Matchers::WithinAbs;

TEST_CASE("ImplicitSurface sphere evaluate", "[implicit]") {
    auto sphere = make_implicit_sphere();
    CHECK_THAT(sphere(1.0, 0.0, 0.0), WithinAbs(0.0, 1e-12));
    CHECK(sphere(2.0, 0.0, 0.0) > 0);
    CHECK(sphere(0.5, 0.0, 0.0) < 0);
}

TEST_CASE("ImplicitSurface sphere project", "[implicit]") {
    auto sphere = make_implicit_sphere();
    Vec3 p{2.0, 0.0, 0.0};
    auto proj = sphere.project(p);
    CHECK_THAT(proj.norm(), WithinAbs(1.0, 1e-8));
}

TEST_CASE("ImplicitSurface sphere normal", "[implicit]") {
    auto sphere = make_implicit_sphere();
    Vec3 p{1.0, 0.0, 0.0};
    auto n = sphere.normal(p);
    CHECK_THAT(n[0], WithinAbs(1.0, 1e-6));
    CHECK_THAT(n.norm(), WithinAbs(1.0, 1e-6));
}

TEST_CASE("ImplicitSurface sphere contains", "[implicit]") {
    auto sphere = make_implicit_sphere();
    CHECK(sphere.contains(Vec3{1.0, 0.0, 0.0}));
    CHECK(sphere.contains(Vec3{0.0, 1.0, 0.0}));
    CHECK_FALSE(sphere.contains(Vec3{2.0, 0.0, 0.0}));
}

TEST_CASE("Marching cubes sphere produces mesh", "[implicit]") {
    auto sphere = make_implicit_sphere();
    auto mesh = marching_cubes(sphere, 10);
    CHECK(mesh.vertices.size() > 0);
    CHECK(mesh.faces.size() > 0);
}

TEST_CASE("Marching cubes sphere vertices near surface", "[implicit]") {
    auto sphere = make_implicit_sphere();
    auto mesh = marching_cubes(sphere, 16);
    for (auto& v : mesh.vertices)
        CHECK_THAT(v.norm(), WithinAbs(1.0, 0.2));  // within grid cell size
}

TEST_CASE("Marching cubes torus produces mesh", "[implicit]") {
    auto torus = make_implicit_torus();
    auto mesh = marching_cubes(torus, 16);
    CHECK(mesh.vertices.size() > 100);
    CHECK(mesh.faces.size() > 100);
}

TEST_CASE("ImplicitSurface exp then project", "[implicit]") {
    auto sphere = make_implicit_sphere();
    Vec3 p{1.0, 0.0, 0.0};
    Vec3 v{0.0, 0.1, 0.0};
    auto q = sphere.exp_map(p, v, 1.0);
    CHECK_THAT(q.norm(), WithinAbs(1.0, 1e-6));
}

TEST_CASE("STL round-trip", "[stl]") {
    auto mesh = mesh::box_mesh();
    auto path = std::filesystem::temp_directory_path() / "spatium_test.stl";

    auto save_r = io::save_stl(mesh, path);
    REQUIRE(save_r.has_value());

    auto load_r = io::load_stl(path);
    REQUIRE(load_r.has_value());

    // STL duplicates vertices per triangle
    CHECK(load_r->faces.size() == mesh.faces.size());
    CHECK(load_r->vertices.size() == mesh.faces.size() * 3);

    std::filesystem::remove(path);
}
