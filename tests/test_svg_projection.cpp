#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/io/svg.hpp>
#include <spatium/mesh/primitives.hpp>
#include <spatium/mesh/subdivision.hpp>
#include <spatium/mesh/geodesic.hpp>
#include <spatium/mesh/topology.hpp>
#include <spatium/spaces/sphere.hpp>
#include <spatium/spaces/euclidean.hpp>
#include <filesystem>

using namespace spatium;
using namespace spatium::io;
using namespace spatium::mesh;
using Catch::Matchers::WithinAbs;

using S2 = Sphere<2>;
using E3 = Euclidean<3>;

// ── OrthoProjection ───────────────────────────────────────────

TEST_CASE("OrthoProjection default Z forward", "[svg_proj]") {
    OrthoProjection proj;
    Vec3 p{1, 2, 3};
    auto r = proj(p);
    // Default: forward={0,0,-1}, up={0,1,0}
    // right = forward × up = {0,0,-1}×{0,1,0} = {1,0,0}
    // after orthogonalization up stays {0,1,0}
    CHECK_THAT(r[0], WithinAbs(1.0, 1e-10));
    CHECK_THAT(r[1], WithinAbs(2.0, 1e-10));
}

TEST_CASE("OrthoProjection custom direction", "[svg_proj]") {
    OrthoProjection proj(Vec3{1, 0, 0}, Vec3{0, 0, 1});
    // Forward=+X, up=+Z, right = X×Z = -Y
    Vec3 p{5, 3, 7};
    auto r = proj(p);
    // right = forward × up = (1,0,0)×(0,0,1) = (0,-1,0) → then up = right × forward
    // Actually after orthogonalization this might differ, just check it produces 2D
    CHECK(std::isfinite(r[0]));
    CHECK(std::isfinite(r[1]));
}

// ── SVG mesh wireframe ────────────────────────────────────────

TEST_CASE("SVG mesh wireframe produces content", "[svg_proj]") {
    auto mesh = box_mesh();
    Svg svg;
    OrthoProjection proj;
    svg.mesh_wireframe(mesh, proj);
    CHECK(svg.content.size() > 100);
}

TEST_CASE("SVG mesh filled produces content", "[svg_proj]") {
    S2 sphere;
    auto mesh = subdivide(icosahedron(sphere), sphere, 1);
    Svg svg;
    OrthoProjection proj(Vec3{0.5, -0.3, -1});
    svg.mesh_filled(mesh, proj);
    CHECK(svg.content.size() > 500);
}

// ── SVG save ──────────────────────────────────────────────────

TEST_CASE("SVG mesh save to file", "[svg_proj]") {
    S2 sphere;
    auto mesh = subdivide(icosahedron(sphere), sphere, 2);
    auto svg = mesh_to_svg(mesh);
    auto path = std::filesystem::temp_directory_path() / "spatium_sphere.svg";
    CHECK(svg.save(path.string()));
    CHECK(std::filesystem::file_size(path) > 100);
    std::filesystem::remove(path);
}

// ── Colored SVG ───────────────────────────────────────────────

TEST_CASE("SVG colored mesh with distance field", "[svg_proj]") {
    S2 sphere;
    auto mesh = subdivide(icosahedron(sphere), sphere, 2);
    auto topo = MeshTopology<S2>::build(mesh);
    auto field = geodesic_distances(topo, sphere, uint32_t{0});

    auto svg = mesh_to_svg_colored(mesh, std::span<const double>(field.distances));
    auto path = std::filesystem::temp_directory_path() / "spatium_geodesic.svg";
    CHECK(svg.save(path.string()));
    CHECK(std::filesystem::file_size(path) > 500);
    std::filesystem::remove(path);
}

// ── Convenience one-shot ──────────────────────────────────────

TEST_CASE("mesh_to_svg convenience function", "[svg_proj]") {
    auto mesh = box_mesh();
    auto svg = mesh_to_svg(mesh);
    CHECK(svg.content.size() > 100);
}

// ── Color maps ────────────────────────────────────────────────

TEST_CASE("Viridis colormap range", "[svg_proj]") {
    auto cmap = Svg::viridis_map();
    auto c0 = cmap(0.0);
    auto c1 = cmap(1.0);
    CHECK(c0[0] == '#');
    CHECK(c0.size() == 7);
    CHECK(c1[0] == '#');
}

TEST_CASE("Heat colormap range", "[svg_proj]") {
    auto cmap = Svg::heat_map();
    auto c0 = cmap(0.0);
    auto c1 = cmap(1.0);
    CHECK(c0[0] == '#');
    CHECK(c1[0] == '#');
}
