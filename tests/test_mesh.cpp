#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/mesh/primitives.hpp>
#include <spatium/mesh/subdivision.hpp>
#include <spatium/mesh/lod.hpp>
#include <numbers>

using namespace spatium;
using namespace spatium::mesh;
using Catch::Matchers::WithinAbs;

TEST_CASE("Icosahedron construction", "[mesh]") {
    S2 sphere;
    auto m = icosahedron(sphere);
    CHECK(m.vertex_count() == 12);
    CHECK(m.face_count() == 20);
}

TEST_CASE("Icosahedron vertices on sphere", "[mesh]") {
    S2 sphere;
    auto m = icosahedron(sphere);
    for (const auto& v : m.vertices) {
        CHECK_THAT(v.norm(), WithinAbs(1.0, 1e-10));
    }
}

TEST_CASE("Tetrahedron construction", "[mesh]") {
    S2 sphere;
    auto m = tetrahedron(sphere);
    CHECK(m.vertex_count() == 4);
    CHECK(m.face_count() == 4);
}

TEST_CASE("Subdivide once", "[mesh]") {
    S2 sphere;
    auto base = icosahedron(sphere);
    auto sub = subdivide_once(base, sphere);

    // 1 subdivision: 20*4=80 faces, ~42 vertices
    CHECK(sub.face_count() == 80);
    CHECK(sub.vertex_count() == 42); // 12 + 30 midpoints

    // All vertices should be on sphere
    for (const auto& v : sub.vertices) {
        CHECK_THAT(v.norm(), WithinAbs(1.0, 1e-10));
    }
}

TEST_CASE("Subdivide twice", "[mesh]") {
    S2 sphere;
    auto base = icosahedron(sphere);
    auto sub2 = subdivide(base, sphere, 2);
    CHECK(sub2.face_count() == 320);
}

TEST_CASE("Subdivision improves sphere approximation", "[mesh]") {
    S2 sphere;
    auto base = icosahedron(sphere);

    // True surface area of unit sphere = 4π
    auto true_area = 4.0 * std::numbers::pi;

    auto area0 = base.area(sphere);
    auto sub1 = subdivide_once(base, sphere);
    auto area1 = sub1.area(sphere);
    auto sub2 = subdivide_once(sub1, sphere);
    auto area2 = sub2.area(sphere);

    // Each level should be closer to true area
    CHECK(std::abs(area1 - true_area) < std::abs(area0 - true_area));
    CHECK(std::abs(area2 - true_area) < std::abs(area1 - true_area));
    // Level 2 should be within 1% of true area
    CHECK_THAT(area2, WithinAbs(true_area, true_area * 0.02));
}

TEST_CASE("LodChain", "[mesh]") {
    S2 sphere;
    auto chain = LodChain<S2>::build(icosahedron(sphere), sphere, 3);

    CHECK(chain.level_count() == 4); // 0, 1, 2, 3
    CHECK(chain.coarsest().face_count() == 20);
    CHECK(chain.at(1).face_count() == 80);
    CHECK(chain.at(2).face_count() == 320);
    CHECK(chain.finest().face_count() == 1280);
}

TEST_CASE("Mesh format", "[mesh]") {
    S2 sphere;
    auto m = icosahedron(sphere);
    auto s = std::format("{}", m);
    CHECK(s.find("V=12") != std::string::npos);
    CHECK(s.find("F=20") != std::string::npos);
}

TEST_CASE("Sphere with custom radius mesh", "[mesh]") {
    Sphere<2> sphere{.radius = 5.0};
    auto m = icosahedron(sphere);
    for (const auto& v : m.vertices)
        CHECK_THAT(v.norm(), WithinAbs(5.0, 1e-8));
}
