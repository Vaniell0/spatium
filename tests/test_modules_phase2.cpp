// Phase 2 verification — consume spatium.spaces, spatium.point and
// spatium.mesh purely through the modular API.
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <numbers>

import spatium.core;
import spatium.algebra;
import spatium.mesh;
import spatium.spaces;
import spatium.point;

using spatium::Vec3;
using spatium::Euclidean;
using spatium::Sphere;
using spatium::E3;
using spatium::Point;
using spatium::Morphism;
using spatium::morph;
using spatium::pt;
using spatium::mesh::Mesh;

TEST_CASE("module spatium.spaces: Euclidean<3> distance", "[modules][phase2][spaces]") {
    E3 e;
    auto d = e.distance(Vec3{0.0, 0.0, 0.0}, Vec3{3.0, 4.0, 0.0});
    REQUIRE(d == 5.0);
}

TEST_CASE("module spatium.spaces: Sphere<2> great-circle distance", "[modules][phase2][spaces]") {
    Sphere<2> s;
    auto north = Vec3{0.0, 0.0, 1.0};
    auto equator = Vec3{1.0, 0.0, 0.0};
    auto d = s.distance(north, equator);
    REQUIRE(std::abs(d - std::numbers::pi / 2.0) < 1e-12);
}

TEST_CASE("module spatium.point: Point + Morphism pipe", "[modules][phase2][point]") {
    auto shift = morph<E3, E3>([](const Vec3& v) { return Vec3{v[0] + 1.0, v[1], v[2]}; });
    auto scale = morph<E3, E3>([](const Vec3& v) { return Vec3{v[0] * 2.0, v[1] * 2.0, v[2] * 2.0}; });
    auto p = pt<E3>(Vec3{1.0, 2.0, 3.0});
    auto q = p | shift | scale;
    REQUIRE(q.raw()[0] == 4.0);
    REQUIRE(q.raw()[1] == 4.0);
    REQUIRE(q.raw()[2] == 6.0);
}

TEST_CASE("module spatium.mesh: Mesh<E3> centroid + bounding box", "[modules][phase2][mesh]") {
    Mesh<E3> m;
    m.vertices = {Vec3{0.0, 0.0, 0.0}, Vec3{2.0, 0.0, 0.0}, Vec3{0.0, 2.0, 0.0}};
    m.faces = {{0, 1, 2}};
    auto c = m.centroid();
    REQUIRE(std::abs(c[0] - (2.0 / 3.0)) < 1e-12);
    REQUIRE(std::abs(c[1] - (2.0 / 3.0)) < 1e-12);
    auto [lo, hi] = m.bounding_box();
    REQUIRE(lo[0] == 0.0);
    REQUIRE(hi[0] == 2.0);
    REQUIRE(hi[1] == 2.0);
}
