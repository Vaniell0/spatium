// Phase 3 verification — consume spatium.geometry through `import`.
#include <catch2/catch_test_macros.hpp>
#include <cmath>

import spatium.core;
import spatium.algebra;
import spatium.spaces;
import spatium.point;
import spatium.geometry;

using spatium::Vec2;
using spatium::Vec3;
using spatium::geometry::Box;
using spatium::geometry::Triangle;
using spatium::geometry::Line;
using spatium::geometry::Hyperplane;
using spatium::geometry::tri;
using spatium::geometry::box;
using spatium::geometry::line;
using spatium::geometry::Quadric;
using spatium::geometry::ray_quadric;

TEST_CASE("module spatium.geometry: Triangle area + centroid", "[modules][phase3][geometry]") {
    auto t = tri(Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0}, Vec3{0.0, 1.0, 0.0});
    REQUIRE(std::abs(t.area() - 0.5) < 1e-12);
    auto c = t.centroid();
    REQUIRE(std::abs(c[0] - (1.0 / 3.0)) < 1e-12);
    REQUIRE(std::abs(c[1] - (1.0 / 3.0)) < 1e-12);
}

TEST_CASE("module spatium.geometry: Box contains + diagonal", "[modules][phase3][geometry]") {
    auto b = box(Vec3{0.0, 0.0, 0.0}, Vec3{2.0, 4.0, 6.0});
    REQUIRE(b.contains(Vec3{1.0, 2.0, 3.0}));
    REQUIRE_FALSE(b.contains(Vec3{3.0, 0.0, 0.0}));
}

TEST_CASE("module spatium.geometry: ray-quadric sphere intersection", "[modules][phase3][geometry]") {
    auto sphere = Quadric<double>::sphere(1.0);
    spatium::geometry::Ray<3, double> r{Vec3{-3.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0}};
    auto hits = ray_quadric(r, sphere);
    REQUIRE(hits.size() == 2);
    REQUIRE(std::abs(hits[0].t - 2.0) < 1e-9);
    REQUIRE(std::abs(hits[1].t - 4.0) < 1e-9);
}
