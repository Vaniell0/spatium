// spatium.render verification — consume it through `import`. write_image.hpp
// stays header-only (see render.cppm) so it isn't exercised here.
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <vector>

import spatium.core;
import spatium.algebra;
import spatium.render;

using spatium::Vec3;
using spatium::render::Camera;
using spatium::render::make_camera_basis;
using spatium::render::hsv_to_rgb255;
using spatium::render::make_starfield;
using spatium::render::sample_sky_color;

TEST_CASE("module spatium.render: camera basis is orthonormal", "[modules][render]") {
    Camera<double> cam{.position = {3.0, 2.0, 1.0}, .target = {0.0, 0.0, 0.0}, .up = {0, 0, 1}};
    auto basis = make_camera_basis(cam);
    REQUIRE(std::abs(basis.fwd.norm() - 1.0) < 1e-12);
    REQUIRE(std::abs(basis.fwd.dot(basis.right)) < 1e-12);
}

TEST_CASE("module spatium.render: hsv_to_rgb255 red hue", "[modules][render]") {
    auto red = hsv_to_rgb255(0.0, 1.0, 1.0);
    REQUIRE(std::abs(red[0] - 255.0) < 1e-9);
    REQUIRE(std::abs(red[1]) < 1e-9);
}

TEST_CASE("module spatium.render: sky sampling finds a placed star", "[modules][render]") {
    auto sky = make_starfield(10, /*seed=*/1);
    REQUIRE(sky.stars.size() == 10);
    auto c = sample_sky_color(sky, sky.stars[0].dir);
    REQUIRE(c[0] >= sky.stars[0].brightness * 255.0 - 1e-6);
}
