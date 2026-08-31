#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <spatium/render/camera.hpp>
#include <spatium/render/color.hpp>
#include <spatium/render/parallel_for_rows.hpp>
#include <spatium/render/sky.hpp>
#include <spatium/render/spectral.hpp>
#include <spatium/render/write_image.hpp>
#include <atomic>
#include <cmath>
#include <filesystem>
#include <vector>

using namespace spatium;
using namespace spatium::render;
using Catch::Matchers::WithinAbs;

TEST_CASE("make_camera_basis produces an orthonormal frame", "[render][camera]") {
    Camera<double> cam{.position = {3.0, 2.0, 1.0}, .target = {0.0, 0.0, 0.0}, .up = {0, 0, 1}};
    auto basis = make_camera_basis(cam);

    CHECK_THAT(basis.fwd.norm(), WithinAbs(1.0, 1e-12));
    CHECK_THAT(basis.right.norm(), WithinAbs(1.0, 1e-12));
    CHECK_THAT(basis.up.norm(), WithinAbs(1.0, 1e-12));
    CHECK_THAT(basis.fwd.dot(basis.right), WithinAbs(0.0, 1e-12));
    CHECK_THAT(basis.fwd.dot(basis.up), WithinAbs(0.0, 1e-12));
    CHECK_THAT(basis.right.dot(basis.up), WithinAbs(0.0, 1e-12));
}

TEST_CASE("camera_pixel_dir at the image center matches the camera's forward direction",
          "[render][camera]") {
    Camera<double> cam{.position = {5.0, 0.0, 0.0}, .target = {0.0, 0.0, 0.0}, .up = {0, 0, 1}};
    auto basis = make_camera_basis(cam);
    // Odd width/height so an exact center pixel exists: for width=2k+1,
    // px=k gives (px+0.5)/width == 0.5 exactly -> nx=ny=0 -> pure forward
    // direction. (An even width/height has no pixel landing exactly on
    // center -- the two nearest are each half a pixel off.)
    auto dir = camera_pixel_dir(cam, basis, 100, 50, 201, 101);
    CHECK_THAT(dir[0], WithinAbs(basis.fwd[0], 1e-9));
    CHECK_THAT(dir[1], WithinAbs(basis.fwd[1], 1e-9));
    CHECK_THAT(dir[2], WithinAbs(basis.fwd[2], 1e-9));
}

TEST_CASE("camera_pixel_dir's two known screen-coordinate conventions agree", "[render][camera]") {
    // parametric_analytical_demo.cpp used ny = (1 - 2t) * tan_half;
    // blackhole_demo.cpp/wormhole_demo.cpp used sy = -(2t - 1) * tan_fov.
    // Algebraically identical -- confirm camera_pixel_dir matches a
    // hand-rolled version of each at an off-center pixel.
    Camera<double> cam{.position = {0.0, 0.0, 10.0}, .target = {0.0, 0.0, 0.0}, .up = {0, 1, 0}};
    auto basis = make_camera_basis(cam);
    int width = 300, height = 150, px = 40, py = 110;
    double aspect = static_cast<double>(width) / height;

    double nx_a = (2.0 * (px + 0.5) / width - 1.0) * aspect * basis.tan_half;
    double ny_a = (1.0 - 2.0 * (py + 0.5) / height) * basis.tan_half;
    double ny_b = -(2.0 * (py + 0.5) / height - 1.0) * basis.tan_half;
    CHECK_THAT(ny_a, WithinAbs(ny_b, 1e-15));

    auto expected = camera_ray_dir(basis, nx_a, ny_a);
    auto actual = camera_pixel_dir(cam, basis, px, py, width, height);
    CHECK_THAT(actual[0], WithinAbs(expected[0], 1e-12));
    CHECK_THAT(actual[1], WithinAbs(expected[1], 1e-12));
    CHECK_THAT(actual[2], WithinAbs(expected[2], 1e-12));
}

TEST_CASE("parallel_for_rows visits every row exactly once", "[render][parallel]") {
    constexpr int kHeight = 733;  // odd, deliberately not a multiple of typical core counts
    std::vector<std::atomic<int>> visits(kHeight);
    for (auto& v : visits) v.store(0);

    parallel_for_rows(kHeight, [&](int y) { visits[static_cast<std::size_t>(y)].fetch_add(1); });

    for (int y = 0; y < kHeight; ++y) CHECK(visits[static_cast<std::size_t>(y)].load() == 1);
}

TEST_CASE("parallel_for_rows handles height smaller than hardware_concurrency",
          "[render][parallel]") {
    std::atomic<int> total{0};
    parallel_for_rows(1, [&](int y) {
        CHECK(y == 0);
        total.fetch_add(1);
    });
    CHECK(total.load() == 1);
}

TEST_CASE("write_png_rgb writes a real, non-empty file", "[render][write_image]") {
    auto path = std::filesystem::temp_directory_path() / "spatium_test_render.png";
    std::filesystem::remove(path);

    constexpr int w = 4, h = 4;
    std::vector<std::uint8_t> rgb(static_cast<std::size_t>(w) * h * 3, 128);
    REQUIRE(write_png_rgb(path, w, h, rgb));

    REQUIRE(std::filesystem::exists(path));
    CHECK(std::filesystem::file_size(path) > 0);

    std::filesystem::remove(path);
}

TEST_CASE("blackbody_to_rgb255 reproduces warm-to-cool ordering", "[render][spectral]") {
    auto warm = blackbody_to_rgb255(1500.0);   // ember/candlelight: strongly red, no blue
    auto neutral = blackbody_to_rgb255(6600.0); // the fit's own red/blue crossover
    auto cool = blackbody_to_rgb255(20000.0);  // O-type star: blue-white

    CHECK_THAT(warm[0], WithinAbs(255.0, 1e-9));
    CHECK(warm[2] == 0.0);                     // temp/100=15 <= 19 branch
    CHECK(warm[0] > warm[2]);

    CHECK_THAT(neutral[0], WithinAbs(255.0, 1.0));
    CHECK_THAT(neutral[2], WithinAbs(255.0, 1.0));

    CHECK_THAT(cool[2], WithinAbs(255.0, 1e-9));
    CHECK(cool[0] < cool[2]);
    CHECK(cool[0] < warm[0]);  // hotter -> less red, the physically-expected direction

    // Clamping holds outside the fit's nominal [1000K, 40000K] domain.
    for (double t : {200.0, 100000.0}) {
        auto c = blackbody_to_rgb255(t);
        for (int i = 0; i < 3; ++i) CHECK((c[i] >= 0.0 && c[i] <= 255.0));
    }
}

TEST_CASE("hsv_to_rgb255 matches known primary/secondary hues", "[render][color]") {
    auto red = hsv_to_rgb255(0.0, 1.0, 1.0);
    CHECK_THAT(red[0], WithinAbs(255.0, 1e-9));
    CHECK_THAT(red[1], WithinAbs(0.0, 1e-9));
    CHECK_THAT(red[2], WithinAbs(0.0, 1e-9));

    auto green = hsv_to_rgb255(1.0 / 3.0, 1.0, 1.0);
    CHECK_THAT(green[0], WithinAbs(0.0, 1e-9));
    CHECK_THAT(green[1], WithinAbs(255.0, 1e-9));
    CHECK_THAT(green[2], WithinAbs(0.0, 1e-9));

    auto white = hsv_to_rgb255(0.5, 0.0, 1.0);
    CHECK_THAT(white[0], WithinAbs(255.0, 1e-9));
    CHECK_THAT(white[1], WithinAbs(255.0, 1e-9));
    CHECK_THAT(white[2], WithinAbs(255.0, 1e-9));
}

TEST_CASE("sample_sky_color's wide_sky gradient is flat at wide_sky=false", "[render][sky]") {
    // Two directions at very different latitudes (theta) would differ
    // sharply under the wide_sky gradient -- with it off, only the fixed
    // tint (plus whatever star/cloud/spiral falls exactly there) should
    // show, so two star/structure-free directions must match exactly.
    Sky sky = make_starfield(0, /*seed=*/1, /*tint=*/{10.0, 10.0, 10.0}, /*wide_sky=*/false);
    sky.spirals.clear();
    sky.clouds.clear();

    auto pole = sample_sky_color(sky, Vec<double, 3>{0.0, 1.0, 0.0});
    auto equator = sample_sky_color(sky, Vec<double, 3>{1.0, 0.0, 0.0});
    CHECK_THAT(pole[0], WithinAbs(equator[0], 1e-12));
    CHECK_THAT(pole[1], WithinAbs(equator[1], 1e-12));
    CHECK_THAT(pole[2], WithinAbs(equator[2], 1e-12));
    CHECK_THAT(pole[0], WithinAbs(10.0, 1e-12));
}

TEST_CASE("sample_sky_color's wide_sky gradient varies with latitude when true", "[render][sky]") {
    Sky sky = make_starfield(0, /*seed=*/1, /*tint=*/{10.0, 10.0, 10.0}, /*wide_sky=*/true);
    sky.spirals.clear();
    sky.clouds.clear();

    auto pole = sample_sky_color(sky, Vec<double, 3>{0.0, 1.0, 0.0});
    auto equator = sample_sky_color(sky, Vec<double, 3>{1.0, 0.0, 0.0});
    CHECK(equator[2] > pole[2]); // brighter near mid-latitude than at the pole
}

TEST_CASE("make_starfield places every star findable via sample_sky_color", "[render][sky]") {
    Sky sky = make_starfield(500, /*seed=*/3);
    REQUIRE(sky.stars.size() == 500);
    for (auto& s : sky.stars) {
        auto c = sample_sky_color(sky, s.dir);
        // A star's own exact direction always falls within its own
        // angular radius, so it must be at least as bright as its
        // brightness*255 floor (clouds/spirals only add, never subtract).
        double expected = s.brightness * 255.0;
        CHECK(c[0] >= expected - 1e-6);
    }
}
