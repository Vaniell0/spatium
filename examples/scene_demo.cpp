// Scene description demo — loads a JSON scene file (io/scene.hpp) and
// renders it with the shared CPU raytracer engine (render::Camera,
// parallel_for_rows, supersample_pixel, write_png_rgb) by calling each
// object's resolved shape's ray_hits() closure -- geometry::ray_quadric/
// ray_torus/the ray-box slab test under the hood, see io/scene.hpp's
// built-in factories.
//
// The render loop below never branches on shape kind: every object is
// resolve_shape()'d once up front and rendered through the exact same
// "call ray_hits(), keep the nearest t" loop, whether the shape came
// from a built-in factory (sphere/box/torus) or one a caller registered
// itself. That's the point of the registry in io/scene.hpp -- adding a
// new shape kind to a scene never touches this file.
//
// Run from the repository's examples/ directory (or pass --scene) so
// the default --scene path resolves to the sample scene shipped
// alongside this file.
//
// The shipped sample scene sticks to sphere/box objects -- "torus" is
// registered and works (see io/scene.hpp's own tests), but a pass
// building this demo found geometry::ray_torus()'s quartic solver
// unreliable for anything but a narrow, exactly-symmetric viewing
// angle, so a torus placed for a normal camera shot renders sparse or
// invisible. See io/scene.hpp's make_torus() comment for the repro; a
// real fix belongs in solve_quartic()/ray_torus(), out of scope here.

// Only `#define` here -- render/write_image.hpp is the one place that
// includes <spatium/vendor/stb_image_write.h> in this file (see
// tumbling_body_demo.cpp's matching comment for why not doing this
// causes a redefinition error).
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "io_helpers.hpp"

#include <spatium/io/scene.hpp>
#include <spatium/render/camera.hpp>
#include <spatium/render/parallel_for_rows.hpp>
#include <spatium/render/supersample.hpp>
#include <spatium/render/write_image.hpp>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <print>
#include <string>
#include <string_view>
#include <vector>

using spatium::Vec;
using spatium::geometry::Ray;
using spatium::geometry::RayHit;
using spatium::io::load_scene;
using spatium::io::resolve_shape;
using spatium::io::ResolvedShape;
using spatium::io::Scene;
using spatium::render::Camera;
using spatium::render::camera_ray_dir;
using spatium::render::make_camera_basis;
using spatium::render::parallel_for_rows;
using spatium::render::supersample_pixel;
using spatium::render::write_png_rgb;

namespace {

// One object, already resolved to a ray-testable shape plus the
// material info the render loop needs; scene.hpp's SceneObject/params
// aren't consulted again after this point.
struct PlacedShape {
    ResolvedShape<double> shape;
    Vec<double, 3> base_color; // linear RGB in [0,1] -- see io/scene.hpp's Material
    std::string name;
};

void print_usage() {
    std::print(
        "Usage: scene_demo [--scene PATH] [--output PATH] [--width N] [--height N] "
        "[--force] [--help]\n"
        "  Loads a JSON scene file (io/scene.hpp) and renders it to a PNG via the\n"
        "  shared CPU raytracer engine, resolving each object's shape through\n"
        "  io/scene.hpp's open shape-kind registry (sphere/box/torus built in).\n"
        "  --scene PATH   scene file to load (default: scene_demo_sample.json)\n"
        "  --output PATH  PNG output path (default: scene_demo.png)\n"
        "  --width N      image width (default 960)\n"
        "  --height N     image height (default 540)\n"
        "  --force        overwrite an existing output file\n"
        "  --help         show this message\n");
}

// Flat two-tone vertical gradient (world is Z-up, matching the default
// camera.up) -- just enough context to read object silhouettes against,
// not a physical sky model like render/sky.hpp's starfield.
Vec<double, 3> background_color(const Vec<double, 3>& dir) {
    double t = 0.5 * (dir[2] + 1.0);
    Vec<double, 3> horizon{225.0, 225.0, 232.0};
    Vec<double, 3> zenith{70.0, 110.0, 195.0};
    return Vec<double, 3>{horizon * (1.0 - t) + zenith * t};
}

} // namespace

int main(int argc, char** argv) {
    std::string scene_path = "scene_demo_sample.json";
    std::string output_path = "scene_demo.png";
    int width = 960, height = 540;
    bool force = false;

    for (int i = 1; i < argc; ++i) {
        std::string_view a = argv[i];
        if (a == "--help" || a == "-h") { print_usage(); return 0; }
        if (a == "--force") { force = true; continue; }
        if (a == "--scene" && i + 1 < argc) { scene_path = argv[++i]; continue; }
        if (a == "--output" && i + 1 < argc) { output_path = argv[++i]; continue; }
        if (a == "--width" && i + 1 < argc) { width = std::atoi(argv[++i]); continue; }
        if (a == "--height" && i + 1 < argc) { height = std::atoi(argv[++i]); continue; }
        std::print(stderr, "unknown option: {}\n", a);
        print_usage();
        return 1;
    }

    auto loaded = load_scene<double>(scene_path);
    if (!loaded) {
        std::print(stderr, "failed to load scene '{}': {}\n", scene_path,
                   loaded.error().to_string());
        return 1;
    }
    const Scene<double>& scene = *loaded;

    std::vector<PlacedShape> placed;
    for (auto& obj : scene.objects) {
        auto resolved = resolve_shape(obj);
        if (!resolved) {
            std::print(stderr, "warning: skipping object '{}': {}\n", obj.name,
                       resolved.error().to_string());
            continue;
        }
        placed.push_back({std::move(*resolved), obj.material.base_color, obj.name});
    }
    if (placed.empty()) {
        std::print(stderr, "scene '{}' has no resolvable objects, nothing to render\n",
                   scene_path);
        return 1;
    }

    Camera<double> cam = scene.camera.value_or(Camera<double>{
        .position = {0.0, -9.0, 3.5}, .target = {0.0, 0.0, 0.0},
        .up = {0.0, 0.0, 1.0}, .fov_deg = 40.0});
    const auto basis = make_camera_basis(cam);
    const Vec<double, 3> light = Vec<double, 3>{Vec<double, 3>{0.4, -0.5, 0.8}.normalized()};

    std::vector<std::uint8_t> img(3 * static_cast<std::size_t>(width) * height, 0);

    parallel_for_rows(height, [&](int y) {
        for (int x = 0; x < width; ++x) {
            std::uint8_t* px = &img[3 * (static_cast<std::size_t>(y) * width + x)];

            auto ray_color = [&](double sx, double sy) -> Vec<double, 3> {
                Vec<double, 3> dir = camera_ray_dir(basis, sx, sy);
                Ray<3, double> ray{cam.position, dir};

                bool found = false;
                double best_t = std::numeric_limits<double>::max();
                RayHit<double> best_hit{};
                const PlacedShape* best_shape = nullptr;

                for (auto& ps : placed) {
                    for (auto& h : ps.shape.ray_hits(ray)) {
                        if (h.t < best_t) {
                            best_t = h.t;
                            best_hit = h;
                            best_shape = &ps;
                            found = true;
                        }
                    }
                }
                if (!found) return background_color(dir);

                double diff = std::max(0.0, best_hit.normal.dot(light));
                double shaded = 0.18 + 0.82 * diff;
                return Vec<double, 3>{best_shape->base_color * 255.0 * shaded};
            };

            supersample_pixel(x, y, width, height, basis.tan_half,
                               static_cast<double>(width) / height, ray_color, px);
        }
    });

    if (spatium::examples::confirm_overwrite(output_path, force))
        write_png_rgb(output_path, width, height, img);

    std::print("scene_demo: rendered {} object(s) from '{}' to '{}' ({}x{})\n",
               placed.size(), scene_path, output_path, width, height);
    return 0;
}
