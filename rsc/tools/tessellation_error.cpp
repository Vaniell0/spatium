// How fine does a tessellation have to be before the eye stops seeing it?
//
// The second dispatch point in the scene pipeline, and it is currently
// decided by magic numbers at the call site: the donut's dough is 240x120,
// its dust specks 96x48, its sprinkles 8x2, and the library defaults are
// 48x24 for a torus and 12x4 for a cylinder. None of those came from a
// measurement, and none of them looks at how large the object actually
// appears.
//
// Which is the whole point: a tessellation is wrong when its deviation
// from the true surface is visible, and visibility is measured in pixels,
// not in triangles. So the quantity to measure is the largest distance
// between the tessellated surface and the real one -- taken at triangle
// midpoints, where a flat facet departs furthest from a curved surface --
// against the world size of one pixel.
//
// The deviation is exact rather than estimated: every surface here can
// `project` a point onto itself, so the error at a midpoint is the length
// of the correction that projection applies. No sampling, no bound, no
// approximation of an approximation.
//
// Build: part of the rsc tools target.

#include <spatium/mesh/mesh.hpp>
#include <spatium/mesh/primitives.hpp>
#include <spatium/spaces/chart.hpp>
#include <spatium/spaces/parametric.hpp>
#include <spatium/spaces/sphere.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <print>
#include <string>
#include <vector>

namespace {

using T = double;
using Surf = spatium::ParametricSurface<T>;
using Vec3 = spatium::Vec<T, 3>;

// Largest distance from the flat tessellation to the surface it
// approximates. Measured at triangle midpoints because that is where a
// chord departs furthest from an arc; the vertices sit on the surface by
// construction and would report zero.
T max_deviation(const Surf& s, std::size_t u_steps, std::size_t v_steps) {
    const auto m = spatium::mesh::parametric_mesh(s, u_steps, v_steps);
    T worst = 0;
    for (const auto& f : m.faces) {
        const Vec3 mid{(m.vertices[f[0]] + m.vertices[f[1]] + m.vertices[f[2]]) * (1.0 / 3.0)};
        const Vec3 on = s.project(mid);
        worst = std::max(worst, (on - mid).norm());
    }
    return worst;
}

std::size_t triangles(std::size_t u, std::size_t v) { return 2 * u * v; }

void sweep(const char* label, const Surf& s, T radius, T pixel) {
    std::println("{} — bounding radius {:.4f}, that is {:.2f} pixels across",
                 label, radius, 2 * radius / pixel);
    std::println("  {:>9} | {:>10} | {:>12} | {:>10}",
                 "steps", "triangles", "deviation", "in pixels");

    for (auto [u, v] : std::vector<std::pair<std::size_t, std::size_t>>{
             {4, 2}, {8, 4}, {12, 6}, {24, 12}, {48, 24}, {96, 48}, {240, 120}}) {
        const T dev = max_deviation(s, u, v);
        std::println("  {:>4}x{:<4} | {:>10} | {:>12.6f} | {:>10.3f}",
                     u, v, triangles(u, v), dev, dev / pixel);
    }
    std::println("");
}

}  // namespace

int main() {
    // One pixel in world units, taken from the donut scene's own camera —
    // the same 5.4e-3 the quaternion round-trip was measured against.
    constexpr T kPixel = 5.4e-3;

    std::println("Tessellation error against the size of a pixel");
    std::println("one pixel = {} world units, from the donut scene's camera", kPixel);
    std::println("A deviation under half a pixel cannot be seen, so every row below");
    std::println("that line is paid for and not delivered.");
    std::println("");

    // The dough, as the demo builds it: a torus at 240x120.
    sweep("torus 2.0 / 1.0 — the dough", spatium::make_torus<T>(2.0, 1.0), 3.0, kPixel);

    // A dust speck, as the demo builds it: a sphere of radius 0.001 at
    // 96x48. Its whole diameter is a fraction of a pixel.
    sweep("sphere r=0.001 — a dust speck",
          chart_of(spatium::Sphere<2, T>{0.001}), 0.001, kPixel);

    // A sprinkle: a cylinder at 8x2, the one place the demo tessellates
    // coarsely.
    sweep("cylinder r=0.011 h=0.062 — a sprinkle",
          spatium::make_cylinder<T>(0.011, 0.062), 0.035, kPixel);

    std::println("The rule this suggests is not a number per shape but a target: the");
    std::println("smallest step count whose deviation falls under half a pixel. What it");
    std::println("costs to get that wrong is visible in the triangle column.");
    return 0;
}
