// Hyperbolic ray marching -- the third and smallest of the three
// manifold-backlog ideas (Riemannian optimization and geodesic procgen are
// done, see docs/ROADMAP.md), and the most visceral proof of the "any
// Riemannian manifold" thesis: Hyperbolic<N> has zero human visual
// intuition behind it.
//
// This is neither of the two rendering techniques already in the tree:
//   - not a mesh + BVH::ray_cast() (that stack is Euclidean ray-primitive
//     intersection -- there is no way to tessellate hyperbolic SPACE
//     itself into Euclidean triangles);
//   - not the RK4 geodesic integration physics/relativity/geodesic.hpp
//     uses for Schwarzschild/Kerr (those metrics have no closed-form
//     geodesic solution). Hyperbolic<N> does: exp_map() is exact.
//
// So this is sphere-tracing (the classic Euclidean SDF-marching technique)
// carried out entirely in the hyperboloid model's own metric: march
// exp_map(origin, dir, t) forward with growing t, using space.distance()
// to the nearest marker as a safe step size, exactly as an ordinary
// raymarcher uses a Euclidean SDF -- just with a non-Euclidean distance
// function.
//
// The camera sits at Hyperbolic<3>::origin() = (1,0,0,0) in the
// hyperboloid embedding. Its tangent space there is genuinely Euclidean:
// metric_at(origin, u, v) = minkowski(u, v), and every tangent vector at
// the origin has u[0] = 0 (orthogonal to the normal (1,0,0,0) under the
// Minkowski form), so minkowski(u, v) = -u[0]*v[0] + sum(u_i*v_i) reduces
// to the ordinary dot product on (u_1, u_2, u_3). That means the existing
// EUCLIDEAN render::Camera/camera_ray_dir pinhole formula gives the
// initial ray direction directly -- it just gets embedded into the
// hyperboloid's ambient Vec<T,4> with a leading 0 before marching. No new
// camera-basis math needed, no Riemannian Gram-Schmidt: reusing the flat
// formula IS correct here, not an approximation, because the tangent space
// at this one fixed point truly is flat.
//
// Markers sit at three shells of hyperbolic distance from the origin
// along the 12 icosahedron vertex directions (mesh/primitives.hpp) --
// reused as unit directions, not as a mesh. This is where hyperbolic
// space's exponential volume growth actually becomes visible: the same 12
// directions occupy far more of the visual field at the outer shell than
// linear (Euclidean) spacing would predict.
//
// Deliberately NOT attempted: Lambertian shading against a surface
// normal. These markers are distance-field balls (points + a hyperbolic-
// distance radius), not literal geometric primitives with a normal
// vector -- shading here is plain distance fog (closer = brighter), an
// honest simplification, not a hidden one.
//
// Output: hyperbolic_world.png

#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "io_helpers.hpp"

#include <spatium/mesh/primitives.hpp>
#include <spatium/mesh/subdivision.hpp>
#include <spatium/render/camera.hpp>
#include <spatium/render/color.hpp>
#include <spatium/render/parallel_for_rows.hpp>
#include <spatium/render/sky.hpp>
#include <spatium/render/supersample.hpp>
#include <spatium/render/write_image.hpp>
#include <spatium/spaces/hyperbolic.hpp>

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <print>
#include <string>
#include <string_view>
#include <vector>

using spatium::H3;
using spatium::Vec;
using spatium::mesh::icosahedron;
using spatium::mesh::subdivide_once;
using spatium::Sphere;
using spatium::render::hsv_to_rgb255;
using spatium::render::make_starfield;
using spatium::render::sample_sky_color;
using spatium::render::Sky;
using spatium::render::Camera;
using spatium::render::make_camera_basis;
using spatium::render::camera_ray_dir;
using spatium::render::parallel_for_rows;
using spatium::render::supersample_pixel;
using spatium::render::write_png_rgb;

namespace {

using Point4 = Vec<double, 4>;

struct Marker {
    Point4 pos;
    double radius;
    Vec<double, 3> color;
};

std::vector<Marker> build_markers(const H3& space) {
    // One subdivision: 12 -> 42 directions, still evenly spread (each new
    // vertex is an edge midpoint re-projected onto the unit sphere) --
    // the bare 12 icosahedron vertices left most of a wide-FOV frame
    // empty since they're spread across the WHOLE sphere of directions.
    Sphere<2, double> unit_sphere{};
    auto ico = subdivide_once(icosahedron<double>(), unit_sphere);
    Point4 origin = H3::origin();

    std::vector<double> shells{1.0, 1.8, 2.6};
    // Small relative to the shells: apparent angular radius of a
    // hyperbolic ball follows sin(alpha) = sinh(radius)/sinh(distance),
    // not radius/distance -- at radius 0.25 even the outer shell's
    // markers filled a third of the frame (hyperbolic distances make
    // things shrink far slower with depth than Euclidean intuition
    // expects). 0.08 keeps the near shell as small readable dots while
    // still letting the far shell's much smaller apparent size make the
    // exponential-volume-growth point.
    constexpr double kRadius = 0.08;

    std::vector<Marker> markers;
    markers.reserve(shells.size() * ico.vertex_count());
    for (std::size_t s = 0; s < shells.size(); ++s) {
        for (std::size_t i = 0; i < ico.vertex_count(); ++i) {
            Vec<double, 3> dir3 = ico.vertices[i]; // already unit (Sphere<2> radius 1)
            Point4 dir4{0.0, dir3[0], dir3[1], dir3[2]};
            Point4 p = space.exp_map(origin, dir4, shells[s]);
            // Golden-ratio-conjugate decorrelation: consecutive icosahedron
            // vertex indices are often spatially adjacent (see the vertex
            // list above), so hue=i/12 clustered same-ish colors together
            // whenever a camera view happened to catch neighboring
            // vertices -- this scatters hues across the wheel instead.
            double hue = std::fmod(static_cast<double>(i) * 0.6180339887498949, 1.0);
            Vec<double, 3> color = hsv_to_rgb255(hue, 0.6, 0.75 + 0.08 * static_cast<double>(s));
            markers.push_back({p, kRadius, color});
        }
    }
    return markers;
}

// Sphere-traces one ray from the origin. Returns the hit marker's color,
// fogged by traveled hyperbolic distance, or sky if nothing was reached
// within T_MAX.
Vec<double, 3> march(const H3& space, const std::vector<Marker>& markers,
                     const Point4& origin, const Vec<double, 3>& dir3, const Sky& sky) {
    constexpr double kTMax = 4.2;
    constexpr double kMinStep = 1e-4;
    constexpr double kHitEps = 1e-3;
    constexpr int kMaxSteps = 200;

    Point4 dir4{0.0, dir3[0], dir3[1], dir3[2]};

    double t = 0.0;
    for (int step = 0; step < kMaxSteps && t < kTMax; ++step) {
        Point4 pos = space.exp_map(origin, dir4, t);

        double d_min = std::numeric_limits<double>::infinity();
        std::size_t nearest = 0;
        for (std::size_t i = 0; i < markers.size(); ++i) {
            double d = space.distance(pos, markers[i].pos) - markers[i].radius;
            if (d < d_min) { d_min = d; nearest = i; }
        }

        if (d_min < kHitEps) {
            double fog = std::exp(-t / 2.2);
            double shade = 0.25 + 0.75 * fog;
            return Vec<double, 3>{markers[nearest].color * shade};
        }
        t += std::max(d_min, kMinStep);
    }
    return sample_sky_color(sky, dir3);
}

} // namespace

int main(int argc, char* argv[]) {
    bool force = false;
    std::string out_path = "hyperbolic_world.png";

    for (int i = 1; i < argc; ++i) {
        std::string_view a = argv[i];
        if (a == "--force") { force = true; continue; }
        if (a == "--out" && i + 1 < argc) { out_path = argv[++i]; continue; }
        if (a == "--help") {
            std::print("hyperbolic_tessellation_demo [--out path] [--force]\n"
                       "  Sphere-traces Hyperbolic<3>'s own hyperboloid-model metric --\n"
                       "  markers at 3 shells of hyperbolic distance along 12 icosahedron\n"
                       "  directions, camera fixed at the space's own origin().\n");
            return 0;
        }
        std::print(stderr, "unknown option: {}\n", a);
        return 1;
    }

    auto t0 = std::chrono::steady_clock::now();

    H3 space;
    Point4 origin = H3::origin();
    auto markers = build_markers(space);

    constexpr int W = 960, H = 720;
    // No spirals/clouds: this is an abstract mathematical space, not a
    // literal starfield -- make_starfield's nebulae are sized for the GR
    // raytracers' wide establishing shots (memory: 0.08-0.30pi angular
    // radius) and read as an oversized, out-of-place blob dominating this
    // demo's much narrower, sparser framing. Plain point stars only.
    Sky sky = make_starfield(1200, /*seed=*/11, /*tint=*/{4.0, 4.0, 10.0}, /*wide_sky=*/false);
    sky.spirals.clear();
    sky.clouds.clear();

    // Dummy Euclidean camera solely to reuse camera_ray_dir()'s pinhole
    // formula -- position/target/up here have no hyperbolic meaning, only
    // the resulting orthonormal (fwd,right,up) frame and tan_half matter.
    // Wide (near-fisheye) FOV deliberately: all 36 markers sit at just 12
    // fixed angular directions (three shells reuse the same 12 icosahedron
    // directions at growing hyperbolic distance) spread over the WHOLE
    // sphere of view, not clustered in a forward hemisphere -- a narrow,
    // "portrait" FOV only ever catches 1-3 of them by chance.
    const Camera<double> cam{
        .position = {0.0, 0.0, -1.0}, .target = {0.0, 0.0, 0.0}, .up = {0.0, 1.0, 0.0},
        .fov_deg = 130.0};
    const auto basis = make_camera_basis(cam);

    std::vector<std::uint8_t> img(3 * static_cast<std::size_t>(W) * H, 0);

    parallel_for_rows(H, [&](int y) {
        for (int x = 0; x < W; ++x) {
            std::uint8_t* px = &img[3 * (static_cast<std::size_t>(y) * W + x)];

            auto ray_color = [&](double sx, double sy) -> Vec<double, 3> {
                Vec<double, 3> dir3 = camera_ray_dir(basis, sx, sy);
                return march(space, markers, origin, dir3, sky);
            };

            supersample_pixel(x, y, W, H, basis.tan_half,
                              static_cast<double>(W) / H, ray_color, px);
        }
    });

    if (spatium::examples::confirm_overwrite(out_path, force))
        write_png_rgb(out_path, W, H, img);

    auto t1 = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::print("hyperbolic_tessellation_demo: {} markers, {:.0f} ms -> {}\n",
               markers.size(), ms, out_path);
    return 0;
}
