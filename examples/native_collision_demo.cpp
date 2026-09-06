// Native rigid-body collision demo -- several spheres fall under
// gravity, collide with each other, and settle on a ground plane,
// entirely through Spatium's own collision path. Zero ipc-toolkit,
// zero SPATIUM_IPC_TOOLKIT: this is the thing that was actually
// missing before this session -- narrow_phase.hpp/xpbd.hpp already
// covered a particle against a STATIC analytical surface (cloth
// draping on a fixed sphere); nothing existed for two moving bodies
// colliding with each other. That gap is `physics/mechanics/
// rigid_contact.hpp` (sphere_sphere_contact / sphere_aabb_contact +
// XpbdSphereCollisionConstraint), and this demo is the visual proof
// it works end to end.
//
// Two different, deliberately different, collision mechanisms are on
// screen at once, matching the plan this code was built from:
//
//   - Sphere vs. sphere (dynamic vs. dynamic): DETECTION is continuous
//     (`sweep_sphere_sphere`, rigid_contact.hpp) -- the earliest time
//     of impact along each substep's actual swept motion, reusing
//     `geometry::ray_quadric`'s existing closed-form root solve rather
//     than only checking the predicted end-of-step positions, which a
//     fast enough pair could tunnel straight through. RESOLUTION is
//     still the hard XPBD position constraint (`XpbdSphereCollisionConstraint`
//     / `xpbd_solve_sphere_collision`, rigid_contact.hpp) -- same
//     compliance/Lagrange-multiplier math as xpbd.hpp's own distance
//     constraint, gated to only push apart, never pull together --
//     applied after `step_physics` rolls a detected pair back to its
//     time of impact.
//   - Sphere vs. ground (dynamic vs. STATIC obstacle): the existing
//     soft IPC barrier force from contact.hpp/narrow_phase.hpp
//     (`ContactQuery` + `ipc_contact_force`), applied as an ordinary
//     external acceleration during the XPBD predict step. The ground
//     plane isn't one of narrow_phase.hpp's built-in shapes (Sphere/
//     Torus/ParametricSurface), so `ground_accel` below constructs the
//     `ContactQuery` directly from `Hyperplane`'s own signed_distance/
//     project -- one line of geometry, zero new IPC math, in the same
//     spirit as narrow_phase.hpp's documented "add your own point_to
//     overload" ADL-extension pattern (see docs/architecture.md's
//     ContactSurface section), just inlined here since this is the
//     only caller.
//
// Broad phase is the plain O(n^2) AABB sweep from rigid_contact.hpp
// (`broad_phase_aabb_pairs`), rebuilt every physics substep -- exactly
// right at "several spheres," see that header's own comment for why
// spatial/bvh.hpp would be the wrong tool here.
//
// A small, deliberate, non-physical addition: `VELOCITY_DAMPING`
// bleeds a fraction of a percent of velocity per substep. The ground
// contact is a pure barrier *force* (not a position projection), so
// nothing in this file's physics dissipates energy -- an undamped
// sphere would bounce on the barrier forever, which looks like a bug
// on screen even though it is the honest behavior of a lossless
// contact potential. This is a demo convenience, not a claim about
// the underlying model; disclosed here rather than left to look like
// an oversight.
//
// Rendering reuses the render/ engine exactly as tumbling_body_demo.cpp
// established it: `render::Camera`/`make_camera_basis`/`camera_ray_dir`,
// `render::parallel_for_rows`, `render::supersample_pixel`,
// `render::write_png_rgb`, and `render::make_starfield`/
// `sample_sky_color` for the background. Spheres are ray-traced via
// `geometry::Quadric::sphere()` + `ray_quadric()` (translate the ray
// into each sphere's local frame -- cheap, and normals come out
// already in world space since translation doesn't rotate anything);
// the ground is `geometry::Hyperplane`/`intersect(Ray, Hyperplane)`
// with a checkerboard tint for depth cues.
//
// Output: collision_frames/frame_%04d.png -- assemble separately:
//   ffmpeg -framerate 60 -i collision_frames/frame_%04d.png -pix_fmt yuv420p native_collision.mp4

#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "io_helpers.hpp"

#include <spatium/algebra/vector.hpp>
#include <spatium/geometry/hyperplane.hpp>
#include <spatium/geometry/intersection.hpp>
#include <spatium/geometry/ray_surface.hpp>
#include <spatium/physics/mechanics/rigid_contact.hpp>
#include <spatium/render/camera.hpp>
#include <spatium/render/color.hpp>
#include <spatium/render/parallel_for_rows.hpp>
#include <spatium/render/sky.hpp>
#include <spatium/render/supersample.hpp>
#include <spatium/render/write_image.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <limits>
#include <print>
#include <random>
#include <string_view>
#include <vector>

using spatium::Vec;
using spatium::geometry::Hyperplane;
using spatium::geometry::intersect;
using spatium::geometry::Quadric;
using spatium::geometry::Ray;
using spatium::geometry::ray_quadric;
using namespace spatium::physics::mechanics;
using spatium::render::Camera;
using spatium::render::camera_ray_dir;
using spatium::render::hsv_to_rgb255;
using spatium::render::make_camera_basis;
using spatium::render::make_starfield;
using spatium::render::parallel_for_rows;
using spatium::render::sample_sky_color;
using spatium::render::Sky;
using spatium::render::supersample_pixel;
using spatium::render::write_png_rgb;

namespace {

constexpr int W = 960;
constexpr int H = 540;

constexpr int    N_SPHERES           = 7;
constexpr double DT_PHYSICS          = 1e-3;
constexpr int    SUBSTEPS_PER_FRAME  = 16;
constexpr int    N_FRAMES            = 240;
constexpr int    XPBD_ITERS          = 6;

constexpr double GRAVITY_Z           = -9.8;
constexpr double GROUND_D_HAT        = 0.05;   // IPC activation band, ground contact
constexpr double GROUND_KAPPA        = 200.0;  // tuned so equilibrium sinkage stays
                                                // well inside the band (~0.53 x d_hat
                                                // for these constants, checked
                                                // numerically against
                                                // ipc_barrier_grad -- not guessed)
constexpr double VELOCITY_DAMPING    = 0.999;  // see file header

constexpr double FOV_DEG = 42.0;

// Falling sphere -- position/velocity live in the reused XpbdParticle,
// radius and render color live alongside it (rigid_contact.hpp keeps
// radius off XpbdParticle deliberately, see XpbdSphereCollisionConstraint's
// doc comment; the demo needs it too, for gravity's ground term and for
// the AABB broad phase).
struct SphereBody {
    XpbdParticle<3> particle;
    double radius;
    Vec<double, 3> color;
};

// Ground contact acceleration for one sphere: builds the ContactQuery
// this file's header comment describes (Hyperplane's own signed_distance
// stands in for a `point_to` overload narrow_phase.hpp doesn't ship) and
// runs it through the existing, unmodified `ipc_contact_force`.
Vec<double, 3> ground_accel(const Vec<double, 3>& x, double radius, double inv_mass) {
    if (inv_mass <= 0.0) return Vec<double, 3>{};
    double signed_d = x[2] - radius;   // ground is the z = 0 plane
    ContactQuery<double> q{std::abs(signed_d), Vec<double, 3>{x[0], x[1], 0.0},
                           Vec<double, 3>{0.0, 0.0, 1.0}, signed_d < 0.0};
    Vec<double, 3> force = ipc_contact_force(q, GROUND_D_HAT, GROUND_KAPPA);
    return Vec<double, 3>{force * inv_mass};
}

// One physics substep for the whole scene: predict (gravity + ground
// barrier force), roll fast-crossing pairs back to their time of
// impact (continuous detection), then resolve remaining overlap by
// Gauss-Seidel XPBD -- broad-phased fresh every substep since every
// body has moved.
void step_physics(std::vector<SphereBody>& bodies, double dt, int n_iter) {
    const std::size_t n = bodies.size();

    // Predict -- same shape as xpbd_step's own predict loop (xpbd.hpp),
    // inlined here because this scene's external acceleration needs
    // each body's radius (for the ground term), which xpbd_step's
    // `external_accel(particle, dt)` callback signature doesn't expose.
    for (std::size_t i = 0; i < n; ++i) {
        auto& p = bodies[i].particle;
        if (p.w <= 0.0) { p.x_prev = p.x; continue; }
        Vec<double, 3> a = Vec<double, 3>{
            Vec<double, 3>{0.0, 0.0, GRAVITY_Z} +
            ground_accel(p.x, bodies[i].radius, p.w)};
        Vec<double, 3> v = Vec<double, 3>{(p.x - p.x_prev) / dt};
        Vec<double, 3> v_pred = Vec<double, 3>{(v + a * dt) * VELOCITY_DAMPING};
        p.x_prev = p.x;
        p.x = Vec<double, 3>{p.x + v_pred * dt};
    }

    // Broad phase over the SWEPT path (x_prev -> predicted x), not
    // just the predicted end position -- a fast mover's end-of-step
    // box alone can miss a body it actually crossed mid-step (see
    // rigid_contact.hpp's `swept_sphere_aabb`).
    std::vector<spatium::geometry::Box<3, double>> aabbs;
    aabbs.reserve(n);
    for (auto& b : bodies)
        aabbs.push_back(swept_sphere_aabb(
            b.particle.x_prev, b.radius,
            Vec<double, 3>{b.particle.x - b.particle.x_prev}));
    auto pairs = broad_phase_aabb_pairs(aabbs);

    // Continuous detection: for each broad-phase candidate, solve the
    // earliest time of impact along this substep's actual swept
    // motion (rigid_contact.hpp's `sweep_sphere_sphere`, reusing
    // `geometry::ray_quadric`'s closed-form root solve). Where a TOI
    // exists strictly before the end of the step, roll the predicted
    // positions back to that moment. This is the actual tunnelling
    // fix, and it is entirely a DETECTION-side change: the resolution
    // step right below (`xpbd_solve_sphere_collision`) is untouched --
    // without the rollback, a pair fast enough to cross paths mid-step
    // can end the substep on opposite, non-overlapping sides of each
    // other, and that unchanged constraint (which only ever looks at
    // CURRENT positions) would then see no overlap and do nothing.
    for (auto& pr : pairs) {
        auto& bi = bodies[pr.first];
        auto& bj = bodies[pr.second];
        Vec<double, 3> disp_i = bi.particle.x - bi.particle.x_prev;
        Vec<double, 3> disp_j = bj.particle.x - bj.particle.x_prev;
        auto sweep = sweep_sphere_sphere(bi.particle.x_prev, bi.radius, disp_i,
                                         bj.particle.x_prev, bj.radius, disp_j);
        if (sweep.hit && sweep.toi < 1.0) {
            bi.particle.x = Vec<double, 3>{bi.particle.x_prev + disp_i * sweep.toi};
            bj.particle.x = Vec<double, 3>{bj.particle.x_prev + disp_j * sweep.toi};
        }
    }

    std::vector<XpbdSphereCollisionConstraint<3>> cons;
    cons.reserve(pairs.size());
    for (auto& pr : pairs)
        cons.push_back({pr.first, pr.second, bodies[pr.first].radius,
                        bodies[pr.second].radius, 0.0, 0.0});

    std::vector<XpbdParticle<3>> parts;
    parts.reserve(n);
    for (auto& b : bodies) parts.push_back(b.particle);

    for (auto& c : cons) c.reset();
    for (int iter = 0; iter < n_iter; ++iter)
        for (auto& c : cons) xpbd_solve_sphere_collision(c, parts, dt);

    for (std::size_t i = 0; i < n; ++i) bodies[i].particle = parts[i];
}

std::vector<SphereBody> make_scene() {
    // A loose vertical column, small horizontal jitter, staggered
    // heights: spheres land almost on top of each other and push apart
    // sideways as they settle, so sphere-sphere collision (not just
    // sphere-ground) is actually visible on screen -- a wide horizontal
    // spread would mostly show independent drops onto the ground plane.
    std::mt19937 rng(2026);
    std::uniform_real_distribution<double> jitter(-0.25, 0.25);

    std::vector<SphereBody> bodies;
    bodies.reserve(N_SPHERES);
    for (int i = 0; i < N_SPHERES; ++i) {
        double radius = 0.30 + 0.10 * (static_cast<double>(i % 3));
        Vec<double, 3> pos{jitter(rng), jitter(rng), 1.0 + 0.9 * i};

        SphereBody b;
        b.particle.x = b.particle.x_prev = pos;
        b.particle.w = 1.0;   // equal mass for every sphere, kept simple on purpose
        b.radius = radius;
        b.color = hsv_to_rgb255(static_cast<double>(i) / N_SPHERES, 0.6, 0.95);
        bodies.push_back(b);
    }
    return bodies;
}

void print_usage() {
    std::print(
        "Usage: native_collision_demo [--frames N] [--force] [--help]\n"
        "  Native XPBD rigid-body collision -- spheres falling, colliding\n"
        "  with each other and a ground plane, zero ipc-toolkit dependency.\n"
        "  --frames N   frame count (default {})\n"
        "  --force      overwrite existing output files\n"
        "  --help       show this message\n"
        "  Output:      collision_frames/frame_%04d.png ({}x{} RGB)\n",
        N_FRAMES, W, H);
}

}  // namespace

int main(int argc, char** argv) {
    bool force = false;
    int n_frames = N_FRAMES;
    for (int i = 1; i < argc; ++i) {
        std::string_view a = argv[i];
        if (a == "--help" || a == "-h") { print_usage(); return 0; }
        if (a == "--force") { force = true; continue; }
        if (a == "--frames" && i + 1 < argc) { n_frames = std::atoi(argv[++i]); continue; }
        std::print(stderr, "unknown option: {}\n", a);
        return 1;
    }

    auto bodies = make_scene();
    Hyperplane<3, double> ground{Vec<double, 3>{0.0, 0.0, 1.0}, 0.0};
    Sky sky = make_starfield(1500, /*seed=*/7, /*tint=*/{4.0, 4.0, 8.0}, /*wide_sky=*/false);

    const Camera<double> cam{
        .position = {3.2, -3.6, 1.7}, .target = {0.0, 0.0, 0.5}, .up = {0.0, 0.0, 1.0},
        .fov_deg = FOV_DEG};
    const auto basis = make_camera_basis(cam);
    const Vec<double, 3> light = Vec<double, 3>{Vec<double, 3>{0.5, -0.6, 1.0}.normalized()};

    std::error_code ec;
    std::filesystem::create_directories("collision_frames", ec);

    auto t0 = std::chrono::steady_clock::now();

    for (int frame = 0; frame < n_frames; ++frame) {
        for (int s = 0; s < SUBSTEPS_PER_FRAME; ++s)
            step_physics(bodies, DT_PHYSICS, XPBD_ITERS);

        std::vector<std::uint8_t> img(3 * static_cast<std::size_t>(W) * H, 0);

        parallel_for_rows(H, [&](int y) {
            for (int x = 0; x < W; ++x) {
                std::uint8_t* px = &img[3 * (static_cast<std::size_t>(y) * W + x)];

                auto ray_color = [&](double sx, double sy) -> Vec<double, 3> {
                    Vec<double, 3> dir = camera_ray_dir(basis, sx, sy);
                    Ray<3, double> ray{cam.position, dir};

                    double best_t = std::numeric_limits<double>::infinity();
                    Vec<double, 3> hit_point{}, hit_normal{};
                    Vec<double, 3> base_color{200.0, 200.0, 210.0};
                    bool hit_anything = false;

                    // Ground plane.
                    if (auto p = intersect(ray, ground)) {
                        double t = (*p - ray.origin).dot(ray.direction);
                        if (t >= 0.0 && t < best_t) {
                            best_t = t;
                            hit_point = *p;
                            hit_normal = ground.normal;
                            bool checker = (static_cast<long long>(std::floor(hit_point[0])) +
                                           static_cast<long long>(std::floor(hit_point[1]))) % 2 == 0;
                            base_color = checker ? Vec<double, 3>{150.0, 150.0, 160.0}
                                                 : Vec<double, 3>{95.0, 95.0, 105.0};
                            hit_anything = true;
                        }
                    }

                    // Spheres -- ray translated into each sphere's local
                    // frame so the shared analytical Quadric::sphere()
                    // path applies unmodified; the hit normal comes back
                    // in world space unchanged since translation carries
                    // no rotation.
                    for (auto& b : bodies) {
                        Ray<3, double> local{Vec<double, 3>{ray.origin - b.particle.x}, dir};
                        auto hits = ray_quadric(local, Quadric<double>::sphere(b.radius));
                        if (hits.empty()) continue;
                        const auto& hit = hits.front();
                        if (hit.t >= 0.0 && hit.t < best_t) {
                            best_t = hit.t;
                            hit_point = ray.origin + dir * hit.t;
                            hit_normal = hit.normal;
                            base_color = b.color;
                            hit_anything = true;
                        }
                    }

                    if (!hit_anything) return sample_sky_color(sky, dir);

                    if (hit_normal.dot(dir) > 0.0) hit_normal = Vec<double, 3>{-hit_normal};
                    double diff = std::max(0.0, hit_normal.dot(light));
                    double shaded = 0.2 + 0.8 * diff;
                    return Vec<double, 3>{base_color * shaded};
                };

                supersample_pixel(x, y, W, H, basis.tan_half,
                                  static_cast<double>(W) / H, ray_color, px);
            }
        });

        char path[64];
        std::snprintf(path, sizeof(path), "collision_frames/frame_%04d.png", frame);
        if (spatium::examples::confirm_overwrite(path, force)) write_png_rgb(path, W, H, img);

        std::print("\r  frame {}/{}", frame + 1, n_frames);
    }

    auto t1 = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::print("\rnative_collision_demo: {} frames at {}x{}, {:.0f} ms ({:.0f} ms/frame)\n",
               n_frames, W, H, ms, ms / n_frames);
    std::print("Assemble with:\n"
               "  ffmpeg -framerate 60 -i collision_frames/frame_%04d.png "
               "-pix_fmt yuv420p native_collision.mp4\n");
    return 0;
}
