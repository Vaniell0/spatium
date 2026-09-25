// Ball pit -- a gallery/promotional showcase for native rigid-body
// collision (`physics/mechanics/rigid_contact.hpp`, see
// native_collision_demo.cpp for the introduction of that header). This
// demo pushes the same machinery further on purpose: instead of one
// ground plane and a handful of falling spheres, a full open-top vat
// (five flat walls) holds several dozen vividly-colored spheres that
// pour in staggered, jostle against each other and the walls, and
// settle into a dense pile -- the "many bodies, real container" scene
// the collision work was ultimately for.
//
// Physics, generalized from native_collision_demo.cpp:
//   - The vat is five `geometry::Hyperplane`s (floor + four side
//     walls) sharing ONE contact function, `wall_accel` below --
//     native_collision_demo.cpp's `ground_accel` handled only the
//     floor because that demo only had a floor; here the floor and
//     every wall are the same kind of query (a point against an
//     infinite plane), so one function generalizes over all five
//     instead of special-casing the floor. Same IPC log-barrier force
//     (`ipc_contact_force`, contact.hpp/narrow_phase.hpp) as before,
//     applied as an external acceleration during the XPBD predict
//     step -- zero changes to that machinery.
//   - Sphere vs. sphere collision is untouched from native_collision_
//     demo.cpp: continuous detection via `sweep_sphere_sphere`
//     (reduces to `geometry::ray_quadric`'s closed-form root solve),
//     resolved by the hard unilateral `XpbdSphereCollisionConstraint`.
//     `step_physics` below is line-for-line the same shape as that
//     demo's, with the single ground plane replaced by a loop over
//     the vat's five walls.
//   - Side walls are infinite planes for physics (an infinite wall
//     can't be tunnelled past sideways, which is one less failure
//     mode to worry about with dozens of bodies in a confined space)
//     but bounded to a finite rectangular patch for rendering --
//     see `on_wall_patch` below.
//   - Every sphere starts PARKED well outside the scene (pinned,
//     inverse mass zero) and is released into a real drop -- at a
//     random point over the vat -- one at a time, spread across most
//     of the run (`make_scene`/`activate_due_bodies` below). This is
//     a fix for a real failure mode a spawn-everything-at-once version
//     of this file hit: the vat is deliberately too small to hold
//     every sphere in one floor layer (see VAT_HALF_EXTENT's comment),
//     which means enough spheres spawning at once end up overlapping
//     each other badly enough that XPBD's correction locks several of
//     them into a permanently wedged, gravity-proof configuration --
//     a frozen pile, not a settled one. Releasing spheres sequentially
//     means at most a handful are ever airborne together, so the pour
//     always reads as a real drop-and-collide, and by the end of the
//     run the vat holds a genuine multi-layer pile rather than a
//     single scattered layer.
//
// Rendering, on top of what native_collision_demo.cpp established
// (render::Camera/make_camera_basis/camera_ray_dir, parallel_for_rows,
// supersample_pixel, write_png_rgb):
//   - The floor is an infinite checkerboard plane exactly like native_
//     collision_demo.cpp's ground; the four side walls reuse the same
//     `intersect(Ray, Hyperplane)` path but additionally reject any
//     intersection point outside their finite patch (lateral extent
//     +/- VAT_HALF_EXTENT, height 0..VAT_WALL_HEIGHT) -- an infinite
//     wall plane would otherwise occlude everything behind it,
//     including the open top a camera looking down into the vat needs
//     to see through.
//   - Deliberately NOT `render::Sky`/`make_starfield()`: that engine
//     piece is a procedural STARFIELD (point stars plus a few spiral-
//     galaxy/gas-cloud structures, see sky.hpp's own header comment)
//     and always adds a few of those regardless of the `wide_sky`
//     flag -- exactly the wrong backdrop for a bright, colorful ball
//     pit, where cosmic structure in the sky would read as a tonal
//     mismatch rather than a stylistic choice. `backdrop_color` below
//     is a small local two-tone gradient instead, sized to this one
//     demo the same way this file's own checkerboard/lighting are.
//   - Each sphere gets Blinn-Phong specular plus a fresnel-style rim
//     term (tinted from the backdrop's cool zenith color, complementing
//     the warm key light) on top of the diffuse term native_collision_
//     demo.cpp already had -- the glossy-plastic look real pit balls
//     have, and worth the extra per-hit cost precisely because spheres
//     are what this demo is showcasing.
//   - A hard shadow ray from every hit point toward the key light,
//     tested against every OTHER sphere (not the floor/walls -- see
//     that code's own comment), grounds the pile visually: contact
//     shadows between touching balls and the shadow a ball casts on
//     the floor are what read as "these are really resting on each
//     other," not just diffuse shading alone.
//   - The camera slowly orbits and dollies in over the run (see
//     `camera_for_frame`) so the frame sequence reads as a short
//     camera move across the settle-and-collide action, not a single
//     static shot repeated -- meant to be assembled into a slideshow
//     or GIF, not just eyeballed one frame at a time.
//
// Output: ball_pit_frames/frame_%04d.png -- assemble separately:
//   ffmpeg -framerate 30 -i ball_pit_frames/frame_%04d.png -pix_fmt yuv420p ball_pit.mp4

#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION

#include "io_helpers.hpp"

#include <spatium/algebra/vector.hpp>
#include <spatium/geometry/hyperplane.hpp>
#include <spatium/geometry/intersection.hpp>
#include <spatium/geometry/ray_surface.hpp>
#include <spatium/physics/mechanics/rigid_contact.hpp>
#include <spatium/render/camera.hpp>
#include <spatium/render/color.hpp>
#include <spatium/render/parallel_for_rows.hpp>
#include <spatium/render/supersample.hpp>
#include <spatium/render/texture.hpp>
#include <spatium/render/write_image.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <limits>
#include <numbers>
#include <numeric>
#include <print>
#include <random>
#include <string>
#include <string_view>
#include <vector>

using spatium::Vec;
using spatium::geometry::Box;
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
using spatium::render::parallel_for_rows;
using spatium::render::supersample_pixel;
using spatium::render::Texture;
using spatium::render::write_png_rgb;

namespace {

constexpr int W = 960;
constexpr int H = 540;

constexpr int    N_SPHERES          = 120;
constexpr double DT_PHYSICS         = 5e-4;
constexpr int    SUBSTEPS_PER_FRAME = 66;   // 66 x DT_PHYSICS = one 30 fps frame: real time
constexpr int    N_FRAMES           = 280;
constexpr int    XPBD_ITERS         = 8;    // more than native_collision_demo's 6 -- up
                                             // to N_SPHERES bodies plus 5 walls means more
                                             // simultaneous contacts per substep, which
                                             // needs a couple more Gauss-Seidel sweeps to
                                             // converge without visible jitter in the pile.

constexpr double GRAVITY_Z        = -9.8;
constexpr double WALL_D_HAT       = 0.08;   // IPC activation band -- wider than native_
                                             // collision_demo.cpp's 0.05 ground tuning.
                                             // At this body count a fast-falling sphere
                                             // can cross a narrow band almost entirely
                                             // within one substep, so the barrier force
                                             // only gets evaluated once at (or near) full
                                             // penetration depth instead of ramping up
                                             // smoothly across several substeps -- a
                                             // stiff-explicit-integration divergence, not
                                             // an occasional visual glitch: confirmed by
                                             // instrumenting a run at the narrower band,
                                             // where scene-wide max speed grew into the
                                             // hundreds and never decayed. Shared by the
                                             // floor and all four side walls, same as
                                             // WALL_KAPPA below -- all five are the same
                                             // kind of contact (a point against an
                                             // infinite plane) at the same sphere-radius/
                                             // mass scale.
constexpr double WALL_KAPPA       = 200.0;
constexpr double VELOCITY_DAMPING = 0.997;  // a touch more than native_collision_demo's
                                             // 0.999 -- far more simultaneous contacts in
                                             // a confined vat means more small-scale
                                             // jitter energy to bleed off for a clean
                                             // settle within N_FRAMES.
constexpr double MAX_SPEED        = 15.0;   // per-substep speed clamp -- see step_physics's
                                             // own comment for why this is needed at this
                                             // body count.

constexpr double VAT_HALF_EXTENT = 2.2;    // vat interior spans [-R, R] in x and y. Big
                                            // enough that the pour spreads into a heap
                                            // rather than filling wall to wall; the pile
                                            // still builds upward because every sphere
                                            // falls from one point -- see SPAWN_XY.
constexpr double VAT_WALL_HEIGHT = 1.3;    // rendered wall height above the floor

constexpr double MIN_RADIUS   = 0.14;
constexpr double MAX_RADIUS   = 0.22;
// Every sphere falls from one point, a stream rather than a scatter, and
// the next one appears only when the last has cleared it -- see
// activate_due_bodies. Off centre, so the heap it builds leans toward one
// wall and the pour reads as a pour.
const Vec<double, 3> SPAWN_POINT{0.7, -0.5, 2.6};
constexpr double SPAWN_JITTER = 0.35;      // the spawn is a small disc, not a point: every
                                            // sphere whose spot in it is clear comes out,
                                            // so the pour keeps pace instead of waiting
                                            // for one sphere at a time to fall clear
constexpr double SPAWN_SPEED  = 3.0;       // downward speed at release
constexpr double SPAWN_SPRAY  = 0.9;       // sideways speed at release, in a random direction.
                                            // Needed, not decoration: dropped straight down
                                            // from one point, spheres stacked into a column
                                            // that stood -- XPBD pushes along the line of
                                            // centres, and nothing tipped it -- and the top
                                            // sphere sat in the spawn point so no more came.
constexpr double RELEASE_WINDOW_FRACTION = 0.7;   // fraction of the run's substeps over which
                                                  // the stream runs

constexpr double FOV_DEG = 36.0;

// Camera slowly orbits and dollies in across the run -- see camera_for_frame.
// Height/radius chosen so the ray to the vat's far interior (aimed at the
// pile's expected surface height, not the bare floor -- this vat is sized
// to pile up, see VAT_HALF_EXTENT's comment) clears the near wall's rim
// with margin (the open-top framing this file's header comment describes)
// while the near strip of interior floor, right behind that wall, still
// resolves to the wall's own inner face -- checked numerically against
// VAT_HALF_EXTENT/VAT_WALL_HEIGHT above, not guessed: at azimuth 0 the ray
// to CAM_TARGET clears the rim by ~0.18-0.2 units at both ends of the
// dolly, and the ray to a point just inside the near wall lands ~0.86
// units below the rim.
constexpr double CAM_HEIGHT       = 6.6;
constexpr double CAM_RADIUS_START = 10.6;
constexpr double CAM_RADIUS_END   = 9.6;
constexpr double CAM_AZ_START_DEG = -18.0;
constexpr double CAM_AZ_END_DEG   = 16.0;
const Vec<double, 3> CAM_TARGET{0.0, 0.0, 0.55};

// "Parking" spot for a not-yet-released sphere -- stacked far above and
// well clear of both the vat and each other (spaced more than any
// possible diameter apart) so pinned, inactive bodies can never be
// mistaken for a broad-phase collision candidate with anything real.
// See make_scene's release-schedule comment for why spheres are held
// here instead of all spawning at once.
constexpr double PARK_BASE_Z = 40.0;
constexpr double PARK_SPACING = 1.0;

// Lighting.
const Vec<double, 3> LIGHT_DIR   = Vec<double, 3>{Vec<double, 3>{0.45, -0.5, 0.85}.normalized()};
const Vec<double, 3> SPEC_TINT{255.0, 248.0, 232.0};   // warm, matches the key light
const Vec<double, 3> RIM_TINT{150.0, 176.0, 208.0};    // cool, matches the backdrop zenith
constexpr double AMBIENT        = 0.22;
constexpr double SHADOW_DARKEN  = 0.55;   // fraction of direct light left in a shadowed
                                           // spot -- softer than a hard 1.0-vs-0.0 cutoff
                                           // so a ball's cast shadow on the floor reads as
                                           // shading, not a flat black cutout
constexpr double SPEC_POWER      = 42.0;
constexpr double SPEC_STRENGTH   = 0.6;    // sphere glossy-plastic highlight
constexpr double RIM_STRENGTH    = 0.28;   // sphere fresnel rim
constexpr double WALL_SPEC_POWER    = 18.0;
constexpr double WALL_SPEC_STRENGTH = 0.12; // faint ceramic sheen on the vat walls --
                                             // enough to separate them from the matte
                                             // floor, not enough to compete with the balls

// One vat sphere -- position/velocity live in the reused XpbdParticle
// (see rigid_contact.hpp), radius and render color live alongside it,
// same split native_collision_demo.cpp's SphereBody uses and for the
// same reason (XpbdSphereCollisionConstraint keeps radius off the
// particle deliberately). `active`/`release_substep`/`drop_pos` are
// this file's own addition for the sequential pour -- see make_scene's
// comment.
struct SphereBody {
    XpbdParticle<3> particle;
    double radius;
    Vec<double, 3> color;
    bool active;
    long long release_substep;
    Vec<double, 3> drop_pos;
    std::array<double, 2> spray{};   // sideways velocity at release -- see SPAWN_SPRAY
};

// Activate every sphere whose release time has arrived: move it from
// its far-away parking spot to its real drop point and switch its
// inverse mass on, so the very next predict step picks it up under
// gravity. See make_scene's comment for why release is staggered by
// TIME rather than by spawning every sphere at once at different
// heights.
void activate_due_bodies(std::vector<SphereBody>& bodies, long long substep_index) {
    // Every parked sphere whose turn has come and whose spot in the spawn
    // disc is clear of every sphere already out -- two spheres born
    // overlapping would start inside each other, and XPBD would pry them
    // apart with a kick no drop produces. Checked against the ones released
    // in this same pass too.
    for (auto& next : bodies) {
        if (next.active || substep_index < next.release_substep) continue;
        bool clear = true;
        for (const auto& b : bodies) {
            if (!b.active) continue;
            if (Vec<double, 3>{b.particle.x - next.drop_pos}.norm() < b.radius + next.radius + 0.02) {
                clear = false;
                break;
            }
        }
        if (!clear) continue;
        next.active = true;
        next.particle.w = 1.0;
        next.particle.x = next.drop_pos;
        // Released moving, as from a chute: x_prev behind x by one substep.
        const Vec<double, 3> v{next.spray[0], next.spray[1], -SPAWN_SPEED};
        next.particle.x_prev = Vec<double, 3>{next.drop_pos - v * DT_PHYSICS};
    }
}

// One side wall: an infinite plane for physics, a bounded rectangular
// patch for rendering. `lateral_axis` says which world axis (0 = x,
// 1 = y) spans the patch's finite width -- the other horizontal axis
// is pinned to the wall's plane, and z always runs 0..VAT_WALL_HEIGHT.
struct SideWall {
    Hyperplane<3, double> plane;
    int lateral_axis;
};

// Contact acceleration for one sphere against one flat plane (floor or
// side wall alike) -- generalizes native_collision_demo.cpp's
// `ground_accel` from "the one plane this scene has" to "any plane":
// `Hyperplane::signed_distance`/`project` stand in for a `point_to`
// overload narrow_phase.hpp doesn't ship for planes (the same
// substitution that file's header comment documents), and the result
// still runs through the existing, unmodified `ipc_contact_force`.
Vec<double, 3> wall_accel(const Hyperplane<3, double>& wall, const Vec<double, 3>& x,
                          double radius, double inv_mass) {
    if (inv_mass <= 0.0) return Vec<double, 3>{};
    double signed_d = wall.signed_distance(x) - radius;
    ContactQuery<double> q{std::abs(signed_d), wall.project(x), wall.normal, signed_d < 0.0};
    Vec<double, 3> force = ipc_contact_force(q, WALL_D_HAT, WALL_KAPPA);
    return Vec<double, 3>{force * inv_mass};
}

// The vat: floor + four side walls, all sharing `wall_accel`. Every
// side wall's inward normal points toward the vat's vertical axis and
// all four share offset = -VAT_HALF_EXTENT -- a direct consequence of
// the vat being a square prism centred on that axis.
std::array<Hyperplane<3, double>, 5> make_physics_walls() {
    constexpr double R = VAT_HALF_EXTENT;
    return {{
        Hyperplane<3, double>{Vec<double, 3>{0.0, 0.0, 1.0}, 0.0},    // floor, z = 0
        Hyperplane<3, double>{Vec<double, 3>{1.0, 0.0, 0.0}, -R},     // x = -R
        Hyperplane<3, double>{Vec<double, 3>{-1.0, 0.0, 0.0}, -R},    // x = +R
        Hyperplane<3, double>{Vec<double, 3>{0.0, 1.0, 0.0}, -R},     // y = -R
        Hyperplane<3, double>{Vec<double, 3>{0.0, -1.0, 0.0}, -R},    // y = +R
    }};
}

std::array<SideWall, 4> make_render_walls() {
    constexpr double R = VAT_HALF_EXTENT;
    return {{
        {Hyperplane<3, double>{Vec<double, 3>{1.0, 0.0, 0.0}, -R}, 1},   // x = -R, patch spans y
        {Hyperplane<3, double>{Vec<double, 3>{-1.0, 0.0, 0.0}, -R}, 1},  // x = +R, patch spans y
        {Hyperplane<3, double>{Vec<double, 3>{0.0, 1.0, 0.0}, -R}, 0},   // y = -R, patch spans x
        {Hyperplane<3, double>{Vec<double, 3>{0.0, -1.0, 0.0}, -R}, 0},  // y = +R, patch spans x
    }};
}

// Is `p` (already known to lie on `wall`'s infinite plane) inside the
// finite rectangle this demo actually renders as the wall? Rejecting
// points outside is what keeps the vat open-topped on screen -- see
// this file's header comment.
bool on_wall_patch(int lateral_axis, const Vec<double, 3>& p) {
    double lateral = p[static_cast<std::size_t>(lateral_axis)];
    return lateral >= -VAT_HALF_EXTENT && lateral <= VAT_HALF_EXTENT &&
           p[2] >= 0.0 && p[2] <= VAT_WALL_HEIGHT;
}

// One physics substep for the whole scene -- identical shape to
// native_collision_demo.cpp's `step_physics`, with that demo's single
// ground plane generalized to a loop over the vat's five walls.
void step_physics(std::vector<SphereBody>& bodies,
                  const std::array<Hyperplane<3, double>, 5>& walls, double dt, int n_iter) {
    const std::size_t n = bodies.size();

    for (std::size_t i = 0; i < n; ++i) {
        auto& p = bodies[i].particle;
        if (p.w <= 0.0) { p.x_prev = p.x; continue; }
        Vec<double, 3> a{0.0, 0.0, GRAVITY_Z};
        // The floor everywhere; a side wall only where the wall is. Above
        // the rim there is no wall to push, so a sphere flung over the edge
        // falls out instead of bouncing off something no one can see.
        a = Vec<double, 3>{a + wall_accel(walls[0], p.x, bodies[i].radius, p.w)};
        for (std::size_t k = 1; k < walls.size(); ++k) {
            const std::size_t lateral = k <= 2 ? 1 : 0;   // x-walls span y, y-walls span x
            const double r = bodies[i].radius;
            if (p.x[2] - r > VAT_WALL_HEIGHT) continue;
            if (std::abs(p.x[lateral]) > VAT_HALF_EXTENT + r) continue;
            a = Vec<double, 3>{a + wall_accel(walls[k], p.x, r, p.w)};
        }
        Vec<double, 3> v = Vec<double, 3>{(p.x - p.x_prev) / dt};
        Vec<double, 3> v_pred = Vec<double, 3>{(v + a * dt) * VELOCITY_DAMPING};

        // Speed clamp -- a deliberate robustness measure, same spirit
        // as VELOCITY_DAMPING above (see this function's own non-
        // physical additions): when two spheres collide hard enough
        // (or one lands deep in an already-dense part of the pile),
        // the XPBD correction for that substep can be large, and the
        // (x - x_prev)/dt velocity that implies for the NEXT substep's
        // predict can be large enough to tunnel a body across a full
        // side wall in one step -- side walls only get the discrete,
        // per-substep barrier force (see this file's header comment;
        // unlike sphere-sphere, there's no swept/continuous check for
        // sphere-wall), so a tunnelled body reappears outside the vat,
        // where the barrier's hard "already inside" push (contact.hpp/
        // narrow_phase.hpp's `ipc_contact_force`) fires at full
        // strength regardless of how deep outside it is and can
        // overshoot back across the whole vat, repeating in the other
        // direction -- an explicit-integration divergence, confirmed
        // by instrumenting a run without this clamp: max|v| across the
        // scene grew from single digits into the hundreds within ~30
        // frames of an otherwise unremarkable settle. MAX_SPEED is set
        // well above any legitimate speed in this scene (free-falling
        // the full DROP_HEIGHT reaches ~6.3 units/s) and well below the
        // observed divergence, so it never touches ordinary pour/settle
        // motion.
        double speed = v_pred.norm();
        if (speed > MAX_SPEED) v_pred = Vec<double, 3>{v_pred * (MAX_SPEED / speed)};

        p.x_prev = p.x;
        p.x = Vec<double, 3>{p.x + v_pred * dt};
    }

    // Broad phase over the SWEPT path -- see native_collision_demo.cpp's
    // matching comment for why (a fast mover's end-of-step box alone
    // can miss a body it actually crossed mid-step).
    std::vector<Box<3, double>> aabbs;
    aabbs.reserve(n);
    for (auto& b : bodies)
        aabbs.push_back(swept_sphere_aabb(
            b.particle.x_prev, b.radius, Vec<double, 3>{b.particle.x - b.particle.x_prev}));
    auto pairs = broad_phase_aabb_pairs(aabbs);

    // Continuous detection: roll fast-crossing pairs back to their
    // time of impact before the (unchanged) discrete XPBD resolution
    // below ever sees them -- see native_collision_demo.cpp's header
    // comment for the full reasoning; unchanged here.
    for (auto& pr : pairs) {
        auto& bi = bodies[pr.first];
        auto& bj = bodies[pr.second];
        Vec<double, 3> disp_i = bi.particle.x - bi.particle.x_prev;
        Vec<double, 3> disp_j = bj.particle.x - bj.particle.x_prev;
        // Only a pair that was apart at the start of the substep is swept.
        // The rewind below sends both spheres back to the moment of touch,
        // tangential motion included, so applied to a pair already resting
        // on each other it forbids them to slide -- and spheres poured onto
        // spheres stood in a column instead of rolling off. Touching pairs
        // are XPBD's to keep apart; CCD is for the ones arriving.
        const double gap0 = Vec<double, 3>{bi.particle.x_prev - bj.particle.x_prev}.norm() -
                            (bi.radius + bj.radius);
        if (gap0 <= 0.01) continue;
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

// Every sphere starts PARKED (pinned, w = 0, stacked far above and
// wide of the vat -- see PARK_BASE_Z/PARK_SPACING) and is released
// into a real drop one at a time, spread across RELEASE_WINDOW_FRACTION
// of the run's total substeps (`activate_due_bodies` does the actual
// hand-off). This is a direct fix for a real, confirmed failure mode
// of the more obvious alternative (spawn every sphere at once, spread
// across the vat footprint and staggered only by starting HEIGHT, an
// earlier version of this file did exactly that): to get a visible
// PILE rather than a single scattered layer, VAT_HALF_EXTENT is small
// enough that the floor can't hold every sphere at once (see that
// constant's own comment) -- which means spawn positions close enough
// together to eventually pile up are also close enough to spawn
// already overlapping if placed all at once. With sixty-plus bodies
// overlapping simultaneously, XPBD's corrections compound: several
// spheres locked into a permanently mutually-wedged configuration that
// gravity could never separate (confirmed by instrumenting a single
// stuck body's per-substep position -- it moved a few thousandths of a
// unit resolving its spawn overlap and then sat there, exactly
// net-zero velocity, indefinitely). Releasing spheres one at a time
// means at most a handful are ever in flight together, so ordinary
// sphere-sphere collision (continuous detection + hard XPBD
// resolution, unchanged from native_collision_demo.cpp) always has a
// realistic, non-degenerate configuration to resolve -- exactly the
// scenario that machinery was built for.
std::vector<SphereBody> make_scene(int n_frames) {
    std::mt19937 rng(20260906);
    std::uniform_real_distribution<double> u01(0.0, 1.0);

    long long total_substeps =
        static_cast<long long>(n_frames) * static_cast<long long>(SUBSTEPS_PER_FRAME);
    long long release_span =
        static_cast<long long>(static_cast<double>(total_substeps) * RELEASE_WINDOW_FRACTION);

    std::vector<int> order(static_cast<std::size_t>(N_SPHERES));
    std::iota(order.begin(), order.end(), 0);
    std::shuffle(order.begin(), order.end(), rng);


    std::vector<SphereBody> bodies;
    bodies.reserve(N_SPHERES);

    for (int i = 0; i < N_SPHERES; ++i) {
        SphereBody b;
        b.radius = MIN_RADIUS + (MAX_RADIUS - MIN_RADIUS) * u01(rng);
        b.color = hsv_to_rgb255(static_cast<double>(i) / N_SPHERES, 0.82, 0.98);
        b.active = false;

        double release_frac = N_SPHERES > 1
            ? static_cast<double>(order[static_cast<std::size_t>(i)]) / (N_SPHERES - 1)
            : 0.0;
        long long jitter = static_cast<long long>(
            (u01(rng) - 0.5) * static_cast<double>(release_span) / N_SPHERES);
        b.release_substep = static_cast<long long>(release_frac *
            static_cast<double>(release_span)) + jitter;
        if (b.release_substep < 0) b.release_substep = 0;

        const double spot_r = SPAWN_JITTER * std::sqrt(u01(rng)), spot_a = u01(rng) * 2.0 * std::numbers::pi;
        b.drop_pos = Vec<double, 3>{SPAWN_POINT[0] + spot_r * std::cos(spot_a),
                                    SPAWN_POINT[1] + spot_r * std::sin(spot_a), SPAWN_POINT[2]};
        const double heading = u01(rng) * 2.0 * std::numbers::pi;
        const double spray = SPAWN_SPRAY * (0.3 + 0.7 * u01(rng));
        b.spray = {spray * std::cos(heading), spray * std::sin(heading)};

        Vec<double, 3> park{0.0, 0.0, PARK_BASE_Z + static_cast<double>(i) * PARK_SPACING};
        b.particle.x = park;
        b.particle.x_prev = park;
        b.particle.w = 0.0;   // pinned until release_substep -- activate_due_bodies flips this

        bodies.push_back(b);
    }
    return bodies;
}

// Slow orbit-and-dolly camera move across the run: azimuth eases from
// CAM_AZ_START_DEG to CAM_AZ_END_DEG and radius eases from
// CAM_RADIUS_START down to CAM_RADIUS_END, both via a smoothstep so
// the move starts and ends at rest instead of cutting in/out at a
// constant angular velocity. Height and target stay fixed -- see this
// file's header comment for why (chosen so the near wall's rendered
// patch stays below the sightline into the vat across the whole
// sweep; a bigger swing risks dipping the camera behind that rim at
// one end).
// ── The vat's UV unwrap ─────────────────────────────────────────────
//
// The floor and the four inner walls laid out flat as one net, the way a
// box is cut from card: the floor in the middle, each wall folded out from
// the floor edge it stands on. A point on a wall at height z lands z away
// from that edge, so the net is continuous across every floor-to-wall
// crease and a texture painted on it wraps the vat without a seam there.
// The four corners of the square the net sits in belong to nothing.
//
// Net coordinates are in world units, the whole net a square of side
// NET_SIZE; (u, v) = net / NET_SIZE indexes a Texture.
constexpr double NET_SIZE = 2.0 * VAT_HALF_EXTENT + 2.0 * VAT_WALL_HEIGHT;

enum class VatFace { Floor, WallXNeg, WallXPos, WallYNeg, WallYPos, Outside };

// Where on the net a point on `face` goes.
Vec<double, 2> vat_net(const Vec<double, 3>& p, VatFace face) {
    constexpr double R = VAT_HALF_EXTENT, Hh = VAT_WALL_HEIGHT;
    switch (face) {
        case VatFace::Floor:    return {Hh + p[0] + R, Hh + p[1] + R};
        case VatFace::WallXNeg: return {Hh - p[2], Hh + p[1] + R};
        case VatFace::WallXPos: return {Hh + 2.0 * R + p[2], Hh + p[1] + R};
        case VatFace::WallYNeg: return {Hh + p[0] + R, Hh - p[2]};
        case VatFace::WallYPos: return {Hh + p[0] + R, Hh + 2.0 * R + p[2]};
        case VatFace::Outside:  break;
    }
    return {};
}

// Which face a net point belongs to, and its height above the floor when
// it is on a wall -- the inverse the default texture is drawn with.
VatFace net_face(double nx, double ny, double& height) {
    constexpr double R = VAT_HALF_EXTENT, Hh = VAT_WALL_HEIGHT;
    const bool mid_x = nx >= Hh && nx <= Hh + 2.0 * R, mid_y = ny >= Hh && ny <= Hh + 2.0 * R;
    height = 0.0;
    if (mid_x && mid_y) return VatFace::Floor;
    if (mid_y && nx < Hh)            { height = Hh - nx;            return VatFace::WallXNeg; }
    if (mid_y && nx > Hh + 2.0 * R)  { height = nx - Hh - 2.0 * R;  return VatFace::WallXPos; }
    if (mid_x && ny < Hh)            { height = Hh - ny;            return VatFace::WallYNeg; }
    if (mid_x && ny > Hh + 2.0 * R)  { height = ny - Hh - 2.0 * R;  return VatFace::WallYPos; }
    return VatFace::Outside;
}

// The default texture: pool tiles with grout, a mosaic band just under the
// rim, one lane line across the floor. `outline` draws the net's cut lines
// on top -- the template to paint a texture of one's own over.
Texture pool_texture(int size, bool outline) {
    constexpr double TILE = 0.24, GROUT = 0.018;
    const double texel = NET_SIZE / size;
    return Texture::generate(size, size, [&](int x, int y) {
        const double nx = (x + 0.5) * texel, ny = (y + 0.5) * texel;
        double height;
        const VatFace face = net_face(nx, ny, height);
        if (face == VatFace::Outside) return Vec<double, 3>{96.0, 96.0, 100.0};
        if (outline) {
            // One texel either side of every edge between faces.
            double h2;
            for (int dx = -1; dx <= 1; ++dx)
                for (int dy = -1; dy <= 1; ++dy)
                    if (net_face(nx + dx * texel, ny + dy * texel, h2) != face)
                        return Vec<double, 3>{255.0, 0.0, 200.0};
        }
        const double gx = std::fmod(nx, TILE), gy = std::fmod(ny, TILE);
        if (gx < GROUT || gy < GROUT) return Vec<double, 3>{226.0, 230.0, 226.0};
        // A tile's own shade, from which tile it is, so the grid does not
        // read as a printed pattern.
        const auto ti = static_cast<std::uint32_t>(nx / TILE), tj = static_cast<std::uint32_t>(ny / TILE);
        std::uint32_t h = (ti * 73856093u) ^ (tj * 19349663u);
        h ^= h >> 13; h *= 0x5bd1e995u; h ^= h >> 15;
        const double jitter = (static_cast<double>(h % 1000u) / 1000.0 - 0.5) * 18.0;
        Vec<double, 3> tile{84.0 + jitter, 174.0 + jitter, 196.0 + jitter * 0.6};
        const bool band = face != VatFace::Floor && height > VAT_WALL_HEIGHT - 0.20 &&
                          height < VAT_WALL_HEIGHT - 0.08;
        const bool lane = face == VatFace::Floor &&
                          std::abs(ny - (VAT_WALL_HEIGHT + VAT_HALF_EXTENT)) < 0.07;
        if (band || lane) tile = Vec<double, 3>{28.0 + jitter, 64.0 + jitter, 138.0 + jitter};
        return tile;
    });
}

Camera<double> camera_for_frame(int frame, int n_frames) {
    double t = n_frames > 1 ? static_cast<double>(frame) / (n_frames - 1) : 0.0;
    double ease = t * t * (3.0 - 2.0 * t);

    double az_deg = CAM_AZ_START_DEG + (CAM_AZ_END_DEG - CAM_AZ_START_DEG) * ease;
    double radius = CAM_RADIUS_START + (CAM_RADIUS_END - CAM_RADIUS_START) * ease;
    double az = az_deg * std::numbers::pi / 180.0;

    Vec<double, 3> pos{CAM_TARGET[0] + radius * std::sin(az),
                       CAM_TARGET[1] - radius * std::cos(az), CAM_HEIGHT};
    return Camera<double>{.position = pos, .target = CAM_TARGET, .up = {0.0, 0.0, 1.0},
                          .fov_deg = FOV_DEG};
}

// Soft two-tone studio backdrop: warm and pale near the horizon,
// easing into a cooler, deeper tone overhead -- a seamless-paper-
// backdrop look, not an outdoor sky. See this file's header comment
// for why `render::Sky`/`make_starfield()` (a procedural STARFIELD,
// always adding a few spiral-galaxy/gas-cloud structures) is the
// wrong tool for this particular demo.
Vec<double, 3> backdrop_color(const Vec<double, 3>& dir) {
    double t = std::clamp(dir[2] * 0.5 + 0.5, 0.0, 1.0);
    double ease = t * t * (3.0 - 2.0 * t);
    const Vec<double, 3> horizon{236.0, 226.0, 213.0};
    const Vec<double, 3> zenith{148.0, 170.0, 200.0};
    return Vec<double, 3>{horizon + (zenith - horizon) * ease};
}

void print_usage() {
    std::print(
        "Usage: ball_pit_demo [--frames N] [--texture PNG] [--unwrap-template PNG] [--force] [--help]\n"
        "  Ball pit -- native XPBD rigid-body collision (rigid_contact.hpp),\n"
        "  {} vividly-colored spheres poured into a five-wall vat, jostling\n"
        "  against each other and the walls before settling into a pile.\n"
        "  --frames N   frame count (default {})\n"
        "  --texture PNG           paint the vat with an image laid on its UV net\n"
        "  --unwrap-template PNG   write the default texture with the net's cut\n"
        "                          lines drawn on it, to paint over, and exit\n"
        "  --force      overwrite existing output files\n"
        "  --help       show this message\n"
        "  Output:      ball_pit_frames/frame_%04d.png ({}x{} RGB)\n",
        N_SPHERES, N_FRAMES, W, H);
}

}  // namespace

int main(int argc, char** argv) {
    bool force = false;
    int n_frames = N_FRAMES;
    std::string texture_path, template_path;
    for (int i = 1; i < argc; ++i) {
        std::string_view a = argv[i];
        if (a == "--help" || a == "-h") { print_usage(); return 0; }
        if (a == "--force") { force = true; continue; }
        if (a == "--frames" && i + 1 < argc) { n_frames = std::atoi(argv[++i]); continue; }
        if (a == "--texture" && i + 1 < argc) { texture_path = argv[++i]; continue; }
        if (a == "--unwrap-template" && i + 1 < argc) { template_path = argv[++i]; continue; }
        std::print(stderr, "unknown option: {}\n", a);
        return 1;
    }

    if (!template_path.empty()) {
        const auto tpl = pool_texture(1024, true);
        write_png_rgb(template_path, tpl.width, tpl.height, tpl.rgb);
        std::print("unwrap template -> {} ({} x {} world units, floor in the middle, walls "
                   "folded out)\n", template_path, NET_SIZE, NET_SIZE);
        return 0;
    }
    Texture vat_texture = pool_texture(1024, false);
    if (!texture_path.empty()) {
        auto t = Texture::load(texture_path);
        if (!t) { std::print(stderr, "{}\n", t.error().message); return 1; }
        vat_texture = std::move(*t);
    }

    auto bodies = make_scene(n_frames);
    const auto physics_walls = make_physics_walls();
    const auto render_walls = make_render_walls();
    const Hyperplane<3, double> floor_plane = physics_walls[0];

    // A warm terracotta vat against a cool neutral floor -- distinct
    // enough in both hue and value that the container reads as a
    // deliberate object, not just undifferentiated background geometry,
    // while staying muted enough that the rainbow-hued spheres (full
    // saturated hue sweep, see make_scene) remain the visual focus.
    const Vec<double, 3> floor_checker_a{188.0, 187.0, 182.0};
    const Vec<double, 3> floor_checker_b{128.0, 127.0, 122.0};

    std::error_code ec;
    std::filesystem::create_directories("ball_pit_frames", ec);

    auto t0 = std::chrono::steady_clock::now();
    long long substep_count = 0;

    for (int frame = 0; frame < n_frames; ++frame) {
        for (int s = 0; s < SUBSTEPS_PER_FRAME; ++s) {
            activate_due_bodies(bodies, substep_count);
            step_physics(bodies, physics_walls, DT_PHYSICS, XPBD_ITERS);
            ++substep_count;
        }

        const Camera<double> cam = camera_for_frame(frame, n_frames);
        const auto basis = make_camera_basis(cam);

        std::vector<std::uint8_t> img(3 * static_cast<std::size_t>(W) * H, 0);

        parallel_for_rows(H, [&](int y) {
            for (int x = 0; x < W; ++x) {
                std::uint8_t* px = &img[3 * (static_cast<std::size_t>(y) * W + x)];

                auto ray_color = [&](double sx, double sy) -> Vec<double, 3> {
                    Vec<double, 3> dir = camera_ray_dir(basis, sx, sy);
                    Ray<3, double> ray{cam.position, dir};

                    double best_t = std::numeric_limits<double>::infinity();
                    Vec<double, 3> hit_point{}, hit_normal{}, base_color{};
                    bool hit_anything = false;
                    bool hit_is_sphere = false;
                    bool hit_is_wall = false;
                    std::size_t hit_sphere_idx = 0;

                    // Floor -- infinite checkerboard plane, exactly as
                    // in native_collision_demo.cpp.
                    if (auto p = intersect(ray, floor_plane)) {
                        double t = (*p - ray.origin).dot(ray.direction);
                        if (t >= 0.0 && t < best_t) {
                            best_t = t;
                            hit_point = *p;
                            hit_normal = floor_plane.normal;
                            if (std::abs(hit_point[0]) <= VAT_HALF_EXTENT &&
                                std::abs(hit_point[1]) <= VAT_HALF_EXTENT) {
                                const auto n = vat_net(hit_point, VatFace::Floor);
                                base_color = vat_texture.sample(n[0] / NET_SIZE, n[1] / NET_SIZE);
                            } else {
                                bool checker = (static_cast<long long>(std::floor(hit_point[0])) +
                                               static_cast<long long>(std::floor(hit_point[1]))) % 2 == 0;
                                base_color = checker ? floor_checker_a : floor_checker_b;
                            }
                            hit_is_sphere = false;
                            hit_is_wall = false;
                            hit_anything = true;
                        }
                    }

                    // Side walls -- same plane intersection, bounded to
                    // each wall's finite rendered patch (see
                    // on_wall_patch's comment for why).
                    for (std::size_t wi = 0; wi < render_walls.size(); ++wi) {
                        const auto& sw = render_walls[wi];
                        if (auto p = intersect(ray, sw.plane)) {
                            double t = (*p - ray.origin).dot(ray.direction);
                            if (t >= 0.0 && t < best_t && on_wall_patch(sw.lateral_axis, *p)) {
                                best_t = t;
                                hit_point = *p;
                                hit_normal = sw.plane.normal;
                                const auto face = static_cast<VatFace>(1 + wi);   // same order as make_render_walls
                                const auto n = vat_net(*p, face);
                                base_color = vat_texture.sample(n[0] / NET_SIZE, n[1] / NET_SIZE);
                                hit_is_sphere = false;
                                hit_is_wall = true;
                                hit_anything = true;
                            }
                        }
                    }

                    // Spheres -- ray translated into each sphere's
                    // local frame, same as native_collision_demo.cpp.
                    for (std::size_t k = 0; k < bodies.size(); ++k) {
                        auto& b = bodies[k];
                        if (!b.active) continue;   // parked, not yet released -- see make_scene
                        Ray<3, double> local{Vec<double, 3>{ray.origin - b.particle.x}, dir};
                        auto hits = ray_quadric(local, Quadric<double>::sphere(b.radius));
                        if (hits.empty()) continue;
                        const auto& hit = hits.front();
                        if (hit.t >= 0.0 && hit.t < best_t) {
                            best_t = hit.t;
                            hit_point = ray.origin + dir * hit.t;
                            hit_normal = hit.normal;
                            base_color = b.color;
                            hit_is_sphere = true;
                            hit_is_wall = false;
                            hit_sphere_idx = k;
                            hit_anything = true;
                        }
                    }

                    if (!hit_anything) return backdrop_color(dir);
                    if (hit_normal.dot(dir) > 0.0) hit_normal = Vec<double, 3>{-hit_normal};

                    // Hard shadow ray toward the key light, tested
                    // against every OTHER sphere -- see this file's
                    // header comment for why the floor/walls are
                    // skipped (contact shadows between balls, and
                    // balls onto the floor, are the dominant visual
                    // cue in a pile this dense; skipping the planes
                    // keeps the cost to one extra ray_quadric sweep
                    // per primary hit instead of two).
                    Vec<double, 3> shadow_origin = Vec<double, 3>{hit_point + hit_normal * 1e-4};
                    bool occluded = false;
                    for (std::size_t k = 0; k < bodies.size(); ++k) {
                        if (hit_is_sphere && k == hit_sphere_idx) continue;
                        auto& b = bodies[k];
                        if (!b.active) continue;   // parked, not yet released -- see make_scene
                        Ray<3, double> local{Vec<double, 3>{shadow_origin - b.particle.x}, LIGHT_DIR};
                        auto hits = ray_quadric(local, Quadric<double>::sphere(b.radius));
                        if (!hits.empty() && hits.front().t > 1e-4) { occluded = true; break; }
                    }
                    double shadow = occluded ? SHADOW_DARKEN : 1.0;

                    double diff = std::max(0.0, hit_normal.dot(LIGHT_DIR));
                    double lit = AMBIENT + (1.0 - AMBIENT) * diff * shadow;
                    Vec<double, 3> color = Vec<double, 3>{base_color * lit};

                    if (hit_is_sphere) {
                        // Blinn-Phong specular plus a fresnel-style rim
                        // term -- the glossy-plastic look real pit
                        // balls have; the full-strength version is
                        // spent on spheres since they're what this
                        // demo showcases (walls get a much fainter
                        // version below).
                        Vec<double, 3> view = Vec<double, 3>{-dir};
                        Vec<double, 3> half = Vec<double, 3>{(LIGHT_DIR + view).normalized()};
                        double spec = std::pow(std::max(0.0, hit_normal.dot(half)), SPEC_POWER);
                        color = Vec<double, 3>{color + SPEC_TINT * (spec * SPEC_STRENGTH * shadow)};

                        double fres = 1.0 - std::max(0.0, hit_normal.dot(view));
                        double rim = std::pow(fres, 3.0) * RIM_STRENGTH;
                        color = Vec<double, 3>{color + RIM_TINT * rim};
                    } else if (hit_is_wall) {
                        // Faint ceramic sheen -- just enough specular to
                        // read as a glazed vat surface, not a matte one,
                        // without competing with the spheres' highlights.
                        Vec<double, 3> view = Vec<double, 3>{-dir};
                        Vec<double, 3> half = Vec<double, 3>{(LIGHT_DIR + view).normalized()};
                        double spec = std::pow(std::max(0.0, hit_normal.dot(half)), WALL_SPEC_POWER);
                        color = Vec<double, 3>{color + SPEC_TINT * (spec * WALL_SPEC_STRENGTH * shadow)};
                    }

                    return color;
                };

                supersample_pixel(x, y, W, H, basis.tan_half,
                                  static_cast<double>(W) / H, ray_color, px);
            }
        });

        char path[64];
        std::snprintf(path, sizeof(path), "ball_pit_frames/frame_%04d.png", frame);
        if (spatium::examples::confirm_overwrite(path, force)) write_png_rgb(path, W, H, img);

        std::print("\r  frame {}/{}", frame + 1, n_frames);
    }

    auto t1 = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::print("\rball_pit_demo: {} frames at {}x{}, {} spheres, {:.0f} ms ({:.0f} ms/frame)\n",
               n_frames, W, H, N_SPHERES, ms, ms / n_frames);
    std::print("Assemble with:\n"
               "  ffmpeg -framerate 30 -i ball_pit_frames/frame_%04d.png "
               "-pix_fmt yuv420p ball_pit.mp4\n");
    return 0;
}
