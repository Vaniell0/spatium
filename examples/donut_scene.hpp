#pragma once

// The donut scene itself, apart from the program that renders it.
//
// Split out of donut_demo.cpp so that something other than the renderer
// can build the *same* scene: tests/test_field_pod.cpp lowers every field
// in it and checks the plain-data interpreter against `eval_into`, and a
// check run over a copy of the scene would be a check of the copy. The
// demo's own header comment still describes what is built and why; this
// file only holds it.

#include <spatium/algebra/noise.hpp>
#include <spatium/io/build.hpp>
#include <spatium/render/camera.hpp>
#include <spatium/spaces/parametric.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <memory>
#include <numbers>
#include <random>
#include <string>
#include <utility>
#include <vector>

namespace donut {

using namespace spatium;
namespace bd = spatium::io::build;
using spatium::io::Material;
using spatium::render::Camera;
using spatium::render::make_camera_basis;

// Build-up timeline, in the same "t" materialize() already takes --
// shared by the dust particles' motion and the donut's grow-in, so
// there's one place that defines "when does what happen."
constexpr double T_EXPLODE = 0.9;   // cube -> dust, flying outward
constexpr double T_CONVERGE = 2.0;  // dust flies from its burst positions to the BOOM letterforms
constexpr double T_HOLD = 2.5;      // BOOM holds
constexpr double T_DISSOLVE = 3.2;  // dust shrinks to nothing
constexpr double T_DONUT_START = 2.9, T_DONUT_END = 3.9; // donut grows in, overlapping the dissolve
constexpr double T_BUILD_END = T_DONUT_END;

inline double smoothstep01(double t) {
    double e = std::clamp(t, 0.0, 1.0);
    return e * e * (3.0 - 2.0 * e);
}

// One dust particle's *core position* (before adding its own tiny local
// offset -- see particle_motion() below) at time t: annihilation, not a
// flight to a fixed point. Every particle swirls the whole time (three
// independent PerlinNoise channels sampled along its own path, not a
// literal orbit formula, but reads as turbulent wandering) growing out
// of the burst, and every particle is pulled toward its own nearby
// point on the BOOM letterforms as t moves through the converge window
// -- but by how much varies per particle (`pull_strength`, in [0,1]):
// near 1 and it lands and holds there; near 0 and the pull is just
// enough to bend its swirl in close past the letterform on the way by,
// not enough to stop it -- flying past *because* of the shape, not
// independently of it, which is also why they end up close enough to
// pick up color (dust_color below) before drifting off again. Only
// strongly-pulled particles get the extra "hits the surface and
// spreads" splat jitter near arrival -- a near-miss brushing past
// doesn't spread on impact, mid particles get a partial jitter.
inline Vec<double, 3> dust_core(double t, const Vec<double, 3>& burst_dir, double burst_dist,
                          const Vec<double, 3>& target, double pull_strength,
                          double swirl_seed, const algebra::PerlinNoise& swirl) {
    double explode_e = smoothstep01(t / T_EXPLODE);
    // A floor on the radius, not 0 at explode_e=0: 220 specks literally
    // coincident at the origin for the first several frames reproducibly
    // stalled BVH::build() for whole minutes (independent of the dust
    // shape -- tried both a UV-sphere and an icosahedron, same stall) --
    // this keeps them spread apart from frame 0 on, avoiding whatever
    // degenerate case that many near-identical bounding boxes hits.
    double radius = burst_dist * (0.08 + 0.92 * explode_e);
    double freq = 0.6 + 0.5 * std::fmod(swirl_seed, 1.0);
    double phase = swirl_seed * 11.0;
    Vec<double, 3> turb{
        swirl(burst_dir[0] * 3.0 + phase, burst_dir[1] * 3.0, t * freq),
        swirl(burst_dir[0] * 3.0 + phase + 7.0, burst_dir[1] * 3.0, t * freq),
        swirl(burst_dir[0] * 3.0 + phase, burst_dir[1] * 3.0 + 7.0, t * freq)};
    Vec<double, 3> swirl_pos = Vec<double, 3>{burst_dir * radius + turb * (0.5 + 0.6 * radius)};

    if (t < T_EXPLODE) return swirl_pos;

    double pull = smoothstep01((t - T_EXPLODE) / (T_CONVERGE - T_EXPLODE)) * pull_strength;
    Vec<double, 3> approach = Vec<double, 3>{swirl_pos * (1.0 - pull) + target * pull};

    // Splat: a jitter burst that peaks as the particle arrives (pull
    // near 1) and settles back down through the hold phase, not before.
    // Scales with pull_strength too -- a brush-past shouldn't spread.
    double arriving = pull * pull;
    double settle = t < T_HOLD ? 1.0 : 1.0 - smoothstep01((t - T_HOLD) / (T_DISSOLVE - T_HOLD + 0.001));
    Vec<double, 3> splat{
        swirl(target[0] * 5.0 + phase, target[1] * 5.0, t * 4.0),
        swirl(target[0] * 5.0 + phase + 3.0, target[1] * 5.0, t * 4.0),
        swirl(target[0] * 5.0 + phase, target[1] * 5.0 + 3.0, t * 4.0)};
    return Vec<double, 3>{approach + splat * (0.22 * arriving * settle * pull_strength)};
}

// How big a speck is at time t: full size until the hold ends, then
// shrinking to nothing through the dissolve. Split out from the motion so
// that "where the particle is" and "how big it is" are separately
// expressible -- which is what lets the motion be written as a placement
// rather than as an opaque map over vertices.
inline double dust_shrink(double t) {
    if (t < T_HOLD) return 1.0;
    return 1.0 - smoothstep01((t - T_HOLD) / (T_DISSOLVE - T_HOLD));
}

// The two halves together, as one point map. Kept because it is the
// clearest statement of what the motion *is* -- and it is exactly the
// form that cannot be instanced, since a reader cannot tell a rigid
// placement from a deformation through a closure. The scene builds the
// split form instead; this stays as the reference the split is checked
// against.
inline Vec<double, 3> particle_motion(const Vec<double, 3>& local_p, double t,
                                const Vec<double, 3>& burst_dir, double burst_dist,
                                const Vec<double, 3>& target, double pull_strength, double swirl_seed,
                                const algebra::PerlinNoise& swirl) {
    Vec<double, 3> pos = dust_core(t, burst_dir, burst_dist, target, pull_strength, swirl_seed, swirl);
    return Vec<double, 3>{pos + local_p * dust_shrink(t)};
}

inline std::vector<std::pair<double, double>> load_boom_points(const std::string& path) {
    std::vector<std::pair<double, double>> pts;
    std::ifstream in(path);
    double u, v;
    while (in >> u >> v) pts.emplace_back(u, v);
    return pts;
}

// Icing guide: the same torus formula make_torus() uses, restricted to a
// v-band around the tube's "top" (v=pi/2 sits opposite the hole -- see
// spaces/parametric.hpp's make_torus) -- so offset() only builds icing
// over the dough's upper cap, not the whole ring. Not promoted to the
// library: the exact band is donut-specific set dressing, not a general
// primitive (unlike offset_surface() itself, which this still goes
// through directly -- this is the base being offset, not a second
// surface projected onto the first).
inline ParametricSurface<double> torus_cap(double major_r, double minor_r,
                                     double v0 = std::numbers::pi * 0.02,
                                     double v1 = std::numbers::pi * 0.98) {
    constexpr double pi = std::numbers::pi;
    return ParametricSurface<double>(
        [=](double u, double v) -> Vec<double, 3> {
            return {
                (major_r + minor_r * std::cos(v)) * std::cos(u),
                (major_r + minor_r * std::cos(v)) * std::sin(u),
                minor_r * std::sin(v)};
        },
        {0.0, 2.0 * pi, v0, v1},               // wide band; the noisy thickness
                                               // falloff below does the actual
                                               // edge shaping, not this domain.
                                               // Being a band, this is an open
                                               // surface -- is_closed() says so,
                                               // and that is why the icing is
                                               // built with offset_shell()
        /*periodic_u=*/true, /*periodic_v=*/false);
}

inline Camera<double> hero_camera() {
    return {.position = {5.0, -4.2, 3.6}, .target = {0.0, 0.0, 0.0}, .up = {0.0, 0.0, 1.0}, .fov_deg = 38.0};
}

// Builds the whole scene into `scene` and returns the index of its root.
//
// Into a trace the caller owns, rather than returning one, because a
// Handle points at its trace: a trace moved out of here would leave every
// handle made while building it pointing at the old address.
inline std::size_t build_scene(bd::Trace<double>& scene, std::size_t sprinkle_count,
                               std::size_t dust_particles,
                               const std::vector<std::pair<double, double>>& boom_uv) {
    // Shared by dough/icing/sprinkles below: grows from nothing to full
    // size during [T_DONUT_START, T_DONUT_END] (overlapping the dust's
    // dissolve), so the donut visibly *replaces* the dust rather than
    // simply being present the whole time -- still one closed-form
    // function of (point, t), no simulation loop.
    // Written as an expression, not as a lambda, and the difference is not
    // style: `p * e(t)` is a uniform scale either way, but only the
    // expression form can be *seen* to be one. An opaque leaf that touches
    // the point sets reads_point, is_placement() answers "deformation",
    // and the node loses instancing -- for how it was spelled rather than
    // for what it does. This one line was costing the dough, the icing and
    // eighteen scatter nodes their instancing.
    //
    // A field is move-only, so each node gets its own rather than sharing
    // one; the closure is stateless, so that costs nothing worth naming.
    auto grow_scale = [] {
        // Written as an expression rather than as a lambda, which is what
        // makes these twenty nodes structural: the same smoothstep, in the
        // field vocabulary, with time as the first parameter.
        using F = bd::ScalarField<double>;
        return scaled(bd::VecField<double>::point(),
                      smoothstep((F::t() - F{T_DONUT_START}) /
                                 F{T_DONUT_END - T_DONUT_START}));
    };

    // Step 0 -- the default cube. Visible briefly, static, then it's
    // gone -- what actually reads as "the cube" from here on is the
    // dust field below (`dust`), not this node's own motion.
    auto cube = scene.cube({0.9, 0.9, 0.9})
                    .colored(Material<double>{.base_color = {0.55, 0.55, 0.58}})
                    // Same reason as grow_scale: a uniform scale, spelled so
                    // that it can be recognised as one. The factor is a hard
                    // step -- 1 until t = 0.12, then 0 -- and `less` is
                    // exactly that, both values bit for bit.
                    .moving(scaled(bd::VecField<double>::point(),
                                   less(bd::ScalarField<double>::t(),
                                        bd::ScalarField<double>{0.12})));

    // Step 0.5 -- delete the cube by *exploding* it: not a shrink this
    // time, real dust -- ~220 tiny cubes flying from the cube's own
    // volume, converging into the shape of the word "BOOM" (points
    // rasterized from a real font via ImageMagick, not hand-placed --
    // see examples/data/boom_points.txt and the shell one-liner that
    // made it), holding briefly, then dissolving before the donut grows
    // in. Each particle is its own tiny Literal node with its own
    // `.moving()` closure -- the DSL's per-node motion slot, just used
    // 220 times instead of once; still no simulation state anywhere.
    Camera<double> text_cam = hero_camera();
    auto text_basis = make_camera_basis(text_cam);
    constexpr double text_scale = 2.6;
    // Positioned relative to the camera's own basis (a point along its
    // forward axis, offset along its screen-up), not a guessed world
    // coordinate -- guarantees it lands on-screen regardless of exactly
    // where the camera sits, rather than trial-and-error against one
    // specific camera pose.
    const Vec<double, 3> text_center{
        Vec<double, 3>{text_cam.position + text_basis.fwd * 7.0 + text_basis.up * 0.85}};

    std::mt19937 burst_rng(11);
    std::uniform_real_distribution<double> unit(-1.0, 1.0);
    std::uniform_real_distribution<double> dist_amt(1.6, 2.6);
    std::uniform_real_distribution<double> unit01(0.0, 1.0);
    // Shared, not copied into every closure: each particle's motion hook
    // needs the same field, and capturing it by value duplicated a
    // 512-byte permutation table 19 800 times (~9.7 MB of identical
    // tables against a 12 MB L3). std::function forced that shape by
    // requiring a copyable callable; the DSL's slots are
    // move_only_function now, so a shared_ptr capture just works.
    auto swirl_noise = std::make_shared<const algebra::PerlinNoise>(4);

    // Grey far from any target, shifting toward red-yellow the closer a
    // particle's *current* position (recomputed live at whatever t this
    // gets asked about, via color_fn -- not baked in once) is to the
    // letterform point it's headed for -- unconditionally: even a
    // near-miss that only brushes past a letter on a weak pull picks up
    // color as it passes close, then fades back toward grey again as it
    // moves off. There's no separate "these ones never get colored"
    // case anymore; distance alone decides it.
    auto dust_color = [](const Vec<double, 3>& p, const Vec<double, 3>& target) {
        Vec<double, 3> grey{0.45, 0.44, 0.42};
        double d = std::clamp((p - target).norm() / 1.4, 0.0, 1.0);
        double hot = 1.0 - d;
        return Vec<double, 3>{grey * (1.0 - hot) + Vec<double, 3>{0.95, 0.35 + 0.35 * hot, 0.06} * hot};
    };

    // Each BOOM-letterform point gets several independent dust specks
    // converging on it (own burst direction/seed/pull each) rather than
    // regenerating a denser point cloud from the font raster -- denser
    // text *and* a denser swirling cloud together, same 220-point
    // letterform data. pull_strength is continuous, not a landed/missed
    // coin flip: most particles are weakly pulled (curve in close past a
    // letter, swept along by the shape rather than drifting
    // independently, then carried on past by the swirl) and a smaller
    // fraction pull strongly enough to actually land and hold.
    // ── The dust, as one node ─────────────────────────────────
    //
    // This used to be 35 200 handles: a loop calling flake() once per
    // particle, each with its own captured burst direction, target
    // letterform, pull, swirl seed and spin. It worked and it does not
    // scale -- a TraceNode is 456 bytes, so a million particles is
    // 456 MB of trace before a single one is drawn, and two million is
    // not reachable at all.
    //
    // One Scatter over an invisible unit sphere replaces all of it. Every
    // particle's parameters now come from *where it started*: the sphere
    // site is its burst direction, and everything else is a hash of that
    // same point. Nothing is stored per particle and the trace is three
    // nodes regardless of the count.
    //
    // What makes that expressible is MotionEnv::origin. A Scatter has one
    // motion field for all its instances, so before the origin existed
    // the only thing a field could tell them apart by was the vertex
    // position -- which makes the motion a deformation and refuses
    // instancing outright. An origin is one value per object, so a field
    // reading it stays affine in the point and stays shareable.
    const std::size_t dust_count = dust_particles;

    // Seeded from the particle's own starting point, so it is a pure
    // function of the site and survives being recomputed anywhere. FNV
    // over the bit patterns, then a final avalanche.
    auto dust_hash = [](const Vec<double, 3>& o, std::uint32_t salt) {
        std::uint32_t h = 2166136261u ^ salt;
        for (int k = 0; k < 3; ++k) {
            std::uint64_t bits = 0;
            double v = o[k];
            std::memcpy(&bits, &v, sizeof(bits));
            h ^= static_cast<std::uint32_t>(bits ^ (bits >> 32));
            h *= 16777619u;
        }
        h ^= h >> 13; h *= 0x85ebca6bu; h ^= h >> 16;
        return h;
    };
    auto unit_from = [dust_hash](const Vec<double, 3>& o, std::uint32_t salt) {
        return static_cast<double>(dust_hash(o, salt) % 1000000u) / 1000000.0;
    };

    auto boom_points = std::make_shared<std::vector<Vec<double, 3>>>();
    boom_points->reserve(boom_uv.size());
    for (auto [u, v] : boom_uv)
        boom_points->push_back(Vec<double, 3>{
            text_center + text_basis.right * (u * text_scale) + text_basis.up * (v * text_scale)});

    // Every per-particle quantity the old loop captured, rebuilt from the
    // origin instead. Same distributions, same meanings -- see the
    // comments on dust_core for why pull_strength is continuous rather
    // than a landed/missed coin flip.
    auto target_of = [unit_from, boom_points](const Vec<double, 3>& o) {
        return (*boom_points)[static_cast<std::size_t>(unit_from(o, 7u) *
                                                       static_cast<double>(boom_points->size())) %
                              boom_points->size()];
    };
    // More of them land, and the ones that do not are pulled harder than
    // before. At 30% landing and no pull on the rest the cloud spread
    // across the whole frame and the letterforms read as a faint tint
    // inside it rather than as writing; the burst is supposed to *become*
    // the word, not drift past it.
    auto pull_of = [unit_from](const Vec<double, 3>& o) {
        double roll = unit_from(o, 11u);
        if (roll < 0.55) return 0.85 + 0.15 * unit_from(o, 13u);   // land and hold
        return 0.25 * unit_from(o, 13u);                            // bent in close on the way past
    };

    auto dust = scene.scatter(scene.flake(Vec<double, 3>{0.010, 0.010, 0.003}),
                              // A *tiny* sphere, and the size is the point.
                              // A Scatter seats each instance on its site,
                              // so the seat is added to wherever the motion
                              // sends it -- which is exactly right when the
                              // scatter means "put these on that surface"
                              // and exactly wrong here, where the surface
                              // is only a way of handing every particle a
                              // distinct starting point. At radius 1 every
                              // flake was displaced by a unit vector from
                              // its own trajectory, which smeared the
                              // letterforms into a haze: 7 393 particles
                              // near their target and 52 actually on it.
                              // At 1e-3 the sites stay distinct for the
                              // hash and the displacement is nothing.
                              scene.sphere(0.001, 96, 48), dust_count, 4242,
                              /*seat=*/0.0)
                    // Slightly see-through, because dust is: thousands of
                    // opaque flecks stack into a wall, and at 0.72 the
                    // cloud has depth.
                    .colored(Material<double>{.roughness = 0.85, .opacity = 0.72})
                    .colored(bd::PointField<double>{[target_of, dust_color](const bd::MotionEnv<double>& e) {
                        return dust_color(e.p, target_of(e.origin));
                    }})
                    // The letterforms light up, and only while a speck is
                    // near one. Squared, so the glow arrives late and
                    // sharply rather than as a haze over the whole cloud.
                    .glowing(bd::PointField<double>{[target_of](const bd::MotionEnv<double>& e) {
                        double d = std::clamp((e.p - target_of(e.origin)).norm() / 0.75, 0.0, 1.0);
                        double hot = (1.0 - d) * (1.0 - d);
                        return Vec<double, 3>{Vec<double, 3>{1.00, 0.52, 0.12} * (1.35 * hot)};
                    }})
                    .moving(bd::VecField<double>::opaque_per_instance(
                                [unit_from, target_of, pull_of, swirl_noise](
                                    const Vec<double, 3>& origin, double time) {
                                    // The site direction is the burst
                                    // direction; its length is an artefact
                                    // of the carrier sphere and is
                                    // normalised away.
                                    Vec<double, 3> dir{Vec<double, 3>{origin}.normalized()};
                                    return dust_core(time, dir, 1.3 + 0.9 * unit_from(origin, 3u),
                                                     target_of(origin), pull_of(origin),
                                                     unit_from(origin, 5u), *swirl_noise);
                                })
                            // The turn wraps the point, not the whole
                            // motion: it orients the flake about its own
                            // centre rather than swinging it round the
                            // origin. Per-instance, or every flake in the
                            // cloud would tumble in lockstep.
                            // The shrink is `dust_shrink` as an expression.
                            // Its `if (t < T_HOLD) return 1.0` is dropped
                            // rather than modelled: below T_HOLD the
                            // argument is negative, the clamp takes it to
                            // zero and the smoothstep with it, so the
                            // branch and the expression agree bit for bit
                            // -- checked over [-1, 5], not assumed.
                            + rotated(scaled(bd::VecField<double>::point(),
                                             bd::ScalarField<double>{1.0} -
                                                 smoothstep((bd::ScalarField<double>::t() -
                                                             bd::ScalarField<double>{T_HOLD}) /
                                                            bd::ScalarField<double>{T_DISSOLVE -
                                                                                    T_HOLD})),
                                      [unit_from](const Vec<double, 3>& o, double time) {
                                          Vec<double, 3> axis{
                                              Vec<double, 3>{unit_from(o, 17u) * 2.0 - 1.0,
                                                             unit_from(o, 19u) * 2.0 - 1.0,
                                                             unit_from(o, 23u) * 2.0 - 1.0}
                                                  .normalized()};
                                          double phase = unit_from(o, 29u) * 6.283185307179586;
                                          double rate = 0.5 + 1.5 * unit_from(o, 31u);
                                          return Vec<double, 3>{axis * (phase + rate * time)};
                                      }));

    // Step 1 -- the dough is a torus, offset by a fine noise bump so it
    // actually has bready surface texture (not just a rough *shading*
    // response on a mathematically perfect torus -- real geometry, so
    // the bumps show up in the silhouette and catch light like bumps
    // should, not a flat-shaded illusion of them).
    // Calibrated against a real reference photo (Evan-Amos, "a pink,
    // frosted doughnut... Dunkin' Donuts", public domain, Wikimedia
    // Commons File:Pink-Frosted-Donut.jpg) rather than guessed -- two
    // things it corrected: real dough is close to smooth (a matte
    // *shading* response, not fine bumpy geometry -- the fine-noise
    // texture from an earlier pass was reading as noise, not "bread");
    // real icing covers almost the entire top uniformly, edge mostly
    // clean along the torus's natural equator, with a *few* localized
    // drips past it, not an all-over ragged boundary.
    //
    // Every surface field below is an expression rather than a lambda, and
    // each reads exactly as the lambda it replaced did, operation for
    // operation -- which is what keeps the frame byte-identical, and why
    // `0.0 - n` stands where `-n` did: the two differ only in the sign of
    // a zero, and subtracting 0.18 next erases it.
    using F = bd::ScalarField<double>;
    algebra::PerlinNoise surface_noise(2);
    auto dough_bump = F{0.020} * noise(surface_noise, F::u() * F{3.0}, F::v() * F{3.0},
                                       F{0.0}); // visible but still broad, not fine-grain
    // The crumb: sparse pits, and the shape of the field is the whole
    // point. The note above records that an earlier pass added symmetric
    // fine-grain noise and it read as sandpaper rather than as bread --
    // that finding stands, and this is built to avoid it rather than to
    // ignore it. Only the deep negative lobe of a second field survives
    // the threshold, so most of the surface stays smooth and a minority
    // of it is pitted, which is what a crumb actually looks like: voids,
    // not roughness. Squared, so a pit is shallow at its rim and deepens
    // toward the middle instead of being a flat-bottomed dent.
    //
    // Applied to the dough only, not to the shared `dough_bump`, because
    // that field is also the icing's base -- glaze pools over the crumb
    // and hides it, so pores pushing through the icing would be wrong.
    algebra::PerlinNoise pore_noise(11);
    auto dough_surface = [&] {
        F broad = F{0.020} * noise(surface_noise, F::u() * F{3.0}, F::v() * F{3.0}, F{0.0});
        // u runs around the major circle and v around the tube, so a
        // radian of u covers ~2.5x the arc a radian of v does; the
        // frequencies are in that ratio so the pits come out round
        // rather than smeared along the ring.
        F n = noise(pore_noise, F::u() * F{11.0}, F::v() * F{4.5}, F{0.0});
        F pit = max(F{0.0}, F{0.0} - n - F{0.18});
        return broad - F{0.075} * pit * pit;
    }();
    auto dough_base = scene.torus(2.0, 1.0, 240, 120);
    auto dough = scene.offset(dough_base, dough_surface)
                     .colored(Material<double>{.base_color = {0.87, 0.58, 0.27}, .roughness = 0.92})
                     .moving(grow_scale());

    // Step 2 -- the icing is the dough's own surface, offset outward --
    // over almost the whole top (torus_cap), clean-edged except for a
    // handful of localized drips, gently undulating (glaze pools and
    // settles in broad waves, not fine ripples) rather than a mirror.
    algebra::PerlinNoise icing_noise(3);
    // offset_shell, not offset: torus_cap is a band with two rims, so
    // the result has an edge and the rule for it has to be stated. The
    // rule is the one this demo always used -- the thickness field below
    // falls to zero before the rim, so the icing meets the dough there
    // and closes against it -- it just has a name now. (The dough above
    // is a plain offset(): a whole torus is closed, so there is no rim
    // to rule on and nothing to state.)
    // The band stops at 0.88pi, not 0.98pi: real glaze stops short of the
    // hole, and the tightest curvature on the whole torus is right there,
    // which is also where an offset of any thickness comes closest to
    // overrunning its own centre of curvature -- the pink fins.
    //
    // It has to be the *domain* that shrinks, not the thickness. Driving
    // the thickness to zero over a band while the surface still exists
    // there lays a whole strip of icing exactly on the dough underneath,
    // and two coincident surfaces read as stripes, which is worse than
    // the fins were. Zero thickness at a single rim is what
    // EdgeRule::ZeroThickness means; zero thickness across a region is
    // just a surface with nothing to do.
    constexpr double icing_rim = std::numbers::pi * 0.88;
    auto icing_base = scene.offset_shell(
        scene.space(torus_cap(2.0, 1.0, std::numbers::pi * 0.02, icing_rim), 160, 64),
        dough_bump, bd::EdgeRule::ZeroThickness);
    auto icing = scene.offset_shell(icing_base, [&] {
                        constexpr double pi = std::numbers::pi;
                        const F u = F::u(), v = F::v();
                        // The two rims are not interchangeable. v -> 0 is
                        // the outer equator and v -> pi is the wall of the
                        // hole, and glaze runs off the outside; it does not
                        // climb into the middle. Adding the drip to both --
                        // which the symmetric `min` used to do -- pushed
                        // icing down the inside of the hole, where the
                        // surface turns sharply and the offset overruns its
                        // own centre of curvature. That is what the pink
                        // fins around the hole were: an offset surface
                        // self-intersecting, not a shading artefact. So the
                        // drip is added to the outer distance only, which
                        // is both the cheap fix and the correct one.
                        F outer = v - F{pi * 0.02};        // toward the outer equator
                        F inner = F{icing_rim} - v;        // toward the hole
                        F wobble = noise(icing_noise, cos(u) * F{2.0}, sin(u) * F{2.0}, F{0.0}) *
                                   F{pi * 0.03};
                        F drip = max(F{0.0}, noise(icing_noise, cos(u) * F{1.3}, sin(u) * F{1.3},
                                                   F{8.0}) - F{0.35}) *
                                 F{pi * 0.35};
                        F edge_dist = min(outer + drip, inner);
                        F falloff = smoothstep((edge_dist + wobble) / F{pi * 0.09}); // soft, not torn
                        // 0.045 until 2026-09-17, tuned when offset()
                        // silently rendered this at 48x24 and the waves
                        // were under-sampled into near-smoothness. With
                        // the base's own 160x64 actually reaching the
                        // mesh they resolve fully, and +-45% of a 0.10
                        // thickness reads as lumps rather than as glaze
                        // settling. The field did not change; what
                        // changed is that it is now being listened to.
                        F pooling = noise(icing_noise, u * F{3.0}, v * F{3.0}, F{4.0}) *
                                    F{0.016}; // broad, gentle waves
                        return (F{0.10} + pooling) * falloff;
                    }(), bd::EdgeRule::ZeroThickness)
                     .colored(Material<double>{.base_color = {0.98, 0.55, 0.68}, .roughness = 0.40})
                     .moving(grow_scale());

    // Step 3 -- sprinkles scatter across the icing, area-weighted, in
    // several colors -- one scatter() call per color (each its own
    // area-weighted placement, not a single-color scatter recolored
    // after the fact) so the mix reads as actual rainbow sprinkles.
    // Flat slivers, not standing cylinders (scatter() maps an item's
    // local z to the target's normal, so a shape *long* in z stands up
    // like a quill; long in x/y and thin in z lies flat instead) --
    // scatter() also gives each instance its own random spin about the
    // normal now, so they no longer all face the same way.
    // solid_cylinder()'s length runs along local z -- exactly the axis
    // scatter() maps onto the target's *normal*, so used as-is it stands
    // the sprinkle straight up again (the quill problem, back). Permute
    // (x,y,z) -> (z,x,y) so the length axis becomes x (a tangent
    // direction instead), then nudge down slightly so it sits sunk into
    // the icing rather than merely resting exactly half-in.
    // A real cylinder now, not a hand-built mesh with its axes permuted.
    //
    // The permutation was there because scatter() could only stand an item
    // *up* on the normal, and a sprinkle lies down -- so the only way to
    // say so was to rotate the vertices, which a closed form cannot
    // follow. cylinder() records a BoundedQuadric about z; permute the
    // mesh and the exact form no longer describes it, so the node had to
    // be a Literal, and a Literal cannot be instanced. SeatAxis::X says
    // "this item's local x points along the normal" instead, which leaves
    // its z -- the cylinder's own axis -- lying in the surface, and the
    // closed form intact.
    //
    // Much smaller than they were and many more of them: a real sprinkle
    // is a fleck, and 600 fat ones read as gravel.
    auto sprinkle = scene.cylinder(0.011, 0.062, 8, 2);
    std::vector<Vec<double, 3>> sprinkle_colors{
        {0.95, 0.20, 0.25}, {0.98, 0.75, 0.15}, {0.25, 0.65, 0.35},
        {0.30, 0.45, 0.90}, {0.85, 0.30, 0.75}, {0.98, 0.98, 0.95}};
    // Two bands rather than one, which is how the density gradient is
    // got without a weight in the sampler: sample_surface_uniform is
    // uniform *by area* and has no notion of "more here than there". A
    // dense scatter over the crown and a sparse one over the whole cap
    // add up to a falloff -- most of them on top, a scattering of
    // stragglers running down the sides, which is what a real donut
    // looks like and is composition rather than new API.
    // Three nested bands, not one scatter and not two, and the reason is
    // that `sample_surface_uniform` is uniform *by area* -- which is the
    // right default and the wrong distribution here. The torus's own area
    // element is (R + r cos v) r, so the outer equator carries 3/2 the
    // area of the crown and therefore 3/2 the sprinkles, pushing them
    // toward the rim exactly where a real donut has fewest.
    //
    // Nesting is how a density gradient is built out of a uniform
    // sampler: each band adds on top of the wider ones under it, so the
    // crown gets all three and the skirt only the widest. That is
    // composition rather than a weight in the sampler, and it keeps the
    // one thing the sampler is good at.
    //
    // It also survives projection, which is the reason the two-band
    // version still read as edge-heavy: near the silhouette the surface
    // turns away from the camera, so the same density lands in fewer
    // pixels and looks denser. The cure is to actually have more on top,
    // not to argue with the perspective.
    constexpr double pi_ = std::numbers::pi;
    auto band = [&](double v0, double v1, std::size_t u_steps, std::size_t v_steps) {
        return scene.offset_shell(scene.space(torus_cap(2.0, 1.0, v0, v1), u_steps, v_steps),
                                  dough_bump, bd::EdgeRule::ZeroThickness);
    };
    auto crown = band(pi_ * 0.40, pi_ * 0.60, 120, 32);
    auto upper = band(pi_ * 0.25, pi_ * 0.75, 120, 40);

    std::vector<bd::Handle<double>> sprinkle_groups;
    std::size_t per_group = sprinkle_count / sprinkle_colors.size();
    for (std::size_t g = 0; g < sprinkle_colors.size(); ++g) {
        auto mat = Material<double>{.base_color = sprinkle_colors[g], .roughness = 0.35};
        auto sow = [&](bd::Handle<double> target, std::size_t count, std::uint32_t salt) {
            // seat 0.45: pressed into the glaze rather than perched on
            // it. Still not the same thing as the glaze closing around
            // them -- that needs a field that can read the scatter's own
            // sites, and it is deliberately not in this release.
            sprinkle_groups.push_back(
                scene.scatter(sprinkle, target, count,
                              static_cast<std::uint32_t>(g * 97 + salt), 0.45,
                              bd::SeatAxis::X)
                    .colored(mat)
                    .moving(grow_scale())); // scatter() places against the icing's *analytic*,
                                          // always-full-size surface (resolve_surface() does not
                                          // see .moving()) -- which is what keeps sprinkles in
                                          // sync with the growing donut
        };
        sow(crown, (per_group * 2) / 5, 11);
        sow(upper, (per_group * 2) / 5, 53);
        sow(icing, per_group / 5, 29);
    }

    // Something for the shadows to land on. Added after measuring what
    // shadows were worth without it: 0.3% of the frame's pixels, for 58%
    // more time. A shadow needs a receiver, and with the donut floating
    // in front of a flat sky the only receivers were the donut itself
    // and the dust. This is the cheapest object in the scene -- twelve
    // triangles -- and it is what makes the rest of the lighting legible.
    //
    // Through the DSL like everything else here, which is the point of
    // the file: `cube()` is a Literal, `.moving()` a translation, and the
    // table is a node in the same trace the donut is.
    // Turned to face the camera edge-on rather than corner-on. Square by
    // construction and seen from an angle, its far *corner* otherwise
    // rises into the middle of the frame as a brown peak -- harmless
    // behind the donut, and directly behind the BOOM letterforms in the
    // build-up, where it costs the text its sky. Yawing the table to the
    // camera's own azimuth puts a straight horizon there instead.
    //
    // Written with rotated(), which is also the first use of the new
    // operation in this file: a rotation is affine in the point, so the
    // table stays a placement and keeps whatever a renderer can do with
    // one.
    constexpr double table_yaw = 0.873;   // ~50 deg; camera sits at atan2(-4.2, 5)
    auto table = scene.cube({6.5, 6.5, 0.05})
                     .colored(Material<double>{.base_color = {0.80, 0.76, 0.70},
                                               .roughness = 0.88})
                     .moving(rotated(bd::VecField<double>::point(),
                                     bd::ScalarField<double>{0.0},
                                     bd::ScalarField<double>{0.0},
                                     bd::ScalarField<double>{table_yaw})
                             + bd::VecField<double>::constant({0.0, 0.0, -1.10}));

    std::vector<bd::Handle<double>> scene_children{table, cube, dough, icing};
    scene_children.insert(scene_children.end(), sprinkle_groups.begin(), sprinkle_groups.end());
    scene_children.push_back(dust);   // one node now, not 35 200
    return scene.compose(scene_children).index;
}

} // namespace donut
