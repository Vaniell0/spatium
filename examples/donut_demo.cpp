// Getting started with Spatium, in the tradition of every 3D beginner's
// first project: a donut. Blender's own famously-loved intro tutorial
// opens with the default cube -- delete it, then build the donut from
// nothing. This is the same idea for spatium::io::build's declarative
// scene DSL -- each step below is one line, and there is no for-loop
// anywhere in this file. See docs/getting-started-dsl.md for the same
// steps as prose.
//
// Step 0: the default cube, then delete it -- here, "delete" means
//         explode: every vertex's displacement is a closed-form function
//         of (its own position, time), built from spatium::algebra::
//         PerlinNoise -- no per-frame simulation state, no velocity
//         accumulated across steps. At t=0 it's a stationary cube; ask
//         for a later t and the same formula places the fragments
//         further out. That's the point of keeping motion declarative:
//         "explode outward with turbulence" is one expression, evaluated
//         wherever you need it, not a loop you have to run.
// Step 1: the dough is a torus -- a real spatium::ParametricSurface,
//         not a mesh.
// Step 2: the icing is an *offset* of the dough's own surface --
//         icing(u,v) = dough(u,v) + thickness * dough.normal_at(u,v) --
//         a new analytic surface, composed from a function, not a
//         separate shape draped/projected onto the first (that was
//         tried and was redundant: a surface built to be projected onto
//         another almost-identical surface teaches nothing offset()
//         doesn't already do directly). It goes through offset_shell()
//         rather than offset(), because real icing covers the top and
//         stops: its base is a band cut out of the torus, so it has a
//         rim, and EdgeRule says what happens there. The dough, being a
//         whole torus with no edge anywhere, uses plain offset().
// Step 3: the sprinkles are small cylinders, scattered across the icing
//         by spatium::sample_surface_uniform -- rejection sampling
//         weighted by the icing's own first-fundamental-form area
//         element, so coverage is even *by surface area*, not by
//         parameter or texture noise. No mesh exists anywhere in this
//         pipeline until the very last step, and only because a
//         renderer needs triangles -- the scene itself (`scene`, a
//         build::Trace) stays a flat, inspectable record of these three
//         operations the whole time; the loop over `scene.node(i)`
//         below is walking that record, not building geometry.
//
// Two output modes:
//   (default)     print the trace and materialized stats to the console
//   --photo PATH  also CPU-raytrace it to a PNG (BVH<Triangle3>, same
//                 engine as every other offline demo in this tree --
//                 tessellation happens here, only for display, not as
//                 part of the scene's own representation)

#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "io_helpers.hpp"

#include <spatium/algebra/noise.hpp>
#include <spatium/geometry/triangle.hpp>
#include <spatium/io/build.hpp>
#include <spatium/render/camera.hpp>
#include <spatium/render/parallel_for_rows.hpp>
#include <spatium/render/supersample.hpp>
#include <spatium/render/write_image.hpp>
#include <spatium/spatial/bvh.hpp>

#include <spatium/mesh/primitives.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <fstream>
#include <functional>
#include <numbers>
#include <optional>
#include <cstring>
#include <map>
#include <set>
#include <memory>
#include <random>
#include <utility>
#include <print>
#include <string>
#include <string_view>
#include <vector>

using namespace spatium;
namespace bd = spatium::io::build;
using spatium::io::Material;
using spatium::geometry::Ray;
using spatium::geometry::Triangle3;
using spatium::geometry::Box;
using spatium::geometry::Instanced;
using spatium::spatial::BVH;
using spatium::render::Camera;
using spatium::render::make_camera_basis;
using spatium::render::camera_ray_dir;
using spatium::render::parallel_for_rows;
using spatium::render::supersample_pixel_hdr;
using spatium::render::write_png_rgb;

namespace {

// Build-up timeline, in the same "t" materialize() already takes --
// shared by the dust particles' motion and the donut's grow-in, so
// there's one place that defines "when does what happen."
constexpr double T_EXPLODE = 0.9;   // cube -> dust, flying outward
constexpr double T_CONVERGE = 2.0;  // dust flies from its burst positions to the BOOM letterforms
constexpr double T_HOLD = 2.5;      // BOOM holds
constexpr double T_DISSOLVE = 3.2;  // dust shrinks to nothing
constexpr double T_DONUT_START = 2.9, T_DONUT_END = 3.9; // donut grows in, overlapping the dissolve
constexpr double T_BUILD_END = T_DONUT_END;

double smoothstep01(double t) {
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
Vec<double, 3> dust_core(double t, const Vec<double, 3>& burst_dir, double burst_dist,
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
double dust_shrink(double t) {
    if (t < T_HOLD) return 1.0;
    return 1.0 - smoothstep01((t - T_HOLD) / (T_DISSOLVE - T_HOLD));
}

// The two halves together, as one point map. Kept because it is the
// clearest statement of what the motion *is* -- and it is exactly the
// form that cannot be instanced, since a reader cannot tell a rigid
// placement from a deformation through a closure. The scene builds the
// split form instead; this stays as the reference the split is checked
// against.
Vec<double, 3> particle_motion(const Vec<double, 3>& local_p, double t,
                                const Vec<double, 3>& burst_dir, double burst_dist,
                                const Vec<double, 3>& target, double pull_strength, double swirl_seed,
                                const algebra::PerlinNoise& swirl) {
    Vec<double, 3> pos = dust_core(t, burst_dir, burst_dist, target, pull_strength, swirl_seed, swirl);
    return Vec<double, 3>{pos + local_p * dust_shrink(t)};
}

std::vector<std::pair<double, double>> load_boom_points(const std::string& path) {
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
ParametricSurface<double> torus_cap(double major_r, double minor_r,
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

// A tiny round dust speck. mesh::uv_sphere_mesh() was the first attempt
// here and reproducibly stalled BVH::build() for whole minutes on
// specific frames of the dust cloud -- its poles are `slices` triangles
// all sharing one exact vertex (near-zero-area at low slice counts),
// and enough of those degenerate triangles clustered near the origin
// early in the explosion looks like what a naive SAH split can pathologically
// stall on. mesh::icosahedron() has no poles at all (12 uniform
// vertices, 20 uniform-ish triangles, a regular solid) -- same round
// silhouette at this size, none of the degenerate geometry.
mesh::Mesh<Euclidean<3, double>> dust_speck(double radius) {
    auto ico = mesh::icosahedron(Sphere<2, double>{radius});
    mesh::Mesh<Euclidean<3, double>> m;
    m.vertices.assign(ico.vertices.begin(), ico.vertices.end());
    m.faces = ico.faces;
    return m;
}

// A solid, capped cylinder -- spaces/parametric.hpp's make_cylinder() is
// the *lateral surface only* (an open tube, no end disks -- fine for a
// ray-traced pipe, wrong for a sprinkle that should read as solid).
// Built directly as a mesh, not a ParametricSurface: a cylinder's two
// flat end caps don't fit one (u,v) formula seamlessly. Centered at the
// origin along z -- scatter() maps local z to the target's normal, so
// half the cylinder sits below the placement point (embedded in the
// surface) and half above, by construction.
mesh::Mesh<Euclidean<3, double>> solid_cylinder(double radius, double height, int sides = 8) {
    mesh::Mesh<Euclidean<3, double>> m;
    constexpr double pi = std::numbers::pi;
    double half = height * 0.5;
    std::uint32_t top_center = 0, bottom_center = 1;
    m.vertices.push_back({0.0, 0.0, half});
    m.vertices.push_back({0.0, 0.0, -half});
    std::vector<std::uint32_t> top_ring, bottom_ring;
    for (int i = 0; i < sides; ++i) {
        double a = 2.0 * pi * static_cast<double>(i) / static_cast<double>(sides);
        double x = radius * std::cos(a), y = radius * std::sin(a);
        top_ring.push_back(static_cast<std::uint32_t>(m.vertices.size()));
        m.vertices.push_back({x, y, half});
        bottom_ring.push_back(static_cast<std::uint32_t>(m.vertices.size()));
        m.vertices.push_back({x, y, -half});
    }
    for (int i = 0; i < sides; ++i) {
        int j = (i + 1) % sides;
        m.faces.push_back({top_center, top_ring[static_cast<std::size_t>(j)], top_ring[static_cast<std::size_t>(i)]});
        m.faces.push_back({bottom_center, bottom_ring[static_cast<std::size_t>(i)], bottom_ring[static_cast<std::size_t>(j)]});
        m.faces.push_back({top_ring[static_cast<std::size_t>(i)], top_ring[static_cast<std::size_t>(j)], bottom_ring[static_cast<std::size_t>(i)]});
        m.faces.push_back({top_ring[static_cast<std::size_t>(j)], bottom_ring[static_cast<std::size_t>(j)], bottom_ring[static_cast<std::size_t>(i)]});
    }
    return m;
}

// Blender-viewport-style coordinate axes (X red, Y green, Z blue, the
// same convention Blender's own gizmo uses) -- three thin solid
// cylinders through the origin. Shown during the build-up (cube ->
// donut), hidden once the donut is finished and the camera starts to
// orbit, so the final hero shot stays clean.
// The arms live in their own small Trace rather than being assembled as
// loose meshes: a Placed is a view onto a node, so anything that wants to
// be a scene object has to be one. Held in a function-local static so the
// trace outlives every view handed out, and because the axes are constant
// geometry there is nothing to rebuild per frame.
// The viewport: coloured axes and a floor grid, shown while the scene is
// being built and gone once it is finished.
//
// The grid is the part that makes the joke land. This demo's whole
// premise is Blender's donut tutorial without Blender, and what a reader
// recognises from that tutorial is not a table -- it is the default cube
// sitting on a grey grid, and then the switch to a render where the grid
// is gone. The build-up frames are the viewport; `--photo` is the render,
// and it never calls this function. Nothing had to be written to say so.
//
// Lives in its own little Trace, outside the scene's, because these are
// not objects in the scene: they are the editor, and a scene that had to
// carry its own gizmos would be describing the tool rather than the
// subject.
const std::vector<bd::Placed<double>>& axes_placed(double length = 3.5, double radius = 0.012) {
    static bd::Trace<double> axes;
    static std::vector<bd::Placed<double>> out = [&] {
        auto arm = [&](std::size_t a, std::size_t b, std::size_t c,
                       const Vec<double, 3>& color) {
            auto m = solid_cylinder(radius, length, 10);
            // The cylinder is built along z; permute its axes to aim it.
            for (auto& v : m.vertices) v = Vec<double, 3>{v[a], v[b], v[c]};
            axes.literal(m).colored(
                Material<double>{.base_color = color, .roughness = 0.6});
        };
        arm(2, 0, 1, {0.85, 0.15, 0.15}); // z-cylinder -> x
        arm(0, 2, 1, {0.15, 0.75, 0.15}); // z-cylinder -> y
        arm(0, 1, 2, {0.15, 0.25, 0.85}); // already along z

        // Unit grid on the floor, axis-aligned in world space -- which is
        // the point of a grid, and the reason it is not turned to match
        // the table under it. Held to +-4.5 so it stays inside the yawed
        // table's inscribed square (6.5 / sqrt(2) is about 4.6) instead of
        // ending in mid-air.
        constexpr double half = 4.5, step = 1.0, w = 0.008;
        constexpr double z = -1.043;   // just clear of the table top at -1.05
        for (int i = -static_cast<int>(half / step); i <= static_cast<int>(half / step); ++i) {
            double c = i * step;
            bool axis = (i == 0);   // the two centre lines read brighter, as in a viewport
            auto line = Material<double>{.base_color = axis ? Vec<double, 3>{0.42, 0.40, 0.38}
                                                            : Vec<double, 3>{0.62, 0.59, 0.55},
                                         .roughness = 1.0};
            auto bar = [&](Vec<double, 3> half_extents, Vec<double, 3> at) {
                auto m = mesh::box_mesh<double>(half_extents);
                for (auto& v : m.vertices) v = Vec<double, 3>{v + at};
                axes.literal(m).colored(line);
            };
            bar({half, w, 0.004}, {0.0, c, z});
            bar({w, half, 0.004}, {c, 0.0, z});
        }

        std::vector<bd::Placed<double>> v;
        for (std::size_t i = 0; i < axes.size(); ++i)
            v.push_back(bd::Placed<double>{&axes, i, 0.0});
        return v;
    }();
    return out;
}

// Smooth per-vertex normals (area-weighted face-normal average) -- flat
// per-triangle normals on a coarse tessellation read as faceted plastic;
// this is the standard fix, needs nothing beyond the mesh's own topology.
std::vector<Vec<double, 3>> smooth_normals(const mesh::Mesh<Euclidean<3, double>>& m) {
    std::vector<Vec<double, 3>> n(m.vertex_count(), Vec<double, 3>{0, 0, 0});
    for (const auto& f : m.faces) {
        auto& a = m.vertices[f[0]];
        auto& b = m.vertices[f[1]];
        auto& c = m.vertices[f[2]];
        Vec<double, 3> face_n{Vec<double, 3>{b - a}.cross(Vec<double, 3>{c - a})}; // area-weighted (unnormalized)
        n[f[0]] = Vec<double, 3>{n[f[0]] + face_n};
        n[f[1]] = Vec<double, 3>{n[f[1]] + face_n};
        n[f[2]] = Vec<double, 3>{n[f[2]] + face_n};
    }
    for (auto& v : n) {
        double len = v.norm();
        v = len > 1e-12 ? Vec<double, 3>{v / len} : Vec<double, 3>{0, 0, 1};
    }
    return n;
}

struct Prim {
    Triangle3 tri;
    std::array<Vec<double, 3>, 3> vertex_normals;
    Vec<double, 3> color;
    double roughness;
    Vec<double, 3> emissive;
};

// One instanced object: which shared shape, where, and what colour.
struct DustInstance {
    Vec<double, 3> color;
    double roughness;
    Vec<double, 3> emissive;
    double opacity;
};

// The halo around something that emits, which is what actually reads as
// light. Brightness alone does not: a clamped hot core and a merely
// bright one arrive at the same pixel value, and the eye reads both as
// paint. What separates them is that light spreads past the silhouette
// of the thing emitting it, and nothing else in a frame does that.
//
// The source is the emission buffer, never the colour buffer, so this is
// a mask and not a threshold -- see `Traced` in render_frame for the
// measurements that ruled a threshold out. The blur is separable, two
// passes of a 1D Gaussian instead of one 2D pass, which is the
// difference between 2*r and r*r samples per pixel; at sigma 7 that is
// 43 against 441. Even so it costs nothing worth measuring next to the
// ray casting that produced the buffer.
void add_bloom(std::vector<Vec<double, 3>>& hdr, const std::vector<Vec<double, 3>>& emission,
               int W, int H) {
    constexpr double kSigma = 7.0;       // halo width in pixels at 960x720
    constexpr double kStrength = 1.25;   // how much of the halo is added back
    const int radius = static_cast<int>(std::ceil(3.0 * kSigma));

    std::vector<double> kernel(static_cast<std::size_t>(radius) + 1);
    double norm = 0.0;
    for (int i = 0; i <= radius; ++i) {
        kernel[static_cast<std::size_t>(i)] = std::exp(-0.5 * (i * i) / (kSigma * kSigma));
        norm += (i == 0 ? 1.0 : 2.0) * kernel[static_cast<std::size_t>(i)];
    }
    for (auto& k : kernel) k /= norm;

    const auto at = [W](int x, int y) {
        return static_cast<std::size_t>(y) * static_cast<std::size_t>(W) +
               static_cast<std::size_t>(x);
    };
    std::vector<Vec<double, 3>> tmp(emission.size()), blurred(emission.size());

    parallel_for_rows(H, [&](int y) {
        for (int x = 0; x < W; ++x) {
            Vec<double, 3> acc{};
            for (int d = -radius; d <= radius; ++d) {
                const int xx = std::clamp(x + d, 0, W - 1);
                acc = Vec<double, 3>{acc + emission[at(xx, y)] *
                                               kernel[static_cast<std::size_t>(std::abs(d))]};
            }
            tmp[at(x, y)] = acc;
        }
    });
    parallel_for_rows(H, [&](int y) {
        for (int x = 0; x < W; ++x) {
            Vec<double, 3> acc{};
            for (int d = -radius; d <= radius; ++d) {
                const int yy = std::clamp(y + d, 0, H - 1);
                acc = Vec<double, 3>{acc + tmp[at(x, yy)] *
                                               kernel[static_cast<std::size_t>(std::abs(d))]};
            }
            blurred[at(x, y)] = acc;
        }
    });

    for (std::size_t i = 0; i < hdr.size(); ++i)
        hdr[i] = Vec<double, 3>{hdr[i] + blurred[i] * kStrength};
}

std::vector<std::uint8_t> render_frame(const bd::Trace<double>& trace,
                                        const bd::Cooked<double>& cooked,
                                        const std::vector<bd::Placed<double>>& gizmos,
                                        const Camera<double>& cam, int W, int H) {
    // Two trees, because the scene has two kinds of object in it.
    //
    // Most of it is ordinary geometry: a torus, a shell, sprinkles --
    // each a handful of triangles, each different, and a triangle tree is
    // the right home for them.
    //
    // The dust is not that. It is thousands of copies of one little cube,
    // and flattening it puts a quarter of a million identical triangles
    // into the tree, where every one of them costs a leaf. An instance
    // holds a *reference* to the one shape plus where this copy sits, so
    // the tree holds thousands of leaves instead of hundreds of
    // thousands -- and for a cube that never rotates the shape is a `Box`,
    // which is not an approximation of it but literally it, tested in six
    // slab comparisons rather than twelve Möller-Trumbore.
    //
    // What makes an object eligible is not that it looks small. It is
    // that its motion is a *placement* -- affine in the point, so it moves
    // the object without reshaping it (see VecField::is_placement). A
    // motion that deforms per vertex cannot share geometry with anything,
    // and those objects go through the triangle path unchanged.
    std::vector<Triangle3> tris;
    std::vector<Prim> prims;
    std::vector<Instanced<geometry::BoundedQuadric<double>>> insts;
    std::vector<DustInstance> inst_info;

    // One quadric per distinct *shape*, not per object -- which is the
    // whole point of reading a cooked scene rather than a list of
    // materialized nodes. Sized up front and never grown, because every
    // instance below holds a pointer into it.
    std::vector<geometry::BoundedQuadric<double>> quadrics(cooked.shape_count());
    std::vector<char> has_quadric(cooked.shape_count(), 0);
    std::vector<std::vector<Vec<double, 3>>> shape_normals(cooked.shape_count());
    for (std::size_t i = 0; i < cooked.shape_count(); ++i) {
        const auto& sh = cooked.shapes()[i];
        if (const auto* q = std::any_cast<geometry::BoundedQuadric<double>>(&sh.exact)) {
            quadrics[i] = *q;
            has_quadric[i] = 1;
        }
    }

    // A shape's geometry is its *rest* form when the object placing it is
    // instanceable, and its *placed* form when the object deforms -- in
    // which case cook() hands back an identity transform. So one formula
    // covers both, and the deforming case is not a special case here.
    // Takes the rotation already expanded, rather than reaching into the
    // object for it: this runs per *vertex*, and a quaternion unpacked
    // here would spend the arithmetic the compact storage was supposed to
    // be paying for. The expansion happens once per object, below.
    auto to_world = [](const Vec<double, 3>& v, const Matrix<double, 3, 3>& R,
                       const bd::Object<double>& o) {
        return Vec<double, 3>{R * Vec<double, 3>{v * o.scale} + o.translation};
    };

    auto emit_triangles = [&](const mesh::Mesh<Euclidean<3, double>>& m,
                              const std::vector<Vec<double, 3>>& vn,
                              const bd::Object<double>& o, const Matrix<double, 3, 3>& R,
                              const Material<double>& mat) {
        for (const auto& f : m.faces) {
            Triangle3 t(to_world(m.vertices[f[0]], R, o), to_world(m.vertices[f[1]], R, o),
                        to_world(m.vertices[f[2]], R, o));
            tris.push_back(t);
            // Normals turn by the rotation alone: the scale is uniform, so
            // it divides out of the inverse transpose.
            prims.push_back(Prim{t,
                                 {Vec<double, 3>{R * vn[f[0]]},
                                  Vec<double, 3>{R * vn[f[1]]},
                                  Vec<double, 3>{R * vn[f[2]]}},
                                 mat.base_color, mat.roughness, mat.emissive});
        }
    };

    for (const auto& obj : cooked.objects()) {
        const auto& mat = obj.material;
        // Once per object. Everything below uses this, including the
        // instance, whose rotation stays a matrix precisely so the
        // traversal never has to unpack anything.
        const Matrix<double, 3, 3> R = obj.rotation_q.to_matrix();

        // Two conditions, and neither is a guess about what the object
        // looks like: the shape must carry a closed form, and the object's
        // motion must be a placement so copies differ only by where they
        // are.
        if (obj.instanceable && has_quadric[obj.shape]) {
            insts.push_back({&quadrics[obj.shape], obj.translation, obj.scale, R});
            inst_info.push_back({mat.base_color, mat.roughness, mat.emissive, mat.opacity});
            continue;
        }

        // Smooth normals once per shape rather than once per object: a
        // thousand instances of one sprinkle used to recompute the same
        // normals a thousand times, through a merged mesh that had to be
        // built first.
        auto& vn = shape_normals[obj.shape];
        if (vn.empty()) vn = smooth_normals(cooked.shapes()[obj.shape].geometry);
        emit_triangles(cooked.shapes()[obj.shape].geometry, vn, obj, R, mat);
    }

    // The gizmos are not in the scene and must not be: they are the
    // editor, not the subject. They arrive as Placed views onto their own
    // little Trace, and they are drawn here because forgetting them is the
    // easy mistake in this refactor -- the build-up frames simply lose
    // their axes and grid and nothing else changes.
    for (const auto& g : gizmos) {
        auto m = g.mesh();
        auto vn = smooth_normals(m);
        auto mat = g.material();
        for (const auto& f : m.faces) {
            Triangle3 t(m.vertices[f[0]], m.vertices[f[1]], m.vertices[f[2]]);
            tris.push_back(t);
            prims.push_back(Prim{t, {vn[f[0]], vn[f[1]], vn[f[2]]}, mat.base_color, mat.roughness,
                                 mat.emissive});
        }
    }

    // Counted before the move, because a moved-from vector is empty and
    // the report below would have quietly started printing zeros. Found
    // by checking what still reads these after this line rather than by
    // assuming nothing did.
    const std::size_t n_inst = insts.size(), n_tris = tris.size();

    // Moved, not copied. `build` takes its shapes by value and keeps
    // them, so handing it an lvalue leaves two copies of the array alive
    // for the rest of the frame -- 112 bytes per instance, which is 214 MB
    // at two million and was simply being spent.
    auto bvh = BVH<Triangle3>::build(std::move(tris));
    auto inst_bvh = BVH<Instanced<geometry::BoundedQuadric<double>>>::build(std::move(insts));

    // Reported rather than assumed, and the last pair is the whole reason
    // a cooked scene exists: what a renderer would have held if every
    // object carried its own copy of its geometry, against what the shape
    // table actually holds.
    std::println("  scene: {} objects -> {} instances + {} triangles; "
                 "{} refused (motion deforms), {} shared; "
                 "vertices {} without instancing, {} stored ({:.1f}x)",
                 cooked.object_count(), n_inst, n_tris,
                 cooked.opaque_refused(), cooked.shared_objects(),
                 cooked.vertices_without_instancing(), cooked.vertices_stored(),
                 cooked.vertices_stored() == 0
                     ? 0.0
                     : static_cast<double>(cooked.vertices_without_instancing()) /
                           static_cast<double>(cooked.vertices_stored()));
    {
        // Which nodes were refused, by kind -- because "21 deformations"
        // is a number nobody can act on, and because a claim about what
        // they are has already been wrong once in this repository.
        std::map<std::string, int> by_kind;
        for (auto i : cooked.refused_nodes())
            ++by_kind[bd::kind_name(trace.node(i).kind)];
        // The worst-case world radius of any object, which is what a
        // rotation's error actually gets multiplied by. A matrix error is
        // dimensionless; a vertex displacement is not.
        double worst_radius = 0.0;
        for (const auto& o : cooked.objects()) {
            double rest = 0.0;
            for (const auto& v : cooked.shapes()[o.shape].geometry.vertices)
                rest = std::max(rest, Vec<double, 3>{v}.norm());
            worst_radius = std::max(worst_radius, rest * std::abs(o.scale) + o.translation.norm());
        }
        std::println("  worst object radius {:.3f} world units", worst_radius);

        // How many *distinct* materials the cooked scene actually holds.
        // The question a palette answers is whether a per-object Material
        // collapses, and the answer is a property of the scene rather
        // than of the idea: comparison is exact, so two colours a single
        // ulp apart are two entries.
        std::set<std::array<double, 8>> palette;
        for (const auto& o : cooked.objects())
            palette.insert({o.material.base_color[0], o.material.base_color[1],
                            o.material.base_color[2], o.material.roughness,
                            o.material.emissive[0], o.material.emissive[1],
                            o.material.emissive[2], o.material.opacity});
        std::println("  distinct materials {} over {} objects", palette.size(),
                     cooked.object_count());
        std::size_t glowing = 0, near_target = 0;
        double max_glow = 0.0;
        for (const auto& o : cooked.objects()) {
            max_glow = std::max(max_glow, o.material.emissive[0]);
            if (o.material.emissive[0] > 0.3) ++glowing;
            if (o.material.base_color[0] > 0.8) ++near_target;
        }
        std::println("  glowing {} , hot-coloured {} , max emissive {:.3f}", glowing, near_target,
                     max_glow);

        std::print("  refused nodes:");
        for (const auto& [k, n] : by_kind) std::print(" {}x{}", n, k);
        std::println("");
    }

    // The last step stays a clamp, and that is a measured decision rather
    // than the one this started as.
    //
    // The complaint it began with was right: the emission maxes out above
    // 1 and clamping flattens every hot particle onto the same orange, so
    // it reads as pigment. The obvious repair -- a tone curve compressing
    // what is above 1 instead of cutting it -- was written, and then it
    // was checked the only way that means anything: render with the
    // emission switched off, before and after. It moved.
    //
    // The reason it moved is that this frame has no headroom. A curve that
    // asymptotes to 1 must put its knee below 1, and the scene's own
    // shading reaches 1.263 -- the icing's specular once the donut is
    // fully grown. So every knee sits inside the range the scene actually
    // occupies, and there is no neutral tone curve for this image at all,
    // only curves that darken the icing by a little or by a lot.
    //
    // Which turned out not to matter, because brightness was never what
    // was missing. What separates a light from a paint is that light
    // spreads past the silhouette of what emits it, and that is add_bloom
    // above, working off the emission buffer rather than off any pixel
    // value. With the halo there, the flattening of the few brightest
    // cores costs nothing visible, and the rest of the frame is left
    // exactly as it was -- byte for byte, which is asserted by rendering
    // a frame with no emission in it and comparing.
    const Vec<double, 3> background{0.55, 0.75, 0.92}; // plain light blue, no starfield
    const auto basis = make_camera_basis(cam);
    // Two lights, and the split is the whole reason the shadow reads.
    //
    // The key was {0.55, -0.4, 1.0}: dominated by +z and pointing back
    // toward the camera, so every visible surface was lit and the shadow
    // fell straight down, hidden under the donut by the donut. Moving it
    // sideways and lowering it puts the shadow on the table where the
    // camera can see it -- the geometry had been casting one all along.
    //
    // That alone left the near-left side in shade, because a single
    // off-axis light means half the object faces away from the only light
    // there is. The fill is what a photographer would reach for: a weaker
    // source from the camera's own side, **casting no shadow**, which
    // recovers the shaded half without touching the shadow the key throws.
    // Skipping the shadow ray for it is not a shortcut -- a fill that cast
    // its own shadow would put a second, contradictory one on the table.
    const Vec<double, 3> light = Vec<double, 3>{Vec<double, 3>{0.85, 0.45, 0.55}.normalized()};
    const Vec<double, 3> fill = Vec<double, 3>{Vec<double, 3>{0.62, -0.70, 0.35}.normalized()};
    constexpr double FILL_STRENGTH = 0.38;
    constexpr int MAX_DEPTH = 3;

    // Is anything between this point and the light? The light is
    // directional -- a unit vector, the sun rather than a bulb -- so the
    // shadow ray has no far limit and any hit at all occludes.
    //
    // Both trees are asked, for the same reason the primary ray asks
    // both: a shadow that ignored the dust would be a shadow of half the
    // scene. This is the change that stops objects floating. Its most
    // visible consequence is not the big one you would expect -- with no
    // ground plane the donut casts onto itself and onto the dust -- it is
    // the sprinkles, each of which now sits in a small dark patch on the
    // icing instead of appearing painted onto it.
    //
    // The offset along the normal is the classic shadow-acne guard: a
    // point on a surface, tested against that same surface, hits itself
    // at t ~ 0 without it and every lit pixel comes back shadowed.
    auto occluded = [&](const Vec<double, 3>& point, const Vec<double, 3>& n) {
        Ray<3, double> shadow{Vec<double, 3>{point + n * 1e-4}, light};
        if (bvh.ray_cast(shadow)) return true;
        return static_cast<bool>(inst_bvh.ray_cast(shadow));
    };

    // What a ray saw, kept as two values rather than one: the colour to
    // show, and the part of that colour which came from something
    // emitting.
    //
    // The split is what lets the bloom below be selective, and it exists
    // because the obvious alternative was measured and does not work. A
    // brightness threshold cannot separate fire from highlight here: the
    // specular term is bounded by (1 - roughness) * 0.6, so a white
    // sprinkle at roughness 0.35 is bounded at 0.98 + 0.39 = 1.37, and
    // the icing at roughness 0.40 mixes in 0.168 of whatever it reflects,
    // which reaches ~1.48 when what it reflects is the fire. The fire's
    // own range starts around 1.63. Those overlap, so every threshold
    // that catches the fire also catches a highlight -- and even if one
    // did fit, it would be a number that silently stops working the next
    // time a material, a light or the camera moves.
    //
    // Carrying the emission instead makes the question structural rather
    // than numeric: a surface whose material has no emissive term
    // contributes exactly zero to it, at any brightness, under any light,
    // for ever. Nothing has to be re-measured when the scene changes.
    struct Traced {
        Vec<double, 3> color{};
        Vec<double, 3> emission{};
    };

    // trace_ray returns linear [0,1] color throughout -- the single
    // conversion back to [0,255] happens exactly once, at the very end
    // of render_frame below, so reflection blending stays on one scale.
    std::function<Traced(const Ray<3, double>&, int)> trace_ray =
        [&](const Ray<3, double>& ray, int depth) -> Traced {
        // Both trees, nearer hit wins. Two structures rather than one
        // because the two kinds of object want different leaves, not
        // because the renderer is special-casing the dust: each tree is
        // asked the same question and the answers are compared by `t`.
        auto hit = bvh.ray_cast(ray);
        auto dhit = inst_bvh.ray_cast(ray);

        Vec<double, 3> n{};
        Vec<double, 3> base_color{};
        double roughness = 1.0;
        Vec<double, 3> hit_point{};
        Vec<double, 3> emissive{};
        double opacity = 1.0;

        const bool dust_won = dhit && (!hit || dhit->t < hit->t);
        if (!hit && !dhit) return Traced{background, Vec<double, 3>{}};

        if (dust_won) {
            // The quadric's own normal at the hit, exact rather than
            // interpolated across a facet -- and the instance transform
            // is a translation and a uniform scale, so it comes back
            // needing no correction.
            n = dhit->normal;
            const auto& info = inst_info[dhit->index];
            base_color = info.color;
            roughness = info.roughness;
            emissive = info.emissive;
            opacity = info.opacity;
            hit_point = dhit->point;
        } else {
            const auto& prim = prims[hit->index];
            // Interpolate the smooth vertex normals via the hit's own
            // barycentric weights, not the triangle's single flat normal.
            double w0 = double{1} - hit->u - hit->v, w1 = hit->u, w2 = hit->v;
            n = Vec<double, 3>{Vec<double, 3>{prim.vertex_normals[0] * w0 + prim.vertex_normals[1] * w1 +
                                              prim.vertex_normals[2] * w2}
                                   .normalized()};
            base_color = prim.color;
            roughness = prim.roughness;
            emissive = prim.emissive;
            hit_point = hit->point;
        }
        if (n.dot(ray.direction) > 0.0) n = Vec<double, 3>{-n};

        double diff = std::max(0.0, n.dot(light));
        roughness = std::clamp(roughness, 0.0, 1.0);

        // A point facing away from the light is already dark and cannot
        // be shadowed any further, so the shadow ray is skipped there --
        // which is most of the scene's back half, and the reason this
        // costs less than doubling the ray count.
        const bool shadowed = diff > 0.0 && occluded(hit_point, n);
        if (shadowed) diff = 0.0;
        diff += FILL_STRENGTH * std::max(0.0, n.dot(fill));
        diff = std::min(diff, 1.0);

        Vec<double, 3> view = Vec<double, 3>{-ray.direction};
        Vec<double, 3> half = Vec<double, 3>{Vec<double, 3>{view + light}.normalized()};
        double shininess = 4.0 + 90.0 * (1.0 - roughness); // rough: broad/dull, glossy: tight/bright
        double spec = std::pow(std::max(0.0, n.dot(half)), shininess) * (1.0 - roughness) * 0.6;
        if (shadowed) spec = 0.0;   // a highlight is the light source seen in the surface

        // The ambient term survives the shadow deliberately. A shadow
        // that goes to black is a shadow with no bounce light in it,
        // which reads as a hole cut out of the object rather than as a
        // shaded part of it.
        Vec<double, 3> local = Vec<double, 3>{base_color * (0.22 + 0.78 * diff)};
        Vec<double, 3> color = Vec<double, 3>{local + Vec<double, 3>{1.0, 1.0, 1.0} * spec};

        // Emission is added, not blended: a surface that emits does so
        // whether or not anything is lighting it, which is exactly why
        // the BOOM letterforms stay bright while the particles around
        // them fall into the donut's shadow.
        color = Vec<double, 3>{color + emissive};

        // The emitted part travels alongside the colour through every
        // blend below, with the same weights, so that a reflection of
        // the fire stays recognisable as fire and a reflection of the
        // table stays recognisable as not-fire. Blending it any other
        // way would make the two disagree about the same pixel.
        Vec<double, 3> emitted = emissive;

        // A real reflected ray for glossy materials, not just a
        // highlight -- what actually earns the word "raytracing" here.
        if (roughness < 0.6 && depth < MAX_DEPTH) {
            Vec<double, 3> refl_dir = Vec<double, 3>{ray.direction - n * (2.0 * ray.direction.dot(n))};
            Vec<double, 3> refl_origin = Vec<double, 3>{hit_point + n * 1e-4};
            Traced refl = trace_ray(Ray<3, double>{refl_origin, refl_dir}, depth + 1);
            double reflectivity = (1.0 - roughness) * 0.28;
            color = Vec<double, 3>{color * (1.0 - reflectivity) + refl.color * reflectivity};
            emitted = Vec<double, 3>{emitted * (1.0 - reflectivity) + refl.emission * reflectivity};
        }
        // Seen through, not bent through. The continuation carries on in
        // the ray's own direction, so a flake tints what is behind it
        // instead of displacing it -- see Material::opacity for why that
        // is the right model at this scale rather than a shortcut.
        if (opacity < 1.0 && depth < MAX_DEPTH) {
            Ray<3, double> through{Vec<double, 3>{hit_point + ray.direction * 1e-4},
                                   ray.direction};
            Traced behind = trace_ray(through, depth + 1);
            color = Vec<double, 3>{color * opacity + behind.color * (1.0 - opacity)};
            emitted = Vec<double, 3>{emitted * opacity + behind.emission * (1.0 - opacity)};
        }
        return Traced{color, emitted};
    };

    // Two linear buffers rather than bytes straight out of the sampler,
    // because the bloom below is a whole-image operation and cannot be
    // done one pixel at a time.
    const std::size_t npix = static_cast<std::size_t>(W) * static_cast<std::size_t>(H);
    std::vector<Vec<double, 3>> hdr(npix), emission(npix);
    parallel_for_rows(H, [&](int y) {
        for (int x = 0; x < W; ++x) {
            // The emitted part is accumulated here and divided by the
            // number of samples actually taken, rather than by the AA
            // factor: the two averages then cannot disagree even if the
            // sampler's default changes underneath this call.
            Vec<double, 3> emis_accum{};
            int taken = 0;
            auto ray_color = [&](double sx, double sy) -> Vec<double, 3> {
                Vec<double, 3> dir = camera_ray_dir(basis, sx, sy);
                Traced tr = trace_ray(Ray<3, double>{cam.position, dir}, 0);
                emis_accum = Vec<double, 3>{emis_accum + tr.emission};
                ++taken;
                // Clamped per sample, before the average, which is where
                // this has always done it. Averaging first and clamping
                // after is the more defensible order -- it is the one
                // that antialiases an overbright edge correctly -- but it
                // is a different picture, and changing it here would ride
                // in on the back of the bloom and be indistinguishable
                // from it. It can be its own change, with its own before
                // and after.
                return Vec<double, 3>{std::clamp(tr.color[0], 0.0, 1.0),
                                      std::clamp(tr.color[1], 0.0, 1.0),
                                      std::clamp(tr.color[2], 0.0, 1.0)};
            };
            const std::size_t i = static_cast<std::size_t>(y) * static_cast<std::size_t>(W) +
                                  static_cast<std::size_t>(x);
            hdr[i] = supersample_pixel_hdr(x, y, W, H, basis.tan_half,
                                           static_cast<double>(W) / H, ray_color);
            emission[i] = Vec<double, 3>{emis_accum * (1.0 / std::max(taken, 1))};
        }
    });

    add_bloom(hdr, emission, W, H);

    std::vector<std::uint8_t> img(3 * npix, 0);
    for (std::size_t i = 0; i < npix; ++i) {
        for (int ch = 0; ch < 3; ++ch)
            img[3 * i + static_cast<std::size_t>(ch)] =
                static_cast<std::uint8_t>(std::clamp(hdr[i][ch], 0.0, 1.0) * 255.0);
    }
    return img;
}

Camera<double> hero_camera() {
    return {.position = {5.0, -4.2, 3.6}, .target = {0.0, 0.0, 0.0}, .up = {0.0, 0.0, 1.0}, .fov_deg = 38.0};
}

// No gizmos here, and that is the whole difference between the still and
// the sequence: --photo is the render, the build-up frames are the
// viewport. Blender's own arc, and nothing had to be written to say so.
void render_photo(const bd::Trace<double>& trace, const bd::Cooked<double>& cooked,
                   const std::string& out_path, bool force) {
    constexpr int W = 960, H = 720;
    auto img = render_frame(trace, cooked, {}, hero_camera(), W, H);
    if (spatium::examples::confirm_overwrite(out_path, force))
        write_png_rgb(out_path, W, H, img);
    std::println("  -> {}", out_path);
}

// The lesson, as a video: build-up (cube deletes, donut grows in its
// place, axes visible like a viewport) then a clean orbit once it's
// finished (axes hidden). Two phases, same materialize() call at
// different t -- no separate animation system, the trace already
// describes motion as a function of time.
void render_video(const bd::Trace<double>& scene, std::size_t lesson_idx, const std::string& dir,
                   bool force, int build_frames, int orbit_frames) {
    namespace fs = std::filesystem;
    fs::create_directories(dir);
    constexpr int W = 960, H = 720;
    Camera<double> cam = hero_camera();
    double orbit_radius = std::sqrt(cam.position[0] * cam.position[0] + cam.position[1] * cam.position[1]);
    double orbit_start_angle = std::atan2(cam.position[1], cam.position[0]);
    int frame = 0;

    auto write_frame = [&](const bd::Cooked<double>& ck,
                           const std::vector<bd::Placed<double>>& gizmos,
                           const Camera<double>& c) {
        std::string path = std::format("{}/frame_{:04d}.png", dir, frame);
        if (spatium::examples::confirm_overwrite(path, force))
            write_png_rgb(path, W, H,
                          render_frame(scene, ck, gizmos, c, W, H));
        ++frame;
    };

    for (int i = 0; i < build_frames; ++i) {
        double t = T_BUILD_END * static_cast<double>(i) / static_cast<double>(build_frames - 1);
        write_frame(bd::cook(scene, lesson_idx, t), axes_placed(), cam);
    }
    auto full = bd::cook(scene, lesson_idx, T_BUILD_END); // finished donut, viewport hidden from here on
    for (int i = 0; i < orbit_frames; ++i) {
        double frac = static_cast<double>(i) / static_cast<double>(orbit_frames);
        double angle = orbit_start_angle + frac * 2.0 * std::numbers::pi;
        Camera<double> orbit_cam = cam;
        orbit_cam.position = {orbit_radius * std::cos(angle), orbit_radius * std::sin(angle), cam.position[2]};
        write_frame(full, {}, orbit_cam);
    }
    std::println("donut_demo: wrote {} frames -> {}/", frame, dir);
    std::println("  ffmpeg -framerate 30 -i {}/frame_%04d.png -pix_fmt yuv420p donut.mp4", dir);
}

} // namespace

int main(int argc, char* argv[]) {
    bool force = false;
    bool photo = false;
    bool video = false;
    std::string out_path = "donut.png";
    std::string video_dir = "donut_frames";
    std::size_t sprinkle_count = 11000; // small flecks, so many more of them
    std::size_t dust_particles = 35200; // one Scatter now, so this is a number rather than a loop
    double t = T_BUILD_END; // how far into the build-up to render (T_BUILD_END = fully formed)
    int build_frames = 75, orbit_frames = 45; // half the sampling density -- ~9x more
                                              // dust triangles from copies_per_point makes
                                              // full 150/90 too costly for today; render at
                                              // 15fps (still smooth) instead of 30fps

    for (int i = 1; i < argc; ++i) {
        std::string_view a = argv[i];
        if (a == "--force") { force = true; continue; }
        if (a == "--photo") { photo = true; if (i + 1 < argc && argv[i + 1][0] != '-') out_path = argv[++i]; continue; }
        if (a == "--video") { video = true; if (i + 1 < argc && argv[i + 1][0] != '-') video_dir = argv[++i]; continue; }
        if (a == "--sprinkles" && i + 1 < argc) { sprinkle_count = static_cast<std::size_t>(std::atoi(argv[++i])); continue; }
        if (a == "--dust" && i + 1 < argc) { dust_particles = static_cast<std::size_t>(std::atol(argv[++i])); continue; }
        if (a == "--frames" && i + 2 < argc) {   // build, orbit -- for previewing a sequence cheaply
            build_frames = std::atoi(argv[++i]);
            orbit_frames = std::atoi(argv[++i]);
            continue;
        }
        if (a == "--t" && i + 1 < argc) { t = std::atof(argv[++i]); continue; }
        if (a == "--help") {
            std::print("donut_demo [--photo [path]] [--video [dir]] [--sprinkles N] [--dust N] [--frames B O] [--t seconds] [--force]\n"
                       "  Spatium's getting-started demo: delete the default cube (explode it,\n"
                       "  --t controls how far), then torus dough, offset-surface icing,\n"
                       "  area-sampled sprinkles -- declarative steps, see the file's own header\n"
                       "  comment. Prints the trace to the console by default; --photo also\n"
                       "  CPU-raytraces it to a PNG; --video renders the whole build-up + orbit\n"
                       "  as a PNG frame sequence (assemble with ffmpeg, printed at the end).\n");
            return 0;
        }
        std::print(stderr, "unknown option: {}\n", a);
        return 1;
    }

    auto t0 = std::chrono::steady_clock::now();

    bd::Trace<double> scene;

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
        return scaled(bd::VecField<double>::point(), [](double time) {
            return smoothstep01((time - T_DONUT_START) / (T_DONUT_END - T_DONUT_START));
        });
    };

    // Step 0 -- the default cube. Visible briefly, static, then it's
    // gone -- what actually reads as "the cube" from here on is the
    // dust field below (`dust`), not this node's own motion.
    auto cube = scene.cube({0.9, 0.9, 0.9})
                    .colored(Material<double>{.base_color = {0.55, 0.55, 0.58}})
                    // Same reason as grow_scale: a uniform scale, spelled so
                    // that it can be recognised as one.
                    .moving(scaled(bd::VecField<double>::point(),
                                   [](double time) { return time < 0.12 ? 1.0 : 0.0; }));

    // Step 0.5 -- delete the cube by *exploding* it: not a shrink this
    // time, real dust -- ~220 tiny cubes flying from the cube's own
    // volume, converging into the shape of the word "BOOM" (points
    // rasterized from a real font via ImageMagick, not hand-placed --
    // see examples/data/boom_points.txt and the shell one-liner that
    // made it), holding briefly, then dissolving before the donut grows
    // in. Each particle is its own tiny Literal node with its own
    // `.moving()` closure -- the DSL's per-node motion slot, just used
    // 220 times instead of once; still no simulation state anywhere.
    auto boom_uv = load_boom_points("examples/data/boom_points.txt");
    if (boom_uv.empty()) boom_uv = load_boom_points("data/boom_points.txt");
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
                            + rotated(scaled(bd::VecField<double>::point(),
                                             [](double time) { return dust_shrink(time); }),
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
    algebra::PerlinNoise surface_noise(2);
    auto dough_bump = bd::ScalarField<double>{[surface_noise](double u, double v) {
        return 0.020 * surface_noise(u * 3.0, v * 3.0, 0.0); // visible but still broad, not fine-grain
    }};
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
    auto dough_surface = bd::ScalarField<double>{[surface_noise, pore_noise](double u, double v) {
        double broad = 0.020 * surface_noise(u * 3.0, v * 3.0, 0.0);
        // u runs around the major circle and v around the tube, so a
        // radian of u covers ~2.5x the arc a radian of v does; the
        // frequencies are in that ratio so the pits come out round
        // rather than smeared along the ring.
        double n = pore_noise(u * 11.0, v * 4.5, 0.0);
        double pit = std::max(0.0, -n - 0.18);
        return broad - 0.075 * pit * pit;
    }};
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
    auto icing = scene.offset_shell(icing_base, bd::ScalarField<double>{[icing_noise](double u, double v) {
                        constexpr double pi = std::numbers::pi;
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
                        double outer = v - pi * 0.02;      // toward the outer equator
                        double inner = icing_rim - v;      // toward the hole
                        double wobble = icing_noise(std::cos(u) * 2.0, std::sin(u) * 2.0, 0.0) * (pi * 0.03);
                        double drip = std::max(0.0, icing_noise(std::cos(u) * 1.3, std::sin(u) * 1.3, 8.0) - 0.35) * (pi * 0.35);
                        double edge_dist = std::min(outer + drip, inner);
                        double falloff = std::clamp((edge_dist + wobble) / (pi * 0.09), 0.0, 1.0);
                        falloff = falloff * falloff * (3.0 - 2.0 * falloff); // soft, not torn
                        // 0.045 until 2026-09-17, tuned when offset()
                        // silently rendered this at 48x24 and the waves
                        // were under-sampled into near-smoothness. With
                        // the base's own 160x64 actually reaching the
                        // mesh they resolve fully, and +-45% of a 0.10
                        // thickness reads as lumps rather than as glaze
                        // settling. The field did not change; what
                        // changed is that it is now being listened to.
                        double pooling = icing_noise(u * 3.0, v * 3.0, 4.0) * 0.016; // broad, gentle waves
                        return (0.10 + pooling) * falloff;
                    }}, bd::EdgeRule::ZeroThickness)
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
                                     [](double) { return Vec<double, 3>{0.0, 0.0, table_yaw}; })
                             + bd::VecField<double>::constant({0.0, 0.0, -1.10}));

    std::vector<bd::Handle<double>> scene_children{table, cube, dough, icing};
    scene_children.insert(scene_children.end(), sprinkle_groups.begin(), sprinkle_groups.end());
    scene_children.push_back(dust);   // one node now, not 35 200
    auto lesson = scene.compose(scene_children);

    std::println("donut_demo: a Trace is real data -- here it is, {} nodes:", scene.size());
    for (std::size_t i = 0; i < scene.size(); ++i)
        std::println("  [{}] {}", i, bd::kind_name(scene.node(i).kind));

    auto placed = bd::materialize(scene, lesson.index, t);

    // What level each object actually resolves to, reported rather than
    // left to be inferred from a frame time. On this scene the answer is
    // "all of them tessellate", and the reason is worth seeing rather
    // than hiding: the dough is an Offset carrying a noise bump, so it
    // stopped being a torus the moment it got bread texture. An exact
    // form is a promise about the shape, and a bumped torus cannot keep
    // it. The machinery reporting zero here is the machinery working.
    {
        std::size_t exact = 0, tess = 0, newton = 0;
        for (const auto& obj : placed) {
            switch (obj.render_level()) {
                case bd::RenderLevel::Exact:       ++exact;  break;
                case bd::RenderLevel::Tessellated: ++tess;   break;
                case bd::RenderLevel::Newton:      ++newton; break;
            }
        }
        // What cook() would find, asked directly so the scene survives to
        // be rendered. A node whose motion is a placement can share one
        // geometry with every other node of the same shape; one whose
        // motion is an opaque point map cannot, because nothing can tell
        // whether it moves the object or reshapes it.
        std::size_t placeable = 0, deforming = 0;
        for (std::size_t i = 0; i < scene.size(); ++i) {
            if (scene.node(i).kind == bd::Kind::Compose) continue;
            (scene.node(i).transform.is_placement() ? placeable : deforming)++;
        }
        std::println("motions: {} placements (shareable geometry), {} deformations (not)",
                     placeable, deforming);

        auto fs = bd::field_report(scene);
        std::println("fields: {} total, {} structural, {} opaque; leaves {} "
                     "({} recognized, {} unknown), {} distinct types, {} B captured [{}]",
                     fs.fields, fs.structural_fields, fs.opaque_fields,
                     fs.opaque_leaves, fs.recognized_leaves, fs.unknown_leaves,
                     fs.distinct_types, fs.payload_bytes,
                     bd::payload_verdict(fs.payload_bytes));
        std::println("render levels: {} exact, {} tessellated, {} newton",
                     exact, tess, newton);
    }
    // The vertex and triangle totals used to be counted here, by walking
    // `placed` and calling `mesh()` on every object. That loop was the
    // single most expensive thing in the program and it existed to print
    // two integers.
    //
    // `mesh()` on a Scatter flattens the whole cloud into one
    // non-deduplicated mesh, so at two million particles it built
    // 64 351 744 vertices and 96 527 616 faces -- 2.70 GB of them --
    // taking twenty seconds, and then threw them away. Measured peak RSS
    // was 2.88 GB, of which this loop was 2.88 GB minus a rounding error;
    // the render path it was supposedly measuring settles under 1 GB.
    // Because it ran before the `--photo` gate, a console run with no
    // output at all paid the whole cost.
    //
    // What makes it worth a comment rather than a silent deletion is what
    // the number was for. It is the "vertices without instancing" figure,
    // the one quoted to show that instancing works -- and it was obtained
    // by performing, in full, the work instancing exists to avoid. The
    // measurement cost more than the thing it measured, and it was also
    // redundant: `cook()` reports the same total two lines below, summed
    // over the shape table without materialising anything
    // (`Cooked::vertices_without_instancing`).
    //
    // The timing went with it. With no mesh built, `materialize()` only
    // hands back views onto trace nodes, so a duration printed here would
    // be measuring nothing and inviting the reader to compare it against
    // `cook()`'s, which does real work.
    std::println("materialized at t={}: exploded cube + dough + icing + {} sprinkles",
                 t, sprinkle_count);

    if (photo) {
        // Timed on its own, because "cook() got faster" has two possible
        // causes -- fewer nodes, or a hash that stopped touching vertices
        // -- and only a number separates them.
        auto c0 = std::chrono::steady_clock::now();
        auto cooked = bd::cook(scene, lesson.index, t);
        auto c1 = std::chrono::steady_clock::now();
        std::println("  cook() {:.1f} ms over {} trace nodes",
                     std::chrono::duration<double, std::milli>(c1 - c0).count(), scene.size());
        render_photo(scene, cooked, out_path, force);
    }
    if (video) render_video(scene, lesson.index, video_dir, force, build_frames, orbit_frames);
    return 0;
}
