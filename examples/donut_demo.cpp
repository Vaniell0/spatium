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
// The steps themselves live in donut_scene.hpp, apart from the renderer
// below, so that a test can build the very same scene -- see that file.
//
// Two output modes:
//   (default)     print the trace and materialized stats to the console
//   --photo PATH  also CPU-raytrace it to a PNG (BVH<Triangle3>, same
//                 engine as every other offline demo in this tree --
//                 tessellation happens here, only for display, not as
//                 part of the scene's own representation)

#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "donut_scene.hpp"
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
using namespace donut;

namespace {


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
    auto boom_uv = load_boom_points("examples/data/boom_points.txt");
    if (boom_uv.empty()) boom_uv = load_boom_points("data/boom_points.txt");
    const std::size_t lesson = build_scene(scene, sprinkle_count, dust_particles, boom_uv);

    std::println("donut_demo: a Trace is real data -- here it is, {} nodes:", scene.size());
    for (std::size_t i = 0; i < scene.size(); ++i)
        std::println("  [{}] {}", i, bd::kind_name(scene.node(i).kind));

    auto placed = bd::materialize(scene, lesson, t);

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
        auto cooked = bd::cook(scene, lesson, t);
        auto c1 = std::chrono::steady_clock::now();
        std::println("  cook() {:.1f} ms over {} trace nodes",
                     std::chrono::duration<double, std::milli>(c1 - c0).count(), scene.size());
        render_photo(scene, cooked, out_path, force);
    }
    if (video) render_video(scene, lesson, video_dir, force, build_frames, orbit_frames);
    return 0;
}
