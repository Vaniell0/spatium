// The cloth-on-an-obstacle scene, apart from any window: a square XPBD
// sheet dropped on a sphere, a torus or a Klein bottle, the ways its
// contact can be resolved, and the numbers that say whether it held.
//
// It is the scene cloth_visual_demo.cpp drew before it was removed in
// April (commit 1b3ce2c62) for the reason that commit gives: stiff
// structural constraints pull the overhanging corners hard enough to drag
// the centre through the obstacle within one substep, faster than a
// contact band can react, and softening them makes the sheet stretch like
// rubber. Kept here so that the failure is measured, not remembered, and
// so the same scene serves every contact method compared on it:
//
//   Projection  the removed demo's: after each substep, a vertex inside
//               the band is put back on it along the analytic surface's
//               normal (point_to).
//   Triangles   the ordinary game route: the obstacle tessellated, and a
//               vertex pushed out of the nearest triangle, found through a
//               BVH -- discrete, and blind to what happened between two
//               positions.
//
// The measures are independent of either method: penetration and how many
// vertices ended up inside come from the exact distance to the analytic
// obstacle, stretch from the structural edges' lengths, and a run that
// produced a NaN or a velocity past 50 m/s is counted as exploded.

#pragma once

#include <spatium/algebra/vector.hpp>
#include <spatium/geometry/ray_surface.hpp>
#include <spatium/geometry/triangle.hpp>
#include <spatium/mesh/mesh.hpp>
#include <spatium/mesh/primitives.hpp>
#include <spatium/physics/mechanics/narrow_phase.hpp>
#include <spatium/physics/mechanics/xpbd.hpp>
#include <spatium/spaces/euclidean.hpp>
#include <spatium/spaces/parametric.hpp>
#include <spatium/spaces/sphere.hpp>
#include <spatium/spatial/bvh.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <tuple>
#include <vector>

namespace cloth {

using namespace spatium;
using namespace spatium::physics::mechanics;
using V3 = Vec<double, 3>;

enum class Obstacle { Sphere, Torus, Klein };
enum class Contact { Projection, Triangles };

struct Config {
    int grid = 19;
    double spacing = 0.075;
    double start_y = 1.4;
    double dhat = 0.04;
    double substep_dt = 1.0 / 480.0;
    double struct_compl = 5e-6;   // stretch compliance (stiffness^-1)
    double bend_compl = 5e-3;     // bend compliance
    int iterations = 8;           // Gauss-Seidel iterations per substep
    double gravity = -9.81;
};

inline const char* name(Obstacle o) {
    return o == Obstacle::Sphere ? "sphere" : o == Obstacle::Torus ? "torus" : "klein";
}
inline const char* name(Contact c) { return c == Contact::Projection ? "projection" : "triangles"; }

// The obstacles, y up, as the removed demo had them.
struct Obstacles {
    Sphere<2, double> sphere{0.55};
    geometry::Torus<double> torus{V3{}, V3{0.0, 1.0, 0.0}, 0.35, 0.22};
    ParametricSurface<double> sphere_chart{
        [](double u, double v) -> V3 {
            const double r = 0.55;
            return {r * std::sin(v) * std::cos(u), r * std::cos(v), r * std::sin(v) * std::sin(u)};
        },
        {0.0, 2.0 * std::numbers::pi, 0.0, std::numbers::pi}, true, false};
    ParametricSurface<double> torus_chart{
        [](double u, double v) -> V3 {
            const double R = 0.35, r = 0.22;
            return {(R + r * std::cos(v)) * std::cos(u), r * std::sin(v), (R + r * std::cos(v)) * std::sin(u)};
        },
        {0.0, 2.0 * std::numbers::pi, 0.0, 2.0 * std::numbers::pi}, true, true};
    // Bonan-Jennings immersion: one smooth polynomial, no seam for a vertex
    // to catch on (the Lawson form the demo used first split at u = pi).
    ParametricSurface<double> klein{
        [](double u, double v) -> V3 {
            const double cu = std::cos(u), su = std::sin(u), cv = std::cos(v), sv = std::sin(v);
            const double cu2 = cu * cu, cu3 = cu2 * cu, cu4 = cu2 * cu2, cu5 = cu4 * cu, cu6 = cu4 * cu2,
                         cu7 = cu6 * cu;
            const double x = -(2.0 / 15.0) * cu *
                             (3.0 * cv - 30.0 * su + 90.0 * cu4 * su - 60.0 * cu6 * su + 5.0 * cu * cv * su);
            const double y = -(1.0 / 15.0) * su *
                             (3.0 * cv - 3.0 * cu2 * cv - 48.0 * cu4 * cv + 48.0 * cu6 * cv - 60.0 * su +
                              5.0 * cu * cv * su - 5.0 * cu3 * cv * su - 80.0 * cu5 * cv * su + 80.0 * cu7 * cv * su);
            const double z = (2.0 / 15.0) * (3.0 + 5.0 * cu * su) * sv;
            constexpr double s = 0.35;
            return {x * s, y * s, z * s};
        },
        {0.0, 2.0 * std::numbers::pi, 0.0, 2.0 * std::numbers::pi}, true, true};

    const ParametricSurface<double>& chart(Obstacle o) const {
        return o == Obstacle::Sphere ? sphere_chart : o == Obstacle::Torus ? torus_chart : klein;
    }
    // The tessellation the Triangles route pushes vertices out of.
    mesh::Mesh<ParametricSurface<double>> tessellation(Obstacle o, int nu = 48, int nv = 24) const {
        return mesh::parametric_mesh(chart(o), nu, nv);
    }
};

struct Sheet {
    std::vector<XpbdParticle<3, double>> parts;
    std::vector<std::array<std::uint32_t, 3>> faces;
    std::vector<XpbdDistanceConstraint<3, double>> cons;
    std::size_t n_struct = 0;
};

// A flat grid, u/v edges structural and the diagonals only as bending --
// a structural diagonal also resists folding and keeps the sheet a plate.
inline Sheet make_sheet(const Config& c) {
    Sheet s;
    const int n = c.grid;
    const double half = (n - 1) * c.spacing * 0.5;
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
            XpbdParticle<3, double> p;
            p.x = V3{i * c.spacing - half, c.start_y, j * c.spacing - half};
            p.x_prev = p.x;
            p.w = 1.0;
            s.parts.push_back(p);
        }
    auto idx = [n](int i, int j) { return static_cast<std::uint32_t>(j * n + i); };
    for (int j = 0; j + 1 < n; ++j)
        for (int i = 0; i + 1 < n; ++i) {
            s.faces.push_back({idx(i, j), idx(i + 1, j), idx(i, j + 1)});
            s.faces.push_back({idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)});
        }
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
            if (i + 1 < n) s.cons.push_back({idx(i, j), idx(i + 1, j), c.spacing, c.struct_compl, 0.0});
            if (j + 1 < n) s.cons.push_back({idx(i, j), idx(i, j + 1), c.spacing, c.struct_compl, 0.0});
        }
    s.n_struct = s.cons.size();
    auto bend = build_bending_distance_constraints<3, double>(s.parts, s.faces, c.bend_compl);
    s.cons.insert(s.cons.end(), bend.begin(), bend.end());
    return s;
}

// The Triangles route's obstacle: a BVH over the tessellation.
struct TriangleObstacle {
    spatial::BVH<geometry::Triangle3> bvh;
    static TriangleObstacle build(const mesh::Mesh<ParametricSurface<double>>& m) {
        std::vector<geometry::Triangle3> tris;
        for (const auto& f : m.faces) tris.emplace_back(m.vertices[f[0]], m.vertices[f[1]], m.vertices[f[2]]);
        return {spatial::BVH<geometry::Triangle3>::build(std::move(tris))};
    }
};

// The exact signed distance to the analytic obstacle, the measure every
// method is judged by; the Klein bottle has none (it is one-sided).
inline bool exact_distance(const Obstacles& ob, Obstacle o, const V3& x, double& d) {
    if (o == Obstacle::Sphere) { d = x.norm() - ob.sphere.radius; return true; }
    if (o == Obstacle::Torus) { d = point_to(x, ob.torus).signed_distance(); return true; }
    return false;
}

struct Metrics {
    double worst_penetration = 0;   // deepest below the surface after contact, >= 0
    double worst_before = 0;        // deepest the solver pulled a vertex before contact
    long tunnelled = 0;             // substep moves that entered the obstacle and left it
    int inside = 0;                 // vertices more than dhat below it at the end
    double max_stretch = 0;         // over structural edges and all substeps
    bool exploded = false;
    double ms_per_substep = 0;
};

// Whether the straight move a -> b went into the solid obstacle and out
// again -- the tunnel a check at the two ends cannot see -- by the exact
// crossings of the segment: a quadric's two roots, a torus's four.
inline bool passed_through(const Obstacles& ob, Obstacle o, const V3& a, const V3& b) {
    const V3 d{b - a};
    const double len = d.norm();
    if (len <= 0) return false;
    const geometry::Ray<3, double> ray{a, V3{d * (1.0 / len)}};
    int crossings = 0;
    if (o == Obstacle::Sphere) {
        for (const auto& h : geometry::ray_quadric(ray, geometry::Quadric<double>::sphere(ob.sphere.radius)))
            if (h.t > 0 && h.t < len) ++crossings;
    } else if (o == Obstacle::Torus) {
        for (const auto& h : geometry::ray_torus(ray, ob.torus))
            if (h.t > 0 && h.t < len) ++crossings;
    }
    return crossings >= 2;
}

// One substep: XPBD, then contact. `m`, when given, records what the
// solver did before contact had its say: how deep it pulled a vertex and
// how many moves went clean through.
inline void substep(Sheet& s, const Config& c, const Obstacles& ob, Obstacle o, Contact mode,
                    const TriangleObstacle* tri, Metrics* m = nullptr) {
    const V3 g{0.0, c.gravity, 0.0};
    std::vector<V3> before;
    if (m) {
        before.reserve(s.parts.size());
        for (const auto& p : s.parts) before.push_back(p.x);
    }
    xpbd_step(s.parts, s.cons.begin(), s.cons.end(), [&](const auto&, double) { return g; }, c.substep_dt,
              c.iterations);
    if (m)
        for (std::size_t i = 0; i < s.parts.size(); ++i) {
            double d;
            if (exact_distance(ob, o, s.parts[i].x, d)) m->worst_before = std::max(m->worst_before, -d);
            if (passed_through(ob, o, before[i], s.parts[i].x)) ++m->tunnelled;
        }
    for (auto& p : s.parts) {
        if (mode == Contact::Projection) {
            ContactQuery<double> q;
            if (o == Obstacle::Sphere) q = point_to(p.x, ob.sphere);
            else if (o == Obstacle::Torus) q = point_to(p.x, ob.torus);
            else q = point_to(p.x, ob.klein);
            if (o == Obstacle::Klein) {
                // One-sided and self-intersecting: the removed demo's filters
                // -- the top half, upward normals, approaching vertices only.
                if (q.closest_point[1] < 0.0 || q.normal[1] < 0.0) continue;
                if (V3{p.x - p.x_prev}.dot(q.normal) > 0.0) continue;
            }
            if (c.dhat - q.signed_distance() > 0.0) p.x = V3{q.closest_point + q.normal * c.dhat};
        } else {
            const auto near = tri->bvh.nearest(p.x);
            if (!near) continue;
            const auto& t = tri->bvh.shapes()[near->index];
            V3 n = t.normal();
            // Which side of the tessellation is out: away from the centre
            // for the sphere, from the tube's centre circle for the torus,
            // and up for the one-sided Klein bottle, as the removed demo
            // filtered it.
            V3 out{0.0, 1.0, 0.0};
            if (o == Obstacle::Sphere) out = V3{near->point};
            else if (o == Obstacle::Torus) {
                const V3 q{near->point};
                const double rho = std::hypot(q[0], q[2]);
                const V3 ring = rho > 0 ? V3{V3{q[0], 0.0, q[2]} * (ob.torus.major_radius / rho)} : V3{};
                out = V3{q - ring};
            }
            if (n.dot(out) < 0) n = V3{n * -1.0};
            const double sd = V3{p.x - near->point}.dot(n);
            if (sd < c.dhat) p.x = V3{p.x + n * (c.dhat - sd)};
        }
    }
}

inline void measure_into(Metrics& m, const Sheet& s, const Config& c, const Obstacles& ob, Obstacle o) {
    for (std::size_t k = 0; k < s.n_struct; ++k) {
        const auto& e = s.cons[k];
        const double len = V3{s.parts[e.i].x - s.parts[e.j].x}.norm();
        m.max_stretch = std::max(m.max_stretch, len / e.rest - 1.0);
    }
    for (const auto& p : s.parts) {
        const V3 v{(p.x - p.x_prev) * (1.0 / c.substep_dt)};
        if (!std::isfinite(p.x[0] + p.x[1] + p.x[2]) || v.norm() > 50.0) m.exploded = true;
        double d;
        if (exact_distance(ob, o, p.x, d)) m.worst_penetration = std::max(m.worst_penetration, -d);
    }
}

} // namespace cloth
