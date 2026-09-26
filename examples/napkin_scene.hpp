// Napkins: several cloth sheets dropped onto a sphere and onto each other,
// the case the cloth-on-a-ball scene does not have. The ball alone needs no
// continuous collision -- a projection onto an analytic surface holds -- but
// two sheets of cloth meeting at speed pass through each other between two
// substeps, where no discrete test looks.
//
// A substep, in order:
//   1. XPBD: stretch and bending, gravity.
//   2. The sphere by projection, as cloth_scene.hpp does it.
//   3. Cloth against cloth, discrete: a vertex within `thickness` of a
//      triangle of another part of the cloth is pushed out, the triangle
//      the other way, on the side the vertex came from.
//   4. Cloth against cloth, continuous, as the guarantee: every
//      vertex-triangle and edge-edge pair whose moves over the substep
//      could meet is asked, through physics/mechanics/surface_ccd.hpp,
//      whether they cross; a pair that does is stopped short of the moment
//      they would, and the pass repeats until none does, or reverts what is
//      left to the substep's start, where nothing crossed.
//
// The broad phase is a spatial hash over the moves' boxes. The measure is
// independent of all of it: how many edges of the cloth pierce a triangle
// of it that they do not belong to, counted by a segment-triangle test at
// the end of every substep -- zero is the claim.

#pragma once

#include <spatium/algebra/vector.hpp>
#include <spatium/physics/mechanics/surface_ccd.hpp>
#include <spatium/physics/mechanics/xpbd.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <unordered_map>
#include <utility>
#include <vector>

namespace napkins {

using namespace spatium;
using namespace spatium::physics::mechanics;
using V3 = Vec<double, 3>;

struct Config {
    int napkins = 4;
    int grid = 13;                 // vertices a side
    double spacing = 0.05;         // 0.6 m a napkin
    double sphere = 0.45;          // radius, at the origin
    double band = 0.015;           // the sphere's and the floor's contact band
    double floor = -0.45;          // a floor at the sphere's foot
    double friction = 0.4;         // of the tangential move, taken back on contact
    double thickness = 0.008;      // cloth against cloth
    double substep_dt = 1.0 / 240.0;
    double struct_compl = 1e-6;
    double bend_compl = 2e-3;
    int iterations = 8;
    double gravity = -9.81;
    int ccd_rounds = 4;
    // Thrown: each napkin starts moving at this speed towards the stack's
    // axis, alternately from either side, so they meet in the air at a
    // speed a substep's discrete test does not see. Zero: dropped.
    double throw_speed = 0;
    bool ccd = true;               // step 4; off shows what it prevents
};

struct Cloth {
    std::vector<XpbdParticle<3, double>> parts;
    std::vector<std::array<std::uint32_t, 3>> faces;
    std::vector<std::array<std::uint32_t, 2>> edges;
    std::vector<XpbdDistanceConstraint<3, double>> cons;
    std::vector<int> part_of;      // napkin of each vertex
};

// Napkins stacked above the sphere, each turned and shifted so they fall
// across each other rather than flat on flat.
inline Cloth make_cloth(const Config& c) {
    Cloth s;
    const int n = c.grid;
    const double half = (n - 1) * c.spacing * 0.5;
    for (int k = 0; k < c.napkins; ++k) {
        const double a = 0.6 * k, tilt = 0.15 * ((k % 2) ? 1 : -1);
        const double ca = std::cos(a), sa = std::sin(a);
        // Layers far enough apart that neighbours, tilted opposite ways,
        // clear each other at their edges: 2 * tilt * half, and a margin.
        const double gap = 2 * 0.15 * half + 0.05;
        const V3 at{0.12 * std::cos(2.1 * k), 0.75 + gap * k, 0.12 * std::sin(2.1 * k)};
        const auto base = static_cast<std::uint32_t>(s.parts.size());
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i) {
                const double x = i * c.spacing - half, z = j * c.spacing - half;
                XpbdParticle<3, double> p;
                p.x = V3{at[0] + ca * x - sa * z, at[1] + tilt * x, at[2] + sa * x + ca * z};
                p.x_prev = p.x;
                p.w = 1.0;
                s.parts.push_back(p);
                s.part_of.push_back(k);
            }
        const auto idx = [n, base](int i, int j) { return base + static_cast<std::uint32_t>(j * n + i); };
        const auto first_face = s.faces.size();
        for (int j = 0; j + 1 < n; ++j)
            for (int i = 0; i + 1 < n; ++i) {
                s.faces.push_back({idx(i, j), idx(i + 1, j), idx(i, j + 1)});
                s.faces.push_back({idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)});
            }
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i) {
                if (i + 1 < n) s.cons.push_back({idx(i, j), idx(i + 1, j), c.spacing, c.struct_compl, 0.0});
                if (j + 1 < n) s.cons.push_back({idx(i, j), idx(i, j + 1), c.spacing, c.struct_compl, 0.0});
                if (i + 1 < n && j + 1 < n)
                    s.cons.push_back({idx(i + 1, j), idx(i, j + 1), c.spacing * std::sqrt(2.0), c.struct_compl * 10, 0.0});
            }
        std::vector<std::array<std::uint32_t, 3>> own(s.faces.begin() + static_cast<std::ptrdiff_t>(first_face),
                                                      s.faces.end());
        auto bend = build_bending_distance_constraints<3, double>(s.parts, own, c.bend_compl);
        s.cons.insert(s.cons.end(), bend.begin(), bend.end());
    }
    // Thrown: a velocity into x_prev, towards the axis from alternate sides.
    if (c.throw_speed > 0) {
        const int per = c.grid * c.grid;
        for (int k = 0; k < c.napkins; ++k) {
            const double side = (k % 2) ? 1.0 : -1.0;
            for (int i = 0; i < per; ++i) {
                auto& p = s.parts[std::size_t(k * per + i)];
                p.x = V3{p.x + V3{-side * 0.6, 0.0, 0.0}};
                p.x_prev = V3{p.x - V3{side * c.throw_speed, 0.0, 0.0} * c.substep_dt};
            }
        }
    }
    // Every edge once.
    std::unordered_map<std::uint64_t, int> seen;
    for (const auto& f : s.faces)
        for (int e = 0; e < 3; ++e) {
            std::uint32_t a = f[e], b = f[(e + 1) % 3];
            if (a > b) std::swap(a, b);
            if (seen.emplace((std::uint64_t(a) << 32) | b, 1).second) s.edges.push_back({a, b});
        }
    return s;
}

// ── Broad phase: a hash grid over boxes ─────────────────────────

struct Box {
    V3 lo, hi;
};
inline Box box_of(std::initializer_list<V3> pts, double pad) {
    Box b{V3{1e300, 1e300, 1e300}, V3{-1e300, -1e300, -1e300}};
    for (const auto& p : pts)
        for (int k = 0; k < 3; ++k) {
            b.lo[k] = std::min(b.lo[k], p[k] - pad);
            b.hi[k] = std::max(b.hi[k], p[k] + pad);
        }
    return b;
}
inline bool overlap(const Box& a, const Box& b) {
    for (int k = 0; k < 3; ++k)
        if (a.hi[k] < b.lo[k] || b.hi[k] < a.lo[k]) return false;
    return true;
}

// A hash grid as a sorted array of (cell, id): built by one sort, looked
// up by binary search, no allocation per bucket.
struct Grid {
    double cell;
    std::vector<std::pair<std::uint64_t, std::uint32_t>> items;
    static std::uint64_t key(long x, long y, long z) {
        return (std::uint64_t(x & 0x1fffff) << 42) | (std::uint64_t(y & 0x1fffff) << 21) | std::uint64_t(z & 0x1fffff);
    }
    template<typename F> void cells(const Box& b, F&& f) const {
        const long x0 = long(std::floor(b.lo[0] / cell)), x1 = long(std::floor(b.hi[0] / cell));
        const long y0 = long(std::floor(b.lo[1] / cell)), y1 = long(std::floor(b.hi[1] / cell));
        const long z0 = long(std::floor(b.lo[2] / cell)), z1 = long(std::floor(b.hi[2] / cell));
        for (long x = x0; x <= x1; ++x)
            for (long y = y0; y <= y1; ++y)
                for (long z = z0; z <= z1; ++z) f(key(x, y, z));
    }
    void insert(const Box& b, std::uint32_t id) { cells(b, [&](std::uint64_t k) { items.push_back({k, id}); }); }
    void finish() { std::sort(items.begin(), items.end()); }
    template<typename F> void each(std::uint64_t k, F&& f) const {
        auto it = std::lower_bound(items.begin(), items.end(), std::pair<std::uint64_t, std::uint32_t>{k, 0});
        for (; it != items.end() && it->first == k; ++it) f(it->second);
    }
};

// Candidate pairs: vertex-triangle and edge-edge whose swept boxes meet,
// with no vertex in common.
struct Candidates {
    std::vector<std::pair<std::uint32_t, std::uint32_t>> vf, ee;
};
// A spatial hash (Teschner et al. 2003): cells of about an element's
// size, hashed into a fixed table filled by counting -- two passes, no
// sort, no allocation per bucket. A box is entered in every cell it
// covers; a query box reads the cells it covers, and a collision of two
// cells in one bucket only costs a box test.
struct SpatialHash {
    double cell;
    std::uint32_t mask;
    std::vector<std::uint32_t> start, ids;
    std::uint32_t bucket(long x, long y, long z) const {
        return std::uint32_t((x * 73856093L) ^ (y * 19349663L) ^ (z * 83492791L)) & mask;
    }
    template<typename F> void cells(const Box& b, F&& f) const {
        const long x0 = long(std::floor(b.lo[0] / cell)), x1 = long(std::floor(b.hi[0] / cell));
        const long y0 = long(std::floor(b.lo[1] / cell)), y1 = long(std::floor(b.hi[1] / cell));
        const long z0 = long(std::floor(b.lo[2] / cell)), z1 = long(std::floor(b.hi[2] / cell));
        for (long x = x0; x <= x1; ++x)
            for (long y = y0; y <= y1; ++y)
                for (long z = z0; z <= z1; ++z) f(bucket(x, y, z));
    }
    SpatialHash(const std::vector<Box>& boxes, double cell_, std::uint32_t bits)
        : cell(cell_), mask((1u << bits) - 1), start((1u << bits) + 1, 0) {
        for (const auto& b : boxes) cells(b, [&](std::uint32_t k) { ++start[k + 1]; });
        for (std::size_t k = 1; k < start.size(); ++k) start[k] += start[k - 1];
        ids.resize(start.back());
        std::vector<std::uint32_t> fill(start.begin(), start.end() - 1);
        for (std::uint32_t i = 0; i < boxes.size(); ++i) cells(boxes[i], [&](std::uint32_t k) { ids[fill[k]++] = i; });
    }
    template<typename F> void near(const Box& b, F&& f) const {
        cells(b, [&](std::uint32_t k) {
            for (std::uint32_t j = start[k]; j < start[k + 1]; ++j) f(ids[j]);
        });
    }
};

inline Candidates broad_phase(const Cloth& s, const std::vector<V3>& from, double pad, double cell) {
    Candidates out;
    const auto& to = s.parts;
    const std::size_t nv = s.parts.size(), nf = s.faces.size(), ne = s.edges.size();
    std::vector<Box> vb(nv), tb(nf), eb(ne);
    for (std::size_t v = 0; v < nv; ++v) vb[v] = box_of({from[v], to[v].x}, pad);
    for (std::size_t f = 0; f < nf; ++f) {
        const auto& F = s.faces[f];
        tb[f] = box_of({from[F[0]], from[F[1]], from[F[2]], to[F[0]].x, to[F[1]].x, to[F[2]].x}, pad);
    }
    for (std::size_t e = 0; e < ne; ++e) {
        const auto& E = s.edges[e];
        eb[e] = box_of({from[E[0]], from[E[1]], to[E[0]].x, to[E[1]].x}, pad);
    }
    const std::uint32_t bits = std::max(10u, std::uint32_t(std::ceil(std::log2(double(nf + ne) * 2))));
    {
        const SpatialHash h(tb, cell, bits);
        std::vector<std::uint32_t> mark(nf, ~0u);
        for (std::uint32_t v = 0; v < nv; ++v)
            h.near(vb[v], [&](std::uint32_t f) {
                if (mark[f] == v) return;
                mark[f] = v;
                const auto& F = s.faces[f];
                if (F[0] != v && F[1] != v && F[2] != v && overlap(vb[v], tb[f])) out.vf.push_back({v, f});
            });
    }
    {
        const SpatialHash h(eb, cell, bits);
        std::vector<std::uint32_t> mark(ne, ~0u);
        for (std::uint32_t e = 0; e < ne; ++e)
            h.near(eb[e], [&](std::uint32_t o) {
                if (o <= e || mark[o] == e) return;
                mark[o] = e;
                const auto& A = s.edges[e];
                const auto& B = s.edges[o];
                if (A[0] == B[0] || A[0] == B[1] || A[1] == B[0] || A[1] == B[1]) return;
                if (overlap(eb[e], eb[o])) out.ee.push_back({e, o});
            });
    }
    return out;
}

// ── Narrow phase through surface_ccd ────────────────────────────

// A query's eight points on the caller's stack, and charts that point at
// them: a lambda holding a pointer fits std::function's own storage, where
// one holding three or six vectors is sent to the heap -- an allocation a
// chart, four a query.
struct Moves {
    std::array<V3, 4> from, to;
};
inline MovingChart<double> point_chart(const Moves& m, int i) {
    MovingChart<double> c;
    const Moves* q = &m;
    c.p = [q, i](double, double) { return q->from[i]; };
    c.w = [q, i](double, double) { return V3{q->to[i] - q->from[i]}; };
    c.u1 = c.v1 = 0;
    return c;
}
inline MovingChart<double> edge_chart(const Moves& m, int i) {
    MovingChart<double> c;
    const Moves* q = &m;
    c.p = [q, i](double u, double) { return V3{q->from[i] + (q->from[i + 1] - q->from[i]) * u}; };
    c.w = [q, i](double u, double) {
        const V3 d0{q->to[i] - q->from[i]}, d1{q->to[i + 1] - q->from[i + 1]};
        return V3{d0 + (d1 - d0) * u};
    };
    c.bp = {V3{m.from[i + 1] - m.from[i]}.norm(), 0, 0, 0};
    c.bw = {V3{V3{m.to[i + 1] - m.from[i + 1]} - V3{m.to[i] - m.from[i]}}.norm(), 0, 0, 0};
    c.v1 = 0;
    return c;
}
// Points i, i+1, i+2 as x(u, v) = a + u ((1 - v)(b - a) + v (c - a)).
inline MovingChart<double> triangle_chart(const Moves& m, int i) {
    MovingChart<double> c;
    const Moves* q = &m;
    const auto tri = [](const V3& a, const V3& b, const V3& c, double u, double v) {
        return V3{a + (V3{b - a} * (1 - v) + V3{c - a} * v) * u};
    };
    c.p = [q, i, tri](double u, double v) { return tri(q->from[i], q->from[i + 1], q->from[i + 2], u, v); };
    c.w = [q, i, tri](double u, double v) {
        return tri(V3{q->to[i] - q->from[i]}, V3{q->to[i + 1] - q->from[i + 1]}, V3{q->to[i + 2] - q->from[i + 2]}, u, v);
    };
    const V3 da{m.to[i] - m.from[i]}, db{m.to[i + 1] - m.from[i + 1]}, dc{m.to[i + 2] - m.from[i + 2]};
    c.bp = {std::max(V3{m.from[i + 1] - m.from[i]}.norm(), V3{m.from[i + 2] - m.from[i]}.norm()),
            V3{m.from[i + 2] - m.from[i + 1]}.norm(), 0, 0};
    c.bw = {std::max(V3{db - da}.norm(), V3{dc - da}.norm()), V3{dc - db}.norm(), 0, 0};
    return c;
}

// The closest point of a triangle to p, as barycentric weights (Ericson).
inline std::array<double, 3> closest_on_triangle(const V3& p, const V3& a, const V3& b, const V3& c) {
    const V3 ab{b - a}, ac{c - a}, ap{p - a};
    const double d1 = ab.dot(ap), d2 = ac.dot(ap);
    if (d1 <= 0 && d2 <= 0) return {1, 0, 0};
    const V3 bp{p - b};
    const double d3 = ab.dot(bp), d4 = ac.dot(bp);
    if (d3 >= 0 && d4 <= d3) return {0, 1, 0};
    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) { const double v = d1 / (d1 - d3); return {1 - v, v, 0}; }
    const V3 cp{p - c};
    const double d5 = ab.dot(cp), d6 = ac.dot(cp);
    if (d6 >= 0 && d5 <= d6) return {0, 0, 1};
    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) { const double w = d2 / (d2 - d6); return {1 - w, 0, w}; }
    const double va = d3 * d6 - d5 * d4;
    if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
        const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return {0, 1 - w, w};
    }
    const double den = 1 / (va + vb + vc);
    const double v = vb * den, w = vc * den;
    return {1 - v - w, v, w};
}

// Closest points of segments p1-q1 and p2-q2 as parameters s, t (Ericson).
inline std::pair<double, double> closest_segments(const V3& p1, const V3& q1, const V3& p2, const V3& q2) {
    const V3 d1{q1 - p1}, d2{q2 - p2}, r{p1 - p2};
    const double a = d1.dot(d1), e = d2.dot(d2), f = d2.dot(r);
    double s = 0, t = 0;
    if (a <= 1e-24 && e <= 1e-24) return {0, 0};
    if (a <= 1e-24) return {0, std::clamp(f / e, 0.0, 1.0)};
    const double c = d1.dot(r);
    if (e <= 1e-24) return {std::clamp(-c / a, 0.0, 1.0), 0};
    const double b = d1.dot(d2), den = a * e - b * b;
    s = den > 1e-24 ? std::clamp((b * f - c * e) / den, 0.0, 1.0) : 0.0;
    t = (b * s + f) / e;
    if (t < 0) { t = 0; s = std::clamp(-c / a, 0.0, 1.0); }
    else if (t > 1) { t = 1; s = std::clamp((b - c) / a, 0.0, 1.0); }
    return {s, t};
}

// Whether segment p-q pierces triangle a-b-c (Moller-Trumbore, both ends).
inline bool pierces(const V3& p, const V3& q, const V3& a, const V3& b, const V3& c) {
    const V3 d{q - p}, e1{b - a}, e2{c - a};
    const V3 h{d.cross(e2)};
    const double det = e1.dot(h);
    if (std::abs(det) < 1e-18) return false;
    const double inv = 1 / det;
    const V3 s{p - a};
    const double u = inv * s.dot(h);
    if (u < 0 || u > 1) return false;
    const V3 qv{s.cross(e1)};
    const double v = inv * d.dot(qv);
    if (v < 0 || u + v > 1) return false;
    const double t = inv * e2.dot(qv);
    return t > 0 && t < 1;
}

// Whether four points moving linearly could become coplanar over [0, 1]:
// the signed volume det[b - a, c - a, d - a] is a cubic in t, and when its
// four Bernstein coefficients on [0, 1] share a sign, clear of rounding by
// a margin far above it, the cubic has no root there -- the convex hull of
// a Bezier curve's control points holds it -- and neither a vertex through
// a triangle nor an edge through an edge can happen. The standard
// coplanarity filter; it answers "maybe" whenever it cannot prove "no".
inline bool may_be_coplanar(const V3& a0, const V3& b0, const V3& c0, const V3& d0, const V3& a1, const V3& b1,
                            const V3& c1, const V3& d1) {
    const auto vol = [&](double t) {
        const auto at = [t](const V3& x, const V3& y) { return V3{x + (y - x) * t}; };
        const V3 a = at(a0, a1);
        return V3{at(b0, b1) - a}.dot(V3{V3{at(c0, c1) - a}.cross(V3{at(d0, d1) - a})});
    };
    // Values at 0, 1/3, 2/3, 1 determine the cubic; to Bernstein form.
    const double f0 = vol(0), f1 = vol(1.0 / 3), f2 = vol(2.0 / 3), f3 = vol(1);
    const double q0 = f0, q3 = f3;
    const double q1 = (-5 * f0 + 18 * f1 - 9 * f2 + 2 * f3) / 6;
    const double q2 = (2 * f0 - 9 * f1 + 18 * f2 - 5 * f3) / 6;
    double L = 0;
    for (const V3* p : std::initializer_list<const V3*>{&a0, &b0, &c0, &d0, &a1, &b1, &c1, &d1})
        for (int k = 0; k < 3; ++k) L = std::max(L, std::abs((*p)[k]));
    const double margin = 1e-10 * (L * L * L + 1e-30);
    const bool pos = q0 > margin && q1 > margin && q2 > margin && q3 > margin;
    const bool neg = q0 < -margin && q1 < -margin && q2 < -margin && q3 < -margin;
    return !(pos || neg);
}

struct Stats {
    double ms_xpbd = 0, ms_broad = 0, ms_push = 0, ms_ccd = 0;   // predict, broad phase, solve, ccd
    std::size_t vf = 0, ee = 0, hits = 0, reverted = 0, asked = 0, near = 0;
    long piercings = 0;            // the independent measure, summed over substeps
    long substeps = 0;
    bool exploded = false;
};

// Edges piercing triangles they are not part of, now. Through a grid, the
// same kind as the broad phase but over the positions alone.
inline long count_piercings(const Cloth& s, double cell) {
    Grid g{cell, {}};
    for (std::uint32_t f = 0; f < s.faces.size(); ++f) {
        const auto& F = s.faces[f];
        g.insert(box_of({s.parts[F[0]].x, s.parts[F[1]].x, s.parts[F[2]].x}, 0.0), f);
    }
    g.finish();
    long n = 0;
    std::vector<std::uint32_t> mark(s.faces.size(), ~0u);
    for (std::uint32_t e = 0; e < s.edges.size(); ++e) {
        const auto& E = s.edges[e];
        const V3 p = s.parts[E[0]].x, q = s.parts[E[1]].x;
        g.cells(box_of({p, q}, 0.0), [&](std::uint64_t k) {
            g.each(k, [&](std::uint32_t f) {
                if (mark[f] == e) return;
                mark[f] = e;
                const auto& F = s.faces[f];
                if (F[0] == E[0] || F[1] == E[0] || F[2] == E[0] || F[0] == E[1] || F[1] == E[1] || F[2] == E[1]) return;
                n += pierces(p, q, s.parts[F[0]].x, s.parts[F[1]].x, s.parts[F[2]].x);
            });
        });
    }
    return n;
}

inline void substep(Cloth& s, const Config& c, Stats& st) {
    using clk = std::chrono::steady_clock;
    const auto ms = [](auto a, auto b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
    const double dt = c.substep_dt;
    auto t0 = clk::now();
    // Predict.
    std::vector<V3> from(s.parts.size());
    for (std::size_t i = 0; i < s.parts.size(); ++i) {
        auto& p = s.parts[i];
        const V3 v{(p.x - p.x_prev) * (1 / dt)};
        from[i] = p.x;
        p.x_prev = p.x;
        p.x = V3{p.x + v * dt + V3{0.0, c.gravity, 0.0} * (dt * dt)};
    }
    for (auto& k : s.cons) k.reset();
    auto t1 = clk::now();

    // Cloth against cloth as constraints in the same iterations as the
    // cloth's own: every pair near enough to matter, its side read from
    // where it started, kept `thickness` apart on that side.
    const double cell = 1.5 * c.spacing;
    auto cand = broad_phase(s, from, c.thickness, cell);
    st.vf += cand.vf.size();
    st.ee += cand.ee.size();
    struct VF { std::uint32_t v, f; double side; };
    struct EE { std::uint32_t e, o; V3 side; };
    std::vector<VF> vfs;
    std::vector<EE> ees;
    // A pair farther apart at the start than the thickness and both moves
    // cannot come within it this substep: neighbours in one flat sheet are
    // most of what the boxes find, and none of them can.
    const auto move = [&](std::uint32_t i) { return V3{s.parts[i].x - from[i]}.norm(); };
    for (const auto& [v, f] : cand.vf) {
        const auto& F = s.faces[f];
        const auto w = closest_on_triangle(from[v], from[F[0]], from[F[1]], from[F[2]]);
        const V3 q{from[F[0]] * w[0] + from[F[1]] * w[1] + from[F[2]] * w[2]};
        const double reach = c.thickness + move(v) + std::max({move(F[0]), move(F[1]), move(F[2])});
        if (V3{from[v] - q}.norm() > reach) continue;
        const V3 n{V3{from[F[1]] - from[F[0]]}.cross(V3{from[F[2]] - from[F[0]]})};
        const double side = V3{from[v] - q}.dot(n);
        if (side != 0) vfs.push_back({v, f, side > 0 ? 1.0 : -1.0});
    }
    for (const auto& [e, o] : cand.ee) {
        const auto& A = s.edges[e];
        const auto& B = s.edges[o];
        const auto [a, b] = closest_segments(from[A[0]], from[A[1]], from[B[0]], from[B[1]]);
        const V3 d{V3{from[A[0]] + (from[A[1]] - from[A[0]]) * a} - V3{from[B[0]] + (from[B[1]] - from[B[0]]) * b}};
        const double len = d.norm();
        const double reach = c.thickness + std::max(move(A[0]), move(A[1])) + std::max(move(B[0]), move(B[1]));
        if (len > 0 && len <= reach) ees.push_back({e, o, V3{d * (1 / len)}});
    }
    const auto contacts = [&] {
        for (const auto& k : vfs) {
            const auto& F = s.faces[k.f];
            auto& P = s.parts[k.v];
            auto& A = s.parts[F[0]];
            auto& B = s.parts[F[1]];
            auto& C = s.parts[F[2]];
            const auto w = closest_on_triangle(P.x, A.x, B.x, C.x);
            const V3 q{A.x * w[0] + B.x * w[1] + C.x * w[2]};
            V3 n{V3{B.x - A.x}.cross(V3{C.x - A.x})};
            const double nl = n.norm();
            if (nl <= 0) continue;
            n = V3{n * (k.side / nl)};
            const V3 r{P.x - q};
            if (V3{r - n * r.dot(n)}.norm() > c.thickness) continue;   // beside the triangle
            const double C0 = r.dot(n) - c.thickness;
            if (C0 >= 0) continue;
            const double wsum = 1 + w[0] * w[0] + w[1] * w[1] + w[2] * w[2];
            const double l = -C0 / wsum;
            P.x = V3{P.x + n * l};
            A.x = V3{A.x - n * (l * w[0])};
            B.x = V3{B.x - n * (l * w[1])};
            C.x = V3{C.x - n * (l * w[2])};
        }
        for (const auto& k : ees) {
            const auto& A = s.edges[k.e];
            const auto& B = s.edges[k.o];
            auto& a0 = s.parts[A[0]];
            auto& a1 = s.parts[A[1]];
            auto& b0 = s.parts[B[0]];
            auto& b1 = s.parts[B[1]];
            const auto [u, v] = closest_segments(a0.x, a1.x, b0.x, b1.x);
            const V3 d{V3{a0.x + (a1.x - a0.x) * u} - V3{b0.x + (b1.x - b0.x) * v}};
            const double C0 = d.dot(k.side) - c.thickness;
            if (C0 >= 0 || d.norm() > 2 * c.thickness) continue;
            const double wa0 = 1 - u, wa1 = u, wb0 = 1 - v, wb1 = v;
            const double wsum = wa0 * wa0 + wa1 * wa1 + wb0 * wb0 + wb1 * wb1;
            if (wsum <= 0) continue;
            const double l = -C0 / wsum;
            a0.x = V3{a0.x + k.side * (l * wa0)};
            a1.x = V3{a1.x + k.side * (l * wa1)};
            b0.x = V3{b0.x - k.side * (l * wb0)};
            b1.x = V3{b1.x - k.side * (l * wb1)};
        }
    };
    const auto bodies = [&] {
        for (auto& p : s.parts) {
            const double r = p.x.norm();
            if (r < c.sphere + c.band && r > 0) p.x = V3{p.x * ((c.sphere + c.band) / r)};
            if (p.x[1] < c.floor + c.band) p.x[1] = c.floor + c.band;
        }
    };
    st.near += vfs.size() + ees.size();
    auto t2 = clk::now();
    for (int it = 0; it < c.iterations; ++it) {
        for (auto& k : s.cons) xpbd_solve_distance(k, s.parts, dt);
        contacts();
        bodies();
    }
    // Friction on the sphere and the floor: the move along the surface,
    // taken back by `friction`.
    for (auto& p : s.parts) {
        const V3 move{p.x - p.x_prev};
        const double r = p.x.norm();
        V3 n{};
        if (r <= c.sphere + c.band * 1.01) n = V3{p.x * (1 / r)};
        else if (p.x[1] <= c.floor + c.band * 1.01) n = V3{0.0, 1.0, 0.0};
        else continue;
        p.x = V3{p.x - V3{move - n * move.dot(n)} * c.friction};
    }
    auto t3 = clk::now();

    // Continuous: nothing crosses between the substep's start and end.
    if (c.ccd) {
        bool clear = false;
        // The candidates once, after the push moved things: a stop only
        // shortens a move, so they cover every later round too.
        cand = broad_phase(s, from, 0.0, cell);
        for (int round = 0; round < c.ccd_rounds && !clear; ++round) {
            std::size_t hits = 0;
            // Stopped short of the moment they would meet: the pair's four
            // vertices back along their moves to 0.8 of it.
            const auto stop = [&](std::initializer_list<std::uint32_t> ids, double toi) {
                const double k = std::max(0.0, 0.8 * toi);
                for (auto i : ids) s.parts[i].x = V3{from[i] + (s.parts[i].x - from[i]) * k};
            };
            for (const auto& [v, f] : cand.vf) {
                const auto& F = s.faces[f];
                if (!may_be_coplanar(from[v], from[F[0]], from[F[1]], from[F[2]], s.parts[v].x, s.parts[F[0]].x,
                                     s.parts[F[1]].x, s.parts[F[2]].x))
                    continue;
                ++st.asked;
                const Moves m{{from[v], from[F[0]], from[F[1]], from[F[2]]},
                              {s.parts[v].x, s.parts[F[0]].x, s.parts[F[1]].x, s.parts[F[2]].x}};
                const auto r = first_contact(point_chart(m, 0), triangle_chart(m, 1), 1e-6, std::size_t{1} << 16);
                if (r.hit) { ++hits; stop({v, F[0], F[1], F[2]}, r.toi); }
            }
            for (const auto& [e, o] : cand.ee) {
                const auto& A = s.edges[e];
                const auto& B = s.edges[o];
                if (!may_be_coplanar(from[A[0]], from[A[1]], from[B[0]], from[B[1]], s.parts[A[0]].x, s.parts[A[1]].x,
                                     s.parts[B[0]].x, s.parts[B[1]].x))
                    continue;
                ++st.asked;
                const Moves m{{from[A[0]], from[A[1]], from[B[0]], from[B[1]]},
                              {s.parts[A[0]].x, s.parts[A[1]].x, s.parts[B[0]].x, s.parts[B[1]].x}};
                const auto r = first_contact(edge_chart(m, 0), edge_chart(m, 2), 1e-6, std::size_t{1} << 16);
                if (r.hit) { ++hits; stop({A[0], A[1], B[0], B[1]}, r.toi); }
            }
            st.hits += hits;
            clear = hits == 0;
        }
        // Still crossing after the rounds: the whole cloth back to the
        // substep's start, where nothing crossed. Rare, and it costs the
        // substep's motion, never an intersection.
        if (!clear) {
            for (std::size_t i = 0; i < s.parts.size(); ++i) s.parts[i].x = from[i];
            ++st.reverted;
        }
    }
    auto t4 = clk::now();
    st.ms_xpbd += ms(t0, t1);
    st.ms_broad += ms(t1, t2);
    st.ms_push += ms(t2, t3);
    st.ms_ccd += ms(t3, t4);
    ++st.substeps;
    for (const auto& p : s.parts) {
        const V3 v{(p.x - p.x_prev) * (1.0 / c.substep_dt)};
        if (!std::isfinite(p.x[0] + p.x[1] + p.x[2]) || v.norm() > 50.0) st.exploded = true;
    }
}

}  // namespace napkins
