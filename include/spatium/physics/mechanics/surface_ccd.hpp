#pragma once
// First contact between two moving surfaces.
//
// Each surface is a chart whose points move linearly over the step:
// x(u, v, t) = p(u, v) + t w(u, v), t in [0, 1] -- a Bezier patch whose
// control points move linearly is one, and so is any chart paired with a
// chart of velocities. The search runs over pairs of cells, one from each
// surface, and an interval of time, and it never splits time unless the
// cells are already as small as asked: each pair's interval is cut down
// to where contact is still possible, exactly, by two floors that are
// both functions of t alone once the cells are fixed.
//
//   A ball. A cell's image at time t lies within r_p + t r_w of the
//   image of its centre (r_p, r_w: the reach of the position and the
//   velocity charts over the cell), so two cells can meet only where
//   |a + t b| <= r0 + t r1 -- a quadratic in t, solved.
//
//   A slab. Along any unit direction n, a cell's height at time t lies
//   between the least and the greatest of its corners' heights -- each a
//   line in t, since a corner moves linearly -- widened by the bilinear
//   interpolation bound (du^2 |f_uu| + dv^2 |f_vv|) / 8 of the chart at
//   that time, itself linear in t. Two cells are apart along n where one's
//   least line clears the other's greatest: the minimum of a few lines, so
//   the times it holds form one interval, found line by line. Any
//   direction is valid; which are tried is the policy's (SurfaceCcdPolicy):
//   the cells' normals across their diagonals, the line between their
//   centres, the cross products of their sides -- what two edges, which
//   have no normal, are apart along -- and the world axes.
//
// What is left of the interval after both is where the pair may touch.
// The earliest pair is refined: the cell of larger reach is halved, or,
// once both cells are within `width` in parameters, the interval of time.
// A pair whose cells and interval are all within `width` answers with the
// start of its interval: no contact before it is possible, so the answer
// is early or exact, never late -- the criterion of Chen et al.'s
// time-dependent inclusion method (SIGGRAPH Asia 2024), kept so the two
// can be compared on the same queries.
//
// The bounds a surface has to supply are the ones any chart can: first
// and second derivatives bounded over the domain, for position and for
// velocity. A Bezier patch reads them off its control net; a torus or a
// Klein bottle from its formula. Nothing here needs control points.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/concepts.hpp>
#  include <algorithm>
#  include <array>
#  include <cmath>
#  include <cstddef>
#  include <cstdint>
#  include <functional>
#  include <limits>
#  include <queue>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::physics::mechanics {

// Bounds on a chart's derivatives over its domain: |f_u|, |f_v|, |f_uu|,
// |f_vv|.
template<Scalar T>
struct DerivativeBounds {
    T u = T{0}, v = T{0}, uu = T{0}, vv = T{0};
};

// A surface moving linearly over the step: position p and velocity w as
// charts on one parameter rectangle, each with its bounds.
template<Scalar T>
struct MovingChart {
    std::function<Vec<T, 3>(T, T)> p, w;
    DerivativeBounds<T> bp, bw;
    T u0 = T{0}, u1 = T{1}, v0 = T{0}, v1 = T{1};
};

// The decisions of the search that no correctness depends on -- which
// directions a slab is tried along, which cell is halved, when time is
// split -- as data, so they can be chosen per kind of query rather than
// once for all of them. Every choice keeps the rule; only the cost moves.
struct SurfaceCcdPolicy {
    enum Axis : std::uint32_t {
        kNormals = 1,        // each cell's normal across its diagonals
        kCentres = 2,        // the line between the cells' centres
        kTangentCross = 4,   // cross products of one cell's sides with the other's --
                             // the axis two edges are apart along, where neither has a normal
        kWorld = 8,          // x, y, z
    };
    std::uint32_t axes = kNormals | kCentres | kTangentCross;
    bool ball = true;
    // Halve by reach (position plus velocity over the interval) or by the
    // widest parameter side.
    bool split_by_reach = true;
};

template<Scalar T>
struct SurfaceContact {
    bool hit = false;
    T toi = T{1};
    T u_a = T{0}, v_a = T{0}, u_b = T{0}, v_b = T{0};   // centres of the answering cells
    std::size_t evaluations = 0;   // of p and w, corners and centres
    std::size_t pairs = 0;         // pairs of cells bounded
};

namespace detail {

// A cell of a moving chart: its rectangle, position and velocity at its
// centre and corners, and its children once made.
template<Scalar T>
struct MovingCell {
    T u0, u1, v0, v1;
    Vec<T, 3> pc, wc;
    std::array<Vec<T, 3>, 4> pk, wk;   // (u0,v0), (u1,v0), (u0,v1), (u1,v1)
    std::int32_t first_child = -1;
};

template<Scalar T>
class MovingCellTree {
public:
    explicit MovingCellTree(const MovingChart<T>& c) : c_(c) {
        std::array<Vec<T, 3>, 4> pk{p(c.u0, c.v0), p(c.u1, c.v0), p(c.u0, c.v1), p(c.u1, c.v1)};
        std::array<Vec<T, 3>, 4> wk{w(c.u0, c.v0), w(c.u1, c.v0), w(c.u0, c.v1), w(c.u1, c.v1)};
        cells_.push_back(make(c.u0, c.u1, c.v0, c.v1, pk, wk));
    }
    const MovingCell<T>& cell(std::size_t i) const { return cells_[i]; }
    std::size_t evaluations() const { return evaluations_; }

    T du(std::size_t i) const { return cells_[i].u1 - cells_[i].u0; }
    T dv(std::size_t i) const { return cells_[i].v1 - cells_[i].v0; }
    T width(std::size_t i) const { return std::max(du(i), dv(i)); }
    // Reach of position and of velocity from the centre over the cell.
    T reach_p(std::size_t i) const { return (c_.bp.u * du(i) + c_.bp.v * dv(i)) / T{2}; }
    T reach_w(std::size_t i) const { return (c_.bw.u * du(i) + c_.bw.v * dv(i)) / T{2}; }
    // The interpolation widening of position and velocity.
    T bend_p(std::size_t i) const { return (du(i) * du(i) * c_.bp.uu + dv(i) * dv(i) * c_.bp.vv) / T{8}; }
    T bend_w(std::size_t i) const { return (du(i) * du(i) * c_.bw.uu + dv(i) * dv(i) * c_.bw.vv) / T{8}; }

    std::uint32_t children(std::size_t i) {
        if (cells_[i].first_child >= 0) return static_cast<std::uint32_t>(cells_[i].first_child);
        const MovingCell<T> c = cells_[i];
        const auto first = static_cast<std::int32_t>(cells_.size());
        const T um = (c.u0 + c.u1) / T{2}, vm = (c.v0 + c.v1) / T{2};
        // Halve across the side along which the chart reaches farther.
        const bool across_u = c_.bp.u * (c.u1 - c.u0) + c_.bw.u * (c.u1 - c.u0) >=
                              c_.bp.v * (c.v1 - c.v0) + c_.bw.v * (c.v1 - c.v0);
        const auto& P = c.pk;
        const auto& W = c.wk;
        if (across_u) {
            const Vec<T, 3> pb = p(um, c.v0), pt = p(um, c.v1), wb = w(um, c.v0), wt = w(um, c.v1);
            cells_.push_back(make(c.u0, um, c.v0, c.v1, {P[0], pb, P[2], pt}, {W[0], wb, W[2], wt}));
            cells_.push_back(make(um, c.u1, c.v0, c.v1, {pb, P[1], pt, P[3]}, {wb, W[1], wt, W[3]}));
        } else {
            const Vec<T, 3> pl = p(c.u0, vm), pr = p(c.u1, vm), wl = w(c.u0, vm), wr = w(c.u1, vm);
            cells_.push_back(make(c.u0, c.u1, c.v0, vm, {P[0], P[1], pl, pr}, {W[0], W[1], wl, wr}));
            cells_.push_back(make(c.u0, c.u1, vm, c.v1, {pl, pr, P[2], P[3]}, {wl, wr, W[2], W[3]}));
        }
        cells_[i].first_child = first;
        return static_cast<std::uint32_t>(first);
    }

private:
    Vec<T, 3> p(T u, T v) { ++evaluations_; return c_.p(u, v); }
    Vec<T, 3> w(T u, T v) { ++evaluations_; return c_.w(u, v); }
    MovingCell<T> make(T u0, T u1, T v0, T v1, const std::array<Vec<T, 3>, 4>& pk,
                       const std::array<Vec<T, 3>, 4>& wk) {
        const T uc = (u0 + u1) / T{2}, vc = (v0 + v1) / T{2};
        return MovingCell<T>{u0, u1, v0, v1, p(uc, vc), w(uc, vc), pk, wk, -1};
    }

    const MovingChart<T>& c_;
    std::vector<MovingCell<T>> cells_;
    std::size_t evaluations_ = 0;
};

// The part of [t0, t1] where alpha_m + beta_m t > 0 for every m: one
// interval, possibly empty (lo > hi).
template<Scalar T, std::size_t N>
std::pair<T, T> all_positive(const std::array<T, N>& alpha, const std::array<T, N>& beta, T t0, T t1) {
    T lo = t0, hi = t1;
    for (std::size_t m = 0; m < N; ++m) {
        const T a = alpha[m], b = beta[m];
        if (b == T{0}) {
            if (!(a > T{0})) return {T{1}, T{0}};
            continue;
        }
        const T root = -a / b;
        if (b > T{0}) lo = std::max(lo, root);   // positive after the root
        else hi = std::min(hi, root);             // positive before it
        if (lo > hi) return {T{1}, T{0}};
    }
    return {lo, hi};
}

// Cut [t0, t1] by a set of times proved apart: the ends it covers go,
// an interior piece is left in (the hull stays valid).
template<Scalar T>
bool cut(T& t0, T& t1, std::pair<T, T> apart) {
    const auto [lo, hi] = apart;
    if (lo > hi) return true;
    if (lo <= t0 && hi >= t1) return false;       // apart throughout
    if (lo <= t0 && hi > t0) t0 = hi;
    else if (hi >= t1 && lo < t1) t1 = lo;
    return t0 <= t1;
}

}  // namespace detail

// The first time two moving surfaces touch, to `width` in parameters and
// time. `budget` caps the pairs bounded; running out answers the earliest
// interval still open -- early, never late.
template<Scalar T>
SurfaceContact<T> first_contact(const MovingChart<T>& a, const MovingChart<T>& b, T width = T{1e-5},
                                std::size_t budget = std::size_t{1} << 24, const SurfaceCcdPolicy& policy = {}) {
    using std::abs; using std::sqrt;
    detail::MovingCellTree<T> ta(a), tb(b);
    SurfaceContact<T> out;
    const T ulp = std::numeric_limits<T>::epsilon() * T{64};

    // Earliest first; among equal starts -- a contact along a whole line
    // or patch at one moment, resting or sliding -- the smallest pair
    // first, so the search goes down one of them instead of across all.
    struct Item { std::uint32_t ia, ib; T t0, t1, size; };
    const auto later = [](const Item& x, const Item& y) { return x.t0 != y.t0 ? x.t0 > y.t0 : x.size > y.size; };
    std::priority_queue<Item, std::vector<Item>, decltype(later)> open(later);

    // Narrow [t0, t1] to where cells ia, ib may touch; false if nowhere.
    const auto narrow = [&](std::uint32_t ia, std::uint32_t ib, T& t0, T& t1) {
        ++out.pairs;
        const auto& A = ta.cell(ia);
        const auto& B = tb.cell(ib);
        // Floating-point slack, relative to the coordinates involved.
        const T scale = A.pc.norm() + B.pc.norm() + A.wc.norm() + B.wc.norm();
        const T pad = ulp * (scale + T{1});

        // Ball: |d0 + t d1| <= r0 + t r1.
        const Vec<T, 3> d0{A.pc - B.pc}, d1{A.wc - B.wc};
        const T r0 = ta.reach_p(ia) + tb.reach_p(ib) + pad, r1 = ta.reach_w(ia) + tb.reach_w(ib);
        if (policy.ball) {
            // q(t) = |d0 + t d1|^2 - (r0 + t r1)^2 > 0 is apart.
            const T qa = d1.dot(d1) - r1 * r1, qb = T{2} * (d0.dot(d1) - r0 * r1), qc = d0.dot(d0) - r0 * r0;
            const auto q = [&](T t) { return (qa * t + qb) * t + qc; };
            // Apart at an end and up to the nearest root inside.
            const T disc = qb * qb - T{4} * qa * qc;
            std::array<T, 2> roots{};
            int n = 0;
            if (qa != T{0} && disc >= T{0}) {
                const T s = sqrt(disc);
                const T k = -(qb + (qb >= T{0} ? s : -s)) / T{2};
                roots = {k / qa, k != T{0} ? qc / k : k / qa};
                if (roots[0] > roots[1]) std::swap(roots[0], roots[1]);
                n = 2;
            } else if (qa == T{0} && qb != T{0}) {
                roots[0] = -qc / qb;
                n = 1;
            }
            // Between consecutive roots the sign of q is fixed; take the
            // pieces of [t0, t1] where it is positive at both their ends.
            std::array<T, 4> knots{t0, t0, t0, t1};
            int m = 1;
            for (int k = 0; k < n; ++k)
                if (roots[k] > t0 && roots[k] < t1) knots[m++] = roots[k];
            knots[m++] = t1;
            T lo = t1, hi = t0;   // hull of where q <= 0
            for (int k = 0; k + 1 < m; ++k) {
                const T mid = (knots[k] + knots[k + 1]) / T{2};
                if (!(q(mid) > T{0})) { lo = std::min(lo, knots[k]); hi = std::max(hi, knots[k + 1]); }
            }
            if (lo > hi) return false;
            t0 = std::max(t0, lo);
            t1 = std::min(t1, hi);
        }

        // Slabs along the policy's directions.
        const T tc = (t0 + t1) / T{2};
        const auto at = [tc](const Vec<T, 3>& p, const Vec<T, 3>& w) { return Vec<T, 3>{p + w * tc}; };
        const auto normal = [&](const auto& C) {
            const Vec<T, 3> d03{at(C.pk[3], C.wk[3]) - at(C.pk[0], C.wk[0])};
            const Vec<T, 3> d12{at(C.pk[2], C.wk[2]) - at(C.pk[1], C.wk[1])};
            return Vec<T, 3>{d03.cross(d12)};
        };
        std::array<Vec<T, 3>, 12> dirs{};
        int nd = 0;
        if (policy.axes & SurfaceCcdPolicy::kNormals) { dirs[nd++] = normal(A); dirs[nd++] = normal(B); }
        if (policy.axes & SurfaceCcdPolicy::kCentres) dirs[nd++] = Vec<T, 3>{d0 + d1 * tc};
        if (policy.axes & SurfaceCcdPolicy::kTangentCross) {
            const auto side_u = [&](const auto& C) { return Vec<T, 3>{at(C.pk[1], C.wk[1]) - at(C.pk[0], C.wk[0])}; };
            const auto side_v = [&](const auto& C) { return Vec<T, 3>{at(C.pk[2], C.wk[2]) - at(C.pk[0], C.wk[0])}; };
            const std::array<Vec<T, 3>, 2> sa{side_u(A), side_v(A)}, sb{side_u(B), side_v(B)};
            for (const auto& x : sa)
                for (const auto& y : sb) dirs[nd++] = Vec<T, 3>{x.cross(y)};
        }
        if (policy.axes & SurfaceCcdPolicy::kWorld) {
            dirs[nd++] = Vec<T, 3>{T{1}, T{0}, T{0}};
            dirs[nd++] = Vec<T, 3>{T{0}, T{1}, T{0}};
            dirs[nd++] = Vec<T, 3>{T{0}, T{0}, T{1}};
        }
        const T bpa = ta.bend_p(ia), bwa = ta.bend_w(ia), bpb = tb.bend_p(ib), bwb = tb.bend_w(ib);
        for (int di = 0; di < nd; ++di) {
            const auto& raw = dirs[di];
            const T len = raw.norm();
            if (!(len > T{0}) || !std::isfinite(len)) continue;
            const Vec<T, 3> n{raw * (T{1} / len)};
            // A's least minus B's greatest, for every corner pair: 16 lines.
            std::array<T, 16> alpha{}, beta{}, alpha2{}, beta2{};
            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j) {
                    const T ha = n.dot(A.pk[i]), va = n.dot(A.wk[i]);
                    const T hb = n.dot(B.pk[j]), vb = n.dot(B.wk[j]);
                    alpha[4 * i + j] = ha - hb - bpa - bpb - pad;
                    beta[4 * i + j] = va - vb - bwa - bwb;
                    alpha2[4 * i + j] = hb - ha - bpa - bpb - pad;
                    beta2[4 * i + j] = vb - va - bwa - bwb;
                }
            if (!detail::cut(t0, t1, detail::all_positive(alpha, beta, t0, t1))) return false;
            if (!detail::cut(t0, t1, detail::all_positive(alpha2, beta2, t0, t1))) return false;
        }
        return true;
    };

    {
        T t0 = T{0}, t1 = T{1};
        if (narrow(0, 0, t0, t1)) open.push(Item{0, 0, t0, t1, ta.width(0) + tb.width(0)});
    }
    const auto finish = [&](bool hit, T toi, const Item* it) {
        out.hit = hit;
        out.toi = toi;
        if (it) {
            const auto& A = ta.cell(it->ia);
            const auto& B = tb.cell(it->ib);
            out.u_a = (A.u0 + A.u1) / T{2}; out.v_a = (A.v0 + A.v1) / T{2};
            out.u_b = (B.u0 + B.u1) / T{2}; out.v_b = (B.v0 + B.v1) / T{2};
        }
        out.evaluations = ta.evaluations() + tb.evaluations();
        return out;
    };
    while (!open.empty()) {
        const Item it = open.top();
        if (out.pairs >= budget) return finish(true, it.t0, &it);
        open.pop();
        const T wa = ta.width(it.ia), wb = tb.width(it.ib);
        if (wa < width && wb < width && it.t1 - it.t0 < width) return finish(true, it.t0, &it);
        const auto push = [&](std::uint32_t ia, std::uint32_t ib, T t0, T t1) {
            if (narrow(ia, ib, t0, t1)) open.push(Item{ia, ib, t0, t1, ta.width(ia) + tb.width(ib) + (t1 - t0)});
        };
        if (wa < width && wb < width) {
            const T tm = (it.t0 + it.t1) / T{2};
            push(it.ia, it.ib, it.t0, tm);
            push(it.ia, it.ib, tm, it.t1);
            continue;
        }
        // Halve the cell that reaches farther over the interval.
        const T ra = policy.split_by_reach ? ta.reach_p(it.ia) + it.t1 * ta.reach_w(it.ia) : wa;
        const T rb = policy.split_by_reach ? tb.reach_p(it.ib) + it.t1 * tb.reach_w(it.ib) : wb;
        if ((ra >= rb && wa >= width) || wb < width) {
            const auto f = ta.children(it.ia);
            push(f, it.ib, it.t0, it.t1);
            push(f + 1, it.ib, it.t0, it.t1);
        } else {
            const auto f = tb.children(it.ib);
            push(it.ia, f, it.t0, it.t1);
            push(it.ia, f + 1, it.t0, it.t1);
        }
    }
    return finish(false, T{1}, nullptr);
}

}  // namespace spatium::physics::mechanics
