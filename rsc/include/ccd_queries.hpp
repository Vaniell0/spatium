#pragma once
// Queries for physics/mechanics/surface_ccd.hpp, of the kinds a policy is
// chosen between (rsc/tools/ccd_policy.cpp) and checked on
// (tests/test_rsc_ccd.cpp): tensor Bezier patches of order 1-3 with
// moving control points, as Chen et al. test on; a vertex and a triangle
// and two edges, the queries of cloth, at cloth's scale; a torus from its
// formula against a cubic patch. `kind` names where a query came from and
// is for reporting only -- a chooser sees SurfaceCcdFeatures, never it.

#include <spatium/physics/mechanics/surface_ccd.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <random>
#include <vector>

namespace rsc::ccd {

namespace mech = spatium::physics::mechanics;
using SV3 = spatium::Vec<double, 3>;
using Chart = mech::MovingChart<double>;

// A tensor Bezier patch of order n, control points P[i * (n + 1) + j].
struct Bezier {
    int n;
    std::vector<SV3> P;
    SV3 operator()(double u, double v) const {
        std::vector<SV3> row(n + 1);
        std::vector<SV3> tmp(n + 1);
        for (int i = 0; i <= n; ++i) {
            for (int j = 0; j <= n; ++j) tmp[j] = P[i * (n + 1) + j];
            for (int r = n; r > 0; --r)
                for (int j = 0; j < r; ++j) tmp[j] = SV3{tmp[j] * (1 - v) + tmp[j + 1] * v};
            row[i] = tmp[0];
        }
        for (int r = n; r > 0; --r)
            for (int i = 0; i < r; ++i) row[i] = SV3{row[i] * (1 - u) + row[i + 1] * u};
        return row[0];
    }
    mech::DerivativeBounds<double> bounds() const {
        mech::DerivativeBounds<double> b;
        const auto at = [&](int i, int j) { return P[i * (n + 1) + j]; };
        for (int i = 0; i <= n; ++i)
            for (int j = 0; j <= n; ++j) {
                if (i < n) b.u = std::max(b.u, n * SV3{at(i + 1, j) - at(i, j)}.norm());
                if (j < n) b.v = std::max(b.v, n * SV3{at(i, j + 1) - at(i, j)}.norm());
                if (i + 1 < n) b.uu = std::max(b.uu, n * (n - 1) * SV3{at(i + 2, j) - at(i + 1, j) * 2.0 + at(i, j)}.norm());
                if (j + 1 < n) b.vv = std::max(b.vv, n * (n - 1) * SV3{at(i, j + 2) - at(i, j + 1) * 2.0 + at(i, j)}.norm());
            }
        return b;
    }
};

inline Chart moving(const Bezier& p, const Bezier& w) {
    return {[p](double u, double v) { return p(u, v); }, [w](double u, double v) { return w(u, v); }, p.bounds(),
            w.bounds()};
}

struct SurfaceQuery {
    Chart a, b;
    int kind;   // for reporting only
};
inline const char* kKinds[] = {"bilinear", "quadratic", "cubic", "vertex-face", "edge-edge", "torus-patch"};

inline std::vector<SurfaceQuery> generate(std::uint64_t seed, int per_kind) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> uni(-1.0, 1.0);
    const auto rnd = [&] { return SV3{uni(rng), uni(rng), uni(rng)}; };
    const auto unit = [&] { SV3 d = rnd(); return SV3{d * (1.0 / d.norm())}; };
    std::vector<SurfaceQuery> out;
    for (int n = 1; n <= 3; ++n)
        for (int k = 0; k < per_kind; ++k) {
            const SV3 dir = unit();
            Bezier p1{n, {}}, w1{n, {}}, p2{n, {}}, w2{n, {}};
            for (int i = 0; i < (n + 1) * (n + 1); ++i) {
                p1.P.push_back(SV3{rnd() - dir});
                w1.P.push_back(SV3{rnd() + dir});
                p2.P.push_back(SV3{rnd() + dir});
                w2.P.push_back(SV3{rnd() - dir});
            }
            out.push_back({moving(p1, w1), moving(p2, w2), n - 1});
        }
    // Cloth-like: small triangles and edges, moving a fraction of their size.
    for (int k = 0; k < per_kind; ++k) {
        const SV3 A = rnd(), B{A + rnd() * 0.3}, C{A + rnd() * 0.3};
        const SV3 dA = SV3{rnd() * 0.2}, dB = SV3{rnd() * 0.2}, dC = SV3{rnd() * 0.2};
        const SV3 c{(A + B + C) * (1.0 / 3.0)};
        const SV3 x{c + rnd() * 0.3}, dx{SV3{c - x} * (0.5 + 1.5 * (uni(rng) + 1) / 2) + rnd() * 0.05};
        Chart tri;
        const auto f = [](SV3 a, SV3 b, SV3 c, double u, double v) { return SV3{a + (SV3{b - a} * (1 - v) + SV3{c - a} * v) * u}; };
        tri.p = [=](double u, double v) { return f(A, B, C, u, v); };
        tri.w = [=](double u, double v) { return f(dA, dB, dC, u, v); };
        tri.bp = {std::max(SV3{B - A}.norm(), SV3{C - A}.norm()), SV3{C - B}.norm(), 0, 0};
        tri.bw = {std::max(SV3{dB - dA}.norm(), SV3{dC - dA}.norm()), SV3{dC - dB}.norm(), 0, 0};
        Chart pt;
        pt.p = [=](double, double) { return x; };
        pt.w = [=](double, double) { return dx; };
        pt.u1 = pt.v1 = 0;
        out.push_back({pt, tri, 3});
    }
    for (int k = 0; k < per_kind; ++k) {
        const auto edge = [&](SV3 e0, SV3 e1, SV3 d0, SV3 d1) {
            Chart c;
            c.p = [=](double u, double) { return SV3{e0 + (e1 - e0) * u}; };
            c.w = [=](double u, double) { return SV3{d0 + (d1 - d0) * u}; };
            c.bp = {SV3{e1 - e0}.norm(), 0, 0, 0};
            c.bw = {SV3{d1 - d0}.norm(), 0, 0, 0};
            c.v1 = 0;
            return c;
        };
        const SV3 m = rnd(), da = unit(), db = unit(), off = unit();
        const SV3 a0{m - da * 0.2}, a1{m + da * 0.2};
        const SV3 b0{m + off * 0.2 - db * 0.2}, b1{m + off * 0.2 + db * 0.2};
        const SV3 push{off * (-0.2 - 0.3 * (uni(rng) + 1) / 2)};
        out.push_back({edge(a0, a1, rnd() * 0.02, rnd() * 0.02), edge(b0, b1, SV3{push + rnd() * 0.02}, SV3{push + rnd() * 0.02}), 4});
    }
    for (int k = 0; k < per_kind; ++k) {
        const double R = 0.6, r = 0.2, tp = 2 * std::numbers::pi;
        const SV3 c = rnd(), dir = unit();
        Chart torus;
        torus.p = [=](double u, double v) {
            const double a = tp * u, b = tp * v;
            return SV3{SV3{(R + r * std::cos(b)) * std::cos(a), (R + r * std::cos(b)) * std::sin(a), r * std::sin(b)} + c - dir};
        };
        const SV3 tw{dir * (1.0 + (uni(rng) + 1) / 2)};
        torus.w = [=](double, double) { return tw; };
        torus.bp = {tp * (R + r), tp * r, tp * tp * (R + r), tp * tp * r};
        Bezier p{3, {}}, w{3, {}};
        for (int i = 0; i < 16; ++i) { p.P.push_back(SV3{rnd() * 0.8 + c + dir}); w.P.push_back(SV3{rnd() * 0.3 - dir}); }
        out.push_back({torus, moving(p, w), 5});
    }
    return out;
}


}  // namespace rsc::ccd
