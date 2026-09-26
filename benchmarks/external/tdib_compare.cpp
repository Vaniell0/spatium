// Against TDIB-CCD (Chen et al., "A Time-Dependent Inclusion-Based Method
// for Continuous Collision Detection between Parametric Surfaces",
// SIGGRAPH Asia 2024, github.com/xw-c/TDIB-CCD) on its own random test:
// the same generator and seed, its time-dependent solver with OBBs -- its
// fastest -- and physics/mechanics/surface_ccd.hpp, per case, with a
// Newton witness for each contact as an upper bound on the truth.
//
// Not part of the build: TDIB-CCD carries no licence, so it is not copied
// here. Build against a clone (it needs Eigen and fmt):
//
//   git clone https://github.com/xw-c/TDIB-CCD && cd TDIB-CCD
//   g++ -std=c++23 -O3 -march=native -DNDEBUG -Icore -Iscene \
//       -I<eigen3> -I<fmt>/include -I<spatium>/include \
//       <spatium>/benchmarks/external/tdib_compare.cpp core/argsParser.cpp -lfmt -o compare
//   ./compare CASES WIDTH REFERENCE_WIDTH ORDER CHECK    e.g. 1000 1e-5 0 3 1
//
#include "scene_random.h"
#include <spatium/physics/mechanics/surface_ccd.hpp>
#include <print>

namespace mech = spatium::physics::mechanics;
using SV = spatium::Vec<double, 3>;

template<typename Obj>
static mech::DerivativeBounds<double> net_bounds(const Obj& o) {
    constexpr int n = Obj::order, m = n + 1;
    auto at = [&](int i, int j) { return o.ctrlp[i * m + j]; };
    mech::DerivativeBounds<double> b;
    for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= n; ++j) {
            if (i < n) b.u = std::max(b.u, n * (at(i + 1, j) - at(i, j)).norm());
            if (j < n) b.v = std::max(b.v, n * (at(i, j + 1) - at(i, j)).norm());
            if (i + 1 < n) b.uu = std::max(b.uu, n * (n - 1) * (at(i + 2, j) - 2 * at(i + 1, j) + at(i, j)).norm());
            if (j + 1 < n) b.vv = std::max(b.vv, n * (n - 1) * (at(i, j + 2) - 2 * at(i, j + 1) + at(i, j)).norm());
        }
    return b;
}

template<typename Obj>
static mech::MovingChart<double> chart(const Obj& pos, const Obj& vel) {
    mech::MovingChart<double> c;
    c.p = [&pos](double u, double v) { const Vector3d x = pos.evaluatePatchPoint(Array2d(u, v)); return SV{x[0], x[1], x[2]}; };
    c.w = [&vel](double u, double v) { const Vector3d x = vel.evaluatePatchPoint(Array2d(u, v)); return SV{x[0], x[1], x[2]}; };
    c.bp = net_bounds(pos);
    c.bw = net_bounds(vel);
    return c;
}


// A witness: some t and parameters where a point of one surface is on
// the other, found by Newton on (t, u2, v2) from grid points of the first.
// Any converged root is an upper bound on the true first contact.
template<typename Obj>
static double witness(const Obj& pa, const Obj& wa, const Obj& pb, const Obj& wb, double t_guess,
                      double ua_c, double va_c, double ub_c, double vb_c) {
    auto X = [](const Obj& p, const Obj& w, double u, double v, double t) {
        return Vector3d(p.evaluatePatchPoint(Array2d(u, v)) + t * w.evaluatePatchPoint(Array2d(u, v)));
    };
    double best = 2.0;
    auto solve = [&](double u, double v, double t, double u2, double v2) {
        const Vector3d a0 = pa.evaluatePatchPoint(Array2d(u, v)), a1 = wa.evaluatePatchPoint(Array2d(u, v));
        for (int it = 0; it < 40; ++it) {
            const Vector3d F = a0 + t * a1 - X(pb, wb, u2, v2, t);
            if (F.norm() < 1e-13) break;
            const double h = 1e-7;
            Eigen::Matrix3d J;
            J.col(0) = a1 - wb.evaluatePatchPoint(Array2d(u2, v2));
            J.col(1) = -(X(pb, wb, u2 + h, v2, t) - X(pb, wb, u2 - h, v2, t)) / (2 * h);
            J.col(2) = -(X(pb, wb, u2, v2 + h, t) - X(pb, wb, u2, v2 - h, t)) / (2 * h);
            const Vector3d d = J.colPivHouseholderQr().solve(-F);
            t += d[0]; u2 += d[1]; v2 += d[2];
            if (!std::isfinite(t) || std::abs(u2) > 3 || std::abs(v2) > 3) return;
        }
        if (t < 0 || t > 1 || u2 < 0 || u2 > 1 || v2 < 0 || v2 > 1) return;
        if ((a0 + t * a1 - X(pb, wb, u2, v2, t)).norm() < 1e-11) best = std::min(best, t);
    };
    const int n = 24;
    for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= n; ++j) solve(double(i) / n, double(j) / n, t_guess, ub_c, vb_c);
    // Densely around the answer.
    for (int i = -8; i <= 8; ++i)
        for (int j = -8; j <= 8; ++j) {
            const double u = std::clamp(ua_c + i * 2e-4, 0.0, 1.0), v = std::clamp(va_c + j * 2e-4, 0.0, 1.0);
            solve(u, v, t_guess, ub_c, vb_c);
        }
    return best;
}

template<typename Obj>
int run(int argc, char** argv) {
    const int kase = argc > 1 ? std::atoi(argv[1]) : 100;
    const double delta = argc > 2 ? std::atof(argv[2]) : 1e-5;
    const double ref_delta = argc > 3 ? std::atof(argv[3]) : 0.0;
    double late_ours = -1, late_td = -1;
    const bool check = argc > 5 && std::atoi(argv[5]) != 0;
    int witnessed = 0;
    double late_w_ours = -1, late_w_td = -1, gap_ours = 0, gap_td = 0;
    std::srand(0);
    using clk = std::chrono::steady_clock;
    double ms_td = 0, ms_ours = 0, worst_later = -1, worst_earlier = 0;
    int both = 0, disagree = 0;
    std::size_t pairs = 0;
    for (int k = 0; k < kase; ++k) {
        Obj pos1, pos2, vel1, vel2;
        generatePatchPair<Obj>(pos1.ctrlp, vel1.ctrlp, pos2.ctrlp, vel2.ctrlp);
        Array2d uv1, uv2;
        auto t0 = clk::now();
        const double td = SolverTD<Obj, Obj, RecParamBound, RecParamBound>::solveCCD(
            pos1, vel1, pos2, vel2, uv1, uv2, BoundingBoxType::OBB, delta);
        auto t1 = clk::now();
        const auto a = chart(pos1, vel1), b = chart(pos2, vel2);
        const auto r = mech::first_contact(a, b, delta);
        auto t2 = clk::now();
        ms_td += std::chrono::duration<double, std::milli>(t1 - t0).count();
        ms_ours += std::chrono::duration<double, std::milli>(t2 - t1).count();
        pairs += r.pairs;
        if (check) {
            const double wab = witness(pos1, vel1, pos2, vel2, r.toi, r.u_a, r.v_a, r.u_b, r.v_b);
            const double wba = witness(pos2, vel2, pos1, vel1, r.toi, r.u_b, r.v_b, r.u_a, r.v_a);
            const double wt = std::min(wab, wba);
            if (r.hit && wt <= 1) {
                ++witnessed;
                late_w_ours = std::max(late_w_ours, r.toi - wt);
                late_w_td = std::max(late_w_td, td - wt);
                gap_ours = std::max(gap_ours, wt - r.toi);
                gap_td = std::max(gap_td, wt - td);
            }
            if (!r.hit && wt <= 1) std::println("case {}: ours misses, a witness at {}", k, wt);
        }
        if (ref_delta > 0) {
            Array2d w1, w2;
            const double ref = SolverTD<Obj, Obj, RecParamBound, RecParamBound>::solveCCD(
                pos1, vel1, pos2, vel2, w1, w2, BoundingBoxType::OBB, ref_delta);
            if (ref >= 0) {
                // The reference is itself early by at most ref_delta.
                late_ours = std::max(late_ours, r.toi - (ref + ref_delta));
                late_td = std::max(late_td, td - (ref + ref_delta));
            } else if (r.hit) std::println("case {}: reference misses, ours hits at {}", k, r.toi);
        }
        const bool hit_td = td >= 0;
        if (hit_td != r.hit) { ++disagree; std::println("case {} disagree: td {} ours hit={} toi={}", k, td, r.hit, r.toi); continue; }
        if (hit_td) {
            ++both;
            worst_later = std::max(worst_later, r.toi - td);
            worst_earlier = std::max(worst_earlier, td - r.toi);
        }
    }
    std::println("{} cases, {} contacts both, {} disagree", kase, both, disagree);
    std::println("ms a query: TDIB td/obb {:.3f}, ours {:.3f}  (x{:.1f}); pairs a query {:.0f}", ms_td / kase, ms_ours / kase,
                 ms_td / ms_ours, double(pairs) / kase);
    std::println("toi ours - td: latest {:.3g}, earliest {:.3g}", worst_later, -worst_earlier);
    if (check)
        std::println("witnessed {}: past the witness ours {:.3g}, TDIB {:.3g} (<= 0 is never late); widest gap to it ours {:.3g}, TDIB {:.3g}",
                     witnessed, late_w_ours, late_w_td, gap_ours, gap_td);
    if (ref_delta > 0) std::println("against a {:g} reference: ours at most {:.3g} past it, TDIB {:.3g} (<= 0 is never late)", ref_delta, late_ours, late_td);
    return 0;
}

int main(int argc, char** argv) {
    const int order = argc > 4 ? std::atoi(argv[4]) : 3;
    if (order == 1) return run<RecLinearBezier>(argc, argv);
    if (order == 2) return run<RecQuadBezier>(argc, argv);
    return run<RecCubicBezier>(argc, argv);
}
