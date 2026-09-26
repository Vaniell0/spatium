// First contact between two moving surfaces (surface_ccd.hpp), held to the
// one-sided rule: an answer may be early by the width asked, never late,
// and a miss is reported only where no contact exists. Lateness is judged
// against a witness -- a time at which a point of one surface is found,
// by Newton, lying on the other -- so the check does not lean on the
// method it checks.
#include <catch2/catch_test_macros.hpp>

#include <spatium/physics/mechanics/surface_ccd.hpp>

#include <array>
#include <cmath>
#include <format>
#include <numbers>
#include <random>

using namespace spatium;
using namespace spatium::physics::mechanics;
using V3 = Vec<double, 3>;

namespace {

// A tensor-product cubic Bezier patch: 16 control points, P[4 * i + j],
// i along u.
struct Cubic {
    std::array<V3, 16> P;
    static V3 casteljau(V3 a, V3 b, V3 c, V3 d, double t) {
        const auto l = [t](const V3& x, const V3& y) { return V3{x * (1 - t) + y * t}; };
        const V3 ab = l(a, b), bc = l(b, c), cd = l(c, d);
        return l(l(ab, bc), l(bc, cd));
    }
    V3 operator()(double u, double v) const {
        std::array<V3, 4> q;
        for (int i = 0; i < 4; ++i) q[i] = casteljau(P[4 * i], P[4 * i + 1], P[4 * i + 2], P[4 * i + 3], v);
        return casteljau(q[0], q[1], q[2], q[3], u);
    }
    // From the control net: a Bezier's derivative is a Bezier of the
    // differences, and a norm is at most its largest control point's.
    DerivativeBounds<double> bounds() const {
        DerivativeBounds<double> b;
        const auto at = [this](int i, int j) { return P[4 * i + j]; };
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) {
                if (i < 3) b.u = std::max(b.u, 3 * V3{at(i + 1, j) - at(i, j)}.norm());
                if (j < 3) b.v = std::max(b.v, 3 * V3{at(i, j + 1) - at(i, j)}.norm());
                if (i < 2) b.uu = std::max(b.uu, 6 * V3{at(i + 2, j) - at(i + 1, j) * 2.0 + at(i, j)}.norm());
                if (j < 2) b.vv = std::max(b.vv, 6 * V3{at(i, j + 2) - at(i, j + 1) * 2.0 + at(i, j)}.norm());
            }
        return b;
    }
};

MovingChart<double> moving(const Cubic& p, const Cubic& w) {
    return {[p](double u, double v) { return p(u, v); }, [w](double u, double v) { return w(u, v); }, p.bounds(),
            w.bounds()};
}

double det3(const V3& a, const V3& b, const V3& c) { return a.dot(V3{b.cross(c)}); }

// The earliest t at which some grid point of `a` is found on `b`, or 2.
double witness(const MovingChart<double>& a, const MovingChart<double>& b, double t0, double ub, double vb) {
    const auto X = [](const MovingChart<double>& c, double u, double v, double t) {
        return V3{c.p(u, v) + c.w(u, v) * t};
    };
    double best = 2;
    const int n = 16;
    for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= n; ++j) {
            const double u = double(i) / n, v = double(j) / n;
            const V3 a0 = a.p(u, v), a1 = a.w(u, v);
            double t = t0, u2 = ub, v2 = vb;
            bool ok = true;
            for (int it = 0; it < 40 && ok; ++it) {
                const V3 F{a0 + a1 * t - X(b, u2, v2, t)};
                if (F.norm() < 1e-13) break;
                const double h = 1e-7;
                const V3 c0{a1 - b.w(u2, v2)};
                const V3 c1{(X(b, u2 + h, v2, t) - X(b, u2 - h, v2, t)) * (-0.5 / h)};
                const V3 c2{(X(b, u2, v2 + h, t) - X(b, u2, v2 - h, t)) * (-0.5 / h)};
                const double D = det3(c0, c1, c2);
                if (std::abs(D) < 1e-300) { ok = false; break; }
                const V3 r{F * -1.0};
                t += det3(r, c1, c2) / D;
                u2 += det3(c0, r, c2) / D;
                v2 += det3(c0, c1, r) / D;
                ok = std::isfinite(t) && std::abs(u2) < 3 && std::abs(v2) < 3;
            }
            if (!ok || t < 0 || t > 1 || u2 < 0 || u2 > 1 || v2 < 0 || v2 > 1) continue;
            if (V3{a0 + a1 * t - X(b, u2, v2, t)}.norm() < 1e-11) best = std::min(best, t);
        }
    return best;
}

Cubic flat(double z, double s) {
    Cubic c;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) c.P[4 * i + j] = V3{s * i, s * j, z};
    return c;
}
Cubic constant(const V3& v) {
    Cubic c;
    c.P.fill(v);
    return c;
}

}  // namespace

TEST_CASE("Two parallel patches meet when their gap closes", "[physics][ccd][surface]") {
    // A plane at z = 0 at rest and one at z = 1 falling at 2: they meet at 0.5.
    const auto a = moving(flat(0.0, 1.0), constant(V3{}));
    const auto b = moving(flat(1.0, 1.0), constant(V3{0.0, 0.0, -2.0}));
    const auto r = first_contact(a, b, 1e-6);
    INFO(std::format("hit={} toi={} pairs={}", r.hit, r.toi, r.pairs));
    REQUIRE(r.hit);
    CHECK(r.toi <= 0.5);
    CHECK(r.toi > 0.5 - 1e-6);

    // Rising instead, they never meet, and that is proved.
    const auto away = moving(flat(1.0, 1.0), constant(V3{0.0, 0.0, 2.0}));
    CHECK_FALSE(first_contact(a, away, 1e-6).hit);
}

// A torus is no polynomial: its bounds come from its formula. Falling onto
// a plane at 2, its lowest circle reaches it at t = (2 - 0.3) / 2.
TEST_CASE("A falling torus touches a plane when its lowest circle does", "[physics][ccd][surface]") {
    const double R = 1.0, r = 0.3, two_pi = 2.0 * std::numbers::pi;
    MovingChart<double> torus;
    torus.p = [=](double u, double v) {
        const double a = two_pi * u, b = two_pi * v;
        return V3{(R + r * std::cos(b)) * std::cos(a), (R + r * std::cos(b)) * std::sin(a), 2.0 + r * std::sin(b)};
    };
    torus.w = [](double, double) { return V3{0.0, 0.0, -2.0}; };
    // Chain rule through u, v in [0, 1]: each derivative gains 2 pi.
    torus.bp = {two_pi * (R + r), two_pi * r, two_pi * two_pi * (R + r), two_pi * two_pi * r};
    const auto plane = moving(flat(0.0, 4.0 / 3.0), constant(V3{}));
    auto shifted = plane;
    shifted.p = [&plane](double u, double v) { return V3{plane.p(u, v) - V3{2.0, 2.0, 0.0}}; };
    // A whole circle touches at once: a flat minimum, refined across the
    // circle, so this asks the width the patch test asks, not finer.
    const auto c = first_contact(torus, shifted, 1e-5);
    const double truth = (2.0 - r) / 2.0;
    INFO(std::format("hit={} toi={} truth={} pairs={}", c.hit, c.toi, truth, c.pairs));
    REQUIRE(c.hit);
    CHECK(c.toi <= truth);
    CHECK(c.toi > truth - 1e-5);
}

// Random patch pairs of the kind Chen et al. (SIGGRAPH Asia 2024) test on:
// control points random in a unit box, the two patches pushed towards
// each other. Never later than a witness, and no miss where one exists.
TEST_CASE("Random moving patches: never later than a witness", "[physics][ccd][surface][fuzz]") {
    std::mt19937_64 rng(20260926);
    std::uniform_real_distribution<double> uni(-1.0, 1.0);
    const auto rnd = [&] { return V3{uni(rng), uni(rng), uni(rng)}; };
    constexpr int kCases = 120;
    int hits = 0, witnessed = 0;
    double worst_late = -1, widest = 0;
    for (int k = 0; k < kCases; ++k) {
        V3 dir = rnd();
        dir = V3{dir * (1.0 / dir.norm())};
        Cubic p1, w1, p2, w2;
        for (int i = 0; i < 16; ++i) {
            p1.P[i] = V3{rnd() - dir};
            w1.P[i] = V3{rnd() + dir};
            p2.P[i] = V3{rnd() + dir};
            w2.P[i] = V3{rnd() - dir};
        }
        const auto a = moving(p1, w1), b = moving(p2, w2);
        const auto r = first_contact(a, b, 1e-5);
        const double w = std::min(witness(a, b, r.toi, r.u_b, r.v_b), witness(b, a, r.toi, r.u_a, r.v_a));
        INFO(std::format("case {}: hit={} toi={:.12f} witness={:.12f}", k, r.hit, r.toi, w));
        if (w <= 1) REQUIRE(r.hit);
        if (!r.hit) continue;
        ++hits;
        if (w > 1) continue;
        ++witnessed;
        REQUIRE(r.toi <= w + 1e-12);
        worst_late = std::max(worst_late, r.toi - w);
        widest = std::max(widest, w - r.toi);
    }
    const auto line = std::format("{} of {} hit, {} witnessed; latest past a witness {:.3g}, widest before one {:.3g}",
                                  hits, kCases, witnessed, worst_late, widest);
    INFO(line);
    WARN(line);
    CHECK(witnessed * 10 >= hits * 8);
}

// Two edges crossing at right angles, one falling onto the other: they
// have no normal to be apart along, only the cross product of their
// directions (SurfaceCcdPolicy::kTangentCross). Edge a lies along x at
// z = 0; edge b along y at z = 1 falls at 2 and meets it at t = 0.5.
TEST_CASE("Two edges meet where they cross, and pass when they do not", "[physics][ccd][surface]") {
    const auto edge = [](V3 e0, V3 e1, V3 w) {
        MovingChart<double> c;
        c.p = [=](double u, double) { return V3{e0 + (e1 - e0) * u}; };
        c.w = [=](double, double) { return w; };
        c.bp = {V3{e1 - e0}.norm(), 0, 0, 0};
        c.v1 = 0;
        return c;
    };
    const auto a = edge(V3{-1, 0, 0}, V3{1, 0, 0}, V3{});
    const auto b = edge(V3{0, -1, 1}, V3{0, 1, 1}, V3{0, 0, -2});
    const auto r = first_contact(a, b, 1e-7);
    INFO(std::format("hit={} toi={} pairs={}", r.hit, r.toi, r.pairs));
    REQUIRE(r.hit);
    CHECK(r.toi <= 0.5);
    CHECK(r.toi > 0.5 - 1e-6);
    // Shifted past a's end, b falls beside it.
    const auto beside = edge(V3{1.5, -1, 1}, V3{1.5, 1, 1}, V3{0, 0, -2});
    CHECK_FALSE(first_contact(a, beside, 1e-7).hit);
}

// A vertex falling through a triangle, and one falling past its edge.
TEST_CASE("A vertex meets a triangle it falls through", "[physics][ccd][surface]") {
    MovingChart<double> tri;
    const V3 A{0, 0, 0}, B{1, 0, 0}, C{0, 1, 0};
    tri.p = [=](double u, double v) { return V3{A + (V3{B - A} * (1 - v) + V3{C - A} * v) * u}; };
    tri.w = [](double, double) { return V3{}; };
    tri.bp = {1.0, std::sqrt(2.0), 0, 0};
    const auto point = [](V3 x, V3 w) {
        MovingChart<double> c;
        c.p = [=](double, double) { return x; };
        c.w = [=](double, double) { return w; };
        c.u1 = c.v1 = 0;
        return c;
    };
    const auto r = first_contact(point(V3{0.25, 0.25, 1}, V3{0, 0, -4}), tri, 1e-7);
    REQUIRE(r.hit);
    CHECK(r.toi <= 0.25);
    CHECK(r.toi > 0.25 - 1e-6);
    CHECK_FALSE(first_contact(point(V3{0.75, 0.75, 1}, V3{0, 0, -4}), tri, 1e-7).hit);
}

// A ball as a point with a thickness: it touches a plane when its centre
// is a radius away -- a falling ball of radius 0.25 from height 1 at speed
// 2 meets z = 0 at t = (1 - 0.25) / 2.
TEST_CASE("A thickness is a distance of contact", "[physics][ccd][surface]") {
    MovingChart<double> ball;
    ball.p = [](double, double) { return V3{1.5, 1.5, 1.0}; };
    ball.w = [](double, double) { return V3{0.0, 0.0, -2.0}; };
    ball.u1 = ball.v1 = 0;
    ball.thickness = 0.25;
    const auto plane = moving(flat(0.0, 1.0), constant(V3{}));
    const auto r = first_contact(ball, plane, 1e-7);
    INFO(std::format("hit={} toi={} pairs={}", r.hit, r.toi, r.pairs));
    REQUIRE(r.hit);
    CHECK(r.toi <= 0.375);
    CHECK(r.toi > 0.375 - 1e-6);
}

// A point grazing a torus given by its formula, against the quartic's first
// root: never later, and to the width asked. Grazes are where a first-order
// bound spends its cells; the torus's second derivatives make them cheap.
TEST_CASE("A point grazing a torus chart: never later than the quartic", "[physics][ccd][surface][fuzz]") {
    const double R = 1.0, r = 0.3, tp = 2 * std::numbers::pi;
    MovingChart<double> torus;
    torus.p = [=](double u, double v) {
        return V3{(R + r * std::cos(v)) * std::cos(u), (R + r * std::cos(v)) * std::sin(u), r * std::sin(v)};
    };
    torus.w = [](double, double) { return V3{}; };
    torus.bp = {R + r, r, R + r, r};
    torus.u1 = torus.v1 = tp;
    std::mt19937_64 rng(20260927);
    std::uniform_real_distribution<double> ang(0, tp), uni(-1, 1);
    int n = 0;
    double worst_early = 0;
    for (int i = 0; i < 150; ++i) {
        const double u = ang(rng), v = ang(rng);
        const V3 nrm{std::cos(v) * std::cos(u), std::cos(v) * std::sin(u), std::sin(v)};
        const V3 foot{V3{std::cos(u), std::sin(u), 0.0} * R + nrm * r};
        V3 tan{uni(rng), uni(rng), uni(rng)};
        tan = V3{tan - nrm * tan.dot(nrm)};
        tan = V3{tan * (1 / tan.norm())};
        const double depth = std::pow(10.0, -2 - 6 * (uni(rng) + 1) / 2);
        const V3 p0{foot - nrm * depth - tan * 2.0}, disp{tan * 4.0};
        // First root of the quartic |p0 + s tan| on the torus, s in [0, 4].
        const auto F = [&](double s) {
            const V3 x{p0 + tan * s};
            const double q = std::hypot(x[0], x[1]) - R;
            return q * q + x[2] * x[2] - r * r;
        };
        double truth = 2;
        for (int k = 0; k < 4000 && truth > 1; ++k)
            if (F(4.0 * k / 4000) > 0 && F(4.0 * (k + 1) / 4000) <= 0) {
                double lo = 4.0 * k / 4000, hi = 4.0 * (k + 1) / 4000;
                for (int it = 0; it < 80; ++it) (F((lo + hi) / 2) > 0 ? lo : hi) = (lo + hi) / 2;
                truth = hi / 4;
            }
        if (F(0) <= 0 || truth > 1) continue;
        MovingChart<double> pt;
        pt.p = [p0](double, double) { return p0; };
        pt.w = [disp](double, double) { return disp; };
        pt.u1 = pt.v1 = 0;
        const auto c = first_contact(pt, torus, 1e-6);
        INFO(std::format("query {}: truth {:.12f} hit={} toi={:.12f}", i, truth, c.hit, c.toi));
        REQUIRE(c.hit);
        REQUIRE(c.toi <= truth + 1e-12);
        worst_early = std::max(worst_early, truth - c.toi);
        ++n;
    }
    INFO(std::format("{} grazes, worst early {:.3g}", n, worst_early));
    CHECK(n > 100);
    CHECK(worst_early < 1e-3);
}
