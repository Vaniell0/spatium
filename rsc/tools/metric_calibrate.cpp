// metric_calibrate -- a constant of the binary's metric found by the vacuum
// residual alone.
//
// The kind of RSC the catalogue calls calibration, on the one oracle the
// merger scene has for free: the exact Ricci residual (vacuum_residual,
// through nested Dual). A superposition of two Kerr-Schild terms is not a
// solution, and nothing tells the scene how to make it closer -- except
// the residual. So the scene's metric gets a parameter, the share of each
// hole's velocity its term is boosted by (SpacetimeScene::boost), and the
// parameter is chosen by minimising the residual over events near both
// holes and between them.
//
// It is a test of the method as much as of the metric, because here the
// answer is known: physics says 1, each term boosted into the frame its
// hole moves in. If minimising the residual finds 1, calibration by this
// oracle works, and the next parameters -- ones with no known answer --
// can be trusted to it.
//
// The residual is taken relative to the local tidal scale, sum m_i / r_i^3,
// so an event by a hole and one far from both weigh alike. Events where the
// superposition is not Lorentzian are left out and counted.
//
// Build: part of the rsc tools target.

#include <spatium/physics/relativity/metric_field.hpp>
#include <spatium/physics/relativity/spacetime_scene.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <print>
#include <random>
#include <vector>

namespace rel = spatium::physics::relativity;
using V4 = spatium::Vec<double, 4>;

namespace {

struct Setup {
    double m1, m2, d;
};

struct Events {
    std::vector<V4> near, middle;
};

// Events a few M from each hole at random times of one orbit, and events
// between the holes.
Events sample(const Setup& s, std::uint64_t seed, int n) {
    rel::SpacetimeScene<double> sc;
    sc.binary(s.m1, s.m2, s.d, false);
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> uni(-1, 1);
    const double M = s.m1 + s.m2, period = 2 * std::numbers::pi * std::sqrt(s.d * s.d * s.d / M);
    Events e;
    for (int i = 0; i < n; ++i) {
        const double t = (uni(rng) + 1) / 2 * period;
        const double phase = std::sqrt(M / (s.d * s.d * s.d)) * t;
        for (int k = 0; k < 2; ++k) {
            const double sgn = k == 0 ? s.m2 / M : -s.m1 / M, m = k == 0 ? s.m1 : s.m2;
            const double cx = sgn * s.d * std::cos(phase), cy = sgn * s.d * std::sin(phase);
            double dx = uni(rng), dy = uni(rng), dz = uni(rng);
            const double len = std::sqrt(dx * dx + dy * dy + dz * dz) + 1e-12;
            const double r = m * (3 + 5 * (uni(rng) + 1) / 2);
            e.near.push_back(V4{t, cx + dx / len * r, cy + dy / len * r, dz / len * r});
        }
        // Between: along the line joining them, off it by a little.
        const double f = 0.3 + 0.4 * (uni(rng) + 1) / 2;
        const double x1 = s.m2 / M * s.d, x2 = -s.m1 / M * s.d, xm = x1 + f * (x2 - x1);
        e.middle.push_back(V4{t, xm * std::cos(phase), xm * std::sin(phase), 0.5 * uni(rng)});
    }
    return e;
}

double tidal(const Setup& s, const V4& x) {
    const double M = s.m1 + s.m2, phase = std::sqrt(M / (s.d * s.d * s.d)) * x[0];
    double k = 0;
    for (int h = 0; h < 2; ++h) {
        const double sgn = h == 0 ? s.m2 / M : -s.m1 / M, m = h == 0 ? s.m1 : s.m2;
        const double dx = x[1] - sgn * s.d * std::cos(phase), dy = x[2] - sgn * s.d * std::sin(phase), dz = x[3];
        const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
        k += m / (r * r * r);
    }
    return k;
}

struct Score {
    double near = 0, middle = 0;   // median relative residual
    int refused = 0;
};

Score score(const Setup& s, const Events& e, double boost) {
    rel::SpacetimeScene<double> sc;
    sc.binary(s.m1, s.m2, s.d, false);
    sc.boost(boost);
    const auto g = sc.metric();
    Score out;
    const auto median_of = [&](const std::vector<V4>& xs) {
        std::vector<double> v;
        for (const auto& x : xs) {
            if (!rel::lorentzian_at(*g, x)) { ++out.refused; continue; }
            v.push_back(rel::vacuum_residual(*g, x) / tidal(s, x));
        }
        if (v.empty()) return std::nan("");
        std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
        return v[v.size() / 2];
    };
    out.near = median_of(e.near);
    out.middle = median_of(e.middle);
    return out;
}

}  // namespace

int main(int argc, char** argv) {
    const int n = argc > 1 ? std::atoi(argv[1]) : 60;
    const std::vector<Setup> setups{{0.5, 0.5, 10}, {0.5, 0.5, 14}, {0.5, 0.5, 20}, {0.5, 0.5, 40},
                                    {0.25, 0.75, 14}, {0.25, 0.75, 20}};
    std::println("{} events near each hole and {} between, per setup; residual over the tidal scale, median", 2 * n, n);
    std::println("{:<14} | {:>10} {:>10} | {:>7} {:>10} {:>10} | {:>10} {:>10} | {:>7}", "masses, d", "near b=0",
                 "mid b=0", "b*", "near b*", "mid b*", "near b=1", "mid b=1", "refused");
    for (const auto& s : setups) {
        const auto e = sample(s, 42, n);
        const auto obj = [&](double b) { return score(s, e, b).near; };
        // Golden section on [0, 1.6] for the near-hole residual, which the
        // boost is for; the middle is reported, not optimised.
        double lo = 0, hi = 1.6;
        const double gr = (std::sqrt(5.0) - 1) / 2;
        double x1 = hi - gr * (hi - lo), x2 = lo + gr * (hi - lo), f1 = obj(x1), f2 = obj(x2);
        for (int it = 0; it < 30; ++it) {
            if (f1 < f2) { hi = x2; x2 = x1; f2 = f1; x1 = hi - gr * (hi - lo); f1 = obj(x1); }
            else { lo = x1; x1 = x2; f1 = f2; x2 = lo + gr * (hi - lo); f2 = obj(x2); }
        }
        const double best = (lo + hi) / 2;
        const auto s0 = score(s, e, 0), sb = score(s, e, best), s1 = score(s, e, 1);
        std::println("{:.2f}:{:.2f}, {:>4.0f} | {:>10.3g} {:>10.3g} | {:>7.4f} {:>10.3g} {:>10.3g} | {:>10.3g} {:>10.3g} | {:>7}",
                     s.m1, s.m2, s.d, s0.near, s0.middle, best, sb.near, sb.middle, s1.near, s1.middle, s0.refused);
    }
}
