// Gravity and flight on the closed 3-sphere.
//
// The consumer of physics/mechanics/sphere_gravity.hpp, of `flight_step` in
// manifold_body.hpp and of spaces/sphere_charts.hpp. Headless; three parts:
//
//   1. N softened bodies pulling each other on S^3 (the Green's function with a
//      uniform neutralising background), advanced by kick - great circle - kick.
//      It prints the drift of the energy and of the six components of angular
//      momentum, the generators of SO(4) (the system is carried to itself by
//      every rotation of R^4, so they are conserved; the sum of m v is not, as
//      the sphere has no translations), and the milliseconds per step.
//   2. Antipodal focusing: a fan of shots of equal speed from P. With no masses
//      every geodesic arrives at -P after an angle pi, whatever the direction;
//      with fixed masses the arrivals spread, more as the mass grows. A table
//      of the spread.
//   3. An SVG of the two stereographic discs, north | south: the fan crosses
//      the equator (the rim, the same coordinates in both discs) and meets at
//      the antipode, without masses (blue) and with them (orange).
//
//   sphere_nbody_demo [--n 8] [--steps 20000] [--dt 0.001] [--eps 0.15]
//                     [--out sphere_charts.svg] [--force]
//
// Units: radius 1, G = 1.

#include "io_helpers.hpp"
#include <spatium/io/svg.hpp>
#include <spatium/physics/mechanics/sphere_gravity.hpp>
#include <spatium/spaces/sphere_charts.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <numbers>
#include <random>
#include <string>
#include <vector>

using namespace spatium;
using namespace spatium::physics::mechanics;

namespace {

using V4 = Vec<double, 4>;
using Clock = std::chrono::steady_clock;

struct Options {
    int n = 8;
    int steps = 20000;
    double dt = 1e-3;
    double eps = 0.15;
    std::string out = "sphere_charts.svg";
    bool force = false;
};

Options parse(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        const auto is = [&](const char* k) { return std::strcmp(argv[i], k) == 0 && i + 1 < argc; };
        if (is("--n")) o.n = std::atoi(argv[++i]);
        else if (is("--steps")) o.steps = std::atoi(argv[++i]);
        else if (is("--dt")) o.dt = std::atof(argv[++i]);
        else if (is("--eps")) o.eps = std::atof(argv[++i]);
        else if (is("--out")) o.out = argv[++i];
        else if (std::strcmp(argv[i], "--force") == 0) o.force = true;
        else { std::fprintf(stderr, "unknown or incomplete option %s\n", argv[i]); std::exit(2); }
    }
    return o;
}

std::vector<SphereBody3> make_bodies(int n, std::mt19937& rng) {
    std::normal_distribution<double> g;
    std::vector<SphereBody3> b;
    for (int i = 0; i < n; ++i) {
        V4 x{g(rng), g(rng), g(rng), g(rng)};
        x = V4{x / x.norm()};
        V4 t{g(rng), g(rng), g(rng), g(rng)};
        t = V4{t - x * x.dot(t)};
        b.push_back({1.0 / n, x, V4{t * (0.4 / t.norm())}});
    }
    return b;
}

// The flight of one shot from P with tangent t, for an angle pi, through fixed masses.
// Returns the path (every `every` steps) and the arrival.
struct Shot { std::vector<V4> path; V4 arrival; };

Shot fly(const SphereGravity<double>& field, const V4& P, const V4& t, int steps, int every) {
    const Sphere<3, double> S;
    SphereBody3 b{1.0, P, t};
    Shot s;
    const double dt = std::numbers::pi / t.norm() / steps;
    s.path.push_back(b.position);
    for (int i = 1; i <= steps; ++i) {
        flight_step(b, S, dt, [&](const V4& x) { return field.acceleration(x); });
        if (i % every == 0) s.path.push_back(b.position);
    }
    s.arrival = b.position;
    return s;
}

SphereGravity<double> fixed_masses(double total, double eps) {
    SphereGravity<double> f({1.0, eps, 1.0, 0.0});
    std::vector<SphereSource<double>> src;
    if (total > 0) {
        const std::array<V4, 3> where{V4{1, 0, 0, 0}, V4{0, 0.6, 0.8, 0}, V4{0.5, 0.5, 0.5, 0.5}};
        for (const auto& w : where) src.push_back({V4{w / w.norm()}, total / 3});
    }
    f.set_sources(std::move(src));
    return f;
}

// The fan lies in the great 2-sphere spanned by e0, e1, e3, so its charts are planar pictures.
std::vector<V4> fan(const V4& P, int shots) {
    const V4 e0{P[3], 0, 0, -P[0]};   // tangent at P = (sin a, 0, 0, cos a): e0 = (cos a, 0, 0, -sin a)
    const V4 e1{0, 1, 0, 0};
    std::vector<V4> out;
    for (int k = 0; k < shots; ++k) {
        const double phi = 2 * std::numbers::pi * k / shots;
        out.push_back(V4{e0 * std::cos(phi) + e1 * std::sin(phi)});
    }
    return out;
}

} // namespace

int main(int argc, char** argv) {
    const Options o = parse(argc, argv);
    const Sphere<3, double> S;
    constexpr double a0 = 0.9;
    const V4 P{std::sin(a0), 0.0, 0.0, std::cos(a0)};

    // 1. N bodies pulling each other.
    {
        std::mt19937 rng(1);
        auto bodies = make_bodies(o.n, rng);
        SphereGravity<double> field(SphereGravityLaw<double>{1.0, o.eps, 1.0, 0.0});
        const double e0 = sphere_nbody_energy(bodies, field);
        const auto L0 = sphere_angular_momentum(bodies);
        double drift_e = 0, drift_l = 0, ms = 0;
        const int every = std::max(1, o.steps / 50);
        for (int i = 1; i <= o.steps; ++i) {
            const auto t0 = Clock::now();
            sphere_nbody_step(bodies, field, S, o.dt);
            ms += std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
            if (i % every == 0) {
                drift_e = std::max(drift_e, std::abs(sphere_nbody_energy(bodies, field) - e0));
                const auto L = sphere_angular_momentum(bodies);
                for (int k = 0; k < 6; ++k) drift_l = std::max(drift_l, std::abs(L[k] - L0[k]));
            }
        }
        double constraint = 0;
        for (const auto& b : bodies) constraint = std::max(constraint, manifold_constraint_residual(b, S));
        double p_norm = 0;   // the sum of m v is not conserved: no translations on a sphere
        {
            V4 p{};
            for (const auto& b : bodies) p = V4{p + b.velocity * b.mass};
            p_norm = p.norm();
        }
        std::printf("n-body on S^3: %d bodies, %d steps of %g, softening %g\n", o.n, o.steps, o.dt, o.eps);
        std::printf("  energy %.9f, max drift %.3e (relative %.2e)\n", e0, drift_e, drift_e / std::abs(e0));
        std::printf("  angular momentum (6 components of SO(4)), max drift %.3e\n", drift_l);
        std::printf("  |sum m v| at the end %.3e (not an invariant)\n", p_norm);
        std::printf("  constraint residual %.2e, %.4f ms per step\n", constraint, ms / o.steps);
    }

    // 2. Antipodal focusing.
    const int shots = 32;
    const auto dirs = fan(P, shots);
    std::printf("\nantipodal focusing: %d shots from P, speed 1, angle pi, 2000 steps\n", shots);
    std::printf("  %-14s %-14s\n", "total mass", "max |arrival + P|");
    for (double mass : {0.0, 1e-3, 1e-2, 1e-1, 1.0}) {
        const auto field = fixed_masses(mass, 0.1);
        double worst = 0;
        for (const auto& t : dirs) worst = std::max(worst, V4{fly(field, P, t, 2000, 2000).arrival + P}.norm());
        std::printf("  %-14g %-14.3e\n", mass, worst);
    }

    // 3. The two discs.
    SphereCharts<3> C{1.0};
    const double d = 1.15;   // centres of the discs, in units of the rim radius
    io::Svg svg(1000, 500, 215);
    svg.circle(-d, 0, 1.0, "#89b4fa", 1.5).circle(d, 0, 1.0, "#89b4fa", 1.5);
    svg.text(-d - 0.35, 1.05, "north chart", "#cdd6f4", 13).text(d - 0.35, 1.05, "south chart", "#cdd6f4", 13);
    const auto draw = [&](const SphereGravity<double>& field, const char* colour) {
        for (const auto& t : dirs) {
            const auto shot = fly(field, P, t, 2000, 10);
            for (std::size_t i = 1; i < shot.path.size(); ++i) {
                const auto &y0 = shot.path[i - 1], &y1 = shot.path[i];
                if (C.chart_of(y0) != C.chart_of(y1)) continue;
                const bool north = C.chart_of(y0) == SphereCharts<3>::Chart::North;
                const auto x0 = C.to_chart(y0, north ? SphereCharts<3>::Chart::North : SphereCharts<3>::Chart::South);
                const auto x1 = C.to_chart(y1, north ? SphereCharts<3>::Chart::North : SphereCharts<3>::Chart::South);
                const double cx = north ? -d : d;
                svg.line(cx + x0[0], x0[1], cx + x1[0], x1[1], colour, 0.8);
            }
        }
    };
    draw(fixed_masses(0.0, 0.1), "#89b4fa");
    draw(fixed_masses(1.0, 0.1), "#fab387");
    const auto xp = C.to_north(P);
    const auto xa = C.to_south(V4{-P});
    svg.point(-d + xp[0], xp[1], 4, "#a6e3a1").text(-d + xp[0], xp[1] + 0.05, "P", "#a6e3a1", 12);
    svg.point(d + xa[0], xa[1], 4, "#f38ba8").text(d + xa[0], xa[1] + 0.05, "-P", "#f38ba8", 12);
    svg.text(-0.95, -1.12, "blue: no masses, orange: total mass 1; the rim is the equator, the same coordinates in both discs",
             "#a6adc8", 11);
    if (examples::confirm_overwrite(o.out, o.force) && svg.save(o.out)) std::printf("\nwrote %s\n", o.out.c_str());
    return 0;
}
