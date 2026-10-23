// Gravity of many sources on the closed 3-sphere, and flight under it.
//
// The potential of a unit mass at geodesic angle psi on the unit S^3 is
// G(psi) = -(pi - psi) cot(psi) / (4 pi^2), up to a constant, with a uniform
// background of opposite total mass. Two derivations that share no code with
// the formula check it: the zonal eigenfunction expansion (eigenvalues n(n+2),
// degeneracy (n+1)^2, zonal functions sin((n+1) psi) / sin(psi)), which also
// fixes the constant -- the mean of G over the sphere is -1/(8 pi^2) -- and the
// equation itself away from the source, G'' + 2 cot(psi) G' = -1/V, with the
// flux of grad G through a sphere around the source.
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <spatium/physics/mechanics/sphere_gravity.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <random>
#include <vector>

using namespace spatium;
using namespace spatium::physics::mechanics;
using Catch::Approx;

namespace {

constexpr double pi = std::numbers::pi;
constexpr double V3 = 2 * pi * pi;   // volume of the unit S^3

using V4 = Vec<double, 4>;

V4 unit4(std::mt19937& g) {
    std::normal_distribution<double> n;
    V4 v{n(g), n(g), n(g), n(g)};
    return V4{v / v.norm()};
}

// A unit tangent at x, random.
V4 tangent4(std::mt19937& g, const V4& x) {
    V4 t = unit4(g);
    t = V4{t - x * x.dot(t)};
    return V4{t / t.norm()};
}

V4 along(const V4& x, const V4& t, double s) {   // great circle from x with unit tangent t
    return V4{x * std::cos(s) + t * std::sin(s)};
}

// Fejer sum of the zonal expansion: -(1/V) sum_{n>=1} (n+1) / (n (n+2)) sin((n+1) psi) / sin(psi).
double zonal_green(double psi, int terms) {
    double s = 0;
    for (int n = 1; n < terms; ++n)
        s += (n + 1.0) / (n * (n + 2.0)) * std::sin((n + 1) * psi) * (1.0 - double(n) / terms);
    return -s / std::sin(psi) / V3;
}

double wrong_green(double psi) { return -std::cos(psi) / std::sin(psi) / (4 * pi * pi); }   // no (pi - psi)

double ode_residual(double (*G)(double), double psi, double h) {
    const double g0 = G(psi), gp = G(psi + h), gm = G(psi - h);
    const double d1 = (gp - gm) / (2 * h), d2 = (gp - 2 * g0 + gm) / (h * h);
    return d2 + 2 * std::cos(psi) / std::sin(psi) * d1 + 1.0 / V3;   // zero when the equation holds
}

double green(double psi) { return s3_green(psi); }

} // namespace

TEST_CASE("S3 Green's function: zonal eigenfunction expansion", "[gravity][s3]") {
    // The expansion has zero mean; the closed form differs by the constant 1/(8 pi^2).
    for (double psi : {0.3, 1.0, 2.0, 3.0, 3.1}) {
        const double series = zonal_green(psi, 400000);
        REQUIRE(series == Approx(s3_green(psi) + 1 / (8 * pi * pi)).margin(4e-6));
    }
    // And the mean over the sphere (weight (2/pi) sin^2) is -1/(8 pi^2): Simpson.
    const int n = 20000;
    double acc = 0;
    for (int i = 1; i < n; ++i) {
        const double psi = pi * i / n;
        acc += (i % 2 ? 4 : 2) * s3_green(psi) * 2 / pi * std::sin(psi) * std::sin(psi);
    }
    REQUIRE(acc * (pi / n) / 3 == Approx(-1 / (8 * pi * pi)).epsilon(1e-6));
}

TEST_CASE("S3 Green's function: Poisson equation and flux", "[gravity][s3]") {
    for (double psi : {0.5, 1.0, 1.5, 2.0, 2.5, 2.9})
        REQUIRE(std::abs(ode_residual(green, psi, 1e-3)) < 1e-6);
    // The verifier must fail on a wrong candidate (no (pi - psi) factor).
    REQUIRE(std::abs(ode_residual(wrong_green, 1.0, 1e-3)) > 1e-2);

    // Flux of grad G out of the geodesic sphere psi0: 1 minus the share of the volume inside
    // (the background); at a small sphere it is 1.
    for (double psi0 : {0.01, 0.5, 1.0, 2.5}) {
        const double h = 1e-6;
        const double dG = (s3_green(psi0 + h) - s3_green(psi0 - h)) / (2 * h);
        const double flux = dG * 4 * pi * std::sin(psi0) * std::sin(psi0);
        const double volume = 2 * pi * (psi0 - std::sin(psi0) * std::cos(psi0));
        REQUIRE(flux == Approx(1 - volume / V3).epsilon(1e-6));
    }
    // The slope-over-sine helper is that derivative divided by sin.
    for (double psi : {0.2, 1.3, 2.7, 3.1}) {
        const double h = 1e-6;
        const double dG = (s3_green(psi + h) - s3_green(psi - h)) / (2 * h);
        REQUIRE(s3_green_slope_over_sin(psi) == Approx(dG / std::sin(psi)).epsilon(1e-6));
    }
    REQUIRE(s3_green_slope_over_sin(pi) == Approx(1 / (6 * pi * pi)).epsilon(1e-12));
    REQUIRE(s3_green(pi) == Approx(1 / (4 * pi * pi)).epsilon(1e-12));
}

TEST_CASE("SphereGravity: force is minus the gradient of the potential", "[gravity][s3]") {
    std::mt19937 g(5);
    for (double eps : {0.0, 0.3}) {
        SphereGravity<double> field({1.7, eps, 1.0, 0.0});
        std::vector<SphereSource<double>> src;
        for (int i = 0; i < 4; ++i) src.push_back({unit4(g), 0.5 + 0.3 * i});
        field.set_sources(src);
        for (int trial = 0; trial < 20; ++trial) {
            const V4 x = unit4(g), t = tangent4(g, x);
            const double h = 1e-5;
            const double dphi = (field.potential(along(x, t, h)) - field.potential(along(x, t, -h))) / (2 * h);
            const V4 a = field.acceleration(x);
            REQUIRE(a.dot(x) == Approx(0).margin(1e-13));       // tangent
            REQUIRE(a.dot(t) == Approx(-dphi).epsilon(1e-6).margin(1e-8));
        }
    }
}

TEST_CASE("SphereGravity: antipode is a smooth zero, degenerate points are finite", "[gravity][s3]") {
    std::mt19937 g(9);
    for (double eps : {0.0, 0.2}) {
        SphereGravity<double> field({1.0, eps, 1.0, 0.0});
        const V4 x = unit4(g), t = tangent4(g, x);
        field.set_sources({{V4{-x}, 2.0}});
        REQUIRE(field.acceleration(x).norm() == Approx(0).margin(1e-14));
        // Next to the antipode the pull is linear in the offset: delta / (6 pi^2) per unit mass.
        for (double delta : {1e-2, 1e-4, 1e-6}) {
            const V4 s = along(V4{-x}, t, delta);
            field.set_sources({{s, 1.0}});
            const V4 a = field.acceleration(x);
            const double expected = delta / (6 * pi * pi) * (eps == 0 ? 1.0 : 1.0 / ((1 + eps * eps / 16) * (1 + eps * eps / 16)) * 1.0);
            REQUIRE(a.norm() == Approx(expected).epsilon(eps == 0 ? 1e-3 : 0.05));
            REQUIRE(a.dot(t) > 0);   // towards the source, which leans to +t from the antipode
        }
    }
    // A source on the receiving point: no direction, so no force; softened potential is finite.
    SphereGravity<double> hard({1.0, 0.0, 1.0, 0.0}), soft({1.0, 0.3, 1.0, 0.0});
    const V4 x = unit4(g);
    hard.set_sources({{x, 1.0}});
    soft.set_sources({{x, 1.0}});
    REQUIRE(hard.acceleration(x).norm() == 0.0);
    REQUIRE(soft.acceleration(x).norm() == Approx(0).margin(1e-14));
    REQUIRE(std::isfinite(soft.potential(x)));
    REQUIRE(std::isinf(hard.potential(x)));
    // Skipping the body's own source.
    REQUIRE(hard.potential(x, 0) == 0.0);
    // Empty set.
    SphereGravity<double> none;
    REQUIRE(none.acceleration(x).norm() == 0.0);
    REQUIRE(none.potential(x) == 0.0);
    REQUIRE(none.potential_energy() == 0.0);
}

TEST_CASE("SphereGravity: Plummer limit, softening to zero, clamp, coupling, radius", "[gravity][s3]") {
    std::mt19937 g(2);
    const V4 x = unit4(g), t = tangent4(g, x);
    // Small angles are flat: a = psi / (4 pi (psi^2 + eps^2)^(3/2)).
    {
        const double eps = 0.01, psi = 0.01;
        SphereGravity<double> f({1.0, eps, 1.0, 0.0});
        f.set_sources({{along(x, t, psi), 1.0}});
        const double plummer = psi / (4 * pi * std::pow(psi * psi + eps * eps, 1.5));
        REQUIRE(f.acceleration(x).norm() == Approx(plummer).epsilon(2e-3));
    }
    // eps -> 0 is the unsoftened law.
    {
        SphereGravity<double> a({1.0, 0.0, 1.0, 0.0}), b({1.0, 1e-9, 1.0, 0.0});
        const V4 s = unit4(g);
        a.set_sources({{s, 1.0}});
        b.set_sources({{s, 1.0}});
        REQUIRE(V4{a.acceleration(x) - b.acceleration(x)}.norm() < 1e-8);
    }
    // Clamp and coupling.
    {
        SphereGravity<double> f({1.0, 0.0, 1.0, 0.05});
        f.set_sources({{along(x, t, 0.01), 1.0}});
        REQUIRE(f.acceleration(x).norm() == Approx(0.05).epsilon(1e-12));
        REQUIRE(f.acceleration(x, SphereGravity<double>::npos, 0.5).norm() == Approx(0.025).epsilon(1e-12));
    }
    // Radius R: potential G/R, acceleration a/R^2, points of norm R.
    {
        const double R = 2.5;
        SphereGravity<double> f1({1.0, 0.0, 1.0, 0.0}), fR({1.0, 0.0, R, 0.0});
        const V4 s = unit4(g);
        f1.set_sources({{s, 1.0}});
        fR.set_sources({{V4{s * R}, 1.0}});
        REQUIRE(fR.potential(V4{x * R}) == Approx(f1.potential(x) / R).epsilon(1e-13));
        REQUIRE(V4{fR.acceleration(V4{x * R}) * (R * R) - f1.acceleration(x)}.norm() < 1e-12);
        // and it is the gradient along arc length R psi
        const double h = 1e-5;
        const double dphi = (fR.potential(V4{along(x, t, h / R) * R}) - fR.potential(V4{along(x, t, -h / R) * R})) / (2 * h);
        REQUIRE(fR.acceleration(V4{x * R}).dot(t) == Approx(-dphi).epsilon(1e-6));
    }
}

TEST_CASE("S3 flight: no force is the exact great circle", "[gravity][s3][flight]") {
    Sphere<3, double> S;
    std::mt19937 g(4);
    const V4 p = unit4(g);
    const V4 v = V4{tangent4(g, p) * 0.8};
    SphereBody3 free_body{1.0, p, v}, forced{1.0, p, v};
    const double dt = 0.1;
    const int steps = 2000;
    for (int i = 0; i < steps; ++i) {
        flight_step(free_body, S, dt);
        flight_step(forced, S, dt, [](const V4&) { return V4{}; });
    }
    const double th = 0.8 * dt * steps;
    const V4 exact{p * std::cos(th) + V4{v / 0.8} * std::sin(th)};
    REQUIRE(V4{free_body.position - exact}.norm() < 1e-12);
    REQUIRE(V4{forced.position - exact}.norm() < 1e-12);
    REQUIRE(manifold_constraint_residual(forced, S) < 1e-12);
    REQUIRE(speed_on_sphere(forced) == Approx(0.8).epsilon(1e-12));
}

namespace {

std::vector<SphereBody3> small_system(unsigned seed, int n) {
    std::mt19937 g(seed);
    std::vector<SphereBody3> b;
    for (int i = 0; i < n; ++i) {
        const V4 x = unit4(g);
        b.push_back({0.2 + 0.1 * i, x, V4{tangent4(g, x) * 0.5}});
    }
    return b;
}

double run_energy_drift(double dt, double span, double eps) {
    Sphere<3, double> S;
    auto bodies = small_system(11, 6);
    SphereGravity<double> field({1.0, eps, 1.0, 0.0});
    const double e0 = sphere_nbody_energy(bodies, field);
    double worst = 0;
    const int steps = int(span / dt);
    for (int i = 0; i < steps; ++i) {
        sphere_nbody_step(bodies, field, S, dt);
        worst = std::max(worst, std::abs(sphere_nbody_energy(bodies, field) - e0));
    }
    return worst;
}

} // namespace

TEST_CASE("S3 n-body: energy drift, second order, angular momentum exact", "[gravity][s3][flight]") {
    Sphere<3, double> S;
    SphereGravity<double> field({1.0, 0.15, 1.0, 0.0});
    auto bodies = small_system(11, 6);
    const double e0 = sphere_nbody_energy(bodies, field);
    const auto L0 = sphere_angular_momentum(bodies);
    double worst_e = 0, worst_l = 0, worst_c = 0;
    for (int i = 0; i < 5000; ++i) {
        sphere_nbody_step(bodies, field, S, 1e-3);
        worst_e = std::max(worst_e, std::abs(sphere_nbody_energy(bodies, field) - e0));
        const auto L = sphere_angular_momentum(bodies);
        for (int k = 0; k < 6; ++k) worst_l = std::max(worst_l, std::abs(L[k] - L0[k]));
        for (const auto& b : bodies) worst_c = std::max(worst_c, manifold_constraint_residual(b, S));
    }
    INFO("energy drift " << worst_e << ", angular momentum drift " << worst_l);
    REQUIRE(worst_e < 1e-6);
    REQUIRE(worst_l < 1e-12);
    REQUIRE(worst_c < 1e-12);

    // Halving the step cuts the energy error by about four.
    const double e1 = run_energy_drift(4e-3, 2.0, 0.15), e2 = run_energy_drift(2e-3, 2.0, 0.15);
    INFO("e(dt) " << e1 << ", e(dt/2) " << e2);
    REQUIRE(e1 / e2 > 3.0);
    REQUIRE(e1 / e2 < 5.0);
}

TEST_CASE("S3 n-body: a clamp breaks the angular momentum, and the check sees it", "[gravity][s3][flight]") {
    Sphere<3, double> S;
    SphereGravity<double> field({1.0, 0.15, 1.0, 0.05});   // max acceleration 0.05, below the pull
    auto bodies = small_system(11, 6);
    const auto L0 = sphere_angular_momentum(bodies);
    for (int i = 0; i < 2000; ++i) sphere_nbody_step(bodies, field, S, 2e-3);
    const auto L = sphere_angular_momentum(bodies);
    double d = 0;
    for (int k = 0; k < 6; ++k) d = std::max(d, std::abs(L[k] - L0[k]));
    REQUIRE(d > 1e-5);
}

TEST_CASE("S3 flight: reversible under a conservative force", "[gravity][s3][flight]") {
    Sphere<3, double> S;
    SphereGravity<double> field({1.0, 0.15, 1.0, 0.0});
    auto bodies = small_system(3, 5);
    const auto start = bodies;
    const int n = 2000;
    for (int i = 0; i < n; ++i) sphere_nbody_step(bodies, field, S, 1e-3);
    for (auto& b : bodies) b.velocity = V4{-b.velocity};
    for (int i = 0; i < n; ++i) sphere_nbody_step(bodies, field, S, 1e-3);
    double err = 0;
    for (std::size_t i = 0; i < bodies.size(); ++i) {
        err = std::max(err, V4{bodies[i].position - start[i].position}.norm());
        err = std::max(err, V4{bodies[i].velocity + start[i].velocity}.norm());
    }
    INFO("round trip error " << err);
    REQUIRE(err < 1e-9);
}

namespace {

// Spread of the arrivals at -P after an angle pi, for a fan of shots from P with equal speed
// through fixed sources of the given total mass.
double fan_spread(double total_mass, int shots, double eps = 0.1) {
    Sphere<3, double> S;
    std::mt19937 g(77);
    const V4 P{0.0, 0.0, 0.0, 1.0};
    SphereGravity<double> field({1.0, eps, 1.0, 0.0});
    std::vector<SphereSource<double>> src;
    if (total_mass > 0) {
        const std::array<V4, 3> where{V4{1, 0, 0, 0}, V4{0, 0.6, 0.8, 0}, V4{0.5, 0.5, 0.5, 0.5}};
        for (auto w : where) src.push_back({V4{w / w.norm()}, total_mass / 3});
    }
    field.set_sources(src);
    const double speed = 1.0;
    const int steps = 2000;
    const double dt = pi / speed / steps;
    double worst = 0;
    for (int k = 0; k < shots; ++k) {
        V4 t{};
        {   // a direction tangent at P: the first three coordinates
            std::normal_distribution<double> n;
            t = V4{n(g), n(g), n(g), 0.0};
            t = V4{t / t.norm() * speed};
        }
        SphereBody3 b{1.0, P, t};
        for (int i = 0; i < steps; ++i)
            flight_step(b, S, dt, [&](const V4& x) { return field.acceleration(x); });
        worst = std::max(worst, V4{b.position + P}.norm());
    }
    return worst;
}

} // namespace

TEST_CASE("S3 flight: antipodal focusing", "[gravity][s3][flight]") {
    // No masses: every geodesic from P reaches -P after an angle pi.
    REQUIRE(fan_spread(0.0, 32) < 1e-9);
    // With masses the fan is spread, and more as the masses grow.
    const double s1 = fan_spread(1e-3, 32), s2 = fan_spread(1e-2, 32), s3 = fan_spread(1e-1, 32);
    INFO("spread " << s1 << " " << s2 << " " << s3);
    REQUIRE(s1 > 1e-9);
    REQUIRE(s2 > s1);
    REQUIRE(s3 > s2);
}
