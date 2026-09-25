// A black-hole scene described by what is in it. One hole at rest is the
// Kerr-Schild metric; a binary's path is the quadrupole inspiral written
// as fields, held to integrating Peters' equations; and the superposition
// reports the residual it has rather than hiding it.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <spatium/physics/relativity/spacetime_scene.hpp>

#include <cmath>
#include <format>

using namespace spatium;
using namespace spatium::physics::relativity;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using V4 = Vec<double, 4>;

namespace {
double at(const io::build::Field<double>& f, double t) { return f(t, 0.0); }
}  // namespace

TEST_CASE("One hole at rest is the Kerr-Schild metric", "[relativity][scene]") {
    SpacetimeScene<double> scene;
    scene.hole(1.0, 0.6);
    const auto g = scene.metric();
    REQUIRE(g);
    const auto ks = kerr_schild(1.0, 0.6);
    for (const V4 x : {V4{0.0, 5.0, 1.0, 2.0}, V4{3.0, -7.0, 4.0, -1.0}, V4{9.0, 0.5, 0.2, 12.0}}) {
        const auto a = (*g)(x), b = ks(x);
        for (std::size_t i = 0; i < 4; ++i)
            for (std::size_t j = 0; j < 4; ++j) CHECK_THAT(a(i, j), WithinAbs(b(i, j), 1e-13));
    }
}

TEST_CASE("A binary's inspiral is Peters' equations, written as fields", "[relativity][scene]") {
    const double m1 = 1.0, m2 = 0.6, a0 = 30.0, stop = 10.0;
    const double M = m1 + m2, beta = 64.0 / 5.0 * m1 * m2 * M;
    SpacetimeScene<double> scene;
    scene.binary(m1, m2, a0, /*inspiral=*/true, stop);
    const auto& h = scene.holes();
    REQUIRE(h.size() == 2);

    // Integrate da/dt = -beta/a^3, dphi/dt = sqrt(M/a^3) by RK4, stopping
    // the shrink at `stop`, and compare the holes' positions.
    const double t_stop = (std::pow(a0, 4) - std::pow(stop, 4)) / (4 * beta);
    auto rhs = [&](double a) { return std::pair{-beta / (a * a * a), std::sqrt(M / (a * a * a))}; };
    double a = a0, phi = 0.0, t = 0.0;
    const double dt = 0.5;
    for (double check : {100.0, 2000.0, 0.5 * t_stop, 0.9 * t_stop, t_stop + 400.0}) {
        while (t + 1e-12 < check) {
            const double step = std::min(dt, check - t);
            if (t >= t_stop) {   // circling at stop
                phi += std::sqrt(M / (stop * stop * stop)) * step;
                t += step;
                continue;
            }
            const double s = std::min(step, t_stop - t);
            const auto [k1a, k1p] = rhs(a);
            const auto [k2a, k2p] = rhs(a + 0.5 * s * k1a);
            const auto [k3a, k3p] = rhs(a + 0.5 * s * k2a);
            const auto [k4a, k4p] = rhs(a + s * k3a);
            a += s / 6 * (k1a + 2 * k2a + 2 * k3a + k4a);
            phi += s / 6 * (k1p + 2 * k2p + 2 * k3p + k4p);
            t += s;
        }
        INFO(std::format("t = {}", check));
        const double x1 = at(h[0].x, check), y1 = at(h[0].y, check);
        const double x2 = at(h[1].x, check), y2 = at(h[1].y, check);
        CHECK_THAT(std::hypot(x1 - x2, y1 - y2), WithinRel(a, 1e-8));
        CHECK_THAT(std::atan2(y1, x1), WithinAbs(std::remainder(phi, 2 * std::numbers::pi), 1e-6));
        // The centre of mass stays at the origin.
        CHECK_THAT(m1 * x1 + m2 * x2, WithinAbs(0.0, 1e-9));
        CHECK_THAT(m1 * y1 + m2 * y2, WithinAbs(0.0, 1e-9));
    }
}

TEST_CASE("A superposed binary reports the residual it has", "[relativity][scene]") {
    SpacetimeScene<double> one;
    one.hole(1.0, 0.0);
    const V4 beside{0.0, 20.0, 0.0, 0.0};
    const auto r1 = one.residual_at(beside);
    REQUIRE(r1);
    CHECK(*r1 < 1e-9);

    SpacetimeScene<double> pair;
    pair.binary(1.0, 1.0, 40.0, /*inspiral=*/false);
    const auto r2 = pair.residual_at(V4{0.0, 0.0, 0.0, 0.0});   // the midpoint at t = 0
    REQUIRE(r2);
    INFO(std::format("residual at the midpoint of a pair 40 M apart: {:.2e}", *r2));
    CHECK(*r2 > 1e-5);   // not a vacuum solution, and it says so
}
