// Free particle on Sphere<2> follows a great-circle geodesic indefinitely
// without drifting off the manifold.
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <spatium/physics/mechanics/manifold_body.hpp>
#include <cmath>
#include <numbers>

using namespace spatium;
using namespace spatium::physics::mechanics;
using Catch::Approx;

TEST_CASE("Sphere<2> geodesic: equator closes in 2π for unit speed", "[manifold][sphere]") {
    Sphere<2, double> S;          // unit sphere
    PointOnManifold<Sphere<2, double>> body{
        /*mass*/ 1.0,
        /*position*/ Vec<double, 3>{1.0, 0.0, 0.0},      // equator start
        /*velocity*/ Vec<double, 3>{0.0, 1.0, 0.0},      // unit tangent along equator
    };

    REQUIRE(manifold_constraint_residual(body, S) < 1e-14);

    // Total arc = 2π at speed 1 ⇒ full equator, period = 2π.
    constexpr double T_period = 2.0 * std::numbers::pi;
    constexpr double dt = T_period / 1000.0;
    for (int i = 0; i < 1000; ++i)
        geodesic_step(body, S, dt);

    // Back to start to high precision (closed-form exp_map ⇒ round-off only).
    REQUIRE(std::abs(body.position[0] - 1.0) < 1e-12);
    REQUIRE(std::abs(body.position[1])      < 1e-12);
    REQUIRE(std::abs(body.velocity[1] - 1.0) < 1e-12);
    // Constraints still satisfied (|pos| = 1, pos · vel = 0).
    REQUIRE(manifold_constraint_residual(body, S) < 1e-12);
}

TEST_CASE("Sphere<2> geodesic: speed preserved over many periods", "[manifold][sphere]") {
    Sphere<2, double> S;
    PointOnManifold<Sphere<2, double>> body{
        1.0,
        Vec<double, 3>{0.0, 0.0, 1.0},                 // north pole
        Vec<double, 3>{0.7, 0.2, 0.0},                 // non-unit tangent
    };
    double v0 = speed_on_sphere(body);

    constexpr double dt = 0.01;
    for (int i = 0; i < 10000; ++i)                    // 100 time-units of motion
        geodesic_step(body, S, dt);

    double v1 = speed_on_sphere(body);
    // Free-particle geodesic motion preserves kinetic energy exactly.
    REQUIRE(std::abs(v1 - v0) / v0 < 1e-12);
    REQUIRE(manifold_constraint_residual(body, S) < 1e-12);
}

TEST_CASE("Sphere<2> geodesic: radius != 1 scales correctly", "[manifold][sphere]") {
    // Sphere of radius 2. Equator circumference is 4π, period at speed 1 is 4π.
    Sphere<2, double> S;
    S.radius = 2.0;
    PointOnManifold<Sphere<2, double>> body{
        1.0,
        Vec<double, 3>{2.0, 0.0, 0.0},                 // on sphere of radius 2
        Vec<double, 3>{0.0, 1.0, 0.0},
    };

    REQUIRE(manifold_constraint_residual(body, S) < 1e-14);

    constexpr double period = 4.0 * std::numbers::pi;
    constexpr double dt = period / 1000.0;
    for (int i = 0; i < 1000; ++i)
        geodesic_step(body, S, dt);

    REQUIRE(std::abs(body.position[0] - 2.0) < 1e-10);
    REQUIRE(std::abs(body.position[1])       < 1e-10);
}
