#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <spatium/physics/mechanics/contact_force.hpp>
#include <spatium/physics/mechanics/integrator.hpp>
#include <spatium/spaces/chart.hpp>
#include <spatium/spaces/sphere.hpp>

#include <cmath>
#include <vector>

using namespace spatium;
using namespace spatium::physics::mechanics;
using Catch::Matchers::WithinAbs;

// The barrier, the contact query and the symplectic integrator all
// existed; nothing joined them. These check the join rather than any of
// the three, which is why they are about invariants and not trajectories:
// a contact force that merely looks plausible is the failure this is
// guarding against.

TEST_CASE("A barrier force keeps a body out of a surface", "[physics][contact]") {
    // Falling straight at a sphere from outside it.
    const auto ball = chart_of(Sphere<2, double>{1.0});
    PointMass<3, double> body{1.0, {0.0, 0.0, 2.5}, {0.0, 0.0, -3.0}};

    auto gravity = [](const PointMass<3, double>&, double) {
        return Vec<double, 3>{0.0, 0.0, -1.0};
    };
    auto contact = surface_contact_force<double>(ball, /*d_hat=*/0.2, /*kappa=*/40.0);
    auto force = sum_forces(gravity, contact);

    double deepest = 1e9;
    for (int i = 0; i < 4000; ++i) {
        verlet_step(body, force, 1e-3, i * 1e-3);
        deepest = std::min(deepest, body.state.position.norm() - 1.0);
    }

    // Never inside. The barrier diverges at zero distance, so this is a
    // property of the potential rather than of a correction applied after
    // the fact.
    CHECK(deepest > 0.0);

    // And it actually got close -- a test that passes because the body
    // stopped a mile away would be testing nothing.
    CHECK(deepest < 0.2);
    INFO("closest approach: " << deepest);
}

TEST_CASE("Several surfaces sum rather than take turns", "[physics][contact]") {
    // Two coincident spheres, which is not a corner and is not pretending
    // to be one -- `Sphere<2, T>` carries a radius and no centre, so two
    // of them cannot be placed apart. What it does test is the property
    // that matters: a barrier is a potential and potentials add, so two
    // surfaces at the same place give exactly twice the force of one, with
    // no ordering between them and nothing to arbitrate.
    //
    // That is what keeps this symplectic. A sequence of pairwise position
    // corrections, which is the obvious alternative, does not sum and does
    // not commute.
    const auto ball = chart_of(Sphere<2, double>{1.0});
    std::vector<ParametricSurface<double>> both{ball, ball};

    PointMass<3, double> probe{1.0, {0.0, 0.0, 1.05}, {}};
    const auto one = surface_contact_force<double>(ball, 0.2, 40.0)(probe, 0.0);
    const auto two = surfaces_contact_force<double>(both, 0.2, 40.0)(probe, 0.0);

    CHECK_THAT(two.norm(), WithinAbs(2.0 * one.norm(), 1e-12));
    CHECK((two - Vec<double, 3>{one * 2.0}).norm() < 1e-12);

    // And clearance reads the same surface set, so a test can ask about
    // non-penetration without integrating anything.
    CHECK_THAT(clearance(probe, both), WithinAbs(0.05, 1e-12));
}

TEST_CASE("The barrier is silent outside its band", "[physics][contact]") {
    const auto ball = chart_of(Sphere<2, double>{1.0});
    auto contact = surface_contact_force<double>(ball, /*d_hat=*/0.1, /*kappa=*/40.0);

    // Far away: exactly zero, not merely small. A force that is nearly
    // zero everywhere is a force that costs everywhere.
    PointMass<3, double> far{1.0, {0.0, 0.0, 5.0}, {}};
    CHECK(contact(far, 0.0).norm() == 0.0);

    // Inside the band: pushing outward, along the surface normal.
    PointMass<3, double> near{1.0, {0.0, 0.0, 1.05}, {}};
    const auto f = contact(near, 0.0);
    CHECK(f.norm() > 0.0);
    CHECK(f[2] > 0.0);
}
