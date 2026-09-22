#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <spatium/physics/mechanics/contact_force.hpp>
#include <spatium/physics/mechanics/integrator.hpp>
#include <spatium/physics/mechanics/symplectic.hpp>
#include <spatium/spaces/euclidean.hpp>
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

// Whether the assembled force is still symplectic, which is the oracle a
// search over compositions would depend on and which had never been run.
//
// The claim inherited from the parts is that a barrier force composed into
// a symplectic integrator stays symplectic, because a barrier is a
// potential and `verlet_step` is symplectic for any potential. Inherited
// claims are the ones worth measuring: a barrier is a *stiff* potential,
// and stiffness is exactly what breaks explicit integrators at a given
// step size. So the interesting number is not whether it holds but where
// it stops holding.
TEST_CASE("The assembled contact force keeps the symplectic form",
          "[physics][contact][symplectic]") {
    using Bundle = CotangentBundle<Euclidean<3, double>>;
    const auto ball = chart_of(Sphere<2, double>{1.0});

    // Inside the barrier's active band, where the force is not zero and
    // there is therefore something to measure. Outside it the map is
    // trivially symplectic and the check would pass without meaning.
    const Bundle::State s0{{0.0, 0.0, 1.08}, {0.0, 0.0, -0.4}};

    auto drift_at = [&](double kappa, double eps) {
        auto contact = surface_contact_force<double>(ball, /*d_hat=*/0.2, kappa);
        auto step = [&](const Bundle::State& s, double h) {
            PointMass<3, double> body{1.0, s.q, s.p};
            for (int i = 0; i < 100; ++i) verlet_step(body, contact, h, i * h);
            return Bundle::State{body.state.position, body.state.velocity};
        };
        return verify_symplecticity_drift<Bundle>(s0, step, eps, 2e-3);
    };

    // The assertion is about *scaling*, not about an absolute bound, and
    // the difference matters. `verify_symplecticity_drift` estimates the
    // differential by finite difference, so a symplectic map reports O(eps)
    // -- its own truncation error -- and a map that fails to preserve the
    // form reports something independent of eps.
    //
    // A stiff barrier makes that truncation term large: measured here at
    // 3.2 * eps, which at eps = 1e-4 is 3.2e-4 and would fail any
    // plausible absolute threshold while being entirely correct. Halving
    // eps has to halve the drift; that is what symplectic means for this
    // probe, and it is a claim a non-symplectic map cannot satisfy.
    const double d1 = drift_at(200.0, 1e-4);
    const double d2 = drift_at(200.0, 1e-5);
    const double d3 = drift_at(200.0, 1e-6);

    CHECK(d1 / d2 > 8.0);    // an order of magnitude down in eps ...
    CHECK(d1 / d2 < 12.0);   // ... is an order of magnitude down in drift
    CHECK(d2 / d3 > 8.0);
    CHECK(d2 / d3 < 12.0);

    // Free flight, where the force is identically zero, is exactly
    // symplectic and pins the other end of the scale.
    CHECK(drift_at(0.0, 1e-4) < 1e-10);
}
