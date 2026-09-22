// Where an explicit contact chain works, and where an implicit solver
// becomes necessary.
//
// The argument this measures: while a Newton step is in the pipeline the
// IR does not compress and there is no interactivity, because a Newton
// step is solve_linear(H, g) then a line search then an update -- a long
// chain whose control flow depends on its own intermediate values, which
// lowers to no kernel and serialises into nothing.
//
// The alternative already exists as of today and is a composition of
// things the library had separately: point_to, then ipc_contact_force,
// then verlet_step. No Hessian, no linear system, no backtracking. Every
// step is arithmetic, so the whole chain is expressible as operations.
//
// So the question stops being "can Newton be replaced" and becomes "where
// is the boundary", which is what this sweeps. Two failures close in from
// opposite sides: a barrier too weak to stop the body, and a barrier too
// stiff for the step. Between them is a window, and it narrows as the
// approach speed rises.
//
// Which reframes what a dispatcher is for. It does not replace Newton; it
// decides where Newton is not needed -- and inside that region the
// description becomes processable. The features are cheap and all three
// are to hand: approach speed, stiffness, step size.
//
// Build: part of the rsc tools target.

#include <spatium/physics/mechanics/contact_force.hpp>
#include <spatium/physics/mechanics/integrator.hpp>
#include <spatium/spaces/chart.hpp>
#include <spatium/spaces/sphere.hpp>
#include <cmath>
#include <print>
using namespace spatium; using namespace spatium::physics::mechanics;
using V3 = Vec<double,3>;

// Drop a body at a sphere and ask two questions of the whole trajectory:
// did it ever get inside, and did the energy stay bounded.
struct Verdict { bool penetrated; double energy_ratio; };

static Verdict run(double kappa, double dt, double v0, double d_hat = 0.1) {
    auto ball = chart_of(Sphere<2,double>{1.0});
    auto contact = surface_contact_force<double>(ball, d_hat, kappa);
    auto gravity = [](const PointMass<3,double>&, double){ return V3{0,0,-1.0}; };
    auto force = sum_forces(gravity, contact);

    PointMass<3,double> b{1.0, {0,0,2.0}, {0,0,-v0}};
    auto energy = [&](const PointMass<3,double>& m){
        const auto q = point_to(m.state.position, ball);
        const double barrier = (!q.inside && q.distance < d_hat)
            ? kappa * ipc_barrier(q.distance, d_hat) : 0.0;
        return 0.5*m.mass*m.state.velocity.norm_squared() + m.mass*1.0*m.state.position[2] + barrier;
    };
    const double e0 = energy(b);
    bool pen = false; double worst = e0;
    const int steps = static_cast<int>(6.0 / dt);
    for (int i = 0; i < steps; ++i) {
        verlet_step(b, force, dt, i*dt);
        if (!std::isfinite(b.state.position.norm())) return {true, 1e9};
        if (point_to(b.state.position, ball).inside) pen = true;
        worst = std::max(worst, energy(b));
    }
    return {pen, worst / e0};
}

int main() {
    std::println("Explicit barrier + Verlet: where it holds, and where Newton starts");
    std::println("drop onto a unit sphere, 6 s, gravity 1.0, d_hat 0.1");
    std::println("");
    for (double v0 : {0.0, 5.0, 50.0}) {
    std::println("  approach speed {:.0f}", v0);
    std::println("  {:>8} | {:>10} | {:>10} | {:>10} | {:>10}", "kappa", "dt=1e-2", "dt=3e-3", "dt=1e-3", "dt=3e-4");
    for (double kappa : {1.0, 100.0, 10000.0, 1000000.0}) {
        std::print("  {:>8.0f} |", kappa);
        for (double dt : {1e-2, 3e-3, 1e-3, 3e-4}) {
            auto v = run(kappa, dt, v0);
            if (v.penetrated)           std::print(" {:>10} |", "THROUGH");
            else if (v.energy_ratio > 2.0) std::print(" {:>10.1f}x |", v.energy_ratio);
            else                        std::print(" {:>10} |", "ok");
        }
        std::println("");
    }
    std::println("");
    }
    std::println("THROUGH = penetrated. A number = energy grew by that factor, so the");
    std::println("explicit step is unstable there. ok = non-penetrating and bounded.");
}
