#pragma once

// The join between the contact barrier and the integrators, which is the
// one piece that was missing.
//
// Everything either side of it already shipped and had done for a while:
// `contact.hpp` carries the IPC barrier as pure mathematics, energy and
// gradient and Hessian; `narrow_phase.hpp` turns a `ContactQuery` into
// `ipc_contact_energy` and `ipc_contact_force`; `point_to` answers that
// query for a sphere, a torus, or -- through a generic overload over
// `project` and `normal` -- for any `Surface` at all. And on the other
// side `verlet_step` and its relatives take a *force functor*, not a fixed
// force.
//
// What did not exist was anything that handed one to the other. `ipc_`
// appears nowhere in `xpbd.hpp`, `variational.hpp`, `lgvi.hpp` or
// `symplectic.hpp`, so the barrier was reachable and never reached.
//
// Which is why this is a force adapter rather than a contact integrator.
// Writing a new integrator would have been the obvious move and the wrong
// one: `verlet_step` is already symplectic, already tested, and already
// takes exactly what is needed. A second integrator would have had to earn
// that again.

#include <spatium/core/concepts.hpp>
#include <spatium/physics/mechanics/body.hpp>
#include <spatium/physics/mechanics/narrow_phase.hpp>

#include <functional>
#include <utility>
#include <vector>

namespace spatium::physics::mechanics {

// A force functor that keeps a body out of a surface, in the shape
// `verlet_step` and `rk4_step` already consume.
//
// `d_hat` is the barrier's active band: outside it the force is exactly
// zero, so a body far from anything pays one `point_to` and nothing else.
// The barrier diverges as the distance goes to zero, which is what makes
// non-penetration a property of the energy rather than a correction
// applied afterwards -- and why this composes with a symplectic
// integrator instead of fighting it.
template<Scalar T, typename S>
auto surface_contact_force(const S& surface, T d_hat, T kappa = T{1}) {
    return [&surface, d_hat, kappa](const PointMass<3, T>& body, T) -> Vec<T, 3> {
        const auto q = point_to(body.state.position, surface);
        return ipc_contact_force(q, d_hat, kappa);
    };
}

// The same against several surfaces, summed.
//
// Summed rather than resolved pairwise on purpose: a barrier is a
// potential, and potentials add. A body in the corner of two surfaces is
// pushed out of both at once, with no ordering between them and nothing to
// arbitrate -- which is the property that makes this stay symplectic when
// a sequence of pairwise corrections would not.
template<Scalar T, typename S>
auto surfaces_contact_force(const std::vector<S>& surfaces, T d_hat, T kappa = T{1}) {
    return [&surfaces, d_hat, kappa](const PointMass<3, T>& body, T) -> Vec<T, 3> {
        Vec<T, 3> total{};
        for (const auto& s : surfaces)
            total = Vec<T, 3>{total + ipc_contact_force(point_to(body.state.position, s),
                                                        d_hat, kappa)};
        return total;
    };
}

// Adds a contact force to a force that already exists, so gravity or a
// spring keeps working and the barrier is one more term rather than a
// replacement.
template<typename A, typename B>
auto sum_forces(A a, B b) {
    return [a = std::move(a), b = std::move(b)](const auto& body, auto t) {
        return decltype(a(body, t)){a(body, t) + b(body, t)};
    };
}

// The smallest signed distance from a body to any of the surfaces, which
// is the cheap non-penetration check a search or a test can run without
// integrating anything. Negative means inside.
template<Scalar T, typename S>
T clearance(const PointMass<3, T>& body, const std::vector<S>& surfaces) {
    T least = std::numeric_limits<T>::max();
    for (const auto& s : surfaces) {
        const auto q = point_to(body.state.position, s);
        least = std::min(least, q.inside ? -q.distance : q.distance);
    }
    return least;
}

}  // namespace spatium::physics::mechanics
