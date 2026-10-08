#include "../pair.hpp"

#include <spatium/algebra/calculus.hpp>
#include <spatium/algebra/dual.hpp>
#include <spatium/physics/mechanics/narrow_phase.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <type_traits>
#include <vector>

namespace symmetry {
namespace {

using spatium::Vec;

// ── The force of the contact barrier is minus the gradient of its energy ──
//
// A point a distance d in [0.01, 0.3) off a unit sphere and off a torus (R = 2, r = 0.5), inside
// the barrier's band (d_hat = 0.3). One path is `ipc_contact_force`, the closed-form force along the
// normal of the contact query; the other is the energy `ipc_contact_energy` differentiated through
// `Dual` with respect to the point, negated. Two paths that share the barrier's formula and nothing
// else: the query's normal must be the gradient of its distance, and `ipc_barrier_grad` the
// derivative of `ipc_barrier`. The bound is 1e-12 of the force; the worst measured is 4e-16.
//
// First run: `point_to_torus<Dual>` did not compile (a qualified std::sqrt).
struct ContactDraw { Vec<double, 3> on_sphere, on_torus; };

inline ContactDraw contact_draw(std::size_t i) {
    const double pi = std::numbers::pi;
    const double d = 0.01 + 0.29 * halton(i, 2);
    const double u = 2 * pi * halton(i, 0), v = 2 * pi * halton(i, 1);
    const double rho = 2.0 + (0.5 + d) * std::cos(v);
    const double w = 2 * halton(i, 3) - 1, ph = 2 * pi * halton(i, 4), q = std::sqrt(1 - w * w);
    return {Vec<double, 3>{(1 + d) * q * std::cos(ph), (1 + d) * q * std::sin(ph), (1 + d) * w},
            Vec<double, 3>{rho * std::cos(u), rho * std::sin(u), (0.5 + d) * std::sin(v)}};
}

inline std::vector<double> contact_forces(const ContactDraw& s) {
    namespace m = spatium::physics::mechanics;
    const Vec<double, 3> origin{};
    const auto fs = m::ipc_contact_force(m::point_to_sphere(s.on_sphere, origin, 1.0), 0.3, 1.0);
    const auto ft = m::ipc_contact_force(m::point_to_torus(s.on_torus, spatium::geometry::Torus<double>{.major_radius = 2.0, .minor_radius = 0.5}), 0.3, 1.0);
    return {fs[0], fs[1], fs[2], ft[0], ft[1], ft[2]};
}

inline std::vector<double> contact_forces_from_energy(const ContactDraw& s) {
    namespace m = spatium::physics::mechanics;
    using D = spatium::Dual<double>;
    const Vec<D, 3> origin{D(0.0), D(0.0), D(0.0)};
    const spatium::geometry::Torus<D> torus{.major_radius = D(2.0), .minor_radius = D(0.5)};
    const auto gs = spatium::gradient([&](const Vec<D, 3>& p) {
        return m::ipc_contact_energy(m::point_to_sphere(p, origin, D(1.0)), D(0.3), D(1.0)); }, s.on_sphere);
    const auto gt = spatium::gradient([&](const Vec<D, 3>& p) {
        return m::ipc_contact_energy(m::point_to_torus(p, torus), D(0.3), D(1.0)); }, s.on_torus);
    return {-gs[0], -gs[1], -gs[2], -gt[0], -gt[1], -gt[2]};
}

inline Entry contact_force_is_minus_energy_gradient(std::size_t n = 3000) {
    return {"contact: the barrier's force is minus the gradient of its energy", n, [n] {
        return check<ContactDraw>(
            "contact: the barrier's force is minus the gradient of its energy", Kind::Function, n, contact_draw,
            contact_forces, contact_forces_from_energy,
            [](const std::vector<double>& a, const std::vector<double>& b) {
                double worst = 0.0, scale = 1.0;
                for (std::size_t k = 0; k < a.size(); ++k) { worst = std::max(worst, std::abs(a[k] - b[k])); scale = std::max(scale, std::abs(a[k])); }
                return worst / scale;
            },
            [](const ContactDraw&) { return 1e-12; });
    }};
}

}  // namespace

const Registrar registered_contact_force{"contact_force_is_minus_energy_gradient", 3000, [](std::size_t n) { return contact_force_is_minus_energy_gradient(n); }};

}  // namespace symmetry
