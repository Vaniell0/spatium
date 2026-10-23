#pragma once

// Gravity of many sources on the closed 3-sphere.
//
// Poisson's equation on a closed manifold has a solution only for a source of
// zero total mass, so, as on the torus (gravity_field.hpp), the field is that
// of the sources plus a uniform background of opposite total mass: the
// potential of one unit mass is the Green's function of the Laplace-Beltrami
// operator, Delta G = delta - 1/V, V the volume of S^3. For the unit sphere
// (V = 2 pi^2), at geodesic angle psi from the source,
//
//     G(psi) = -(pi - psi) cot(psi) / (4 pi^2)            up to a constant,
//
// which is -1/(4 pi psi) near the source (the flat 1/r law) and has the
// finite value 1/(4 pi^2) at the antipode, where the field is smooth and zero.
// The mean of G over the sphere is -1/(8 pi^2). It was checked against the
// zonal eigenfunction expansion (eigenvalues n(n+2), degeneracy (n+1)^2, zonal
// functions sin((n+1) psi) / sin(psi)), against the equation away from the
// source, and by the flux through a small sphere (tests/test_sphere_gravity.cpp).
// On a sphere of radius R the potential is G(psi) / R and the acceleration
// scales with 1 / R^2.
//
// Positions are points of the sphere in R^4 (norm R); the acceleration is a
// vector of R^4 tangent to the sphere at the point, towards the source along
// the great circle:  a = G m H(psi) (s - cos(psi) x) / R^2  with
// H = G'(psi) / sin(psi) = [(pi - psi) + sin(psi) cos(psi)] / (4 pi^2 sin^3(psi)),
// which is finite at the antipode (1 / 6 pi^2) where s - cos(psi) x vanishes.
//
// Softening is Plummer's in the angle, written so that it stays smooth through
// the antipode: tan^2(psi_s / 2) = tan^2(psi / 2) + eps^2 / 4, so psi_s^2 is
// psi^2 + eps^2 near the source, and the force of the softened potential
// G(psi_s) is  a = G m H(psi_s) (4 / (4 + eps^2 b / 4))^2 (s - cos(psi) x) / R^2,
// b = |x + s|^2 (unit sphere). The pair law is symmetric in x and s, so the
// torque of every pair cancels: the angular momentum of the whole system in
// R^4 (six components, the generators of SO(4)) is conserved by the flow,
// softened or not, and so is the energy.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/physics/mechanics/manifold_body.hpp>
#  include <spatium/spaces/sphere.hpp>
#  include <array>
#  include <cmath>
#  include <cstddef>
#  include <limits>
#  include <numbers>
#  include <stdexcept>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::physics::mechanics {

namespace sphere_gravity_detail {

// (pi - psi) + sin(psi) cos(psi), written in u = pi - psi as u - sin(2u) / 2,
// by its series for small u (the two terms cancel to u^3 there).
template<Scalar T>
T green_numerator(T u) {
    if (u < T{0.5}) {
        const T x = T{2} * u;
        T term = x * x * x / T{6}, sum{};
        for (int k = 1; k <= 12; ++k) {
            sum += term;
            term *= -x * x / (T(2 * k + 2) * T(2 * k + 3));
        }
        return sum / T{2};
    }
    return u - std::sin(T{2} * u) / T{2};
}

} // namespace sphere_gravity_detail

namespace sphere_gravity_detail {

// G and G' / sin from sin(psi) and cos(psi) (sin >= 0), each angle taken from
// atan2 so that neither end of (0, pi) loses digits.
template<Scalar T>
T green_sc(T sn, T cs) {
    constexpr T four_pi2 = T{4} * std::numbers::pi_v<T> * std::numbers::pi_v<T>;
    const T u = std::atan2(sn, -cs);                 // pi - psi
    if (sn <= T{0}) return cs > T{0} ? -std::numeric_limits<T>::infinity() : T{1} / four_pi2;
    return -u * cs / (sn * four_pi2);
}

template<Scalar T>
T slope_sc(T sn, T cs) {
    constexpr T pi2 = std::numbers::pi_v<T> * std::numbers::pi_v<T>;
    const T u = std::atan2(sn, -cs);
    if (u < T{1e-4}) return (T{1} / (T{6} * pi2)) * (T{1} + T{0.3} * u * u);   // antipode
    return green_numerator(u) / (T{4} * pi2 * sn * sn * sn);
}

} // namespace sphere_gravity_detail

// G(psi) on the unit sphere, psi in [0, pi]; minus infinity at 0.
template<Scalar T = double>
T s3_green(T psi) {
    return sphere_gravity_detail::green_sc(std::sin(psi), std::cos(psi));
}

// G'(psi) / sin(psi) on the unit sphere: the field divided by the direction's
// length; finite, 1 / (6 pi^2), at the antipode.
template<Scalar T = double>
T s3_green_slope_over_sin(T psi) {
    return sphere_gravity_detail::slope_sc(std::sin(psi), std::cos(psi));
}

template<Scalar T = double>
struct SphereGravityLaw {
    T G{1};
    T softening{0};            // Plummer epsilon, an angle
    T radius{1};
    T max_acceleration{0};     // 0: no clamp
};

template<Scalar T = double>
struct SphereSource {
    Vec<T, 4> position{};      // a point of the sphere of the law's radius
    T         mass{};          // gravitational charge
};

template<Scalar T = double>
class SphereGravity {
public:
    using ScalarType = T;
    using Law = SphereGravityLaw<T>;
    using Source = SphereSource<T>;
    using Point = Vec<T, 4>;
    static constexpr std::size_t npos = ~std::size_t{0};

    explicit SphereGravity(Law law = Law{}) : law_(law) {
        if (law_.softening < T{0}) throw std::invalid_argument("SphereGravity: negative softening");
        if (!(law_.radius > T{0})) throw std::invalid_argument("SphereGravity: radius must be positive");
        q_ = law_.softening * law_.softening / T{4};
    }

    const Law& law() const { return law_; }
    void set_sources(std::vector<Source> s) { sources_ = std::move(s); }
    const std::vector<Source>& sources() const { return sources_; }

    // Potential of one unit source at s, at x (G excluded): G(psi_s) / R. Minus
    // infinity for coincident points with no softening.
    T pair_potential(const Point& x, const Point& s) const {
        const auto k = angles(x, s);
        return sphere_gravity_detail::green_sc(k.sin_s, k.cos_s) / law_.radius;
    }

    // Acceleration towards one unit source at s, at x (G excluded), tangent at x.
    Point pair_acceleration(const Point& x, const Point& s) const {
        const auto k = angles(x, s);
        if (k.a == T{0}) return Point{};                       // coincident: no direction
        const T h = sphere_gravity_detail::slope_sc(k.sin_s, k.cos_s) * k.f * k.f;
        const T inv_r2 = T{1} / (law_.radius * law_.radius);
        Point out;
        for (std::size_t i = 0; i < 4; ++i)
            out[i] = h * (k.xs[i] - k.cos_psi * k.xx[i]) * inv_r2;
        return out;
    }

    // Potential at x of all sources but `skip`, times G (up to a constant).
    T potential(const Point& x, std::size_t skip = npos) const {
        T phi{};
        for (std::size_t j = 0; j < sources_.size(); ++j) {
            if (j == skip) continue;
            phi += sources_[j].mass * pair_potential(x, sources_[j].position);
        }
        return law_.G * phi;
    }

    // Sum of m_i m_j G(psi_ij) over pairs, G included; a constant differs from the
    // physical energy by a multiple of the total mass squared, which no motion changes.
    T potential_energy() const {
        T u{};
        for (std::size_t i = 0; i < sources_.size(); ++i)
            u += sources_[i].mass * potential(sources_[i].position, i);
        return T{0.5} * u;
    }

    // Acceleration per unit coupling at x: the sum over the sources but `skip`,
    // clamped if the law says so, times the body's coupling (gravitational
    // charge over inertia). Zero for an empty source set.
    Point acceleration(const Point& x, std::size_t skip = npos, T coupling = T{1}) const {
        Point sum{};
        for (std::size_t j = 0; j < sources_.size(); ++j) {
            if (j == skip) continue;
            sum = Point{sum + pair_acceleration(x, sources_[j].position) * sources_[j].mass};
        }
        sum = Point{sum * law_.G};
        if (law_.max_acceleration > T{0}) {
            const T n2 = sum.dot(sum);
            const T m2 = law_.max_acceleration * law_.max_acceleration;
            if (n2 > m2) sum = Point{sum * (law_.max_acceleration / std::sqrt(n2))};
        }
        return Point{sum * coupling};
    }

private:
    Law law_;
    T q_{};
    std::vector<Source> sources_;

    struct Pair {
        Point xs, xx;          // unit vectors: s and x
        T a, b, f;             // |x - s|^2, |x + s|^2, 4 / (4 + q b)
        T cos_psi, sin_s, cos_s;
    };

    // Everything the pair law needs, from |x -+ s|^2 so neither end loses digits.
    Pair angles(const Point& x, const Point& s) const {
        Pair k;
        k.xx = Point{x / law_.radius};
        k.xs = Point{s / law_.radius};
        const Point dm{k.xx - k.xs}, dp{k.xx + k.xs};
        k.a = dm.dot(dm);
        k.b = dp.dot(dp);
        k.f = T{4} / (T{4} + q_ * k.b);
        k.cos_psi = (k.b - k.a) / T{4};
        const T den = T{4} + q_ * k.b;
        k.cos_s = (k.b - k.a - q_ * k.b) / den;
        k.sin_s = T{2} * std::sqrt(k.b * (k.a + q_ * k.b)) / den;
        return k;
    }
};

// ── N bodies on the 3-sphere ──────────────────────────────────

using SphereBody3 = PointOnManifold<spatium::Sphere<3, double>>;

// Advance all bodies one step: kick, great-circle drift, kick, with the force of
// every other body (a body's mass is both its charge and its inertia). The
// drift of every body is exact, the kicks are the same pair law for all, so the
// angular momentum is conserved to rounding and the energy to second order.
inline void sphere_nbody_step(std::vector<SphereBody3>& bodies, SphereGravity<double>& field,
                              const spatium::Sphere<3, double>& S, double dt) {
    using Src = SphereSource<double>;
    const auto publish = [&] {
        std::vector<Src> src;
        src.reserve(bodies.size());
        for (const auto& b : bodies) src.push_back(Src{b.position, b.mass});
        field.set_sources(std::move(src));
    };
    const auto kick_all = [&](double h) {
        for (std::size_t i = 0; i < bodies.size(); ++i)
            kick(bodies[i], S, field.acceleration(bodies[i].position, i), h);
    };
    publish();
    kick_all(0.5 * dt);
    for (auto& b : bodies) { geodesic_step(b, S, dt); renormalise(b, S); }
    publish();
    kick_all(0.5 * dt);
}

// Total energy of the bodies: kinetic plus the pair potential (G(psi) summed over pairs).
inline double sphere_nbody_energy(const std::vector<SphereBody3>& bodies, SphereGravity<double>& field) {
    using Src = SphereSource<double>;
    std::vector<Src> src;
    for (const auto& b : bodies) src.push_back(Src{b.position, b.mass});
    field.set_sources(std::move(src));
    double kin = 0;
    for (const auto& b : bodies) kin += 0.5 * b.mass * b.velocity.dot(b.velocity);
    return kin + field.potential_energy();
}

// The six conserved components L_ab = sum m (x_a v_b - x_b v_a), a < b: the
// generators of SO(4), whose action on the sphere carries the system to itself.
template<Scalar T>
std::array<T, 6> sphere_angular_momentum(const std::vector<PointOnManifold<spatium::Sphere<3, T>>>& bodies) {
    std::array<T, 6> L{};
    for (const auto& b : bodies) {
        std::size_t k = 0;
        for (std::size_t a = 0; a < 4; ++a)
            for (std::size_t c = a + 1; c < 4; ++c)
                L[k++] += b.mass * (b.position[a] * b.velocity[c] - b.position[c] * b.velocity[a]);
    }
    return L;
}

} // namespace spatium::physics::mechanics
