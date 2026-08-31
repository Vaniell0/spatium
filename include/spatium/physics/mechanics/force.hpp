#pragma once

// Force models for PointMass-based flat-space mechanics. Each force is a
// callable
//     Vec<T,N> f(const PointMass<N,T>& body, T time)
// returning the force vector acting on the body in N-D Euclidean space.
//
// This header ships the textbook trio: uniform Gravity, inverse-square
// Gravity (Newton), Hookean Spring, viscous Damper. Custom forces are just
// lambdas with the right signature.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/epsilon.hpp>
#  include <spatium/physics/mechanics/body.hpp>
#  include <cmath>
#  include <cstddef>
#  include <tuple>
#  include <utility>
#endif

SPATIUM_EXPORT namespace spatium::physics::mechanics {

// ── Uniform gravity (constant acceleration) ───────────────────
template<std::size_t N, Scalar T = double>
struct UniformGravity {
    Vec<T, N> g{};  // acceleration vector (m/s²); Earth surface = (0, -9.81, 0) typical.

    constexpr UniformGravity() = default;
    constexpr explicit UniformGravity(Vec<T, N> accel) : g(accel) {}

    constexpr Vec<T, N> operator()(const PointMass<N, T>& body, T /*t*/) const {
        return Vec<T, N>{g * body.mass};   // F = m·g
    }
};

// ── Newtonian inverse-square gravity from a point source ─────
// F = -G · M · m · (r - r₀) / |r - r₀|³
template<std::size_t N, Scalar T = double>
struct PointGravity {
    Vec<T, N> source{};   // source position
    T         GM{T{1}};   // G · M product (subsume the attractor mass)

    constexpr PointGravity() = default;
    constexpr PointGravity(Vec<T, N> source_pos, T G_times_M)
        : source(source_pos), GM(G_times_M) {}

    Vec<T, N> operator()(const PointMass<N, T>& body, T /*t*/) const {
        Vec<T, N> r = Vec<T, N>{body.state.position - source};
        T r2 = r.dot(r);
        if (r2 < epsilon<T>()) return Vec<T, N>{};  // singular — caller should bail
        T inv_r3 = T{1} / (r2 * std::sqrt(r2));
        return Vec<T, N>{r * (-GM * body.mass * inv_r3)};
    }
};

// ── Hookean spring: F = -k(r - rest_anchor) ───────────────────
template<std::size_t N, Scalar T = double>
struct Spring {
    Vec<T, N> anchor{};
    T         k{T{1}};        // stiffness (N/m)
    T         rest_length{T{0}};

    constexpr Spring() = default;
    constexpr Spring(Vec<T, N> a, T stiffness, T rest = T{0})
        : anchor(a), k(stiffness), rest_length(rest) {}

    Vec<T, N> operator()(const PointMass<N, T>& body, T /*t*/) const {
        Vec<T, N> d = Vec<T, N>{body.state.position - anchor};
        T len = std::sqrt(d.dot(d));
        if (len < epsilon<T>()) return Vec<T, N>{};
        T extension = len - rest_length;
        T scale = -k * extension / len;
        return Vec<T, N>{d * scale};
    }
};

// ── Linear damper: F = -c · v ─────────────────────────────────
template<std::size_t N, Scalar T = double>
struct Damper {
    T c{T{1}};

    constexpr Damper() = default;
    constexpr explicit Damper(T coeff) : c(coeff) {}

    constexpr Vec<T, N> operator()(const PointMass<N, T>& body, T /*t*/) const {
        return Vec<T, N>{body.state.velocity * (-c)};
    }
};

// ── Sum of forces ─────────────────────────────────────────────
// Blender for an arbitrary tuple of force callables. Use `composite_force(g,
// spring, damper)` to get a single object you can pass to integrators.
template<typename... Forces>
struct CompositeForce {
    std::tuple<Forces...> parts;

    constexpr explicit CompositeForce(Forces... fs) : parts(std::move(fs)...) {}

    template<std::size_t N, Scalar T>
    Vec<T, N> operator()(const PointMass<N, T>& body, T t) const {
        Vec<T, N> total{};
        std::apply([&](auto const&... f) {
            ((total = Vec<T, N>{total + f(body, t)}), ...);
        }, parts);
        return total;
    }
};

template<typename... Forces>
constexpr auto composite_force(Forces... fs) {
    return CompositeForce<Forces...>{std::move(fs)...};
}

} // namespace spatium::physics::mechanics
