#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/causal.hpp>
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/algebra/matrix.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <cmath>
#  include <limits>
#endif

// Minkowski spacetime in geometric units, c = 1: events are (t, x_1, ..., x_N)
// -- (ct, 3d) for N = 3 -- with the signature (-, +, ..., +), the same Minkowski
// form `Hyperbolic` is written on.
//
// c is not a parameter of the type and not a field: nothing flows through the
// operations, and `Kerr` and `Schwarzschild` are written the same way (G = c = 1).
// Where a dimensional quantity is wanted, c comes back as an ordinary scalar
// argument of an ordinary function (`kinetic_energy` below), and a limit in c --
// the Newtonian one, and the post-Newtonian corrections in 1/c^2 -- is a limit of
// that function, which `Series` reads exactly (tests/test_minkowski.cpp).
//
// A Lorentzian space is a different concept from a Riemannian one, not a weaker
// one (core/concepts.hpp, `LorentzianManifold`): there is no distance. What it has
// is the squared interval, a causal class, a proper time between events a signal
// can join, and the reverse triangle inequality, which `verify_lorentzian` holds.
//
// The velocities of special relativity are a hyperbolic space: the four-velocities
// gamma (1, v) are exactly the points of `Hyperbolic<N>` (the same form, the same
// sheet), a rapidity is a hyperbolic distance, and two boosts compose as isometries
// of it -- so the Thomas-Wigner rotation of two boosts is the angular defect of a
// hyperbolic triangle, two paths that share nothing.

SPATIUM_EXPORT namespace spatium {

template<std::size_t N, Scalar T = double>
struct Minkowski {
    using ScalarType    = T;
    using PointType     = Vec<T, N + 1>;
    using TangentVector = Vec<T, N + 1>;
    using BoostMatrix   = Matrix<T, N + 1, N + 1>;

    // N spatial dimensions, N + 1 for the spacetime.
    static constexpr std::size_t dimension = N + 1;
    static constexpr bool is_complete = true;

    static constexpr T minkowski(const PointType& a, const PointType& b) {
        T sum = -a[0] * b[0];
        for (std::size_t i = 1; i <= N; ++i) sum += a[i] * b[i];
        return sum;
    }

    constexpr bool contains(const PointType&) const { return true; }

    // The squared interval, negative for timelike separation.
    ScalarType interval(const PointType& a, const PointType& b) const {
        const PointType d = b - a;
        return minkowski(d, d);
    }

    // Null within `tolerance` (relative to the larger of the time and space parts).
    Causal causal(const PointType& a, const PointType& b) const { return causal(a, b, T(epsilon<T>())); }
    Causal causal(const PointType& a, const PointType& b, T tolerance) const {
        const PointType d = b - a;
        T space{0};
        for (std::size_t i = 1; i <= N; ++i) space += d[i] * d[i];
        const T time2 = d[0] * d[0];
        const T s = space - time2;
        const T scale = time2 > space ? time2 : space;
        if (s < -tolerance * (scale > T{1} ? scale : T{1})) return Causal::Timelike;
        if (s > tolerance * (scale > T{1} ? scale : T{1})) return Causal::Spacelike;
        return Causal::Null;
    }

    // sqrt(-interval): the time a clock carried from a to b reads, for events a
    // signal can join; NaN for a spacelike pair, not a wrong number.
    ScalarType proper_time(const PointType& a, const PointType& b) const {
        using std::sqrt;
        const T s = interval(a, b);
        if (!(s <= T{0})) return T(std::numeric_limits<double>::quiet_NaN());
        return sqrt(T(-s));
    }

    // sqrt(interval): the length measured in the frame where a and b are simultaneous.
    ScalarType proper_length(const PointType& a, const PointType& b) const {
        using std::sqrt;
        const T s = interval(a, b);
        if (!(s >= T{0})) return T(std::numeric_limits<double>::quiet_NaN());
        return sqrt(s);
    }

    // b is in the causal future of a: timelike or null separated, and later.
    bool precedes(const PointType& a, const PointType& b) const {
        return causal(a, b) != Causal::Spacelike && b[0] > a[0];
    }

    // Geodesics are straight lines.
    PointType exp_map(const PointType& p, const TangentVector& v, ScalarType t) const {
        return PointType{p + v * t};
    }
    TangentVector log_map(const PointType& p, const PointType& q) const { return TangentVector{q - p}; }

    constexpr ScalarType metric_at(const PointType&, const TangentVector& u, const TangentVector& v) const {
        return minkowski(u, v);
    }

    static constexpr PointType origin() { return PointType{}; }

    // ── boosts ────────────────────────────────────────────────

    // gamma = 1 / sqrt(1 - v^2) for a velocity |v| < 1 (c = 1); NaN at or beyond light.
    static T gamma(const Vec<T, N>& v) {
        using std::sqrt;
        const T b2 = v.norm_squared();
        if (!(b2 < T{1})) return T(std::numeric_limits<double>::quiet_NaN());
        return T(T{1} / sqrt(T(T{1} - b2)));
    }

    // gamma (1, v): the four-velocity, a point of Hyperbolic<N>.
    static PointType four_velocity(const Vec<T, N>& v) {
        const T g = gamma(v);
        PointType u{};
        u[0] = g;
        for (std::size_t i = 0; i < N; ++i) u[i + 1] = g * v[i];
        return u;
    }

    // The boost that takes an event at rest to one moving with velocity v:
    //   [ gamma        gamma v^T                    ]
    //   [ gamma v      I + gamma^2/(gamma + 1) v v^T ]
    // gamma^2/(gamma + 1) is (gamma - 1)/v^2 without the 0/0 at v = 0.
    static BoostMatrix boost(const Vec<T, N>& v) {
        const T g = gamma(v);
        const T k = T(g * g / (g + T{1}));
        BoostMatrix m = BoostMatrix::identity();
        m(0, 0) = g;
        for (std::size_t i = 0; i < N; ++i) {
            m(0, i + 1) = g * v[i];
            m(i + 1, 0) = g * v[i];
            for (std::size_t j = 0; j < N; ++j) m(i + 1, j + 1) = m(i + 1, j + 1) + k * v[i] * v[j];
        }
        return m;
    }
};

// The kinetic energy (gamma - 1) m c^2 as an ordinary function of c, generic over
// the scalar: on a `Series` with c = 1/eps it is the post-Newtonian expansion,
// 1/2 m v^2 + 3/8 m v^4 / c^2 + ..., read exactly, where a double at large c is
// 1 - 1 of two nearly equal numbers.
template<Scalar T>
T kinetic_energy(const T& mass, const T& speed, const T& c) {
    using std::sqrt;
    const T beta2 = T(speed * speed / (c * c));
    const T gamma = T(T{1} / sqrt(T(T{1} - beta2)));
    return T((gamma - T{1}) * mass * c * c);
}

static_assert(LorentzianManifold<Minkowski<3>>);
static_assert(!MetricSpace<Minkowski<3>>, "a Lorentzian space has no distance");

} // namespace spatium
