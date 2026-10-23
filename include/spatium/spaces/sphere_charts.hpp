#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/spaces/constant_curvature.hpp>
#  include <spatium/spaces/sphere.hpp>
#  include <cmath>
#  include <cstddef>
#endif

// The two stereographic charts of the sphere S^N of radius R = 1 / sqrt(kappa),
// as the sphere `Sphere<N>` in R^(N+1) sees them and as `ConstantCurvature<N>{kappa}`
// does: the chart's own space, whose exp, log and distance are the model's.
//
// One chart cannot cover the sphere (the model leaves out the point at
// infinity), so a flight that crosses the equator needs two. Convention: the
// polar axis is the last coordinate, y_N; a point is y = (y', y_N), y' in R^N.
//
//   north chart  x = y' / (1 + y_N / R)   projection from the south pole, the
//                                         ball |x| < R is the hemisphere y_N > 0,
//                                         the north pole is x = 0
//   south chart  x = y' / (1 - y_N / R)   projection from the north pole, the
//                                         ball |x| < R is the hemisphere y_N < 0,
//                                         the south pole is x = 0
//
// Both project along the same axes, so a point of the equator (y_N = 0, |x| = R)
// has the SAME coordinates x = y' in the two charts -- no reflection is applied
// -- and the transition is the inversion  x -> x / (kappa |x|^2)  (|x| |x'| = R^2,
// the same direction). The inversion reverses orientation: seen from outside, the
// south chart is the mirror image of the north one, as a map and its reflection
// in the equator would be. Charts are conformal: a tangent vector u in a chart
// is the vector w of R^(N+1) of length lambda(x) |u|, lambda = 2 / (1 + kappa |x|^2),
// the conformal factor of ConstantCurvature; `push_*` and `pull_*` are the two
// ways of that map, the exact derivatives of the point maps. A geodesic of the
// sphere (a great circle) is the geodesic of the chart: the embedded flight and
// `ConstantCurvature::exp_map` with the pushed velocity trace the same curve.
//
// The pole a chart leaves out has no coordinates (infinite ones: the maps
// return inf or NaN there).

SPATIUM_EXPORT namespace spatium {

template<std::size_t N, Scalar T = double>
struct SphereCharts {
    using Embedded = Vec<T, N + 1>;
    using Coords   = Vec<T, N>;
    enum class Chart { North, South };

    T kappa{1};

    constexpr SphereCharts() = default;
    explicit constexpr SphereCharts(T k) : kappa(k) {}

    T radius() const { return T{1} / std::sqrt(kappa); }
    ConstantCurvature<N, T> space() const { return ConstantCurvature<N, T>{kappa}; }
    Sphere<N, T> sphere() const { return Sphere<N, T>{radius()}; }

    // The hemisphere a point belongs to; the equator is the north's.
    Chart chart_of(const Embedded& y) const { return y[N] >= T{0} ? Chart::North : Chart::South; }

    // embedding -> chart
    Coords to_north(const Embedded& y) const { return to_chart(y, T{1}); }
    Coords to_south(const Embedded& y) const { return to_chart(y, T{-1}); }
    Coords to_chart(const Embedded& y, Chart c) const { return to_chart(y, c == Chart::North ? T{1} : T{-1}); }

    // chart -> embedding
    Embedded from_north(const Coords& x) const { return from_chart(x, T{1}); }
    Embedded from_south(const Coords& x) const { return from_chart(x, T{-1}); }
    Embedded from_chart(const Coords& x, Chart c) const { return from_chart(x, c == Chart::North ? T{1} : T{-1}); }

    // The transition between the two charts, an involution.
    Coords transition(const Coords& x) const { return Coords{x / T(kappa * x.norm_squared())}; }

    // Its derivative: the chart vector u at x carried to the other chart at transition(x).
    Coords transition_push(const Coords& x, const Coords& u) const {
        const T x2 = x.norm_squared();
        return Coords{(u - x * T(T{2} * x.dot(u) / x2)) / T(kappa * x2)};
    }

    // tangent vector of the embedding at y -> chart vector at the chart's point of y
    Coords push_north(const Embedded& y, const Embedded& w) const { return push(y, w, T{1}); }
    Coords push_south(const Embedded& y, const Embedded& w) const { return push(y, w, T{-1}); }
    // chart vector u at x -> tangent vector of the embedding at the point of x
    Embedded pull_north(const Coords& x, const Coords& u) const { return pull(x, u, T{1}); }
    Embedded pull_south(const Coords& x, const Coords& u) const { return pull(x, u, T{-1}); }

private:
    // sign +1: north, -1: south.
    Coords to_chart(const Embedded& y, T sign) const {
        const T den = T{1} + sign * y[N] / radius();
        Coords x;
        for (std::size_t i = 0; i < N; ++i) x[i] = y[i] / den;
        return x;
    }

    Embedded from_chart(const Coords& x, T sign) const {
        const T d = T{1} + kappa * x.norm_squared();
        Embedded y;
        for (std::size_t i = 0; i < N; ++i) y[i] = T{2} * x[i] / d;
        y[N] = sign * radius() * (T{1} - kappa * x.norm_squared()) / d;
        return y;
    }

    Coords push(const Embedded& y, const Embedded& w, T sign) const {
        const T R = radius();
        const T den = T{1} + sign * y[N] / R;
        Coords u;
        for (std::size_t i = 0; i < N; ++i)
            u[i] = w[i] / den - y[i] * (sign * w[N] / R) / (den * den);
        return u;
    }

    Embedded pull(const Coords& x, const Coords& u, T sign) const {
        const T d = T{1} + kappa * x.norm_squared();
        const T xu = x.dot(u);
        Embedded w;
        for (std::size_t i = 0; i < N; ++i)
            w[i] = T{2} * u[i] / d - T{4} * kappa * x[i] * xu / (d * d);
        w[N] = -sign * T{4} * radius() * kappa * xu / (d * d);
        return w;
    }
};

} // namespace spatium
