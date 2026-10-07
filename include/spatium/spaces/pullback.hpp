#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/dual.hpp>
#  include <spatium/algebra/matrix.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/concepts.hpp>
#  include <array>
#  include <cstddef>
#  include <utility>
#endif

// The pullback of a metric through a map: a space given by coordinates x in R^K and
// a map f: R^K -> R^M into an ambient space with a metric, is the K-dimensional space
// whose metric is
//
//   g(x)_ij = <d_i f(x), d_j f(x)>_ambient,   d_i f = the partial of f along e_i
//
// with the partials exact, by a Dual seeded along each coordinate. The ambient metric
// is a callable `(point, u, v) -> scalar` generic over the scalar, so a Euclidean one
// (a surface in R^3: the induced metric J^T J), a Lorentzian one (the hyperboloid as a
// ball: the Poincare metric) and a sphere's all go through the same code. The result
// is the callable `spaces::metric_chart` takes, so the pulled-back space has Christoffel
// symbols, geodesics, exp and log with nothing more written.
//
//   pullback_metric<2, 3>(stereographic, EuclideanMetric{})   the unit sphere's metric
//   pullback_metric<2, 3>(poincare_to_hyperboloid, MinkowskiMetric{})   the hyperbolic plane's
//
// `ParametricSurface`'s induced metric is the first case with its own copy of the code
// no longer: it is this, with the Euclidean metric of R^3 (tests/test_pullback.cpp holds
// the two to the same numbers and to the closed forms).

SPATIUM_EXPORT namespace spatium::spaces {

// <a, b> of R^M, at any point.
struct EuclideanMetric {
    template<class S, std::size_t M>
    S operator()(const Vec<S, M>&, const Vec<S, M>& a, const Vec<S, M>& b) const { return a.dot(b); }
};

// -a0 b0 + a1 b1 + ...: the Minkowski form, signature (-, +, ..., +).
struct MinkowskiMetric {
    template<class S, std::size_t M>
    S operator()(const Vec<S, M>&, const Vec<S, M>& a, const Vec<S, M>& b) const {
        S sum = -a[0] * b[0];
        for (std::size_t i = 1; i < M; ++i) sum = sum + a[i] * b[i];
        return sum;
    }
};

// f: Vec<Dual<S>, K> -> Vec<Dual<S>, M>, called on any scalar S (so a Dual of a Dual for
// the second derivatives the Christoffel symbols need).
template<std::size_t K, std::size_t M, class F, class Ambient = EuclideanMetric>
struct PullbackMetric {
    F f;
    Ambient ambient{};

    template<class S>
    Matrix<S, K, K> operator()(const Vec<S, K>& x) const {
        std::array<Vec<S, M>, K> partial{};
        Vec<S, M> at{};
        for (std::size_t i = 0; i < K; ++i) {
            Vec<Dual<S>, K> xi;
            for (std::size_t j = 0; j < K; ++j)
                xi[j] = (i == j) ? Dual<S>::variable(x[j]) : Dual<S>::constant(x[j]);
            const Vec<Dual<S>, M> y = f(xi);
            for (std::size_t a = 0; a < M; ++a) {
                partial[i][a] = y[a].deriv;
                if (i == 0) at[a] = y[a].value;
            }
        }
        Matrix<S, K, K> g{};
        for (std::size_t i = 0; i < K; ++i)
            for (std::size_t j = i; j < K; ++j) {
                const S gij = ambient(at, partial[i], partial[j]);
                g(i, j) = gij;
                g(j, i) = gij;
            }
        return g;
    }
};

template<std::size_t K, std::size_t M, class F, class Ambient = EuclideanMetric>
PullbackMetric<K, M, F, Ambient> pullback_metric(F f, Ambient ambient = {}) {
    return {std::move(f), std::move(ambient)};
}

} // namespace spatium::spaces
