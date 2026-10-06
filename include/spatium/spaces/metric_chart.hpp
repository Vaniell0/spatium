#pragma once

// A space given by nothing but its metric in coordinates.
//
// A Riemannian manifold is a metric, and everything else a space offers --
// geodesics, exp, log, distance -- follows from it. core/access.hpp states
// the basis (exp, log, a metric) and derives distance, midpoint and the rest
// from it; this header closes the other end: given the metric alone, exp and
// log are derived, by the geodesic equation
//
//     x''^k = -Gamma^k_{ij}(x) x'^i x'^j,
//     Gamma^k_{ij} = 1/2 g^{kl} (d_i g_{lj} + d_j g_{li} - d_l g_{ij}),
//
// with the partials of g taken exactly by Dual, so a metric written once as
// a callable on any scalar needs no derivative of its own. It is the
// N-dimensional form of what physics/relativity/geodesic.hpp did for four,
// which now calls these.
//
//   metric:  Vec<S, N> -> Matrix<S, N, N>, templated on the scalar S -- a
//            generic lambda or a functor, never std::function, so that S
//            can be Dual<T> and the derivatives come out of the arithmetic.
//   exp:     the geodesic flow, integrated by Gragg-Bulirsch-Stoer to the
//            precision of the scalar -- the same call is right on double and
//            on Real50.
//   log:     Newton's shooting method on exp, started from the chord.
//
// What a derived space is for: it is the same object as a closed-form one
// reached the long way round. Sphere<2> in (theta, phi), and the hyperbolic
// upper half-plane, have closed forms in the library; the derived exp and
// distance agree with them (tests/test_metric_chart.cpp), which is a
// symmetry that holds the closed forms, the derivation and the integrator
// at once, and a space with no closed form -- Kerr in Boyer-Lindquist
// coordinates, a metric someone measured -- is the same code.
//
// Limits, stated: coordinates are not points of a manifold but of one
// chart, so a geodesic that leaves the chart, or meets a coordinate
// singularity (theta = 0 for the sphere), gives NaN rather than a wrong
// number; and log is the shooting solution nearest the chord, unique only
// inside the injectivity radius.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/core/error.hpp>
#  include <spatium/algebra/dual.hpp>
#  include <spatium/algebra/linear_solve.hpp>
#  include <spatium/algebra/matrix.hpp>
#  include <spatium/algebra/ode.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <array>
#  include <cmath>
#  include <cstddef>
#  include <limits>
#  include <utility>
#endif

SPATIUM_EXPORT namespace spatium::spaces {

// g_{ij}(x) and its exact partials d_k g_{ij}(x) for all (k, i, j) at once:
// one Dual-seeded evaluation of the whole metric per coordinate axis (N
// metric evaluations, not N*N*N) -- the rank-2 generalisation of
// calculus.hpp's per-axis-seeded gradient().
template<Scalar T, std::size_t N>
struct MetricDerivatives {
    Matrix<T, N, N> g{};
    std::array<Matrix<T, N, N>, N> dg{};   // dg[k](i,j) = d g_{ij} / d x^k
};

template<Scalar T, std::size_t N, typename Metric>
MetricDerivatives<T, N> metric_derivatives(const Metric& metric, const Vec<T, N>& x) {
    MetricDerivatives<T, N> md{};
    for (std::size_t k = 0; k < N; ++k) {
        Vec<Dual<T>, N> dx;
        for (std::size_t j = 0; j < N; ++j)
            dx[j] = (j == k) ? Dual<T>::variable(x[j]) : Dual<T>::constant(x[j]);
        const Matrix<Dual<T>, N, N> gd = metric(dx);
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < N; ++j) {
                if (k == 0) md.g(i, j) = gd(i, j).value;
                md.dg[k](i, j) = gd(i, j).deriv;
            }
    }
    return md;
}

namespace metric_chart_detail {

template<Scalar T, std::size_t N>
std::array<Matrix<T, N, N>, N> connection(const MetricDerivatives<T, N>& md, const Matrix<T, N, N>& ginv) {
    std::array<Matrix<T, N, N>, N> Gamma{};
    for (std::size_t l = 0; l < N; ++l)
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < N; ++j) {
                T sum{0};
                for (std::size_t s = 0; s < N; ++s)
                    sum += ginv(l, s) * (md.dg[i](s, j) + md.dg[j](s, i) - md.dg[s](i, j));
                Gamma[l](i, j) = T{0.5} * sum;
            }
    return Gamma;
}

}  // namespace metric_chart_detail

// Christoffel symbols of the second kind, Gamma[l](i, j) = Gamma^l_{ij}.
// The inverse is a general N x N one, so a non-diagonal metric (Kerr's
// g_{t phi}) needs nothing special. Where g does not invert, the inverse is
// taken as zero and so are the symbols -- the behaviour physics/relativity
// has relied on, which does not refuse; `christoffel_checked` does.
template<Scalar T, std::size_t N, typename Metric>
std::array<Matrix<T, N, N>, N> christoffel(const Metric& metric, const Vec<T, N>& x) {
    const auto md = metric_derivatives(metric, x);
    const auto inv = invert(md.g);
    return metric_chart_detail::connection(md, inv ? *inv : Matrix<T, N, N>{});
}

// The same, refusing where the metric does not invert: a flow through
// symbols that silently read 0 is a straight line, and a wrong answer.
template<Scalar T, std::size_t N, typename Metric>
Result<std::array<Matrix<T, N, N>, N>> christoffel_checked(const Metric& metric, const Vec<T, N>& x) {
    const auto md = metric_derivatives(metric, x);
    const auto inv = invert(md.g);
    if (!inv) return std::unexpected(Error{ErrorCode::SingularMatrix, "the metric does not invert here"});
    return metric_chart_detail::connection(md, *inv);
}

// The state of a geodesic is its position and velocity in one vector,
// [0, N) = x, [N, 2N) = dx/dlambda, so the equation is the first-order
// system ode.hpp takes: x' = u, u'^l = -Gamma^l_{ij} u^i u^j. The Christoffel
// call is the argument's own type of scalar -- Dual for derivatives of a flow
// -- and a singular metric gives NaN, not zero.
template<Scalar T, std::size_t M, typename Metric>
    requires (M % 2 == 0)
Vec<T, M> geodesic_rhs(const Metric& metric, T /*lambda*/, const Vec<T, M>& state) {
    constexpr std::size_t N = M / 2;
    Vec<T, N> x, u;
    for (std::size_t i = 0; i < N; ++i) { x[i] = state[i]; u[i] = state[N + i]; }
    const auto G = christoffel_checked(metric, x);

    Vec<T, M> ds{};
    for (std::size_t i = 0; i < N; ++i) ds[i] = u[i];
    for (std::size_t l = 0; l < N; ++l) {
        T acc{0};
        if (G) {
            for (std::size_t i = 0; i < N; ++i)
                for (std::size_t j = 0; j < N; ++j) acc += (*G)[l](i, j) * u[i] * u[j];
        } else {
            acc = T(std::numeric_limits<double>::quiet_NaN());
        }
        ds[N + l] = -acc;
    }
    return ds;
}

// One RK4 step of the geodesic equation, for a caller that wants a fixed
// step (a renderer's march); `MetricChart::exp_map` below takes the
// adaptive one.
template<Scalar T, std::size_t M, typename Metric>
    requires (M % 2 == 0)
Vec<T, M> geodesic_step(const Metric& metric, const Vec<T, M>& state, T dlambda) {
    auto f = [&metric](T lam, const Vec<T, M>& y) { return geodesic_rhs(metric, lam, y); };
    return rk4_step(f, T{0}, state, dlambda);
}

// ── The space ──────────────────────────────────────────────────

template<Scalar T, std::size_t N, typename Metric>
class MetricChart {
public:
    using ScalarType = T;
    using PointType = Vec<T, N>;
    using TangentVector = Vec<T, N>;
    static constexpr std::size_t dimension = N;
    static constexpr bool is_complete = false;

    // `accuracy` is the integrator's tolerance and the shooting's, absolute
    // for coordinates of order one; 0 takes eps^(2/3) of the scalar -- 4e-11
    // for double, 2e-33 for Real50 -- far below the sqrt(eps) an acos-based
    // closed form keeps, so a derived and a closed form compare at theirs.
    explicit MetricChart(Metric metric, double accuracy = 0)
        : metric_(std::move(metric)),
          tol_(accuracy > 0 ? accuracy : std::pow(machine_epsilon<T>(), 2.0 / 3.0)) {}

    Matrix<T, N, N> metric_tensor(const PointType& p) const { return metric_(p); }

    ScalarType metric_at(const PointType& p, const TangentVector& u, const TangentVector& v) const {
        const Matrix<T, N, N> g = metric_(p);
        T s{0};
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < N; ++j) s += u[i] * g(i, j) * v[j];
        return s;
    }

    // The point the geodesic from p with velocity v reaches at parameter t.
    // NaN where the flow cannot be integrated: a singular metric, a
    // coordinate singularity, a geodesic that runs off the chart.
    PointType exp_map(const PointType& p, const TangentVector& v, ScalarType t) const {
        const auto r = flow(p, v, t);
        return r ? *r : nan_point();
    }

    // The tangent at p whose geodesic reaches q at parameter 1: Newton on
    // exp with a central-difference Jacobian, from the chord q - p, which
    // is right to first order. The Jacobian's step only sets the rate; the
    // solution is a fixed point of exp itself, so it carries exp's accuracy.
    TangentVector log_map(const PointType& p, const PointType& q) const {
        TangentVector v{q - p};
        const T one{1};
        auto residual = [&](const TangentVector& w) -> Result<TangentVector> {
            const auto e = flow(p, w, one);
            if (!e) return std::unexpected(e.error());
            return TangentVector{*e - q};
        };
        auto size = [](const TangentVector& r) {
            double m = 0;
            for (std::size_t i = 0; i < N; ++i) {
                const double a = std::abs(primal_double(r[i]));
                if (!(a <= m)) m = a;
            }
            return m;
        };
        const double stop = 8 * tol_;
        const double step = std::cbrt(machine_epsilon<T>());
        auto r = residual(v);

        // Newton with the Jacobian kept while it keeps working: the 2N
        // exps it costs are paid again only when a step fails to cut the
        // residual to a third, which, near the solution, it does not need.
        Matrix<T, N, N> Jinv{};
        bool fresh = false;
        for (int it = 0; it < 60 && r; ++it) {
            const double rs = size(*r);
            if (rs <= stop) return v;

            if (!fresh) {
                Matrix<T, N, N> J{};
                for (std::size_t j = 0; j < N; ++j) {
                    const T h = T(step * (1.0 + std::abs(primal_double(v[j]))));
                    TangentVector vp = v, vm = v;
                    vp[j] += h;
                    vm[j] -= h;
                    const auto rp = residual(vp), rm = residual(vm);
                    if (!rp || !rm) return nan_point();
                    for (std::size_t i = 0; i < N; ++i) J(i, j) = ((*rp)[i] - (*rm)[i]) / (T{2} * h);
                }
                const auto inv = invert(J);
                if (!inv) return nan_point();
                Jinv = *inv;
                fresh = true;
            }

            // The step, halved until the residual falls: a chord far from
            // the geodesic can overshoot where exp bends.
            bool moved = false;
            double scale = 1.0;
            for (int back = 0; back < 10; ++back, scale *= 0.5) {
                TangentVector next = v;
                for (std::size_t i = 0; i < N; ++i) {
                    T d{0};
                    for (std::size_t j = 0; j < N; ++j) d += Jinv(i, j) * (*r)[j];
                    next[i] = v[i] - d * T(scale);
                }
                const auto rn = residual(next);
                if (rn && size(*rn) < rs) {
                    fresh = fresh && size(*rn) <= rs / 3;
                    v = next;
                    r = rn;
                    moved = true;
                    break;
                }
            }
            if (!moved) {
                if (fresh) break;      // a new Jacobian and still no step: stuck
                fresh = false;
            }
        }
        // Stagnated above `stop` but close: the integrator's noise floor is
        // the answer. Further away it is a failure and a NaN says so.
        if (r && size(*r) <= std::sqrt(machine_epsilon<T>())) return v;
        return nan_point();
    }

private:
    Metric metric_;
    double tol_;

    static PointType nan_point() {
        PointType n;
        for (std::size_t i = 0; i < N; ++i) n[i] = T(std::numeric_limits<double>::quiet_NaN());
        return n;
    }

    Result<PointType> flow(const PointType& p, const TangentVector& v, T t) const {
        Vec<T, 2 * N> y0;
        for (std::size_t i = 0; i < N; ++i) { y0[i] = p[i]; y0[N + i] = v[i]; }
        auto f = [this](T lam, const Vec<T, 2 * N>& y) { return geodesic_rhs(metric_, lam, y); };
        const auto y = integrate_extrapolated(f, T{0}, y0, t, tol_);
        if (!y) return std::unexpected(y.error());
        PointType out;
        for (std::size_t i = 0; i < N; ++i) out[i] = (*y)[i];
        return out;
    }
};

// A space from its metric, a callable on any scalar:
//   auto h2 = metric_chart<double, 2>([](const auto& x) {
//       using S = std::remove_cvref_t<decltype(x[0])>;
//       Matrix<S, 2, 2> g{};  g(0, 0) = g(1, 1) = S{1} / (x[1] * x[1]);  return g; });
template<Scalar T, std::size_t N, typename Metric>
MetricChart<T, N, Metric> metric_chart(Metric metric, double accuracy = 0) {
    return MetricChart<T, N, Metric>(std::move(metric), accuracy);
}

} // namespace spatium::spaces
