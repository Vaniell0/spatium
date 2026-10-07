#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/algebra/dual.hpp>
#  include <spatium/algebra/matrix.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <spatium/spaces/metric_chart.hpp>
#  include <spatium/spaces/pullback.hpp>
#  include <array>
#  include <cmath>
#  include <concepts>
#  include <functional>
#  include <limits>
#  include <utility>
#endif

SPATIUM_EXPORT namespace spatium {

// ParametricSurface: f(u,v) → R³ becomes a full Surface.
// Automatically computes exp/log/distance/project/normal from the parameterization.
// Works with all geodesic algorithms.
//
// The function is a type, `F`, and by default the erased one a
// std::function is -- so `ParametricSurface<T>` means what it always did and
// every algorithm that takes one still does. Give it a callable generic over
// the scalar instead (`[](auto u, auto v) { ... }`, made by `make_parametric`)
// and the description stays lazy: it can be called on a `Dual`, so its
// partials are exact where the erased form takes finite differences (about
// 1e-8), and its geodesics -- exp and log -- are those of the induced metric,
// through `spaces::MetricChart`, where the erased form can only offer a step
// in the ambient space projected back. The erased form is an explicit
// boundary, `erased()`, for what takes only that: tessellate, the ray code,
// the scene DSL.

namespace parametric_detail {

// F can be called on a Dual<T>: a generic callable, not a
// std::function<Vec<T,3>(T,T)>.
template<class F, class T>
concept DualEvaluable = requires(const F& f, Dual<T> u, Dual<T> v) {
    { f(u, v) } -> std::convertible_to<Vec<Dual<T>, 3>>;
};

// The first fundamental form of f, g = J^T J, as a metric on the parameters:
// a callable on any scalar S (it evaluates f on Dual<S> for J), which is the
// shape spaces::MetricChart takes -- and, there, differentiated once more for
// the Christoffel symbols.
template<class F>
struct InducedMetric {
    F f;
    template<class S>
    Matrix<S, 2, 2> operator()(const Vec<S, 2>& x) const {
        // The pullback of R^3's metric through (u, v) -> f(u, v), by spaces/pullback.hpp.
        const auto as_map = [this](const auto& uv) { return f(uv[0], uv[1]); };
        return spaces::pullback_metric<2, 3>(as_map)(x);
    }
};

}  // namespace parametric_detail

template<Scalar T = double, class F = std::function<Vec<T, 3>(T, T)>>
class ParametricSurface {
public:
    using ScalarType = T;
    using PointType = Vec<T, 3>;
    using TangentVector = Vec<T, 3>;
    using VectorType = Vec<T, 3>;
    static constexpr std::size_t dimension = 2;
    static constexpr bool is_complete = false;  // not always

    using ParamFn = F;

    // True when F can be called on a Dual<T>: the partials are then exact and
    // the geodesics are the induced metric's.
    static constexpr bool exact = parametric_detail::DualEvaluable<F, T>;

    struct Domain {
        T u_min, u_max, v_min, v_max;
    };

    ParametricSurface(F fn, Domain domain, bool periodic_u = false, bool periodic_v = false)
        : fn_(std::move(fn)), domain_(domain),
          periodic_u_(periodic_u), periodic_v_(periodic_v) {}

    // A chart that knows its own projection in closed form says so, and
    // `project`/`normal` stop searching for what is already known.
    //
    // This is not an optimisation, it is a correctness fix. `find_params`
    // is a grid search plus refinement, and where a chart degenerates it
    // cannot converge in principle: at a sphere's pole every `u` names the
    // same point, so the parameter is invisible to the gradient while
    // still choosing the direction of any later step. Measured on a unit
    // sphere at 240x120, the iterative path was wrong on 12.9% of sampled
    // points by up to 20% of the radius, against 1.5e-4 for the exact
    // projection.
    //
    // The library's own principle applied to itself: where there is a
    // closed form, use it, and iterate only where there is not -- which is
    // the same argument `RenderLevel::Exact` is built on.
    ParametricSurface& with_closed_forms(
        std::function<PointType(const PointType&)> project_fn,
        std::function<TangentVector(const PointType&)> normal_fn = {}) {
        project_fn_ = std::move(project_fn);
        normal_fn_  = std::move(normal_fn);
        return *this;
    }

    // ── Evaluate ──────────────────────────────────────────────

    PointType operator()(T u, T v) const { return fn_(u, v); }

    // ── Surface concept methods ───────────────────────────────

    bool contains(const PointType& p) const {
        auto [u, v] = find_params(p);
        return (fn_(u, v) - p).norm() < epsilon<T>() * T{100};
    }

    // Euclidean distance in ambient R³ — a valid metric but NOT the geodesic
    // (intrinsic) distance along the surface. Acts as a lower bound for geodesic distance.
    // Only for the erased form: a surface that has its geodesics has no member
    // `distance`, so core/access.hpp derives it as |log|_g, the intrinsic one.
    T distance(const PointType& a, const PointType& b) const requires (!exact) {
        return (a - b).norm();
    }

    // The chord, for any surface, under its own name.
    T chord_distance(const PointType& a, const PointType& b) const {
        return (a - b).norm();
    }

    PointType project(const PointType& p) const {
        if (project_fn_) return project_fn_(p);
        auto [u, v] = find_params(p);
        return fn_(u, v);
    }

    TangentVector normal(const PointType& p) const {
        if (normal_fn_) return normal_fn_(p);
        auto [u, v] = find_params(p);
        return normal_at(u, v);
    }

    // The erased form: a retraction, not the exponential map -- a step in the
    // ambient space projected back, first order in t. exp(log) is not the
    // identity on it, and the matrix says so. Make the surface with
    // `make_parametric` for the geodesic.
    PointType exp_map(const PointType& p, const TangentVector& v, T t) const requires (!exact) {
        return project(PointType{p + v * t});
    }

    // The geodesic: the induced metric g = J^T J on (u, v), its Christoffel
    // symbols by Dual, the flow by spaces::MetricChart, and f evaluated at the
    // end. The tangent vector, in R^3, is put into the chart's basis by the
    // normal equations. Where the chart degenerates -- a sphere's pole in
    // latitude-longitude -- g does not invert and the answer is NaN.
    PointType exp_map(const PointType& p, const TangentVector& w, T t) const requires (exact) {
        const auto [u, v] = find_params(p);
        const auto coords = chart_coordinates(u, v, w);
        const auto flow = spaces::metric_chart<T, 2>(parametric_detail::InducedMetric<F>{fn_});
        const Vec<T, 2> end = flow.exp_map(Vec<T, 2>{u, v}, coords, t);
        return fn_(end[0], end[1]);
    }

    // Shooting in the parameters from p to q, and the result carried back into
    // R^3 as the combination of the partials. On a periodic chart q's
    // parameters are the representative nearest p's: u = 0.1 and u = 6.2 on a
    // 2 pi chart are 0.18 apart, not 6.1.
    TangentVector log_map(const PointType& p, const PointType& q) const requires (exact) {
        using std::floor;
        const auto [u, v] = find_params(p);
        auto [u2, v2] = find_params(q);
        if (periodic_u_) {
            const T range = domain_.u_max - domain_.u_min;
            u2 = u2 - range * floor(T((u2 - u) / range + T{0.5}));
        }
        if (periodic_v_) {
            const T range = domain_.v_max - domain_.v_min;
            v2 = v2 - range * floor(T((v2 - v) / range + T{0.5}));
        }
        const auto flow = spaces::metric_chart<T, 2>(parametric_detail::InducedMetric<F>{fn_});
        const Vec<T, 2> d = flow.log_map(Vec<T, 2>{u, v}, Vec<T, 2>{u2, v2});
        return TangentVector{du(u, v) * d[0] + dv(u, v) * d[1]};
    }

    TangentVector log_map(const PointType& p, const PointType& q) const requires (!exact) {
        // Project (q-p) onto tangent plane at p
        // A TangentVector, not `auto`: q - p is an expression template that
        // refers to its operands, and assigning to it below did not compile
        // -- for every scalar, found by the connectivity matrix.
        TangentVector diff{q - p};
        auto n = normal(p);
        auto n2 = n.dot(n);
        if (n2 > epsilon<T>())
            diff = PointType{diff - n * (diff.dot(n) / n2)};
        return diff;
    }

    ScalarType metric_at(const PointType&,
                         const TangentVector& u, const TangentVector& v) const {
        return u.dot(v);
    }

    const Domain& domain() const { return domain_; }
    bool periodic_u() const { return periodic_u_; }
    bool periodic_v() const { return periodic_v_; }
    PointType evaluate(T u, T v) const { return fn_(u, v); }

    // Exact analytic normal at a known (u,v) -- cheaper and more direct
    // than normal(p), which first has to recover (u,v) from a 3D point
    // via find_params()'s Newton search. Useful whenever the caller
    // already has the parameters on hand (tessellation, mesh generation,
    // distortion analysis below).
    TangentVector normal_at(T u, T v) const {
        auto fu = du(u, v);
        auto fv = dv(u, v);
        return fu.cross(fv).normalized();
    }

    // Local anisotropy of the parametrization at (u,v): ratio of the two
    // singular values of the Jacobian [fu fv], via the eigenvalues of the
    // first fundamental form (E=fu.fu, F=fu.fv, G=fv.fv). 1.0 means
    // locally isometric -- a UV-uniform mesh step maps to comparable R^3
    // edge lengths in every direction; large values mean it doesn't, e.g.
    // a thin-ring torus (fu, going around the big loop, grows much faster
    // than fv near the outer rim) or a cone near its apex (fu shrinks to
    // zero while fv doesn't). This is the general, measurable form of what
    // caused the Klein-bottle seam mesh twist (see
    // examples/primitives_demo.cpp's KleinBottle comment): the two
    // branches' fu/fv scale mismatched discontinuously right there.
    T parametrization_anisotropy(T u, T v) const {
        auto fu = du(u, v);
        auto fv = dv(u, v);
        T E = fu.dot(fu), Fuv = fu.dot(fv), G = fv.dot(fv);   // Fuv: F is the function's type
        T tr = E + G;
        T disc_sq = tr * tr - T{4} * (E * G - Fuv * Fuv);
        using std::sqrt;   // unqualified calls from here on: ADL finds Dual's and Real50's
        T disc = disc_sq > T{0} ? sqrt(disc_sq) : T{0};
        T lambda_min = (tr - disc) / T{2};
        T lambda_max = (tr + disc) / T{2};
        if (lambda_min < epsilon<T>()) return std::numeric_limits<T>::max();
        return sqrt(lambda_max / lambda_min);
    }

    // Area element sqrt(EG - F^2) of the first fundamental form at (u,v) --
    // how much a unit (du,dv) patch stretches into R^3 area here. The
    // per-point weight uniform-by-surface-area sampling needs (see
    // spaces/sample.hpp): sampling (u,v) uniformly instead would bunch
    // points wherever the parametrization compresses space, e.g. near a
    // thin-ring torus's inner rim.
    T area_element(T u, T v) const {
        auto fu = du(u, v);
        auto fv = dv(u, v);
        T E = fu.dot(fu), Fuv = fu.dot(fv), G = fv.dot(fv);
        T disc = E * G - Fuv * Fuv;
        using std::sqrt;
        return disc > T{0} ? sqrt(disc) : T{0};
    }

    // The same surface with its function stored as a std::function: the form
    // tessellate, the ray code and the scene DSL take. An explicit boundary,
    // not a conversion -- what is lost crossing it is the exact partials and
    // the geodesics.
    ParametricSurface<T> erased() const {
        if constexpr (std::same_as<F, std::function<PointType(T, T)>>) {
            return *this;
        } else {
            // Domain is a nested type, one per instantiation: its fields cross, not the type.
            ParametricSurface<T> out(std::function<PointType(T, T)>(fn_),
                                     typename ParametricSurface<T>::Domain{domain_.u_min, domain_.u_max, domain_.v_min, domain_.v_max},
                                     periodic_u_, periodic_v_);
            if (project_fn_) out.with_closed_forms(project_fn_, normal_fn_);
            return out;
        }
    }

private:
    F fn_;
    std::function<PointType(const PointType&)> project_fn_;
    std::function<TangentVector(const PointType&)> normal_fn_;
    Domain domain_;
    bool periodic_u_, periodic_v_;

    // Partial derivatives: exact, through a Dual, when F can be called on one;
    // central differences otherwise.
    PointType du(T u, T v) const {
        if constexpr (exact) {
            const Vec<Dual<T>, 3> r = fn_(Dual<T>::variable(u), Dual<T>::constant(v));
            return PointType{r[0].deriv, r[1].deriv, r[2].deriv};
        } else {
            T h = (domain_.u_max - domain_.u_min) * T{1e-6};
            return (fn_(u + h, v) - fn_(u - h, v)) / (T{2} * h);
        }
    }

    PointType dv(T u, T v) const {
        if constexpr (exact) {
            const Vec<Dual<T>, 3> r = fn_(Dual<T>::constant(u), Dual<T>::variable(v));
            return PointType{r[0].deriv, r[1].deriv, r[2].deriv};
        } else {
            T h = (domain_.v_max - domain_.v_min) * T{1e-6};
            return (fn_(u, v + h) - fn_(u, v - h)) / (T{2} * h);
        }
    }

    // A tangent vector of R^3 at f(u, v) as coordinates (u', v') in the
    // chart's basis: the least-squares solution of u' f_u + v' f_v = w, by the
    // normal equations. NaN where f_u and f_v are parallel -- a pole.
    Vec<T, 2> chart_coordinates(T u, T v, const TangentVector& w) const {
        const PointType a = du(u, v), b = dv(u, v);
        const T E = a.dot(a), Fuv = a.dot(b), G = b.dot(b);
        const T det = E * G - Fuv * Fuv;
        const T r0 = a.dot(w), r1 = b.dot(w);
        if (!(det > T{0})) {
            const T nan = T(std::numeric_limits<double>::quiet_NaN());
            return Vec<T, 2>{nan, nan};
        }
        return Vec<T, 2>{T((G * r0 - Fuv * r1) / det), T((E * r1 - Fuv * r0) / det)};
    }

    // Find closest UV parameters for a 3D point (Newton-like search)
    // Nearest (u, v) to a point, by grid search then refinement.
    //
    // Rewritten 2026-09-22 after it was measured wrong on the easiest
    // surface there is. The previous version took an 8x8 grid and five
    // Gauss-Newton steps, broke out of the loop when the 2x2 system went
    // singular, and returned whatever it had. On a unit sphere that is
    // 12.9% of sampled points wrong by up to 20% of the radius, because
    // `fu` and `fv` degenerate at the poles: the determinant collapses,
    // the loop exits, and the grid answer survives -- accurate to the grid
    // spacing, 2*pi/8, about 45 degrees.
    //
    // Four changes, each addressing something measured rather than
    // suspected:
    //
    //   * **Keep the best point seen, not the last one.** A refinement
    //     that wanders can now never return worse than the grid it
    //     started from. This alone bounds the old failure.
    //   * **Degeneracy is a direction, not an exit.** Where the
    //     parametrization collapses -- a pole, where every `u` names the
    //     same point -- Gauss-Newton has nothing to say, but the residual
    //     can still be reduced along the gradient. Stepping there makes
    //     progress where breaking made none.
    //   * **Backtracking.** A step that increases the residual is halved
    //     rather than taken, which stops the overshoot-then-clamp
    //     oscillation the domain clamp used to produce.
    //   * **Iterate until it stops improving**, not a fixed five. The
    //     common case still converges in a handful; the hard case is no
    //     longer cut off mid-descent.
    //
    // Still not a `Result<T>`, and that is the remaining honest gap: this
    // returns parameters whether or not it converged, so a caller cannot
    // tell a projection from a best effort. Everything fallible in this
    // library returns `Result<T>`; that this does not is the same class of
    // defect the polynomial solvers were fixed for, and it wants its own
    // change rather than riding in on this one.
    std::pair<T, T> find_params(const PointType& target) const {
        // A target that is not finite has no nearest parameters. Every
        // comparison with NaN is false, so the search below kept its first
        // grid point and handed back a finite surface point for NaN or an
        // infinite step -- an answer where there is none, found by the
        // connectivity matrix's infinity probe. Non-finite in, non-finite out.
        if (!std::isfinite(primal_double(target.norm_squared()))) {
            const T nan = T(std::numeric_limits<double>::quiet_NaN());
            return {nan, nan};
        }
        // A finer grid than before: 17x17 samples rather than 9x9, which
        // costs nothing measurable and starts the refinement inside a
        // basin rather than up to 45 degrees outside one.
        constexpr int GRID = 16;
        T best_u = domain_.u_min, best_v = domain_.v_min;
        T best_dist = std::numeric_limits<T>::max();

        const T du_step = (domain_.u_max - domain_.u_min) / GRID;
        const T dv_step = (domain_.v_max - domain_.v_min) / GRID;

        for (int j = 0; j <= GRID; ++j) {
            T v = domain_.v_min + static_cast<T>(j) * dv_step;
            for (int i = 0; i <= GRID; ++i) {
                T u = domain_.u_min + static_cast<T>(i) * du_step;
                T d = (fn_(u, v) - target).norm_squared();
                if (d < best_dist) { best_dist = d; best_u = u; best_v = v; }
            }
        }

        // Bring a candidate back into the domain, wrapping where the
        // surface is periodic and clamping where it is not.
        const auto settle = [&](T& u, T& v) {
            using std::fmod;   // not std::fmod(...): qualified, it blocks ADL for Dual and Real50
            if (periodic_u_) {
                const T range = domain_.u_max - domain_.u_min;
                u = domain_.u_min + fmod(T(u - domain_.u_min), range);
                if (u < domain_.u_min) u += range;
            } else {
                u = std::clamp(u, domain_.u_min, domain_.u_max);
            }
            if (periodic_v_) {
                const T range = domain_.v_max - domain_.v_min;
                v = domain_.v_min + fmod(T(v - domain_.v_min), range);
                if (v < domain_.v_min) v += range;
            } else {
                v = std::clamp(v, domain_.v_max < domain_.v_min ? domain_.v_max : domain_.v_min,
                               domain_.v_max);
            }
        };

        T cur_u = best_u, cur_v = best_v;
        for (int iter = 0; iter < 32; ++iter) {
            const auto p = fn_(cur_u, cur_v);
            const auto fu = this->du(cur_u, cur_v);
            const auto fv = this->dv(cur_u, cur_v);
            const auto r = target - p;

            const T a11 = fu.dot(fu), a12 = fu.dot(fv);
            const T a22 = fv.dot(fv);
            const T b1 = fu.dot(r),  b2 = fv.dot(r);
            const T det = a11 * a22 - a12 * a12;

            T step_u, step_v;
            // Scale the singularity test against the matrix itself: an
            // absolute epsilon calls a small surface degenerate
            // everywhere and a large one degenerate nowhere.
            using std::abs;
            if (abs(det) > epsilon<T>() * std::max(T{1}, T(a11 * a22))) {
                step_u = (a22 * b1 - a12 * b2) / det;
                step_v = (a11 * b2 - a12 * b1) / det;
            } else {
                // Rank-deficient: descend the residual instead of solving.
                const T scale = std::max(T(a11 + a22), epsilon<T>());
                step_u = b1 / scale;
                step_v = b2 / scale;
            }

            // Backtrack rather than accept a step that makes it worse.
            bool improved = false;
            T t = T{1};
            for (int back = 0; back < 8; ++back) {
                T try_u = cur_u + t * step_u, try_v = cur_v + t * step_v;
                settle(try_u, try_v);
                const T d = (fn_(try_u, try_v) - target).norm_squared();
                if (d < best_dist) {
                    best_dist = d; best_u = try_u; best_v = try_v;
                    cur_u = try_u; cur_v = try_v;
                    improved = true;
                    break;
                }
                t *= T{0.5};
            }
            if (!improved) break;   // no direction left that helps
        }

        // A derivative-free pass, because at a degenerate point no
        // derivative method can succeed, and this one does not need to.
        //
        // Measured before it was written, and the mechanism is exact
        // rather than suspected: after the refinement above, every
        // remaining failure on a sphere sat within 5.5 degrees of a pole,
        // 961 of 961. At a pole `fu` vanishes, so `b1 = fu·r` is zero and
        // `u` never updates -- while `u` there still chooses which
        // meridian a later step in `v` travels along. The parameter is
        // invisible to the gradient and steers the descent at the same
        // time, which Gauss-Newton cannot recover from by construction.
        //
        // Nor can the initial grid: at a degenerate point every `u` gives
        // the same position, so all of them tie on distance and `best_u`
        // is whichever the comparison happened to see first. A window
        // around an arbitrary `u` searches the wrong neighbourhood however
        // far it shrinks.
        //
        // So the meridian is resolved first, by sweeping the whole `u`
        // range at a `v` nudged off the degeneracy -- one step away from a
        // pole the dependence on `u` is real again and a sweep sees it --
        // and only then does a shrinking local search take over. That
        // search samples `u` rather than differentiating it, so it has no
        // blind spot where the parametrization collapses.
        {
            const T nudge = dv_step * T{0.25};
            for (const T v_probe : {T(best_v - nudge), T(best_v + nudge)}) {   // T(...): an initializer_list of Boost expressions does not deduce
                for (int i = 0; i <= GRID; ++i) {
                    T u = domain_.u_min + static_cast<T>(i) * du_step;
                    T v = v_probe;
                    settle(u, v);
                    const T d = (fn_(u, v) - target).norm_squared();
                    if (d < best_dist) { best_dist = d; best_u = u; best_v = v; }
                }
            }

            T wu = du_step, wv = dv_step;
            for (int round = 0; round < 4; ++round) {
                const T cu = best_u, cv = best_v;
                for (int j = -2; j <= 2; ++j) {
                    for (int i = -2; i <= 2; ++i) {
                        if (!i && !j) continue;
                        T u = cu + static_cast<T>(i) * wu * T{0.5};
                        T v = cv + static_cast<T>(j) * wv * T{0.5};
                        settle(u, v);
                        const T d = (fn_(u, v) - target).norm_squared();
                        if (d < best_dist) { best_dist = d; best_u = u; best_v = v; }
                    }
                }
                wu /= T{3}; wv /= T{3};
            }
        }

        return {best_u, best_v};
    }
};

// A surface whose function stays a type, generic over the scalar: write it
// `[](auto u, auto v) { using std::cos; using std::sin; return Vec<decltype(u), 3>{...}; }`
// and give the scalar T it lives over. Its partials are exact and exp and log
// are its geodesics (see ParametricSurface above); `.erased()` is the way back.
template<Scalar T = double, class F>
ParametricSurface<T, F> make_parametric(F fn, typename ParametricSurface<T, F>::Domain domain,
                                        bool periodic_u = false, bool periodic_v = false) {
    return ParametricSurface<T, F>(std::move(fn), domain, periodic_u, periodic_v);
}

// ── Closure: does the surface end anywhere? ───────────────────

// A parameter direction closes up either because the map is periodic
// there -- the torus, in both u and v -- or because both of its edge
// curves collapse to a single point. The second case is a pole: a place
// the surface passes through, not one it ends at, which is how a sphere
// closes in v while its v domain [0, pi] stays a plain interval. A
// direction that does neither has a real edge, like the open ends of a
// cylinder's tube or the two rims of a band cut out of a torus.
//
// Both directions have to close for the surface to: half a sphere is
// periodic in nothing and poles in nothing, and its u edges are a rim.
//
// Periodicity is a flag the caller sets and we trust -- tessellation
// already trusts it. Poles are not recorded anywhere, so this looks
// instead of asking. Sampling an edge curve can only be fooled by a
// parametrization that returns the same point at all eight samples and
// wanders off between them, which is not a thing anyone writes.
template<Scalar T>
bool is_closed(const ParametricSurface<T>& surf) {
    constexpr std::size_t samples = 8;
    const auto dom = surf.domain();
    const T du = dom.u_max - dom.u_min, dv = dom.v_max - dom.v_min;

    // Tolerance relative to how big the surface is, and measured as a
    // coordinate spread so it does not depend on where the origin sits.
    Vec<T, 3> lo = surf.evaluate(dom.u_min, dom.v_min), hi = lo;
    for (std::size_t i = 0; i <= 4; ++i)
        for (std::size_t j = 0; j <= 4; ++j) {
            auto p = surf.evaluate(dom.u_min + du * static_cast<T>(i) / T{4},
                                   dom.v_min + dv * static_cast<T>(j) / T{4});
            for (std::size_t k = 0; k < 3; ++k) {
                using std::min, std::max;
                lo[k] = min(lo[k], p[k]);
                hi[k] = max(hi[k], p[k]);
            }
        }
    T extent{0};
    for (std::size_t k = 0; k < 3; ++k) { using std::max; extent = max(extent, hi[k] - lo[k]); }
    const T tol = extent * epsilon<T>() * T{64} + epsilon<T>();

    // Is the edge at the fixed end of one axis a single point? `vary_u`
    // says which axis runs along the edge: the v = v_min edge is a curve
    // in u, and vice versa.
    auto edge_is_point = [&](bool vary_u, T fixed) {
        auto at = [&](std::size_t k) {
            T s = static_cast<T>(k) / static_cast<T>(samples - 1);
            return vary_u ? surf.evaluate(dom.u_min + du * s, fixed)
                          : surf.evaluate(fixed, dom.v_min + dv * s);
        };
        auto p0 = at(0);
        for (std::size_t k = 1; k < samples; ++k)
            if (Vec<T, 3>{at(k) - p0}.norm() > tol) return false;
        return true;
    };

    bool u_closed = surf.periodic_u() ||
                    (edge_is_point(false, dom.u_min) && edge_is_point(false, dom.u_max));
    bool v_closed = surf.periodic_v() ||
                    (edge_is_point(true, dom.v_min) && edge_is_point(true, dom.v_max));
    return u_closed && v_closed;
}

// ── Convenience factories ─────────────────────────────────────

template<Scalar T = double>
ParametricSurface<T> make_torus(T major_r = T{2}, T minor_r = T{1}) {
    using std::acos;
    return ParametricSurface<T>(
        [=](T u, T v) -> Vec<T, 3> {
            using std::cos; using std::sin;
            return {
                (major_r + minor_r * cos(v)) * cos(u),
                (major_r + minor_r * cos(v)) * sin(u),
                minor_r * sin(v)
            };
        },
        {T{0}, T{2} * acos(T{-1}), T{0}, T{2} * acos(T{-1})},
        true, true  // periodic in both u and v
    );
}

template<Scalar T = double>
ParametricSurface<T> make_cylinder(T radius = T{1}, T height = T{2}) {
    using std::acos;
    return ParametricSurface<T>(
        [=](T u, T v) -> Vec<T, 3> {
            using std::cos; using std::sin;
            return {radius * cos(u), radius * sin(u), v};
        },
        {T{0}, T{2} * acos(T{-1}), T{0}, height},
        true, false
    );
}

template<Scalar T = double>
ParametricSurface<T> make_cone(T radius = T{1}, T height = T{2}) {
    using std::acos;
    return ParametricSurface<T>(
        [=](T u, T v) -> Vec<T, 3> {
            using std::cos; using std::sin;
            T r = radius * (T{1} - v / height);
            return {r * cos(u), r * sin(u), v};
        },
        {T{0}, T{2} * acos(T{-1}), T{0}, height},
        true, false
    );
}

template<Scalar T = double>
ParametricSurface<T> make_mobius(T radius = T{2}, T half_width = T{0.5}) {
    using std::acos;
    return ParametricSurface<T>(
        [=](T u, T v) -> Vec<T, 3> {
            using std::cos; using std::sin;
            T half_u = u / T{2};
            return {
                (radius + v * cos(half_u)) * cos(u),
                (radius + v * cos(half_u)) * sin(u),
                v * sin(half_u)
            };
        },
        {T{0}, T{2} * acos(T{-1}), -half_width, half_width},
        false, false  // Möbius is not globally periodic in the simple sense
    );
}

// ── DSL factories ────────────────────────────────────────────

// Domain builder: periodic(u_min, u_max, v_min, v_max) marks both axes periodic
template<Scalar T = double>
struct PeriodicDomain {
    typename ParametricSurface<T>::Domain domain;
};

template<Scalar T = double>
PeriodicDomain<T> periodic(T u_min, T u_max, T v_min, T v_max) {
    return {{u_min, u_max, v_min, v_max}};
}

// parametric(fn, domain) — non-periodic
template<Scalar T = double, typename F>
    requires std::invocable<F, T, T>
ParametricSurface<T> parametric(F&& fn, typename ParametricSurface<T>::Domain domain,
                                 bool pu = false, bool pv = false) {
    return ParametricSurface<T>(std::forward<F>(fn), domain, pu, pv);
}

// parametric(fn, periodic(...)) — both axes periodic
template<Scalar T = double, typename F>
    requires std::invocable<F, T, T>
ParametricSurface<T> parametric(F&& fn, PeriodicDomain<T> pd) {
    return ParametricSurface<T>(std::forward<F>(fn), pd.domain, true, true);
}

} // namespace spatium
