#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <spatium/mesh/mesh.hpp>
#  include <array>
#  include <cmath>
#  include <functional>
#endif

SPATIUM_EXPORT namespace spatium {

// ParametricSurface: f(u,v) → R³ becomes a full Surface.
// Automatically computes exp/log/distance/project/normal from the parameterization.
// Works with all geodesic algorithms.

template<Scalar T = double>
class ParametricSurface {
public:
    using ScalarType = T;
    using PointType = Vec<T, 3>;
    using TangentVector = Vec<T, 3>;
    using VectorType = Vec<T, 3>;
    static constexpr std::size_t dimension = 2;
    static constexpr bool is_complete = false;  // not always

    using ParamFn = std::function<PointType(T, T)>;

    struct Domain {
        T u_min, u_max, v_min, v_max;
    };

    ParametricSurface(ParamFn fn, Domain domain, bool periodic_u = false, bool periodic_v = false)
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
    T distance(const PointType& a, const PointType& b) const {
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

    PointType exp_map(const PointType& p, const TangentVector& v, T t) const {
        // First-order: move in ambient space, project back
        return project(PointType{p + v * t});
    }

    TangentVector log_map(const PointType& p, const PointType& q) const {
        // Project (q-p) onto tangent plane at p
        auto diff = q - p;
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
        T E = fu.dot(fu), F = fu.dot(fv), G = fv.dot(fv);
        T tr = E + G;
        T disc_sq = tr * tr - T{4} * (E * G - F * F);
        T disc = disc_sq > T{0} ? std::sqrt(disc_sq) : T{0};
        T lambda_min = (tr - disc) / T{2};
        T lambda_max = (tr + disc) / T{2};
        if (lambda_min < epsilon<T>()) return std::numeric_limits<T>::max();
        using std::sqrt;
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
        T E = fu.dot(fu), F = fu.dot(fv), G = fv.dot(fv);
        T disc = E * G - F * F;
        using std::sqrt;
        return disc > T{0} ? sqrt(disc) : T{0};
    }

private:
    ParamFn fn_;
    std::function<PointType(const PointType&)> project_fn_;
    std::function<TangentVector(const PointType&)> normal_fn_;
    Domain domain_;
    bool periodic_u_, periodic_v_;

    // Partial derivatives (finite differences)
    PointType du(T u, T v) const {
        T h = (domain_.u_max - domain_.u_min) * T{1e-6};
        return (fn_(u + h, v) - fn_(u - h, v)) / (T{2} * h);
    }

    PointType dv(T u, T v) const {
        T h = (domain_.v_max - domain_.v_min) * T{1e-6};
        return (fn_(u, v + h) - fn_(u, v - h)) / (T{2} * h);
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
            if (periodic_u_) {
                const T range = domain_.u_max - domain_.u_min;
                u = domain_.u_min + std::fmod(u - domain_.u_min, range);
                if (u < domain_.u_min) u += range;
            } else {
                u = std::clamp(u, domain_.u_min, domain_.u_max);
            }
            if (periodic_v_) {
                const T range = domain_.v_max - domain_.v_min;
                v = domain_.v_min + std::fmod(v - domain_.v_min, range);
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
            if (std::abs(det) > epsilon<T>() * std::max(T{1}, a11 * a22)) {
                step_u = (a22 * b1 - a12 * b2) / det;
                step_v = (a11 * b2 - a12 * b1) / det;
            } else {
                // Rank-deficient: descend the residual instead of solving.
                const T scale = std::max(a11 + a22, epsilon<T>());
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
            for (const T v_probe : {best_v - nudge, best_v + nudge}) {
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

// ── Tessellation (free function — class must be complete for Mesh<Surface>) ──

template<Scalar T>
mesh::Mesh<ParametricSurface<T>> tessellate(const ParametricSurface<T>& surf,
                                             std::size_t nu, std::size_t nv) {
    mesh::Mesh<ParametricSurface<T>> m;
    auto dom = surf.domain();
    bool pu = surf.periodic_u(), pv = surf.periodic_v();

    T du = (dom.u_max - dom.u_min) / static_cast<T>(nu);
    T dv = (dom.v_max - dom.v_min) / static_cast<T>(nv);

    std::size_t nu_verts = pu ? nu : nu + 1;
    std::size_t nv_verts = pv ? nv : nv + 1;

    m.vertices.reserve(nu_verts * nv_verts);
    for (std::size_t j = 0; j < nv_verts; ++j) {
        T v = dom.v_min + static_cast<T>(j) * dv;
        for (std::size_t i = 0; i < nu_verts; ++i) {
            T u = dom.u_min + static_cast<T>(i) * du;
            m.vertices.push_back(surf.evaluate(u, v));
        }
    }

    m.faces.reserve(nu * nv * 2);
    for (std::size_t j = 0; j < nv; ++j) {
        for (std::size_t i = 0; i < nu; ++i) {
            auto idx = [&](std::size_t ii, std::size_t jj) -> uint32_t {
                if (pu) ii %= nu;
                if (pv) jj %= nv;
                return static_cast<uint32_t>(jj * nu_verts + ii);
            };
            uint32_t a = idx(i, j), b = idx(i + 1, j);
            uint32_t c = idx(i + 1, j + 1), d = idx(i, j + 1);
            m.faces.push_back({a, b, c});
            m.faces.push_back({a, c, d});
        }
    }
    return m;
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
    return ParametricSurface<T>(
        [=](T u, T v) -> Vec<T, 3> {
            return {
                (major_r + minor_r * std::cos(v)) * std::cos(u),
                (major_r + minor_r * std::cos(v)) * std::sin(u),
                minor_r * std::sin(v)
            };
        },
        {T{0}, T{2} * std::acos(T{-1}), T{0}, T{2} * std::acos(T{-1})},
        true, true  // periodic in both u and v
    );
}

template<Scalar T = double>
ParametricSurface<T> make_cylinder(T radius = T{1}, T height = T{2}) {
    return ParametricSurface<T>(
        [=](T u, T v) -> Vec<T, 3> {
            return {radius * std::cos(u), radius * std::sin(u), v};
        },
        {T{0}, T{2} * std::acos(T{-1}), T{0}, height},
        true, false
    );
}

template<Scalar T = double>
ParametricSurface<T> make_cone(T radius = T{1}, T height = T{2}) {
    return ParametricSurface<T>(
        [=](T u, T v) -> Vec<T, 3> {
            T r = radius * (T{1} - v / height);
            return {r * std::cos(u), r * std::sin(u), v};
        },
        {T{0}, T{2} * std::acos(T{-1}), T{0}, height},
        true, false
    );
}

template<Scalar T = double>
ParametricSurface<T> make_mobius(T radius = T{2}, T half_width = T{0.5}) {
    return ParametricSurface<T>(
        [=](T u, T v) -> Vec<T, 3> {
            T half_u = u / T{2};
            return {
                (radius + v * std::cos(half_u)) * std::cos(u),
                (radius + v * std::cos(half_u)) * std::sin(u),
                v * std::sin(half_u)
            };
        },
        {T{0}, T{2} * std::acos(T{-1}), -half_width, half_width},
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
