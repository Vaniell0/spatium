#pragma once

// Ray intersection with a ParametricSurface via Newton UV.
//
// Solves S(u,v) = o + t·d for (u, v, t) — three equations, three unknowns.
// Seeded from a coarse UV grid, with per-seed rejection when the cell
// sample is too far from the ray line. Domain handling uses wrap for
// periodic axes, reject-on-exit for non-periodic ones.
//
// Complementary to ray_quadric / ray_torus: those close analytical roots
// in polynomial form (quadratic / quartic) for specific surface families.
// ray_parametric handles any f : (u, v) → R³ at the cost of iteration.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/geometry/line.hpp>
#  include <spatium/spaces/parametric.hpp>
#  include <algorithm>
#  include <cmath>
#  include <cstddef>
#  include <limits>
#  include <optional>
#  include <utility>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::geometry {

template<Scalar T>
struct RayParametricHit {
    T t;
    T u, v;
    Vec<T, 3> point;
    Vec<T, 3> normal;
};

template<Scalar T = double>
struct RayParametricConfig {
    std::size_t grid_u = 12;
    std::size_t grid_v = 12;
    int max_iter = 24;
    T tol = T{1e-8};        // ||S(u,v) - ray(t)|| threshold
    T dedup_t = T{1e-3};    // merge near-duplicate hits in t
    T cell_slack = T{2.5};  // cell_radius multiplier for seed rejection
};

// Estimates how far a Newton seed may need to travel in ambient space to
// reach its true root — depends only on the surface and grid resolution,
// never on any particular ray. ray_parametric()/ray_parametric_first()
// compute this internally by default when no hint is given, but a caller
// tracing many rays against the same surface (a per-pixel renderer) should
// compute it ONCE up front and pass it back in via cell_radius_hint:
// recomputing it per ray costs 25 extra surface evaluations on every
// single pixel for no reason — measured as part of a real ~1000-1700x
// per-ray slowdown this was contributing to in
// examples/parametric_analytical_demo.cpp (2026-08-25), fixed here.
template<Scalar T>
T estimate_cell_radius(const ParametricSurface<T>& surf, const RayParametricConfig<T>& cfg = {}) {
    const auto dom = surf.domain();
    const T u_range = dom.u_max - dom.u_min;
    const T v_range = dom.v_max - dom.v_min;
    const T du_step = u_range / static_cast<T>(cfg.grid_u);
    const T dv_step = v_range / static_cast<T>(cfg.grid_v);

    // Sample max chord between adjacent grid samples — bounds how far the
    // true hit can be from a seeded surface point.
    T cell_radius = T{0};
    for (int j = 0; j <= 4; ++j) {
        for (int i = 0; i <= 4; ++i) {
            T u = dom.u_min + u_range * static_cast<T>(i) / T{4};
            T v = dom.v_min + v_range * static_cast<T>(j) / T{4};
            Vec<T, 3> p0 = surf(u, v);
            T uu = std::min(u + du_step, dom.u_max);
            T vv = std::min(v + dv_step, dom.v_max);
            Vec<T, 3> p1 = surf(uu, v);
            Vec<T, 3> p2 = surf(u, vv);
            T c1 = Vec<T, 3>{p1 - p0}.norm();
            T c2 = Vec<T, 3>{p2 - p0}.norm();
            if (c1 > cell_radius) cell_radius = c1;
            if (c2 > cell_radius) cell_radius = c2;
        }
    }
    return cell_radius * cfg.cell_slack;
}

namespace detail {

// Everything a Newton solve needs that's shared across every seed tried
// against one (ray, surface) pair — built once per ray_parametric()/
// ray_parametric_first() call, not per seed.
template<Scalar T>
struct ParametricRayContext {
    const ParametricSurface<T>* surf;
    const RayParametricConfig<T>* cfg;
    Vec<T, 3> o, d, neg_d;
    typename ParametricSurface<T>::Domain dom;
    T u_range, v_range, du_step, dv_step, h_u, h_v, cell_radius;
    bool pu, pv;

    void partials(T u, T v, Vec<T, 3>& p, Vec<T, 3>& su, Vec<T, 3>& sv) const {
        p  = (*surf)(u, v);
        su = Vec<T, 3>{((*surf)(u + h_u, v) - (*surf)(u - h_u, v)) / (T{2} * h_u)};
        sv = Vec<T, 3>{((*surf)(u, v + h_v) - (*surf)(u, v - h_v)) / (T{2} * h_v)};
    }

    // Keeps (u, v) in-domain during Newton: wrap for periodic axes,
    // clamp for non-periodic.
    void clamp_or_wrap(T& u, T& v) const {
        if (pu) {
            u = dom.u_min + std::fmod(u - dom.u_min, u_range);
            if (u < dom.u_min) u += u_range;
        } else {
            u = std::clamp(u, dom.u_min, dom.u_max);
        }
        if (pv) {
            v = dom.v_min + std::fmod(v - dom.v_min, v_range);
            if (v < dom.v_min) v += v_range;
        } else {
            v = std::clamp(v, dom.v_min, dom.v_max);
        }
    }
};

template<Scalar T>
ParametricRayContext<T> make_context(const Ray<3, T>& r, const ParametricSurface<T>& surf,
                                      const RayParametricConfig<T>& cfg, T cell_radius) {
    const auto dom = surf.domain();
    ParametricRayContext<T> ctx{
        .surf = &surf, .cfg = &cfg,
        .o = r.origin, .d = r.direction, .neg_d = Vec<T, 3>{-r.direction},
        .dom = dom,
        .u_range = dom.u_max - dom.u_min, .v_range = dom.v_max - dom.v_min,
        .du_step = (dom.u_max - dom.u_min) / static_cast<T>(cfg.grid_u),
        .dv_step = (dom.v_max - dom.v_min) / static_cast<T>(cfg.grid_v),
        .h_u = (dom.u_max - dom.u_min) * T{1e-6}, .h_v = (dom.v_max - dom.v_min) * T{1e-6},
        .cell_radius = cell_radius,
        .pu = surf.periodic_u(), .pv = surf.periodic_v(),
    };
    return ctx;
}

// One Newton solve from a single (u0, v0) seed. Shared by the exhaustive
// grid scan in ray_parametric() and the single warm-start attempt
// ray_parametric_first()'s uv_hint uses.
template<Scalar T>
std::optional<RayParametricHit<T>> solve_from_seed(const ParametricRayContext<T>& ctx, T u0, T v0) {
    Vec<T, 3> p0 = (*ctx.surf)(u0, v0);
    Vec<T, 3> diff{p0 - ctx.o};
    T t0 = diff.dot(ctx.d);
    if (t0 < T{0}) t0 = T{0};

    Vec<T, 3> perp{diff - ctx.d * t0};
    if (perp.norm() > ctx.cell_radius) return std::nullopt;

    T u = u0, v = v0, t = t0;
    Vec<T, 3> p, su, sv;
    bool converged = false;
    T prev_f_norm = std::numeric_limits<T>::max();
    for (int it = 0; it < ctx.cfg->max_iter; ++it) {
        ctx.partials(u, v, p, su, sv);
        Vec<T, 3> F{p - (ctx.o + ctx.d * t)};
        T f_norm = F.norm();
        if (f_norm < ctx.cfg->tol) { converged = true; break; }

        // Stall/divergence guard: a converging Newton solve shrinks its
        // residual fast (quadratically once inside the root's basin), so
        // a seed whose residual isn't dropping by at least 30% a step
        // past the first correction is never going to reach cfg.tol --
        // bail out now instead of burning the rest of max_iter. This is
        // the fix for a real, measured cost: seeds near a silhouette/
        // grazing region pass the coarse cell_radius cull (above) but
        // never actually converge, and without this guard every one of
        // them ran the full 24-iteration budget -- dominating a Klein-
        // bottle scanline's miss pixels at ~568,000 ns/ray vs. ~1,900
        // ns/ray for a genuine hit (2026-08-25). One-iteration grace
        // (checked from it=1, not it=0) covers the usual "coarse first
        // step, then quadratic" convergence pattern; a threshold/grace
        // sweep against a real render (0.5-0.9 x 0-3 iterations) found
        // this the fastest setting that never dropped a real hit.
        if (it >= 1 && f_norm > prev_f_norm * T{0.7}) return std::nullopt;
        prev_f_norm = f_norm;

        // J = [su | sv | -d] ∈ R^{3x3}. Solve J·δ = -F via Cramer's rule.
        Vec<T, 3> sv_x_nd = sv.cross(ctx.neg_d);
        T det = su.dot(sv_x_nd);
        if (std::abs(det) < epsilon<T>()) return std::nullopt;
        Vec<T, 3> rhs = -F;
        T delta_u = rhs.dot(sv_x_nd) / det;
        T delta_v = su.dot(Vec<T, 3>{rhs.cross(ctx.neg_d)}) / det;
        T delta_t = su.dot(Vec<T, 3>{sv.cross(rhs)})        / det;

        // Clamp step to avoid blowing out of the basin of attraction.
        T step = std::sqrt(delta_u * delta_u + delta_v * delta_v);
        T max_step = std::max(ctx.du_step, ctx.dv_step) * T{2};
        if (step > max_step) {
            T s = max_step / step;
            delta_u *= s; delta_v *= s; delta_t *= s;
        }

        u += delta_u;
        v += delta_v;
        t += delta_t;

        ctx.clamp_or_wrap(u, v);
    }

    if (!converged) {
        ctx.partials(u, v, p, su, sv);
        Vec<T, 3> F{p - (ctx.o + ctx.d * t)};
        if (F.norm() > ctx.cfg->tol * T{1e3}) return std::nullopt;
    }
    if (t < T{0}) return std::nullopt;

    // Domain final check with small slack — boundary hits on
    // non-periodic axes remain valid (e.g. Möbius v = ±half_width).
    const T slack_u = ctx.u_range * T{1e-6};
    const T slack_v = ctx.v_range * T{1e-6};
    if (!ctx.pu && (u < ctx.dom.u_min - slack_u || u > ctx.dom.u_max + slack_u)) return std::nullopt;
    if (!ctx.pv && (v < ctx.dom.v_min - slack_v || v > ctx.dom.v_max + slack_v)) return std::nullopt;

    Vec<T, 3> n{su.cross(sv)};
    T nlen = n.norm();
    if (nlen > epsilon<T>()) n = n / nlen;

    return RayParametricHit<T>{t, u, v, Vec<T, 3>{ctx.o + ctx.d * t}, n};
}

} // namespace detail

template<Scalar T>
std::vector<RayParametricHit<T>> ray_parametric(
    const Ray<3, T>& r,
    const ParametricSurface<T>& surf,
    const RayParametricConfig<T>& cfg = {},
    std::optional<T> cell_radius_hint = std::nullopt)
{
    T cell_radius = cell_radius_hint ? *cell_radius_hint : estimate_cell_radius(surf, cfg);
    auto ctx = detail::make_context(r, surf, cfg, cell_radius);

    std::vector<RayParametricHit<T>> out;
    for (std::size_t j = 0; j < cfg.grid_v; ++j) {
        T v0 = ctx.dom.v_min + (static_cast<T>(j) + T{0.5}) * ctx.dv_step;
        for (std::size_t i = 0; i < cfg.grid_u; ++i) {
            T u0 = ctx.dom.u_min + (static_cast<T>(i) + T{0.5}) * ctx.du_step;
            if (auto hit = detail::solve_from_seed(ctx, u0, v0)) out.push_back(*hit);
        }
    }

    std::sort(out.begin(), out.end(),
              [](const auto& a, const auto& b) { return a.t < b.t; });

    std::vector<RayParametricHit<T>> uniq;
    uniq.reserve(out.size());
    for (const auto& h : out) {
        if (!uniq.empty() && std::abs(uniq.back().t - h.t) < cfg.dedup_t) continue;
        uniq.push_back(h);
    }
    return uniq;
}

// Nearest forward hit, or std::nullopt-equivalent via Result.
//
// uv_hint: a warm-start (u, v) — typically the previous scanline pixel's
// converged surface coordinate — tried as a single Newton solve before
// paying for the full grid_u*grid_v scan. Adjacent pixels usually hit
// nearly the same (u, v), so this converges in one solve instead of up
// to grid_u*grid_v; measured to close most of a real ~1000-1700x per-ray
// slowdown on Klein bottle/Möbius/bumpy-sphere renders (2026-08-25).
//
// Trade-off, accepted deliberately: on success this returns the hit in
// the hint's own local Newton basin, not necessarily the true
// globally-nearest hit on a self-occluding, multi-layer surface (the
// full scan below always finds that, and is still what runs whenever the
// hint is absent or fails to converge). Fine for coherent scanline
// rendering; not a drop-in replacement for exact geometric queries that
// need every intersection layer considered.
template<Scalar T>
Result<RayParametricHit<T>> ray_parametric_first(
    const Ray<3, T>& r,
    const ParametricSurface<T>& surf,
    const RayParametricConfig<T>& cfg = {},
    std::optional<std::pair<T, T>> uv_hint = std::nullopt,
    std::optional<T> cell_radius_hint = std::nullopt)
{
    if (uv_hint) {
        T cell_radius = cell_radius_hint ? *cell_radius_hint : estimate_cell_radius(surf, cfg);
        auto ctx = detail::make_context(r, surf, cfg, cell_radius);
        if (auto hit = detail::solve_from_seed(ctx, uv_hint->first, uv_hint->second))
            return *hit;
        // Hint missed (new silhouette/occlusion boundary, first pixel of
        // a scanline, etc.) — fall through to the exhaustive scan, no
        // worse off than not having a hint at all.
    }
    auto hits = ray_parametric(r, surf, cfg, cell_radius_hint);
    if (hits.empty())
        return std::unexpected(Error{ErrorCode::NoIntersection, "no forward hit"});
    return hits.front();
}

} // namespace spatium::geometry
