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
        auto [u, v] = find_params(p);
        return fn_(u, v);
    }

    TangentVector normal(const PointType& p) const {
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

private:
    ParamFn fn_;
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
    std::pair<T, T> find_params(const PointType& target) const {
        // Grid search for initial guess
        constexpr int GRID = 8;
        T best_u = domain_.u_min, best_v = domain_.v_min;
        T best_dist = std::numeric_limits<T>::max();

        T du_step = (domain_.u_max - domain_.u_min) / GRID;
        T dv_step = (domain_.v_max - domain_.v_min) / GRID;

        for (int j = 0; j <= GRID; ++j) {
            T v = domain_.v_min + static_cast<T>(j) * dv_step;
            for (int i = 0; i <= GRID; ++i) {
                T u = domain_.u_min + static_cast<T>(i) * du_step;
                T d = (fn_(u, v) - target).norm_squared();
                if (d < best_dist) { best_dist = d; best_u = u; best_v = v; }
            }
        }

        // Newton refinement (few iterations)
        for (int iter = 0; iter < 5; ++iter) {
            auto p = fn_(best_u, best_v);
            auto fu = this->du(best_u, best_v);
            auto fv = this->dv(best_u, best_v);
            auto r = target - p;

            // Solve 2x2 system: [fu·fu  fu·fv] [du] = [fu·r]
            //                   [fv·fu  fv·fv] [dv]   [fv·r]
            T a11 = fu.dot(fu), a12 = fu.dot(fv);
            T a21 = a12,        a22 = fv.dot(fv);
            T b1 = fu.dot(r),   b2 = fv.dot(r);

            T det = a11 * a22 - a12 * a21;
            if (std::abs(det) < epsilon<T>()) break;

            T delta_u = (a22 * b1 - a12 * b2) / det;
            T delta_v = (a11 * b2 - a21 * b1) / det;

            best_u += delta_u;
            best_v += delta_v;

            // Clamp to domain (or wrap for periodic)
            if (periodic_u_) {
                T range = domain_.u_max - domain_.u_min;
                best_u = domain_.u_min + std::fmod(best_u - domain_.u_min, range);
                if (best_u < domain_.u_min) best_u += range;
            } else {
                best_u = std::clamp(best_u, domain_.u_min, domain_.u_max);
            }
            if (periodic_v_) {
                T range = domain_.v_max - domain_.v_min;
                best_v = domain_.v_min + std::fmod(best_v - domain_.v_min, range);
                if (best_v < domain_.v_min) best_v += range;
            } else {
                best_v = std::clamp(best_v, domain_.v_min, domain_.v_max);
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
