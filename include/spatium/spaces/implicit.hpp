#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/algebra/dual.hpp>
#  include <spatium/algebra/matrix.hpp>
#  include <spatium/algebra/ode.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <spatium/spaces/euclidean.hpp>
#  include <array>
#  include <cmath>
#  include <concepts>
#  include <cstdint>
#  include <functional>
#  include <limits>
#  include <numbers>
#  include <utility>
#  include <unordered_map>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium {

// ImplicitSurface: F(x,y,z) = 0 defines a surface.
// Marching cubes extracts a mesh. Gradient gives normal.
// Satisfies Surface concept for geodesic algorithms.
//
// As ParametricSurface, the function is a type, `F`, by default the erased
// std::function. Make one with `make_implicit` from a generic lambda
// (`[](auto x, auto y, auto z) { ... }`) and it can be called on a Dual and a
// Dual of a Dual: the gradient is exact, the Hessian is exact, the Gaussian
// curvature follows from the two, and exp and log are geodesics -- where the
// erased form has a gradient taken by a difference quotient and a step in the
// ambient space projected back.

namespace implicit_detail {

// F can be called on a Dual and on a Dual of a Dual: a generic callable.
template<class F, class T>
concept DualEvaluable =
    requires(const F& f, Dual<T> a, Dual<T> b, Dual<T> c) {
        { f(a, b, c) } -> std::convertible_to<Dual<T>>;
    } &&
    requires(const F& f, Dual<Dual<T>> a, Dual<Dual<T>> b, Dual<Dual<T>> c) {
        { f(a, b, c) } -> std::convertible_to<Dual<Dual<T>>>;
    };

}  // namespace implicit_detail

template<Scalar T = double, class F = std::function<T(T, T, T)>>
class ImplicitSurface {
public:
    using ScalarType = T;
    using PointType = Vec<T, 3>;
    using TangentVector = Vec<T, 3>;
    static constexpr std::size_t dimension = 2;
    static constexpr bool is_complete = false;

    using ImplicitFn = F;

    // True when F can be called on a Dual and a Dual of a Dual: the gradient
    // and the Hessian are then exact and the geodesics are the surface's.
    static constexpr bool exact = implicit_detail::DualEvaluable<F, T>;

    struct Bounds { T x_min, x_max, y_min, y_max, z_min, z_max; };

    ImplicitSurface(F fn, Bounds bounds)
        : fn_(std::move(fn)), bounds_(bounds) {}

    T operator()(T x, T y, T z) const { return fn_(x, y, z); }
    T operator()(const PointType& p) const { return fn_(p[0], p[1], p[2]); }

    // ── Surface concept ───────────────────────────────────────

    bool contains(const PointType& p) const {
        using std::abs;   // ADL: a qualified std::abs does not compile for Dual or Real50
        return abs(fn_(p[0], p[1], p[2])) < epsilon<T>() * T{100};
    }

    // Euclidean distance in ambient R³ — a valid metric but NOT the geodesic
    // (intrinsic) distance along the surface. Acts as a lower bound for geodesic distance.
    // Only for the erased form: a surface that has its geodesics has no member
    // `distance`, so core/access.hpp derives |log|_g, the intrinsic one.
    T distance(const PointType& a, const PointType& b) const requires (!exact) {
        return (a - b).norm();
    }

    // The chord, for any surface, under its own name.
    T chord_distance(const PointType& a, const PointType& b) const {
        return (a - b).norm();
    }

    PointType project(const PointType& p) const {
        // Newton projection onto F=0
        using std::abs;
        auto q = p;
        for (int i = 0; i < 20; ++i) {
            T val = fn_(q[0], q[1], q[2]);
            if (abs(val) < epsilon<T>()) break;
            auto g = gradient(q);
            T g2 = g.dot(g);
            if (g2 < epsilon<T>()) break;
            q = PointType{q - g * (val / g2)};
        }
        return q;
    }

    TangentVector normal(const PointType& p) const {
        return gradient(p).normalized();
    }

    // The erased form: a retraction -- a step in the ambient space projected
    // back, first order in t. exp(log) is not the identity on it.
    PointType exp_map(const PointType& p, const TangentVector& v, T t) const requires (!exact) {
        return project(PointType{p + v * t});
    }

    TangentVector log_map(const PointType& p, const PointType& q) const requires (!exact) {
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

    // The geodesic, in the ambient space. A curve on F = 0 has F(x(t)) = 0, so
    // grad F . x' = 0, and differentiating again grad F . x'' = -x'^T H x'; a
    // geodesic's acceleration is normal to the surface, x'' = lambda grad F, so
    //   x'' = -(x'^T H x') / |grad F|^2  grad F.
    // Regular wherever grad F is -- there is no chart to degenerate at a pole.
    // Integrated by the extrapolated flow to the scalar's own precision, then
    // put back on the surface, which the ODE alone keeps to its tolerance only.
    // NaN where grad F vanishes or the flow cannot be integrated.
    PointType exp_map(const PointType& p, const TangentVector& w, T t) const requires (exact) {
        const PointType x0 = project(p);
        const TangentVector n = normal(x0);
        const TangentVector v0{w - n * w.dot(n)};
        const Vec<T, 6> y0{x0[0], x0[1], x0[2], v0[0], v0[1], v0[2]};
        const auto flow = integrate_extrapolated(
            [this](T, const Vec<T, 6>& y) { return acceleration(y); }, T{0}, y0, t, tolerance());
        if (!flow) return nan_point();
        return project(PointType{(*flow)[0], (*flow)[1], (*flow)[2]});
    }

    // Gauss-Newton shooting in the tangent plane at p: unknowns (a, b) in an
    // orthonormal basis of it, residual exp(a e1 + b e2) - q in R^3, started
    // from the chord's tangent part. The Jacobian by central differences of exp
    // only sets the rate; the answer is a fixed point of exp itself.
    TangentVector log_map(const PointType& p, const PointType& q) const requires (exact) {
        const PointType x0 = project(p);
        const TangentVector n = normal(x0);
        // the coordinate axis least aligned with n makes the first basis vector
        std::size_t axis = 0;
        for (std::size_t k = 1; k < 3; ++k) {
            using std::abs;
            if (abs(n[k]) < abs(n[axis])) axis = k;
        }
        TangentVector unit{};
        unit[axis] = T{1};
        const TangentVector e1 = n.cross(unit).normalized();
        const TangentVector e2 = n.cross(e1);

        const auto reach = [&](T a, T b) { return exp_map(x0, TangentVector{e1 * a + e2 * b}, T{1}); };
        const auto miss = [&](const PointType& r) {
            double m = 0;
            for (std::size_t i = 0; i < 3; ++i) {
                const double d = std::abs(primal_double(T(r[i] - q[i])));
                if (!(d <= m)) m = d;
            }
            return m;
        };

        const TangentVector chord{q - x0};
        T a = chord.dot(e1), b = chord.dot(e2);
        const double stop = 8 * tolerance();
        const double step = std::cbrt(machine_epsilon<T>());
        double rs = miss(reach(a, b));
        for (int it = 0; it < 40; ++it) {
            if (rs <= stop) break;
            const T ha = T(step * (1.0 + std::abs(primal_double(a))));
            const T hb = T(step * (1.0 + std::abs(primal_double(b))));
            const PointType r0 = reach(a, b);
            const TangentVector ja{(reach(T(a + ha), b) - reach(T(a - ha), b)) * (T{1} / (T{2} * ha))};
            const TangentVector jb{(reach(a, T(b + hb)) - reach(a, T(b - hb))) * (T{1} / (T{2} * hb))};
            const TangentVector res{r0 - q};
            // (J^T J) d = J^T res
            const T saa = ja.dot(ja), sab = ja.dot(jb), sbb = jb.dot(jb);
            const T det = saa * sbb - sab * sab;
            if (!(det > T{0})) return TangentVector{nan_point()};
            const T ra = ja.dot(res), rb = jb.dot(res);
            const T da = (sbb * ra - sab * rb) / det, db = (saa * rb - sab * ra) / det;
            bool moved = false;
            double scale = 1.0;
            for (int back = 0; back < 10; ++back, scale *= 0.5) {
                const T na = T(a - da * T(scale)), nb = T(b - db * T(scale));
                const double rn = miss(reach(na, nb));
                if (rn < rs) { a = na; b = nb; rs = rn; moved = true; break; }
            }
            if (!moved) break;
        }
        if (!(rs <= std::sqrt(machine_epsilon<T>()))) return TangentVector{nan_point()};
        return TangentVector{e1 * a + e2 * b};
    }

    ScalarType metric_at(const PointType&,
                         const TangentVector& u, const TangentVector& v) const {
        return u.dot(v);
    }

    // The gradient: exact, through a Dual, when F can be called on one; a
    // central difference otherwise.
    TangentVector gradient(const PointType& p) const {
        if constexpr (exact) {
            TangentVector g{};
            for (std::size_t k = 0; k < 3; ++k) {
                Dual<T> c[3];
                for (std::size_t j = 0; j < 3; ++j)
                    c[j] = (j == k) ? Dual<T>::variable(p[j]) : Dual<T>::constant(p[j]);
                g[k] = fn_(c[0], c[1], c[2]).deriv;
            }
            return g;
        } else {
            T h = epsilon<T>() * T{1000};
            return {
                (fn_(p[0] + h, p[1], p[2]) - fn_(p[0] - h, p[1], p[2])) / (T{2} * h),
                (fn_(p[0], p[1] + h, p[2]) - fn_(p[0], p[1] - h, p[2])) / (T{2} * h),
                (fn_(p[0], p[1], p[2] + h) - fn_(p[0], p[1], p[2] - h)) / (T{2} * h),
            };
        }
    }

    // The Hessian, exact: F on a Dual of a Dual, one seed in each. Only for the
    // typed form -- a second difference quotient has no accuracy to speak of.
    Matrix<T, 3, 3> hessian(const PointType& p) const requires (exact) {
        using DD = Dual<Dual<T>>;
        Matrix<T, 3, 3> H{};
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = i; j < 3; ++j) {
                DD c[3];
                for (std::size_t k = 0; k < 3; ++k) {
                    const Dual<T> inner = (k == j) ? Dual<T>::variable(p[k]) : Dual<T>::constant(p[k]);
                    c[k] = (k == i) ? DD::variable(inner) : DD::constant(inner);
                }
                H(i, j) = H(j, i) = fn_(c[0], c[1], c[2]).deriv.deriv;
            }
        return H;
    }

    // Gaussian curvature at a point of the surface, from the gradient and the
    // Hessian alone: K = grad F . adj(H) grad F / |grad F|^4 (Goldman 2005),
    // adj the adjugate. Exact -- no tessellation, no fitted patch.
    T gaussian_curvature(const PointType& p) const requires (exact) {
        const TangentVector g = gradient(p);
        const Matrix<T, 3, 3> H = hessian(p);
        const T c00 = H(1, 1) * H(2, 2) - H(1, 2) * H(1, 2);
        const T c11 = H(0, 0) * H(2, 2) - H(0, 2) * H(0, 2);
        const T c22 = H(0, 0) * H(1, 1) - H(0, 1) * H(0, 1);
        const T c01 = H(0, 2) * H(1, 2) - H(0, 1) * H(2, 2);
        const T c02 = H(0, 1) * H(1, 2) - H(1, 1) * H(0, 2);
        const T c12 = H(0, 1) * H(0, 2) - H(0, 0) * H(1, 2);
        const T quad = g[0] * (c00 * g[0] + c01 * g[1] + c02 * g[2]) +
                       g[1] * (c01 * g[0] + c11 * g[1] + c12 * g[2]) +
                       g[2] * (c02 * g[0] + c12 * g[1] + c22 * g[2]);
        const T g2 = g.dot(g);
        return quad / (g2 * g2);
    }

    const Bounds& bounds() const { return bounds_; }

    // The same surface with its function stored as a std::function: the form
    // tessellate and marching cubes take. An explicit boundary -- what is lost
    // crossing it is the exact derivatives and the geodesics.
    ImplicitSurface<T> erased() const {
        if constexpr (std::same_as<F, std::function<T(T, T, T)>>) {
            return *this;
        } else {
            return ImplicitSurface<T>(std::function<T(T, T, T)>(fn_),
                typename ImplicitSurface<T>::Bounds{bounds_.x_min, bounds_.x_max, bounds_.y_min,
                                                    bounds_.y_max, bounds_.z_min, bounds_.z_max});
        }
    }

private:
    F fn_;
    Bounds bounds_;

    // The tolerance of the flow and the shooting: eps^(2/3) of the scalar.
    static double tolerance() { return std::pow(machine_epsilon<T>(), 2.0 / 3.0); }

    static PointType nan_point() {
        const T nan = T(std::numeric_limits<double>::quiet_NaN());
        return PointType{nan, nan, nan};
    }

    // d/dt of (x, x'): x'' = -(x'^T H x') / |grad F|^2 grad F.
    Vec<T, 6> acceleration(const Vec<T, 6>& y) const requires (exact) {
        const PointType x{y[0], y[1], y[2]};
        const TangentVector v{y[3], y[4], y[5]};
        const TangentVector g = gradient(x);
        const Matrix<T, 3, 3> H = hessian(x);
        const T g2 = g.dot(g);
        const T nan = T(std::numeric_limits<double>::quiet_NaN());
        T lambda = (g2 > T{0}) ? T(-(v[0] * (H(0, 0) * v[0] + H(0, 1) * v[1] + H(0, 2) * v[2]) +
                                     v[1] * (H(1, 0) * v[0] + H(1, 1) * v[1] + H(1, 2) * v[2]) +
                                     v[2] * (H(2, 0) * v[0] + H(2, 1) * v[1] + H(2, 2) * v[2])) / g2)
                                : nan;
        return Vec<T, 6>{v[0], v[1], v[2], lambda * g[0], lambda * g[1], lambda * g[2]};
    }
};

// A surface whose function stays a type, generic over the scalar: write it
// `[](auto x, auto y, auto z) { return x * x + y * y + z * z - one; }` and give
// the scalar T it lives over. Its gradient and Hessian are exact, its curvature
// follows, exp and log are its geodesics; `.erased()` is the way back.
template<Scalar T = double, class F>
ImplicitSurface<T, F> make_implicit(F fn, typename ImplicitSurface<T, F>::Bounds bounds) {
    return ImplicitSurface<T, F>(std::move(fn), bounds);
}

// ── Convenience factories ─────────────────────────────────────

template<Scalar T = double>
ImplicitSurface<T> make_implicit_sphere(T radius = T{1}) {
    return ImplicitSurface<T>(
        [=](T x, T y, T z) { return x * x + y * y + z * z - radius * radius; },
        {-radius * T{1.5}, radius * T{1.5},
         -radius * T{1.5}, radius * T{1.5},
         -radius * T{1.5}, radius * T{1.5}}
    );
}

template<Scalar T = double>
ImplicitSurface<T> make_implicit_torus(T major_r = T{2}, T minor_r = T{1}) {
    T bound = (major_r + minor_r) * T{1.5};
    return ImplicitSurface<T>(
        [=](T x, T y, T z) {
            using std::sqrt;
            T d = sqrt(x * x + y * y) - major_r;
            return d * d + z * z - minor_r * minor_r;
        },
        {-bound, bound, -bound, bound, -minor_r * T{1.5}, minor_r * T{1.5}}
    );
}

// A sphere perturbed by a handful of fixed low-frequency directional
// terms -- no RNG/noise-library dependency, just a small fixed sum of
// sines over the surface direction, giving an irregular "little planet"
// silhouette. `amplitude` is deliberately kept well under `radius` (and
// the term weights sum to 1) so the zero-level-set stays a single closed
// topological sphere instead of pinching into multiple components.
template<Scalar T = double>
ImplicitSurface<T> make_bumpy_sphere(T radius = T{1}, T amplitude = T{0.15}) {
    struct Term { T dx, dy, dz, freq, phase, weight; };
    constexpr std::array<Term, 4> terms{{
        {0.36, 0.93, 0.10, 3, 0.7, 0.45},
        {-0.71, 0.20, 0.67, 4, 2.1, 0.30},
        {0.12, -0.58, 0.81, 5, 4.4, 0.15},
        {-0.44, -0.65, -0.62, 2, 5.9, 0.10},
    }};

    return ImplicitSurface<T>(
        [=](T x, T y, T z) {
            using std::sin; using std::sqrt;
            T r = sqrt(x * x + y * y + z * z);
            if (r < radius * T{0.01}) return r - radius; // avoid direction at the origin
            T inv_r = T{1} / r;
            T bump{0};
            for (auto& term : terms) {
                T d = term.dx * x + term.dy * y + term.dz * z; // unnormalized dir . (x,y,z)
                bump += term.weight * sin(term.freq * d * inv_r + term.phase);
            }
            return r - radius * (T{1} + amplitude * bump);
        },
        {-radius * T{1.6}, radius * T{1.6},
         -radius * T{1.6}, radius * T{1.6},
         -radius * T{1.6}, radius * T{1.6}}
    );
}

template<Scalar T = double>
ImplicitSurface<T> make_gyroid() {
    constexpr T pi = std::numbers::pi_v<T>;
    T bound = T{2} * pi;
    return ImplicitSurface<T>(
        [](T x, T y, T z) {
            using std::sin; using std::cos;
            return sin(x) * cos(y) + sin(y) * cos(z) + sin(z) * cos(x);
        },
        {-bound, bound, -bound, bound, -bound, bound}
    );
}

} // namespace spatium
