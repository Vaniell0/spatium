#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/algebra/dual.hpp>
#  include <spatium/algebra/matrix.hpp>
#  include <spatium/algebra/ode.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <spatium/mesh/mesh.hpp>
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

// ── Marching Cubes tessellation ───────────────────────────────

// internal — do not use, no API stability. Holds the Paul Bourke
// MC lookup tables and the edge-interpolation helpers used by
// tessellate(). Public users should call tessellate() or the
// surface wrappers; this namespace's contents may change layout
// or be moved to a .cpp without notice.
namespace detail {

// Edge table and triangle table for marching cubes (compressed).
// Only the 256 cube configurations, 12 edges, up to 5 triangles per cube.
// Using Paul Bourke's tables.

// Edge vertices: edge index → (vertex_a, vertex_b)
inline constexpr std::array<std::array<int, 2>, 12> mc_edge_verts = {{
    {0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4},{0,4},{1,5},{2,6},{3,7}
}};

// Cube vertex offsets (i,j,k) for the 8 corners
inline constexpr std::array<std::array<int, 3>, 8> mc_corner = {{
    {0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}
}};

// Paul Bourke edge table: for each of the 256 cube configurations,
// a 12-bit mask indicating which edges are intersected by the isosurface.
inline constexpr std::array<uint16_t, 256> mc_edge_table = {{
    0x000, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
    0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
    0x190, 0x099, 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
    0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
    0x230, 0x339, 0x033, 0x13a, 0x636, 0x73f, 0x435, 0x53c,
    0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
    0x3a0, 0x2a9, 0x1a3, 0x0aa, 0x7a6, 0x6af, 0x5a5, 0x4ac,
    0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
    0x460, 0x569, 0x663, 0x76a, 0x066, 0x16f, 0x265, 0x36c,
    0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
    0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0x0ff, 0x3f5, 0x2fc,
    0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
    0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x055, 0x15c,
    0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
    0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0x0cc,
    0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
    0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
    0x0cc, 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
    0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
    0x15c, 0x055, 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
    0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
    0x2fc, 0x3f5, 0x0ff, 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
    0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
    0x36c, 0x265, 0x16f, 0x066, 0x76a, 0x663, 0x569, 0x460,
    0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
    0x4ac, 0x5a5, 0x6af, 0x7a6, 0x0aa, 0x1a3, 0x2a9, 0x3a0,
    0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
    0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x033, 0x339, 0x230,
    0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
    0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x099, 0x190,
    0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
    0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x000,
}};

// Paul Bourke triangle table: for each of the 256 cube configurations,
// a list of edge triplets forming triangles. -1 = end of list.
inline constexpr std::array<std::array<int8_t, 16>, 256> mc_tri_table = {{
    {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 0, 8, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 0, 1, 9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 1, 8, 3, 9, 8, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 1, 2,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 0, 8, 3, 1, 2,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 9, 2,10, 0, 2, 9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 2, 8, 3, 2,10, 8,10, 9, 8,-1,-1,-1,-1,-1,-1,-1},
    { 3,11, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 0,11, 2, 8,11, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 1, 9, 0, 2, 3,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 1,11, 2, 1, 9,11, 9, 8,11,-1,-1,-1,-1,-1,-1,-1},
    { 3,10, 1,11,10, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 0,10, 1, 0, 8,10, 8,11,10,-1,-1,-1,-1,-1,-1,-1},
    { 3, 9, 0, 3,11, 9,11,10, 9,-1,-1,-1,-1,-1,-1,-1},
    { 9, 8,10,10, 8,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 4, 7, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 4, 3, 0, 7, 3, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 0, 1, 9, 8, 4, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 4, 1, 9, 4, 7, 1, 7, 3, 1,-1,-1,-1,-1,-1,-1,-1},
    { 1, 2,10, 8, 4, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 3, 4, 7, 3, 0, 4, 1, 2,10,-1,-1,-1,-1,-1,-1,-1},
    { 9, 2,10, 9, 0, 2, 8, 4, 7,-1,-1,-1,-1,-1,-1,-1},
    { 2,10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4,-1,-1,-1,-1},
    { 8, 4, 7, 3,11, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {11, 4, 7,11, 2, 4, 2, 0, 4,-1,-1,-1,-1,-1,-1,-1},
    { 9, 0, 1, 8, 4, 7, 2, 3,11,-1,-1,-1,-1,-1,-1,-1},
    { 4, 7,11, 9, 4,11, 9,11, 2, 9, 2, 1,-1,-1,-1,-1},
    { 3,10, 1, 3,11,10, 7, 8, 4,-1,-1,-1,-1,-1,-1,-1},
    { 1,11,10, 1, 4,11, 1, 0, 4, 7,11, 4,-1,-1,-1,-1},
    { 4, 7, 8, 9, 0,11, 9,11,10,11, 0, 3,-1,-1,-1,-1},
    { 4, 7,11, 4,11, 9, 9,11,10,-1,-1,-1,-1,-1,-1,-1},
    { 9, 5, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 9, 5, 4, 0, 8, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 0, 5, 4, 1, 5, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 8, 5, 4, 8, 3, 5, 3, 1, 5,-1,-1,-1,-1,-1,-1,-1},
    { 1, 2,10, 9, 5, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 3, 0, 8, 1, 2,10, 4, 9, 5,-1,-1,-1,-1,-1,-1,-1},
    { 5, 2,10, 5, 4, 2, 4, 0, 2,-1,-1,-1,-1,-1,-1,-1},
    { 2,10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8,-1,-1,-1,-1},
    { 9, 5, 4, 2, 3,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 0,11, 2, 0, 8,11, 4, 9, 5,-1,-1,-1,-1,-1,-1,-1},
    { 0, 5, 4, 0, 1, 5, 2, 3,11,-1,-1,-1,-1,-1,-1,-1},
    { 2, 1, 5, 2, 5, 8, 2, 8,11, 4, 8, 5,-1,-1,-1,-1},
    {10, 3,11,10, 1, 3, 9, 5, 4,-1,-1,-1,-1,-1,-1,-1},
    { 4, 9, 5, 0, 8, 1, 8,10, 1, 8,11,10,-1,-1,-1,-1},
    { 5, 4, 0, 5, 0,11, 5,11,10,11, 0, 3,-1,-1,-1,-1},
    { 5, 4, 8, 5, 8,10,10, 8,11,-1,-1,-1,-1,-1,-1,-1},
    { 9, 7, 8, 5, 7, 9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 9, 3, 0, 9, 5, 3, 5, 7, 3,-1,-1,-1,-1,-1,-1,-1},
    { 0, 7, 8, 0, 1, 7, 1, 5, 7,-1,-1,-1,-1,-1,-1,-1},
    { 1, 5, 3, 3, 5, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 9, 7, 8, 9, 5, 7,10, 1, 2,-1,-1,-1,-1,-1,-1,-1},
    {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3,-1,-1,-1,-1},
    { 8, 0, 2, 8, 2, 5, 8, 5, 7,10, 5, 2,-1,-1,-1,-1},
    { 2,10, 5, 2, 5, 3, 3, 5, 7,-1,-1,-1,-1,-1,-1,-1},
    { 7, 9, 5, 7, 8, 9, 3,11, 2,-1,-1,-1,-1,-1,-1,-1},
    { 9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7,11,-1,-1,-1,-1},
    { 2, 3,11, 0, 1, 8, 1, 7, 8, 1, 5, 7,-1,-1,-1,-1},
    {11, 2, 1,11, 1, 7, 7, 1, 5,-1,-1,-1,-1,-1,-1,-1},
    { 9, 5, 8, 8, 5, 7,10, 1, 3,10, 3,11,-1,-1,-1,-1},
    { 5, 7, 0, 5, 0, 9, 7,11, 0, 1, 0,10,11,10, 0,-1},
    {11,10, 0,11, 0, 3,10, 5, 0, 8, 0, 7, 5, 7, 0,-1},
    {11,10, 5, 7,11, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {10, 6, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 0, 8, 3, 5,10, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 9, 0, 1, 5,10, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 1, 8, 3, 1, 9, 8, 5,10, 6,-1,-1,-1,-1,-1,-1,-1},
    { 1, 6, 5, 2, 6, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 1, 6, 5, 1, 2, 6, 3, 0, 8,-1,-1,-1,-1,-1,-1,-1},
    { 9, 6, 5, 9, 0, 6, 0, 2, 6,-1,-1,-1,-1,-1,-1,-1},
    { 5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8,-1,-1,-1,-1},
    { 2, 3,11,10, 6, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {11, 0, 8,11, 2, 0,10, 6, 5,-1,-1,-1,-1,-1,-1,-1},
    { 0, 1, 9, 2, 3,11, 5,10, 6,-1,-1,-1,-1,-1,-1,-1},
    { 5,10, 6, 1, 9, 2, 9,11, 2, 9, 8,11,-1,-1,-1,-1},
    { 6, 3,11, 6, 5, 3, 5, 1, 3,-1,-1,-1,-1,-1,-1,-1},
    { 0, 8,11, 0,11, 5, 0, 5, 1, 5,11, 6,-1,-1,-1,-1},
    { 3,11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9,-1,-1,-1,-1},
    { 6, 5, 9, 6, 9,11,11, 9, 8,-1,-1,-1,-1,-1,-1,-1},
    { 5,10, 6, 4, 7, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 4, 3, 0, 4, 7, 3, 6, 5,10,-1,-1,-1,-1,-1,-1,-1},
    { 1, 9, 0, 5,10, 6, 8, 4, 7,-1,-1,-1,-1,-1,-1,-1},
    {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4,-1,-1,-1,-1},
    { 6, 1, 2, 6, 5, 1, 4, 7, 8,-1,-1,-1,-1,-1,-1,-1},
    { 1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7,-1,-1,-1,-1},
    { 8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6,-1,-1,-1,-1},
    { 7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9,-1},
    { 3,11, 2, 7, 8, 4,10, 6, 5,-1,-1,-1,-1,-1,-1,-1},
    { 5,10, 6, 4, 7, 2, 4, 2, 0, 2, 7,11,-1,-1,-1,-1},
    { 0, 1, 9, 4, 7, 8, 2, 3,11, 5,10, 6,-1,-1,-1,-1},
    { 9, 2, 1, 9,11, 2, 9, 4,11, 7,11, 4, 5,10, 6,-1},
    { 8, 4, 7, 3,11, 5, 3, 5, 1, 5,11, 6,-1,-1,-1,-1},
    { 5, 1,11, 5,11, 6, 1, 0,11, 7,11, 4, 0, 4,11,-1},
    { 0, 5, 9, 0, 6, 5, 0, 3, 6,11, 6, 3, 8, 4, 7,-1},
    { 6, 5, 9, 6, 9,11, 4, 7, 9, 7,11, 9,-1,-1,-1,-1},
    {10, 4, 9, 6, 4,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 4,10, 6, 4, 9,10, 0, 8, 3,-1,-1,-1,-1,-1,-1,-1},
    {10, 0, 1,10, 6, 0, 6, 4, 0,-1,-1,-1,-1,-1,-1,-1},
    { 8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1,10,-1,-1,-1,-1},
    { 1, 4, 9, 1, 2, 4, 2, 6, 4,-1,-1,-1,-1,-1,-1,-1},
    { 3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4,-1,-1,-1,-1},
    { 0, 2, 4, 4, 2, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 8, 3, 2, 8, 2, 4, 4, 2, 6,-1,-1,-1,-1,-1,-1,-1},
    {10, 4, 9,10, 6, 4,11, 2, 3,-1,-1,-1,-1,-1,-1,-1},
    { 0, 8, 2, 2, 8,11, 4, 9,10, 4,10, 6,-1,-1,-1,-1},
    { 3,11, 2, 0, 1, 6, 0, 6, 4, 6, 1,10,-1,-1,-1,-1},
    { 6, 4, 1, 6, 1,10, 4, 8, 1, 2, 1,11, 8,11, 1,-1},
    { 9, 6, 4, 9, 3, 6, 9, 1, 3,11, 6, 3,-1,-1,-1,-1},
    { 8,11, 1, 8, 1, 0,11, 6, 1, 9, 1, 4, 6, 4, 1,-1},
    { 3,11, 6, 3, 6, 0, 0, 6, 4,-1,-1,-1,-1,-1,-1,-1},
    { 6, 4, 8,11, 6, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 7,10, 6, 7, 8,10, 8, 9,10,-1,-1,-1,-1,-1,-1,-1},
    { 0, 7, 3, 0,10, 7, 0, 9,10, 6, 7,10,-1,-1,-1,-1},
    {10, 6, 7, 1,10, 7, 1, 7, 8, 1, 8, 0,-1,-1,-1,-1},
    {10, 6, 7,10, 7, 1, 1, 7, 3,-1,-1,-1,-1,-1,-1,-1},
    { 1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7,-1,-1,-1,-1},
    { 2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9,-1},
    { 7, 8, 0, 7, 0, 6, 6, 0, 2,-1,-1,-1,-1,-1,-1,-1},
    { 7, 3, 2, 6, 7, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 2, 3,11,10, 6, 8,10, 8, 9, 8, 6, 7,-1,-1,-1,-1},
    { 2, 0, 7, 2, 7,11, 0, 9, 7, 6, 7,10, 9,10, 7,-1},
    { 1, 8, 0, 1, 7, 8, 1,10, 7, 6, 7,10, 2, 3,11,-1},
    {11, 2, 1,11, 1, 7,10, 6, 1, 6, 7, 1,-1,-1,-1,-1},
    { 8, 9, 6, 8, 6, 7, 9, 1, 6,11, 6, 3, 1, 3, 6,-1},
    { 0, 9, 1,11, 6, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 7, 8, 0, 7, 0, 6, 3,11, 0,11, 6, 0,-1,-1,-1,-1},
    { 7,11, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 7, 6,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 3, 0, 8,11, 7, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 0, 1, 9,11, 7, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 8, 1, 9, 8, 3, 1,11, 7, 6,-1,-1,-1,-1,-1,-1,-1},
    {10, 1, 2, 6,11, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 1, 2,10, 3, 0, 8, 6,11, 7,-1,-1,-1,-1,-1,-1,-1},
    { 2, 9, 0, 2,10, 9, 6,11, 7,-1,-1,-1,-1,-1,-1,-1},
    { 6,11, 7, 2,10, 3,10, 8, 3,10, 9, 8,-1,-1,-1,-1},
    { 7, 2, 3, 6, 2, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 7, 0, 8, 7, 6, 0, 6, 2, 0,-1,-1,-1,-1,-1,-1,-1},
    { 2, 7, 6, 2, 3, 7, 0, 1, 9,-1,-1,-1,-1,-1,-1,-1},
    { 1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6,-1,-1,-1,-1},
    {10, 7, 6,10, 1, 7, 1, 3, 7,-1,-1,-1,-1,-1,-1,-1},
    {10, 7, 6, 1, 7,10, 1, 8, 7, 1, 0, 8,-1,-1,-1,-1},
    { 0, 3, 7, 0, 7,10, 0,10, 9, 6,10, 7,-1,-1,-1,-1},
    { 7, 6,10, 7,10, 8, 8,10, 9,-1,-1,-1,-1,-1,-1,-1},
    { 6, 8, 4,11, 8, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 3, 6,11, 3, 0, 6, 0, 4, 6,-1,-1,-1,-1,-1,-1,-1},
    { 8, 6,11, 8, 4, 6, 9, 0, 1,-1,-1,-1,-1,-1,-1,-1},
    { 9, 4, 6, 9, 6, 3, 9, 3, 1,11, 3, 6,-1,-1,-1,-1},
    { 6, 8, 4, 6,11, 8, 2,10, 1,-1,-1,-1,-1,-1,-1,-1},
    { 1, 2,10, 3, 0,11, 0, 6,11, 0, 4, 6,-1,-1,-1,-1},
    { 4,11, 8, 4, 6,11, 0, 2, 9, 2,10, 9,-1,-1,-1,-1},
    {10, 9, 3,10, 3, 2, 9, 4, 3,11, 3, 6, 4, 6, 3,-1},
    { 8, 2, 3, 8, 4, 2, 4, 6, 2,-1,-1,-1,-1,-1,-1,-1},
    { 0, 4, 2, 4, 6, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8,-1,-1,-1,-1},
    { 1, 9, 4, 1, 4, 2, 2, 4, 6,-1,-1,-1,-1,-1,-1,-1},
    { 8, 1, 3, 8, 6, 1, 8, 4, 6, 6,10, 1,-1,-1,-1,-1},
    {10, 1, 0,10, 0, 6, 6, 0, 4,-1,-1,-1,-1,-1,-1,-1},
    { 4, 6, 3, 4, 3, 8, 6,10, 3, 0, 3, 9,10, 9, 3,-1},
    {10, 9, 4, 6,10, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 4, 9, 5, 7, 6,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 0, 8, 3, 4, 9, 5,11, 7, 6,-1,-1,-1,-1,-1,-1,-1},
    { 5, 0, 1, 5, 4, 0, 7, 6,11,-1,-1,-1,-1,-1,-1,-1},
    {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5,-1,-1,-1,-1},
    { 9, 5, 4,10, 1, 2, 7, 6,11,-1,-1,-1,-1,-1,-1,-1},
    { 6,11, 7, 1, 2,10, 0, 8, 3, 4, 9, 5,-1,-1,-1,-1},
    { 7, 6,11, 5, 4,10, 4, 2,10, 4, 0, 2,-1,-1,-1,-1},
    { 3, 4, 8, 3, 5, 4, 3, 2, 5,10, 5, 2,11, 7, 6,-1},
    { 7, 2, 3, 7, 6, 2, 5, 4, 9,-1,-1,-1,-1,-1,-1,-1},
    { 9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7,-1,-1,-1,-1},
    { 3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0,-1,-1,-1,-1},
    { 6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8,-1},
    { 9, 5, 4,10, 1, 6, 1, 7, 6, 1, 3, 7,-1,-1,-1,-1},
    { 1, 6,10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4,-1},
    { 4, 0,10, 4,10, 5, 0, 3,10, 6,10, 7, 3, 7,10,-1},
    { 7, 6,10, 7,10, 8, 5, 4,10, 4, 8,10,-1,-1,-1,-1},
    { 6, 9, 5, 6,11, 9,11, 8, 9,-1,-1,-1,-1,-1,-1,-1},
    { 3, 6,11, 0, 6, 3, 0, 5, 6, 0, 9, 5,-1,-1,-1,-1},
    { 0,11, 8, 0, 5,11, 0, 1, 5, 5, 6,11,-1,-1,-1,-1},
    { 6,11, 3, 6, 3, 5, 5, 3, 1,-1,-1,-1,-1,-1,-1,-1},
    { 1, 2,10, 9, 5,11, 9,11, 8,11, 5, 6,-1,-1,-1,-1},
    { 0,11, 3, 0, 6,11, 0, 9, 6, 5, 6, 9, 1, 2,10,-1},
    {11, 8, 5,11, 5, 6, 8, 0, 5,10, 5, 2, 0, 2, 5,-1},
    { 6,11, 3, 6, 3, 5, 2,10, 3,10, 5, 3,-1,-1,-1,-1},
    { 5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2,-1,-1,-1,-1},
    { 9, 5, 6, 9, 6, 0, 0, 6, 2,-1,-1,-1,-1,-1,-1,-1},
    { 1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8,-1},
    { 1, 5, 6, 2, 1, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 1, 3, 6, 1, 6,10, 3, 8, 6, 5, 6, 9, 8, 9, 6,-1},
    {10, 1, 0,10, 0, 6, 9, 5, 0, 5, 6, 0,-1,-1,-1,-1},
    { 0, 3, 8, 5, 6,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {10, 5, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {11, 5,10, 7, 5,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {11, 5,10,11, 7, 5, 8, 3, 0,-1,-1,-1,-1,-1,-1,-1},
    { 5,11, 7, 5,10,11, 1, 9, 0,-1,-1,-1,-1,-1,-1,-1},
    {10, 7, 5,10,11, 7, 9, 8, 1, 8, 3, 1,-1,-1,-1,-1},
    {11, 1, 2,11, 7, 1, 7, 5, 1,-1,-1,-1,-1,-1,-1,-1},
    { 0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2,11,-1,-1,-1,-1},
    { 9, 7, 5, 9, 2, 7, 9, 0, 2, 2,11, 7,-1,-1,-1,-1},
    { 7, 5, 2, 7, 2,11, 5, 9, 2, 3, 2, 8, 9, 8, 2,-1},
    { 2, 5,10, 2, 3, 5, 3, 7, 5,-1,-1,-1,-1,-1,-1,-1},
    { 8, 2, 0, 8, 5, 2, 8, 7, 5,10, 2, 5,-1,-1,-1,-1},
    { 9, 0, 1, 5,10, 3, 5, 3, 7, 3,10, 2,-1,-1,-1,-1},
    { 9, 8, 2, 9, 2, 1, 8, 7, 2,10, 2, 5, 7, 5, 2,-1},
    { 1, 3, 5, 3, 7, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 0, 8, 7, 0, 7, 1, 1, 7, 5,-1,-1,-1,-1,-1,-1,-1},
    { 9, 0, 3, 9, 3, 5, 5, 3, 7,-1,-1,-1,-1,-1,-1,-1},
    { 9, 8, 7, 5, 9, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 5, 8, 4, 5,10, 8,10,11, 8,-1,-1,-1,-1,-1,-1,-1},
    { 5, 0, 4, 5,11, 0, 5,10,11,11, 3, 0,-1,-1,-1,-1},
    { 0, 1, 9, 8, 4,10, 8,10,11,10, 4, 5,-1,-1,-1,-1},
    {10,11, 4,10, 4, 5,11, 3, 4, 9, 4, 1, 3, 1, 4,-1},
    { 2, 5, 1, 2, 8, 5, 2,11, 8, 4, 5, 8,-1,-1,-1,-1},
    { 0, 4,11, 0,11, 3, 4, 5,11, 2,11, 1, 5, 1,11,-1},
    { 0, 2, 5, 0, 5, 9, 2,11, 5, 4, 5, 8,11, 8, 5,-1},
    { 9, 4, 5, 2,11, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 2, 5,10, 3, 5, 2, 3, 4, 5, 3, 8, 4,-1,-1,-1,-1},
    { 5,10, 2, 5, 2, 4, 4, 2, 0,-1,-1,-1,-1,-1,-1,-1},
    { 3,10, 2, 3, 5,10, 3, 8, 5, 4, 5, 8, 0, 1, 9,-1},
    { 5,10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2,-1,-1,-1,-1},
    { 8, 4, 5, 8, 5, 3, 3, 5, 1,-1,-1,-1,-1,-1,-1,-1},
    { 0, 4, 5, 1, 0, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5,-1,-1,-1,-1},
    { 9, 4, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 4,11, 7, 4, 9,11, 9,10,11,-1,-1,-1,-1,-1,-1,-1},
    { 0, 8, 3, 4, 9, 7, 9,11, 7, 9,10,11,-1,-1,-1,-1},
    { 1,10,11, 1,11, 4, 1, 4, 0, 7, 4,11,-1,-1,-1,-1},
    { 3, 1, 4, 3, 4, 8, 1,10, 4, 7, 4,11,10,11, 4,-1},
    { 4,11, 7, 9,11, 4, 9, 2,11, 9, 1, 2,-1,-1,-1,-1},
    { 9, 7, 4, 9,11, 7, 9, 1,11, 2,11, 1, 0, 8, 3,-1},
    {11, 7, 4,11, 4, 2, 2, 4, 0,-1,-1,-1,-1,-1,-1,-1},
    {11, 7, 4,11, 4, 2, 8, 3, 4, 3, 2, 4,-1,-1,-1,-1},
    { 2, 9,10, 2, 7, 9, 2, 3, 7, 7, 4, 9,-1,-1,-1,-1},
    { 9,10, 7, 9, 7, 4,10, 2, 7, 8, 7, 0, 2, 0, 7,-1},
    { 3, 7,10, 3,10, 2, 7, 4,10, 1,10, 0, 4, 0,10,-1},
    { 1,10, 2, 8, 7, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 4, 9, 1, 4, 1, 7, 7, 1, 3,-1,-1,-1,-1,-1,-1,-1},
    { 4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1,-1,-1,-1,-1},
    { 4, 0, 3, 7, 4, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 4, 8, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 9,10, 8,10,11, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 3, 0, 9, 3, 9,11,11, 9,10,-1,-1,-1,-1,-1,-1,-1},
    { 0, 1,10, 0,10, 8, 8,10,11,-1,-1,-1,-1,-1,-1,-1},
    { 3, 1,10,11, 3,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 1, 2,11, 1,11, 9, 9,11, 8,-1,-1,-1,-1,-1,-1,-1},
    { 3, 0, 9, 3, 9,11, 1, 2, 9, 2,11, 9,-1,-1,-1,-1},
    { 0, 2,11, 8, 0,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 3, 2,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 2, 3, 8, 2, 8,10,10, 8, 9,-1,-1,-1,-1,-1,-1,-1},
    { 9,10, 2, 0, 9, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 2, 3, 8, 2, 8,10, 0, 1, 8, 1,10, 8,-1,-1,-1,-1},
    { 1,10, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 1, 3, 8, 9, 1, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 0, 9, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 0, 3, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
}};

} // namespace detail

template<Scalar T>
mesh::Mesh<ImplicitSurface<T>> marching_cubes(
    const ImplicitSurface<T>& surf, std::size_t resolution)
{
    // Simplified marching cubes: evaluate grid, find sign changes on edges,
    // interpolate vertices, emit triangles.

    auto& b = surf.bounds();
    auto res = static_cast<T>(resolution);
    T dx = (b.x_max - b.x_min) / res;
    T dy = (b.y_max - b.y_min) / res;
    T dz = (b.z_max - b.z_min) / res;

    // Evaluate scalar field on grid
    auto idx = [&](std::size_t i, std::size_t j, std::size_t k) {
        return k * (resolution + 1) * (resolution + 1) + j * (resolution + 1) + i;
    };

    std::size_t n = resolution + 1;
    std::vector<T> values(n * n * n);
    for (std::size_t k = 0; k < n; ++k)
        for (std::size_t j = 0; j < n; ++j)
            for (std::size_t i = 0; i < n; ++i) {
                T x = b.x_min + static_cast<T>(i) * dx;
                T y = b.y_min + static_cast<T>(j) * dy;
                T z = b.z_min + static_cast<T>(k) * dz;
                values[idx(i, j, k)] = surf(x, y, z);
            }

    mesh::Mesh<ImplicitSurface<T>> m;

    // Edge vertex cache: edge key → vertex index
    struct PairHash {
        std::size_t operator()(std::pair<std::size_t, std::size_t> p) const {
            return std::hash<std::size_t>{}(p.first * 2654435761u ^ p.second);
        }
    };
    std::unordered_map<std::pair<std::size_t, std::size_t>, uint32_t, PairHash> edge_cache;

    auto interp_edge = [&](std::size_t ci, std::size_t cj, std::size_t ck, int edge) -> uint32_t {
        auto [ea, eb] = detail::mc_edge_verts[edge];
        auto [ai, aj, ak] = detail::mc_corner[ea];
        auto [bi, bj, bk] = detail::mc_corner[eb];

        std::size_t ia = idx(ci + ai, cj + aj, ck + ak);
        std::size_t ib = idx(ci + bi, cj + bj, ck + bk);

        auto key = std::minmax(ia, ib);
        if (auto it = edge_cache.find(key); it != edge_cache.end())
            return it->second;

        T va = values[ia], vb = values[ib];
        using std::abs;
        T t = (abs(T(va - vb)) > epsilon<T>()) ? va / (va - vb) : T{0.5};
        t = std::clamp(t, T{0}, T{1});

        Vec<T, 3> pa{
            b.x_min + static_cast<T>(ci + ai) * dx,
            b.y_min + static_cast<T>(cj + aj) * dy,
            b.z_min + static_cast<T>(ck + ak) * dz
        };
        Vec<T, 3> pb{
            b.x_min + static_cast<T>(ci + bi) * dx,
            b.y_min + static_cast<T>(cj + bj) * dy,
            b.z_min + static_cast<T>(ck + bk) * dz
        };

        auto vidx = static_cast<uint32_t>(m.vertices.size());
        m.vertices.push_back(Vec<T,3>{pa + (pb - pa) * t});
        edge_cache[key] = vidx;
        return vidx;
    };

    // Process each cube cell
    for (std::size_t k = 0; k < resolution; ++k)
        for (std::size_t j = 0; j < resolution; ++j)
            for (std::size_t i = 0; i < resolution; ++i) {
                // Classify corners
                uint8_t cube_idx = 0;
                std::array<T, 8> cv;
                for (int c = 0; c < 8; ++c) {
                    auto [ci, cj, ck] = detail::mc_corner[c];
                    cv[c] = values[idx(i + ci, j + cj, k + ck)];
                    if (cv[c] < T{0}) cube_idx |= (1 << c);
                }

                if (cube_idx == 0 || cube_idx == 255) continue;

                uint16_t edge_mask = detail::mc_edge_table[cube_idx];
                if (edge_mask == 0) continue;

                // Interpolate vertices on intersected edges
                std::array<uint32_t, 12> edge_verts{};
                for (int e = 0; e < 12; ++e)
                    if (edge_mask & (1 << e))
                        edge_verts[e] = interp_edge(i, j, k, e);

                // Emit triangles from lookup table
                auto& tri = detail::mc_tri_table[cube_idx];
                for (int t = 0; t < 16 && tri[t] >= 0; t += 3)
                    m.faces.push_back({edge_verts[tri[t]],
                                       edge_verts[tri[t + 1]],
                                       edge_verts[tri[t + 2]]});
            }

    return m;
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
