#pragma once
// The connectivity matrix's cells: a scalar, a space, a probe.
//
// A cell instantiates one generic algorithm of the library on one space
// over one scalar and reports how far it got:
//
//   L0  it compiles              -- found only by building the cell alone,
//                                   since a concept checks declarations,
//                                   not bodies (scripts/gen_connectivity.py)
//   L1  its results are finite
//   L2  an axiom or a symmetry holds within the scalar's own tolerance --
//       two paths to one answer agreeing: a closed-form distance and the
//       one derived from log and the metric, exp after log, a midpoint
//       equidistant, a derivative through Dual and by finite differences
//   L3  its numbers agree with the same cell over another scalar: double
//       against Real50, every other scalar against double
//
// L3 is decided across cells from the signatures they print, by the
// script and by tests/test_connectivity.cpp alike; a cell alone reports
// L1 or L2. Everything a cell computes goes through the customization
// points of core/access.hpp where they exist, so a space reached by ADL or
// by derivation is probed the same way as one with members.

#include <spatium/algebra/calculus.hpp>
#include <spatium/algebra/dual.hpp>
#include <spatium/core/access.hpp>
#include <spatium/core/verify.hpp>
#include <spatium/spaces/euclidean.hpp>
#include <spatium/spaces/hyperbolic.hpp>
#include <spatium/spaces/implicit.hpp>
#include <spatium/spaces/parametric.hpp>
#include <spatium/spaces/product.hpp>
#include <spatium/spaces/spd.hpp>
#include <spatium/spaces/sphere.hpp>
#if SPATIUM_HAS_BOOST_MULTIPRECISION
#  include <spatium/core/precision.hpp>
#endif

#include <cmath>
#include <cstdio>
#include <limits>
#include <numbers>
#include <random>
#include <string>
#include <type_traits>
#include <vector>

namespace connectivity {

using namespace spatium;

// ── Scalars ─────────────────────────────────────────────────────

template<class T> struct is_dual : std::false_type {};
template<class U> struct is_dual<Dual<U>> : std::true_type {};

template<class T>
T from_double(double x) {
    if constexpr (is_dual<T>::value) return T(from_double<decltype(T{}.value)>(x));
    else return T(x);
}
template<class T>
double to_double(const T& x) {
    if constexpr (is_dual<T>::value) return to_double(x.value);
    else return static_cast<double>(x);
}
// The machine epsilon of the arithmetic underneath: a Dual carries its
// base type's rounding.
template<class T>
double base_epsilon() {
    if constexpr (is_dual<T>::value) return base_epsilon<decltype(T{}.value)>();
    else return static_cast<double>(std::numeric_limits<T>::epsilon());
}
// A tolerance a correct result meets: the square root of the arithmetic's
// epsilon, which is what acos- and acosh-based closed forms keep near
// coincident points, times the problem's size.
template<class T>
double tolerance(double scale = 1.0) {
    return std::sqrt(base_epsilon<T>()) * 64 * (scale + 1);
}

template<class T, std::size_t N>
Vec<T, N> cast(const Vec<double, N>& v) {
    Vec<T, N> r;
    for (std::size_t i = 0; i < N; ++i) r[i] = from_double<T>(v[i]);
    return r;
}
template<class T, std::size_t R, std::size_t C>
Matrix<T, R, C> cast(const Matrix<double, R, C>& m) {
    Matrix<T, R, C> r;
    for (std::size_t i = 0; i < R; ++i)
        for (std::size_t j = 0; j < C; ++j) r(i, j) = from_double<T>(m(i, j));
    return r;
}

// A product's points when a factor's are not Vecs.
template<class T, class A, class B>
auto cast(const spatium::product_detail::Pair<A, B>& p) {
    using CA = decltype(cast<T>(p.first));
    using CB = decltype(cast<T>(p.second));
    return spatium::product_detail::Pair<CA, CB>{cast<T>(p.first), cast<T>(p.second)};
}

// ── Spaces ──────────────────────────────────────────────────────
// Each gives the space over any scalar and sample points and tangents,
// generated once in double and cast, so every scalar sees the same inputs.

struct Rand {
    std::mt19937_64 rng{20260928};
    double operator()(double a = -1, double b = 1) { return std::uniform_real_distribution<double>(a, b)(rng); }
};

struct E3 {
    static constexpr const char* name = "E3";
    template<class T> static Euclidean<3, T> make() { return {}; }
    static std::vector<Vec<double, 3>> points(Rand& r) {
        std::vector<Vec<double, 3>> p;
        for (int i = 0; i < 6; ++i) p.push_back({r(), r(), r()});
        return p;
    }
    static Vec<double, 3> tangent(Rand& r, const Vec<double, 3>&) { return {r(), r(), r()}; }
};

struct S2 {
    static constexpr const char* name = "S2";
    template<class T> static Sphere<2, T> make() { return {}; }
    static std::vector<Vec<double, 3>> points(Rand& r) {
        std::vector<Vec<double, 3>> p;
        for (int i = 0; i < 6; ++i) {
            Vec<double, 3> x{r(), r(), r() + 1.5};   // one hemisphere and a margin: log is unique
            p.push_back(Vec<double, 3>{x * (1 / x.norm())});
        }
        return p;
    }
    static Vec<double, 3> tangent(Rand& r, const Vec<double, 3>& at) {
        const Vec<double, 3> v{r(), r(), r()};
        return Vec<double, 3>{(v - at * v.dot(at)) * 0.5};
    }
};

struct H2 {
    static constexpr const char* name = "H2";
    template<class T> static Hyperbolic<2, T> make() { return {}; }
    static std::vector<Vec<double, 3>> points(Rand& r) {
        const Hyperbolic<2> h;
        std::vector<Vec<double, 3>> p;
        for (int i = 0; i < 6; ++i) p.push_back(h.exp_map(Hyperbolic<2>::origin(), Vec<double, 3>{0, r(), r()}, 1.0));
        return p;
    }
    static Vec<double, 3> tangent(Rand& r, const Vec<double, 3>& at) {
        // Minkowski-orthogonal to `at`: <v, at>_L = -v0 a0 + v1 a1 + v2 a2 = 0.
        const double v1 = r(), v2 = r();
        return {(v1 * at[1] + v2 * at[2]) / at[0], v1, v2};
    }
};

struct S2xE1 {
    static constexpr const char* name = "S2xE1";
    template<class T> static ProductSpace<Sphere<2, T>, Euclidean<1, T>> make() { return {}; }
    static std::vector<Vec<double, 4>> points(Rand& r) {
        std::vector<Vec<double, 4>> p;
        for (const auto& s : S2::points(r)) p.push_back({s[0], s[1], s[2], r(-2, 2)});
        return p;
    }
    static Vec<double, 4> tangent(Rand& r, const Vec<double, 4>& at) {
        const auto s = S2::tangent(r, Vec<double, 3>{at[0], at[1], at[2]});
        return {s[0], s[1], s[2], r()};
    }
};

struct SPDLogE {
    static constexpr const char* name = "SPD2-logE";
    template<class T> static SPD<2, T> make() { return {}; }
    static std::vector<Vec<double, 3>> points(Rand& r) { return E3::points(r); }
    static Vec<double, 3> tangent(Rand& r, const Vec<double, 3>&) { return {r(), r(), r()}; }
};

struct SPDAff {
    static constexpr const char* name = "SPD2-affine";
    template<class T> static SPDAffineInvariant<2, T> make() { return {}; }
    static std::vector<Matrix<double, 2, 2>> points(Rand& r) {
        std::vector<Matrix<double, 2, 2>> p;
        for (int i = 0; i < 6; ++i) {
            Matrix<double, 2, 2> m;
            const double a = r(0.3, 2), c = r(0.3, 2), b = r(-0.25, 0.25);
            m(0, 0) = a; m(0, 1) = m(1, 0) = b; m(1, 1) = c;
            p.push_back(m);
        }
        return p;
    }
    static Matrix<double, 2, 2> tangent(Rand& r, const Matrix<double, 2, 2>&) {
        Matrix<double, 2, 2> m;
        m(0, 0) = r(); m(0, 1) = m(1, 0) = r(); m(1, 1) = r();
        return m;
    }
};

// Surfaces of the library's two erased kinds, over any scalar. A torus as
// a chart and a sphere as a level set.
struct TorusChart {
    static constexpr const char* name = "Torus-chart";
    template<class T> static ParametricSurface<T> make() {
        return ParametricSurface<T>(
            [](T u, T v) {
                using std::cos; using std::sin;
                const T R = from_double<T>(1.0), q = from_double<T>(0.3);
                return Vec<T, 3>{(R + q * cos(v)) * cos(u), (R + q * cos(v)) * sin(u), q * sin(v)};
            },
            {from_double<T>(0), from_double<T>(2 * std::numbers::pi), from_double<T>(0), from_double<T>(2 * std::numbers::pi)},
            true, true);
    }
    static std::vector<Vec<double, 3>> points(Rand& r) {
        std::vector<Vec<double, 3>> p;
        for (int i = 0; i < 6; ++i) {
            const double u = r(0, 0.5), v = r(0, 0.5);   // one patch, where log is local
            p.push_back({(1 + 0.3 * std::cos(v)) * std::cos(u), (1 + 0.3 * std::cos(v)) * std::sin(u), 0.3 * std::sin(v)});
        }
        return p;
    }
    static Vec<double, 3> tangent(Rand& r, const Vec<double, 3>& at) {
        const double u = std::atan2(at[1], at[0]);
        return Vec<double, 3>{Vec<double, 3>{-std::sin(u), std::cos(u), 0.0} * (0.1 * r())};
    }
};

struct SphereLevelSet {
    static constexpr const char* name = "Sphere-implicit";
    template<class T> static ImplicitSurface<T> make() {
        using std::sqrt;
        const T lo = from_double<T>(-2), hi = from_double<T>(2);
        return ImplicitSurface<T>([](T x, T y, T z) { return sqrt(x * x + y * y + z * z) - from_double<T>(1.0); },
                                  {lo, hi, lo, hi, lo, hi});
    }
    static std::vector<Vec<double, 3>> points(Rand& r) { return S2::points(r); }
    static Vec<double, 3> tangent(Rand& r, const Vec<double, 3>& at) { return S2::tangent(r, at); }
};

// ── A space with no members, over any scalar ────────────────────
// The flat cylinder R x S^1 of test_space_access.cpp, templated: exp, log
// and the metric as free functions found by ADL, the types through
// space_traits. Every cell of it goes through the customization points'
// ADL route; a cell that needs a member (the library's own verifier, a
// product, which asks for MetricSpace) shows where that route stops.

}  // namespace connectivity

namespace third_party {
template<class T> struct Cyl { T radius{2}; };
template<class T> T wrap(T a) {
    using std::fmod;
    const T pi = connectivity::from_double<T>(std::numbers::pi);
    a = fmod(a + pi, pi + pi);
    if (a < T{0}) a = a + pi + pi;
    return a - pi;
}
template<class T>
spatium::Vec<T, 2> exp_map(const Cyl<T>&, const spatium::Vec<T, 2>& p, const spatium::Vec<T, 2>& v, T t) {
    return spatium::Vec<T, 2>{p[0] + t * v[0], wrap(p[1] + t * v[1])};
}
template<class T>
spatium::Vec<T, 2> log_map(const Cyl<T>&, const spatium::Vec<T, 2>& p, const spatium::Vec<T, 2>& q) {
    return spatium::Vec<T, 2>{q[0] - p[0], wrap(q[1] - p[1])};
}
template<class T>
T metric_at(const Cyl<T>& c, const spatium::Vec<T, 2>&, const spatium::Vec<T, 2>& u, const spatium::Vec<T, 2>& v) {
    return u[0] * v[0] + c.radius * c.radius * u[1] * v[1];
}
}  // namespace third_party

template<class T>
struct spatium::spaces::space_traits<third_party::Cyl<T>> {
    using point_type = spatium::Vec<T, 2>;
    using tangent_type = spatium::Vec<T, 2>;
    using scalar_type = T;
};

namespace connectivity {

struct CylinderADL {
    static constexpr const char* name = "Cylinder-ADL";
    template<class T> static third_party::Cyl<T> make() { return {from_double<T>(2)}; }
    static std::vector<Vec<double, 2>> points(Rand& r) {
        std::vector<Vec<double, 2>> p;
        for (int i = 0; i < 6; ++i) p.push_back({r(-2, 2), r(-1, 1)});
        return p;
    }
    static Vec<double, 2> tangent(Rand& r, const Vec<double, 2>&) { return {r(), r() * 0.5}; }
};

// ── Products of spaces ──────────────────────────────────────────
// Two spaces the matrix finds green one by one must be green together: a
// product is a symmetry of the matrix. ProductSpace joins points through
// the double product itself, so a factor it cannot take fails in the
// library, not here.
template<class A, class B>
struct Prod {
    template<class T> static auto make() {
        using SA = decltype(A::template make<T>());
        using SB = decltype(B::template make<T>());
        return ProductSpace<SA, SB>{A::template make<T>(), B::template make<T>()};
    }
    static auto points(Rand& r) {
        const auto s = make<double>();
        const auto pa = A::points(r), pb = B::points(r);
        std::vector<decltype(s.join(pa[0], pb[0]))> p;
        for (std::size_t i = 0; i < pa.size() && i < pb.size(); ++i) p.push_back(s.join(pa[i], pb[i]));
        return p;
    }
    template<class Pt>
    static auto tangent(Rand& r, const Pt& at) {
        const auto s = make<double>();
        return s.join(A::tangent(r, s.first(at)), B::tangent(r, s.second(at)));
    }
};
using S2xH2 = Prod<S2, H2>;
using H2xE3 = Prod<H2, E3>;
using S2xS2 = Prod<S2, S2>;
using SPDLogExS2 = Prod<SPDLogE, S2>;
using SPDAffxS2 = Prod<SPDAff, S2>;
using CylxE3 = Prod<CylinderADL, E3>;

// ── Probes ──────────────────────────────────────────────────────
// run(space, points, tangents, signature) returns 1 or 2 and appends the
// numbers the level-3 comparison reads.

template<class S> using P = spaces::point_t<S>;
template<class S> using V = spaces::tangent_t<S>;
template<class S> using T_ = spaces::scalar_t<S>;

struct MetricAxioms {
    static constexpr const char* name = "metric-axioms";
    template<class S>
    static int run(const S& s, const std::vector<P<S>>& p, const std::vector<V<S>>&, std::vector<double>& sig) {
        using T = T_<S>;
        const double tol = tolerance<T>();
        bool ok = true;
        for (std::size_t i = 0; i < p.size(); ++i)
            for (std::size_t j = 0; j < p.size(); ++j) {
                const double dij = to_double(spaces::distance(s, p[i], p[j]));
                const double dji = to_double(spaces::distance(s, p[j], p[i]));
                if (!std::isfinite(dij)) return 0;
                if (j > i) sig.push_back(dij);
                ok = ok && dij >= -tol && std::abs(dij - dji) <= tol * (1 + dij);
                if (i == j) ok = ok && std::abs(dij) <= tol;
                for (std::size_t k = 0; k < p.size(); ++k)
                    ok = ok && dij <= to_double(spaces::distance(s, p[i], p[k])) +
                                          to_double(spaces::distance(s, p[k], p[j])) + tol;
            }
        return ok ? 2 : 1;
    }
};

// The closed form a space writes against |log_p q|_g from its own log and
// metric -- the symmetry that found the acos imprecision.
struct DerivedDistance {
    static constexpr const char* name = "distance=|log|";
    template<class S>
    static int run(const S& s, const std::vector<P<S>>& p, const std::vector<V<S>>&, std::vector<double>& sig) {
        using T = T_<S>;
        bool ok = true;
        for (std::size_t i = 0; i < p.size(); ++i)
            for (std::size_t j = i + 1; j < p.size(); ++j) {
                const auto v = V<S>{spaces::log_map(s, p[i], p[j])};
                const double derived = to_double(spaces::norm_at(s, p[i], v));
                const double closed = to_double(spaces::distance(s, p[i], p[j]));
                if (!std::isfinite(derived) || !std::isfinite(closed)) return 0;
                sig.push_back(derived);
                ok = ok && std::abs(derived - closed) <= tolerance<T>(closed);
            }
        return ok ? 2 : 1;
    }
};

struct ExpLog {
    static constexpr const char* name = "exp(log)=id";
    template<class S>
    static int run(const S& s, const std::vector<P<S>>& p, const std::vector<V<S>>&, std::vector<double>& sig) {
        using T = T_<S>;
        bool ok = true;
        for (std::size_t i = 0; i < p.size(); ++i)
            for (std::size_t j = 0; j < p.size(); ++j) {
                if (i == j) continue;
                const auto v = V<S>{spaces::log_map(s, p[i], p[j])};
                const auto q = P<S>{spaces::exp_map(s, p[i], v, from_double<T>(1))};
                const double err = to_double(spaces::distance(s, q, p[j]));
                if (!std::isfinite(err)) return 0;
                sig.push_back(to_double(spaces::distance(s, p[i], q)));
                ok = ok && err <= tolerance<T>();
            }
        return ok ? 2 : 1;
    }
};

struct Midpoint {
    static constexpr const char* name = "midpoint";
    template<class S>
    static int run(const S& s, const std::vector<P<S>>& p, const std::vector<V<S>>&, std::vector<double>& sig) {
        using T = T_<S>;
        bool ok = true;
        for (std::size_t i = 0; i + 1 < p.size(); ++i) {
            const auto m = spaces::midpoint(s, p[i], p[i + 1]);
            const double a = to_double(spaces::distance(s, p[i], m)), b = to_double(spaces::distance(s, m, p[i + 1]));
            const double d = to_double(spaces::distance(s, p[i], p[i + 1]));
            if (!std::isfinite(a) || !std::isfinite(b)) return 0;
            sig.push_back(a);
            ok = ok && std::abs(a - b) <= tolerance<T>(d) && std::abs(a - d / 2) <= tolerance<T>(d);
        }
        return ok ? 2 : 1;
    }
};

// The mean of a point and its reflections through a centre is the centre.
struct FrechetMean {
    static constexpr const char* name = "frechet-mean";
    template<class S>
    static int run(const S& s, const std::vector<P<S>>& p, const std::vector<V<S>>&, std::vector<double>& sig) {
        using T = T_<S>;
        const auto& c = p[0];
        // Tangents at c, whatever the space: log from c to other samples.
        std::vector<P<S>> pts;
        for (std::size_t k = 1; k < 4 && k < p.size(); ++k) {
            const auto w = V<S>{spaces::log_map(s, c, p[k])};
            pts.push_back(P<S>{spaces::exp_map(s, c, w, from_double<T>(1))});
            pts.push_back(P<S>{spaces::exp_map(s, c, w, from_double<T>(-1))});
        }
        const auto m = algebra::frechet_mean(s, pts, pts.front(), from_double<T>(tolerance<T>() * 1e-3), 200);
        const double err = to_double(spaces::distance(s, m, c));
        if (!std::isfinite(err)) return 0;
        sig.push_back(err);
        return err <= tolerance<T>() * 10 ? 2 : 1;
    }
};

// The library's own verifier, which asks for the member-based concepts.
struct VerifyExpLog {
    static constexpr const char* name = "verify_exp_log";
    template<class S>
    static int run(const S& s, const std::vector<P<S>>& p, const std::vector<V<S>>&, std::vector<double>& sig) {
        using T = T_<S>;
        const auto r = verify_exp_log(s, std::span<const P<S>>(p), from_double<T>(tolerance<T>()));
        sig.push_back(r.passed ? 1.0 : 0.0);
        return r.passed ? 2 : 1;
    }
};

// d/dt of the distance from p to exp_q(t v) at t = 0: through the scalar's
// own derivative when it is a Dual, against a central difference in double.
// For any other scalar the cell reports the difference quotient alone.
struct Derivative {
    static constexpr const char* name = "d/dt distance";
    template<class S>
    static int run(const S& s, const std::vector<P<S>>& p, const std::vector<V<S>>& v, std::vector<double>& sig) {
        using T = T_<S>;
        const auto dist_at = [&](const T& t, std::size_t i) {
            return spaces::distance(s, p[0], P<S>{spaces::exp_map(s, p[i], v[i], t)});
        };
        bool ok = true;
        for (std::size_t i = 1; i < p.size(); ++i) {
            // The step that balances truncation against rounding.
            const double h = std::cbrt(base_epsilon<T>());
            // The difference taken in the cell's own arithmetic: over Real50
            // it is below double's resolution of the distances themselves.
            const double fd = to_double(T{(dist_at(from_double<T>(h), i) - dist_at(from_double<T>(-h), i)) / from_double<T>(2 * h)});
            if (!std::isfinite(fd)) return 0;
            if constexpr (is_dual<T>::value) {
                T t = from_double<T>(0);
                t.deriv = from_double<decltype(T{}.deriv)>(1);
                const double exact = to_double(dist_at(t, i).deriv);
                if (!std::isfinite(exact)) return 0;
                sig.push_back(exact);
                ok = ok && std::abs(exact - fd) <= 1e-5 * (1 + std::abs(fd));
            } else {
                sig.push_back(fd);
            }
        }
        return ok ? 2 : 1;
    }
};

// ── A cell ──────────────────────────────────────────────────────

template<class T, class Space, class Probe>
int run_cell(std::vector<double>& sig) {
    const auto s = Space::template make<T>();
    Rand r;
    const auto pd = Space::points(r);
    std::vector<P<decltype(s)>> p;
    std::vector<V<decltype(s)>> v;
    // Points made in double lie on a curved space only to double's
    // precision; in a finer arithmetic they are put back on it first, or
    // every cell over Real50 would be measuring double's rounding.
    for (const auto& x : pd) {
        auto q = cast<T>(x);
        if constexpr (requires { spaces::project(s, q); }) q = P<decltype(s)>{spaces::project(s, q)};
        p.push_back(q);
    }
    for (const auto& x : pd) v.push_back(cast<T>(Space::tangent(r, x)));
    return Probe::run(s, p, v, sig);
}

// One line a cell prints and the script reads: level, then the signature.
inline void print_cell(const char* scalar, const char* space, const char* probe, int level,
                       const std::vector<double>& sig) {
    std::printf("%s|%s|%s|%d|", scalar, space, probe, level);
    for (double x : sig) std::printf("%.17g ", x);
    std::printf("\n");
}

}  // namespace connectivity
