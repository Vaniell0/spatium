#pragma once
// A metric written in the scene's own field language: ten scalar fields of
// the spacetime point (`Field::coord(0..3)`), one per independent entry of
// the symmetric g_{mu nu}, each lowered to plain data (io/field_pod.hpp).
//
// What that buys is that a metric stops being code. `KerrMetric` in
// kerr.hpp is a C++ callable, so changing coordinates -- Boyer-Lindquist,
// singular on the spin axis, to Kerr-Schild, which is not -- means
// rewriting it and every kernel transcribed from it. Here it means another
// ten fields. The same data lowers to GLSL for the device, serialises
// into a scene description, and is differentiated by evaluating it on
// `Dual<T>`: `operator()` is templated on the scalar exactly as
// `KerrMetric`'s is, so geodesic.hpp's Christoffel symbols and
// `geodesic_step` take a MetricField unchanged.
//
// On the CPU the interpreted metric costs more than a hand-written one:
// frame 40 of the since-removed blackhole_gr_demo at 160x90 took 5.8 s through the Boyer-
// Lindquist fields against 3.1 s through KerrMetric, pixel for pixel the
// same image; before the ten entries shared one pool it was 16 s. On the
// device the same data is compiled, not interpreted.
//
// A metric reads only arithmetic, Min/Max, Sin/Cos/Sqrt, Less and Coord.
// Noise, gathers and the per-instance inputs have no meaning at a point of
// spacetime and are refused when the field is made, with the op named.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/matrix.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/error.hpp>
#  include <spatium/io/field.hpp>
#  include <spatium/io/field_pod.hpp>
#  include <spatium/physics/relativity/geodesic.hpp>
#  include <array>
#  include <cmath>
#  include <cstring>
#  include <format>
#  include <limits>
#  include <map>
#  include <tuple>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::physics::relativity {

namespace detail {

inline bool metric_op(io::build::Op o) {
    using io::build::Op;
    switch (o) {
        case Op::Const: case Op::Add: case Op::Sub: case Op::Mul: case Op::Div:
        case Op::Min: case Op::Max: case Op::Sin: case Op::Cos: case Op::Sqrt:
        case Op::Less: case Op::Coord:
            return true;
        default:
            return false;
    }
}

// The POD interpreter's switch, on any scalar S: what makes a metric's
// derivatives exact is running this with S = Dual<T>.
template<typename S, Scalar T>
void evaluate_as(const std::vector<io::build::PodOp<T>>& ops, const Vec<S, 4>& x, S* s) {
    using io::build::Op;
    const auto& f = ops;
    const auto count = static_cast<std::uint32_t>(f.size());
    for (std::uint32_t i = 0; i < count; ++i) {
        const auto& n = f[i];
        switch (static_cast<Op>(n.code)) {
            case Op::Const: s[i] = S(n.value); break;
            case Op::Add:   s[i] = s[n.a] + s[n.b]; break;
            case Op::Sub:   s[i] = s[n.a] - s[n.b]; break;
            case Op::Mul:   s[i] = s[n.a] * s[n.b]; break;
            case Op::Div:   s[i] = s[n.a] / s[n.b]; break;
            case Op::Min:   s[i] = s[n.a] < s[n.b] ? s[n.a] : s[n.b]; break;
            case Op::Max:   s[i] = s[n.a] < s[n.b] ? s[n.b] : s[n.a]; break;
            case Op::Sin:   { using std::sin; s[i] = sin(s[n.a]); break; }
            case Op::Cos:   { using std::cos; s[i] = cos(s[n.a]); break; }
            case Op::Sqrt:  { using std::sqrt; s[i] = sqrt(s[n.a]); break; }
            case Op::Less:  s[i] = s[n.a] < s[n.b] ? S(T{1}) : S(T{0}); break;
            case Op::Coord: s[i] = x[n.k]; break;
            default:        s[i] = S(T{0}); break;   // refused by MetricField::make
        }
    }
}

}  // namespace detail

template<Scalar T = double>
class MetricField {
public:
    // The independent entries, in this order.
    static constexpr std::array<std::pair<int, int>, 10> kEntries{{
        {0, 0}, {0, 1}, {0, 2}, {0, 3}, {1, 1}, {1, 2}, {1, 3}, {2, 2}, {2, 3}, {3, 3}}};

    static Result<MetricField> make(const std::array<io::build::Field<T>, 10>& g) {
        MetricField m;
        for (std::size_t e = 0; e < 10; ++e) {
            auto pod = io::build::lower(g[e]);
            if (!pod) return std::unexpected(std::move(pod.error()));
            for (std::size_t i = 0; i < pod->ops.size(); ++i) {
                const auto op = static_cast<io::build::Op>(pod->ops[i].code);
                if (!detail::metric_op(op))
                    return std::unexpected(Error(ErrorCode::InvalidArgument,
                        std::format("MetricField: g{}{} op {} is {}, which a metric cannot read",
                                    kEntries[e].first, kEntries[e].second, i, io::build::op_name(op))));
            }
            m.g_[e] = std::move(*pod);
        }
        m.build_pool();
        return m;
    }

    template<Scalar S>
    Matrix<S, 4, 4> operator()(const Vec<S, 4>& x) const {
        // Scratch on the stack for every metric in this file; a larger one
        // takes the heap.
        std::array<S, 256> local;
        std::vector<S> heap;
        S* s = local.data();
        if (pool_.size() > local.size()) { heap.resize(pool_.size()); s = heap.data(); }
        detail::evaluate_as(pool_, x, s);
        Matrix<S, 4, 4> out{};
        for (std::size_t e = 0; e < 10; ++e) {
            const auto [i, j] = kEntries[e];
            out(i, j) = s[root_[e]];
            out(j, i) = s[root_[e]];
        }
        return out;
    }

    // Each entry lowered on its own, as written.
    const std::array<io::build::PodField<T>, 10>& components() const { return g_; }
    // All ten in one pool with every repeated op kept once, and where each
    // entry's value lands in it. What operator() evaluates.
    const std::vector<io::build::PodOp<T>>& pool() const { return pool_; }
    // Where entry e's value lands in the pool.
    std::uint32_t root(std::size_t e) const { return root_[e]; }

private:
    // The ten entries share most of their arithmetic -- Kerr's Sigma, sin
    // and cos of theta appear in several -- and each field owns its own
    // copy of it, since a Field is a value. Merging them into one pool,
    // keeping an op once when its code, operands and immediate are the same
    // as one already there, evaluates each shared subexpression once; it
    // is the value numbering the GLSL emitter does, done on the data.
    void build_pool() {
        using Key = std::tuple<std::uint8_t, std::uint32_t, std::uint32_t, std::uint32_t, std::uint32_t,
                               std::uint64_t>;
        std::map<Key, std::uint32_t> seen;
        for (std::size_t e = 0; e < 10; ++e) {
            const auto& ops = g_[e].ops;
            std::vector<std::uint32_t> at(ops.size());
            for (std::size_t i = 0; i < ops.size(); ++i) {
                auto n = ops[i];
                const auto k = io::build::arity(static_cast<io::build::Op>(n.code));
                if (k >= 1) n.a = at[n.a];
                if (k >= 2) n.b = at[n.b];
                if (k >= 3) n.c = at[n.c];
                std::uint64_t bits = 0;
                if constexpr (sizeof(T) == sizeof(std::uint64_t)) std::memcpy(&bits, &n.value, sizeof bits);
                else { const double v = static_cast<double>(n.value); std::memcpy(&bits, &v, sizeof bits); }
                const Key key{n.code, k >= 1 ? n.a : 0u, k >= 2 ? n.b : 0u, k >= 3 ? n.c : 0u, n.k, bits};
                auto [it, fresh] = seen.try_emplace(key, static_cast<std::uint32_t>(pool_.size()));
                if (fresh) pool_.push_back(n);
                at[i] = it->second;
            }
            root_[e] = at.back();
        }
    }

    std::array<io::build::PodField<T>, 10> g_{};
    std::vector<io::build::PodOp<T>> pool_;
    std::array<std::uint32_t, 10> root_{};
};

// ── Two forms of the Kerr metric, written as fields ─────────────

// Boyer-Lindquist (t, r, theta, phi): the form of kerr.hpp's KerrMetric,
// entry for entry. Singular on the spin axis, where 1/sin(theta) terms in
// the Christoffel symbols blow up.
template<Scalar T>
MetricField<T> kerr_boyer_lindquist(T mass, T spin) {
    using F = io::build::Field<T>;
    const F r = F::coord(1), theta = F::coord(2);
    const F M(mass), a(spin);
    const F c = cos(theta), s = sin(theta);
    const F sigma = r * r + a * a * c * c;
    const F delta = r * r - F(T{2}) * M * r + a * a;
    const F zero(T{0});
    std::array<F, 10> g{
        F(T{0}) - (F(T{1}) - F(T{2}) * M * r / sigma),               // tt
        zero, zero,                                                  // tr, ttheta
        F(T{0}) - F(T{2}) * M * r * a * s * s / sigma,               // tphi
        sigma / delta,                                               // rr
        zero, zero,                                                  // rtheta, rphi
        sigma,                                                       // theta theta
        zero,                                                        // theta phi
        (r * r + a * a + F(T{2}) * M * r * a * a * s * s / sigma) * s * s};   // phi phi
    return *MetricField<T>::make(g);
}

// Kerr-Schild, Cartesian (t, x, y, z): g = eta + f l (x) l, with eta the
// Minkowski metric (-,+,+,+), l the null covector
//   l = (1, (r x + a y)/(r^2 + a^2), (r y - a x)/(r^2 + a^2), z / r),
//   f = 2 M r^3 / (r^4 + a^2 z^2),
// and r the Boyer-Lindquist radius, the root of
//   r^4 - (rho^2 - a^2) r^2 - a^2 z^2 = 0,  rho^2 = x^2 + y^2 + z^2.
// Regular on the spin axis and across the horizon. The form in which two
// holes can be added -- eta plus one such term per hole -- which is what
// the superposed binary is built on.
template<Scalar T>
MetricField<T> kerr_schild(T mass, T spin) {
    using F = io::build::Field<T>;
    const F x = F::coord(1), y = F::coord(2), z = F::coord(3);
    const F M(mass), a(spin), one(T{1}), two(T{2}), half(T{0.5}), quarter(T{0.25});
    const F rho2 = x * x + y * y + z * z;
    const F w = rho2 - a * a;
    const F r2 = half * w + sqrt(quarter * w * w + a * a * z * z);
    const F r = sqrt(r2);
    const F den = r2 + a * a;
    const std::array<F, 4> l{one, (r * x + a * y) / den, (r * y - a * x) / den, z / r};
    const F f = two * M * r2 * r / (r2 * r2 + a * a * z * z);
    std::array<F, 10> g;
    for (std::size_t e = 0; e < 10; ++e) {
        const auto [i, j] = MetricField<T>::kEntries[e];
        const T eta = i != j ? T{0} : (i == 0 ? T{-1} : T{1});
        g[e] = F(eta) + f * l[static_cast<std::size_t>(i)] * l[static_cast<std::size_t>(j)];
    }
    return *MetricField<T>::make(g);
}

// ── Is the metric a spacetime here ───────────────────────────────
//
// det g < 0 with g invertible: the one-negative-eigenvalue signature a
// 4x4 symmetric metric needs to be Lorentzian (with det < 0 there is an
// odd number of negative eigenvalues; three negative would need g_00 > 0
// on every timelike direction, which no metric here has). A superposition
// of Kerr-Schild terms loses it between two holes: measured for equal
// masses, det g >= 0 at d <= 4 M, det within [-0.39, -0.04] at 5 M. Checked
// before anything inverts g, since christoffel() does not refuse.
template<Scalar T, typename Metric>
bool lorentzian_at(const Metric& metric, const Vec<T, 4>& x) {
    const Matrix<T, 4, 4> g = metric(x);
    const auto inv = invert(g);
    return inv.has_value() && g.determinant() < T{0};
}

// ── How far a metric is from solving the vacuum equations ────────
//
// The largest |R_{mu nu}| at an event, the Ricci tensor built from the
// Christoffel symbols and their derivatives by central differences of
// step h:
//   R_{mu nu} = d_l G^l_{mu nu} - d_nu G^l_{mu l} + G^l_{l s} G^s_{mu nu} - G^l_{nu s} G^s_{mu l}.
// Zero for a vacuum solution up to O(h^2) of truncation: 1.2e-6 for
// Boyer-Lindquist Kerr and 4.4e-10 for Kerr-Schild at h = 1e-4 over
// tests/test_metric_field.cpp's points, against 6.1e-2 for a metric with a
// source. The measure of a superposition's error, which is not zero.
template<Scalar T, typename Metric>
T vacuum_residual_fd(const Metric& metric, const Vec<T, 4>& x, T h = T{1e-4}) {
    using std::abs;
    const auto G = christoffel(metric, x);
    std::array<std::array<Matrix<T, 4, 4>, 4>, 4> dG{};   // dG[k][l](mu,nu) = d_k G^l_{mu nu}
    for (std::size_t k = 0; k < 4; ++k) {
        Vec<T, 4> xp = x, xm = x;
        xp[k] += h;
        xm[k] -= h;
        const auto Gp = christoffel(metric, xp), Gm = christoffel(metric, xm);
        for (std::size_t l = 0; l < 4; ++l)
            for (std::size_t mu = 0; mu < 4; ++mu)
                for (std::size_t nu = 0; nu < 4; ++nu)
                    dG[k][l](mu, nu) = (Gp[l](mu, nu) - Gm[l](mu, nu)) / (T{2} * h);
    }
    T worst{0};
    for (std::size_t mu = 0; mu < 4; ++mu)
        for (std::size_t nu = 0; nu < 4; ++nu) {
            T v{0};
            for (std::size_t l = 0; l < 4; ++l) {
                v += dG[l][l](mu, nu) - dG[nu][l](mu, l);
                for (std::size_t s = 0; s < 4; ++s)
                    v += G[l](l, s) * G[s](mu, nu) - G[l](nu, s) * G[s](mu, l);
            }
            worst = std::max(worst, abs(v));
        }
    return worst;
}

// The same residual with the derivatives of the Christoffel symbols taken
// exactly: christoffel() run on Dual<T>, which metric_derivatives() wraps
// once more, so the metric is evaluated on Dual<Dual<T>> and the second
// derivatives of g come out of the arithmetic rather than a step size.
// Nothing in the metric or in geodesic.hpp changes for it. Measured when it
// was written: Kerr (Kerr-Schild and Boyer-Lindquist) 1e-16 to 1e-18 against
// vacuum_residual_fd's 1e-9 to 1e-12, at the same cost per event (~29 us on
// the binary's 155-op pool), and the Schwarzschild Kretschmann scalar to 12
// digits. This is the oracle a regime choice is graded against.
//
// At a point where the metric does not invert -- a degenerate signature,
// which a superposition of two holes reaches between them at d <= ~4.5 M --
// christoffel() substitutes a zero inverse and the curvature would read as
// 0. So the inverse is checked first and a NaN is returned instead: an
// answer nobody can mistake for "vacuum".
template<Scalar T, typename Metric>
T vacuum_residual(const Metric& metric, const Vec<T, 4>& x) {
    using std::abs;
    if (!lorentzian_at(metric, x)) return std::numeric_limits<T>::quiet_NaN();
    std::array<Matrix<T, 4, 4>, 4> G{};
    std::array<std::array<Matrix<T, 4, 4>, 4>, 4> dG{};   // dG[k][l](mu,nu) = d_k G^l_{mu nu}
    for (std::size_t k = 0; k < 4; ++k) {
        Vec<Dual<T>, 4> xd;
        for (std::size_t j = 0; j < 4; ++j)
            xd[j] = (j == k) ? Dual<T>::variable(x[j]) : Dual<T>::constant(x[j]);
        const auto Gd = christoffel(metric, xd);
        for (std::size_t l = 0; l < 4; ++l)
            for (std::size_t mu = 0; mu < 4; ++mu)
                for (std::size_t nu = 0; nu < 4; ++nu) {
                    if (k == 0) G[l](mu, nu) = Gd[l](mu, nu).value;
                    dG[k][l](mu, nu) = Gd[l](mu, nu).deriv;
                }
    }
    T worst{0};
    for (std::size_t mu = 0; mu < 4; ++mu)
        for (std::size_t nu = 0; nu < 4; ++nu) {
            T v{0};
            for (std::size_t l = 0; l < 4; ++l) {
                v += dG[l][l](mu, nu) - dG[nu][l](mu, l);
                for (std::size_t s = 0; s < 4; ++s)
                    v += G[l](l, s) * G[s](mu, nu) - G[l](nu, s) * G[s](mu, l);
            }
            worst = std::max(worst, abs(v));
        }
    return worst;
}

} // namespace spatium::physics::relativity
