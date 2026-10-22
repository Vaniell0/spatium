#pragma once

// Gravity of many sources, in N dimensions, free or on a torus.
//
// `PointGravity` (force.hpp) is one source, no softening, a fixed law. A
// field of hundreds of bodies pulling each other needs more, and the
// choices are made here, once:
//
//   law         a source of mass m pulls with acceleration G m r / (r^2 + eps^2)^((p+1)/2),
//               p the integer `exponent`. p = N-1 is the field of a point
//               source in N dimensions (Poisson's equation: 1/r in 2D,
//               1/r^2 in 3D); p = 1 in 3D is the softer 1/r law. eps is
//               the Plummer softening, so two coincident bodies do not blow up.
//   clamp       an optional cap on the length of the total acceleration. A
//               clamp breaks conservation of energy and momentum: it is a
//               safety net for close encounters, not part of the physics.
//   coupling    the ratio of a body's gravitational charge to its inertial
//               mass, per body. A source's `mass` is its charge; a body with
//               ratio c and inertia M is pulled by c * g(x) and should be
//               entered as a source of mass c * M, which keeps the third law
//               (sum of M v) exact. `GravityPull` applies one ratio to the
//               existing integrators without changing them.
//   periodic    a torus: a box that wraps. Poisson's equation on a closed
//               space has no solution unless the net source is zero, so the
//               field is that of the sources plus a uniform background of
//               opposite total mass. It is summed by Ewald's method: the
//               kernel is split by a Gaussian at scale 1/alpha into a short
//               part summed over the nearest images in real space and a long
//               part summed over reciprocal vectors, both converged to about
//               exp(-eta^2), with k = 0 (the background) left out. The split
//               is exact for every integer exponent 1 <= p <= N-1, softened
//               or not: the Plummer factor exp(-t eps^2) goes into the
//               Gaussian representation of the kernel and the reciprocal
//               coefficient is a one-dimensional integral per vector,
//               computed once. The result does not depend on alpha.
//
// `acceleration(x)` is the direct sum over all sources, O(sources) per
// point; `gravity_tree.hpp` is the O(log) per point version on a BVH. Both
// share `kernel_real` and `reciprocal_field`, so they differ only in how
// the short-range part is summed.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/physics/mechanics/body.hpp>
#  include <algorithm>
#  include <array>
#  include <cmath>
#  include <complex>
#  include <cstddef>
#  include <cstdint>
#  include <limits>
#  include <numbers>
#  include <stdexcept>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::physics::mechanics {

// One source: where it is and its gravitational charge.
template<std::size_t N, Scalar T = double>
struct GravitySource {
    Vec<T, N> position{};
    T         mass{};
};

template<std::size_t N, Scalar T = double>
struct GravityLaw {
    unsigned exponent = N - 1;   // acceleration ~ 1 / r^exponent far from the source
    T G{1};
    T softening{0};              // Plummer epsilon
    T max_acceleration{0};       // 0: no clamp

    static constexpr GravityLaw newtonian(T G = T{1}, T eps = T{0}) {
        return GravityLaw{N - 1, G, eps, T{0}};
    }
    // 1/r acceleration in any dimension: the softer law (N = 3: "inverse distance").
    static constexpr GravityLaw inverse_distance(T G = T{1}, T eps = T{0}) {
        return GravityLaw{1, G, eps, T{0}};
    }
};

// Ewald tuning for a periodic field. `alpha_L` is alpha times the shortest box
// side; `eta` is the number of Gaussian widths both sums run to (the dropped
// tail is about exp(-eta^2)). The real-space cutoff eta / alpha must not
// exceed the shortest side, so only the nearest images are summed.
struct EwaldOptions {
    double alpha_L = 6.0;
    double eta     = 4.4;
};

namespace gravity_detail {

// Regularized upper incomplete gamma Q(two_a / 2, x), two_a >= 1 an integer.
template<Scalar T>
T gamma_q(unsigned two_a, T x) {
    using std::exp; using std::erfc; using std::sqrt; using std::pow; using std::tgamma;
    T q;
    T a;
    if (two_a % 2 == 1) { q = erfc(sqrt(x)); a = T{0.5}; }
    else                { q = exp(-x);       a = T{1}; }
    const T ex = exp(-x);
    for (unsigned two = (two_a % 2 == 1) ? 1 : 2; two < two_a; two += 2) {
        q += pow(x, a) * ex / tgamma(a + T{1});
        a += T{1};
    }
    return q;
}

// Exponential integral E1(x), x > 0.
template<Scalar T>
T expint_e1(T x) {
    using std::exp; using std::log;
    if (x < T{1}) {
        T sum{}, term{1};
        for (int k = 1; k < 60; ++k) {
            term *= -x / T(k);
            sum -= term / T(k);
        }
        return -T{0.5772156649015329} - log(x) + sum;
    }
    // Continued fraction (modified Lentz).
    T b = x + T{1}, c = T{1} / T{1e-300}, d = T{1} / b, h = d;
    for (int i = 1; i < 200; ++i) {
        const T an = -T(i) * T(i);
        b += T{2};
        d = T{1} / (an * d + b);
        c = b + an / c;
        const T del = c * d;
        h *= del;
        if (std::abs(del - T{1}) < T{1e-16}) break;
    }
    return h * exp(-x);
}

// 12-point Gauss-Legendre nodes and weights on [-1, 1].
inline const std::array<std::array<double, 2>, 12>& gauss_legendre_12() {
    static const auto table = [] {
        std::array<std::array<double, 2>, 12> t{};
        constexpr int n = 12;
        for (int i = 0; i < n / 2; ++i) {
            double x = std::cos(std::numbers::pi * (i + 0.75) / (n + 0.5));
            double dp = 1;
            for (int it = 0; it < 100; ++it) {
                double p0 = 1, p1 = x;
                for (int k = 2; k <= n; ++k) {
                    const double p2 = ((2 * k - 1) * x * p1 - (k - 1) * p0) / k;
                    p0 = p1; p1 = p2;
                }
                dp = n * (x * p1 - p0) / (x * x - 1);
                const double dx = p1 / dp;
                x -= dx;
                if (std::abs(dx) < 1e-16) break;
            }
            const double w = 2 / ((1 - x * x) * dp * dp);
            t[i]         = {-x, w};
            t[n - 1 - i] = { x, w};
        }
        return t;
    }();
    return table;
}

} // namespace gravity_detail

template<std::size_t N, Scalar T = double>
class GravityField {
    static_assert(N >= 2, "GravityField needs at least two dimensions");
public:
    using ScalarType = T;
    static constexpr std::size_t dimension = N;
    using Law = GravityLaw<N, T>;
    using Source = GravitySource<N, T>;
    static constexpr std::size_t npos = ~std::size_t{0};

    // Free space.
    explicit GravityField(Law law = Law{}) : law_(law) { init_law(); }

    // A torus with sides `box`.
    GravityField(Law law, Vec<T, N> box, EwaldOptions ewald = {})
        : law_(law), box_(box), periodic_(true), ewald_(ewald) {
        init_law();
        init_periodic();
    }

    const Law& law() const { return law_; }
    bool periodic() const { return periodic_; }
    const Vec<T, N>& box() const { return box_; }
    const EwaldOptions& ewald() const { return ewald_; }
    // Real-space cutoff of the short-range part (periodic); infinite otherwise.
    T cutoff() const { return rc_; }
    T alpha() const { return alpha_; }
    std::size_t reciprocal_vectors() const { return kvecs_.size(); }

    // Replace the sources. On a torus, positions are wrapped into [0, L) and
    // the structure factor of the reciprocal sum is rebuilt: O(sources x vectors).
    void set_sources(std::vector<Source> sources) {
        sources_ = std::move(sources);
        if (periodic_) {
            for (auto& s : sources_) s.position = wrap(s.position);
            build_structure_factors();
        }
    }
    const std::vector<Source>& sources() const { return sources_; }

    Vec<T, N> wrap(Vec<T, N> x) const {
        if (!periodic_) return x;
        for (std::size_t i = 0; i < N; ++i) {
            x[i] -= box_[i] * std::floor(x[i] / box_[i]);
            if (x[i] >= box_[i]) x[i] = T{0};
        }
        return x;
    }

    // ── Kernel pieces (shared with the tree) ──────────────────────

    // The short-range part of the field kernel at displacement d = x - source:
    // d (d^2 + eps^2)^(-(p+1)/2), times Q(q, alpha^2 (d^2 + eps^2)) on a torus,
    // and zero beyond the cutoff. The pull is -G m times this.
    Vec<T, N> kernel_real(const Vec<T, N>& d) const {
        const T r2 = d.dot(d);
        if (periodic_ && r2 >= rc2_) return Vec<T, N>{};
        const T s = r2 + eps2_;
        if (s < T{1e-300}) return Vec<T, N>{};
        T g = inv_power(s);
        if (periodic_) g *= gravity_detail::gamma_q<T>(law_.exponent + 1, alpha2_ * s);
        return Vec<T, N>{d * g};
    }

    // Potential counterpart of kernel_real: Phi(d) with grad Phi = kernel.
    T potential_real(const Vec<T, N>& d) const {
        const T r2 = d.dot(d);
        if (periodic_ && r2 >= rc2_) return T{0};
        const T s = r2 + eps2_;
        if (s < T{1e-300}) return T{0};
        if (!periodic_) return free_potential(s);
        const T x = alpha2_ * s;
        if (law_.exponent == 1) return -gravity_detail::expint_e1<T>(x) / T{2};
        const T q = T(law_.exponent + 1) / T{2};
        const T gq = std::tgamma(q);
        const T inc = gq * gravity_detail::gamma_q<T>(law_.exponent + 1, x);
        const T I = (std::pow(s, T{1} - q) * inc - std::pow(alpha2_, q - T{1}) * std::exp(-x))
                  / (q - T{1});
        return -I / (T{2} * gq);
    }

    // The long-range part of  sum_j m_j kernel(x - s_j)  on a torus (0 in free space).
    Vec<T, N> reciprocal_field(const Vec<T, N>& x) const {
        Vec<T, N> out{};
        if (!periodic_ || kvecs_.empty()) return out;
        thread_local Tables tab;
        phase_tables(x, tab);
        const T scale = T{2} / volume_;
        for (std::size_t v = 0; v < kvecs_.size(); ++v) {
            const auto& kv = kvecs_[v];
            const std::complex<T> ph = phase(tab, kv);
            const T f = kv.coef * (ph * structure_[v]).imag() * scale;
            for (std::size_t i = 0; i < N; ++i) out[i] += kv.k[i] * f;
        }
        return out;
    }

    // Sum over sources of m_j Phi_per(x - s_j), up to an additive constant that
    // does not depend on any position. Potential energy of an isolated pair set
    // is G/2 sum_i m_i potential(x_i, i).
    T potential(const Vec<T, N>& x, std::size_t skip = npos) const {
        T phi{};
        for (std::size_t j = 0; j < sources_.size(); ++j) {
            if (j == skip) continue;
            phi += sources_[j].mass * pair_potential_real(x - sources_[j].position);
        }
        if (periodic_ && !kvecs_.empty()) {
            thread_local Tables tab;
            phase_tables(x, tab);
            T rec{};
            for (std::size_t v = 0; v < kvecs_.size(); ++v) {
                const auto& kv = kvecs_[v];
                rec += kv.coef * (phase(tab, kv) * structure_[v]).real();
            }
            phi -= T{2} / volume_ * rec;
        }
        return phi;
    }

    // Total potential energy, up to a constant independent of positions.
    T potential_energy() const {
        T u{};
        for (std::size_t i = 0; i < sources_.size(); ++i)
            u += sources_[i].mass * potential(sources_[i].position, i);
        return T{0.5} * law_.G * u;
    }

    // ── The field ─────────────────────────────────────────────────

    // Acceleration per unit coupling at x: direct sum over all sources, the
    // source `skip` (the body asking, if it is one) excluded. Clamped if the
    // law says so.
    Vec<T, N> acceleration(const Vec<T, N>& x, std::size_t skip = npos) const {
        return clamped(acceleration_unclamped(x, skip));
    }

    Vec<T, N> acceleration_unclamped(const Vec<T, N>& x, std::size_t skip = npos) const {
        Vec<T, N> sum{};
        for (std::size_t j = 0; j < sources_.size(); ++j) {
            if (j == skip) continue;
            add_images(sum, x - sources_[j].position, sources_[j].mass);
        }
        sum = Vec<T, N>{sum + reciprocal_field(x)};
        return Vec<T, N>{sum * (-law_.G)};
    }

    // Apply the law's clamp to a total acceleration.
    Vec<T, N> clamped(Vec<T, N> a) const {
        if (law_.max_acceleration > T{0}) {
            const T n2 = a.dot(a);
            const T m2 = law_.max_acceleration * law_.max_acceleration;
            if (n2 > m2) a = Vec<T, N>{a * (law_.max_acceleration / std::sqrt(n2))};
        }
        return a;
    }

    // Short-range kernel summed over the images of one source that reach x:
    // adds m * kernel to `sum` (minimum image first, then the neighbours within the cutoff).
    void add_images(Vec<T, N>& sum, Vec<T, N> d, T m) const {
        if (!periodic_) {
            sum = Vec<T, N>{sum + kernel_real(d) * m};
            return;
        }
        std::array<std::array<T, 2>, N> opt{};
        std::array<int, N> cnt{};
        for (std::size_t i = 0; i < N; ++i) {
            d[i] -= box_[i] * std::round(d[i] / box_[i]);
            opt[i][0] = d[i];
            cnt[i] = 1;
            const T other = d[i] > T{0} ? d[i] - box_[i] : d[i] + box_[i];
            if (std::abs(other) < rc_) { opt[i][1] = other; cnt[i] = 2; }
        }
        std::array<int, N> pick{};
        for (;;) {
            Vec<T, N> dd;
            for (std::size_t i = 0; i < N; ++i) dd[i] = opt[i][pick[i]];
            sum = Vec<T, N>{sum + kernel_real(dd) * m};
            std::size_t i = 0;
            for (; i < N; ++i) {
                if (++pick[i] < cnt[i]) break;
                pick[i] = 0;
            }
            if (i == N) break;
        }
    }

private:
    struct KVec {
        std::array<int, N> n{};   // integer index per axis
        std::array<T, N> k{};     // wave vector
        T coef{};                 // C(k)
    };
    using Tables = std::array<std::vector<std::complex<T>>, N>;

    Law law_;
    Vec<T, N> box_{};
    bool periodic_ = false;
    EwaldOptions ewald_{};
    std::vector<Source> sources_;

    T eps2_{}, alpha_{}, alpha2_{}, rc_ = std::numeric_limits<T>::infinity(), rc2_{}, volume_{};
    std::array<int, N> nmax_{};
    std::vector<KVec> kvecs_;                    // half space: one of each +-k pair
    std::vector<std::complex<T>> structure_;     // S(k) = sum m_j exp(-i k . s_j)

    void init_law() {
        if (law_.exponent < 1)
            throw std::invalid_argument("GravityField: exponent must be at least 1");
        if (law_.softening < T{0})
            throw std::invalid_argument("GravityField: negative softening");
        eps2_ = law_.softening * law_.softening;
    }

    void init_periodic() {
        if (law_.exponent > N - 1)
            throw std::invalid_argument(
                "GravityField: a periodic field needs exponent <= N-1 (the background sum diverges otherwise)");
        T lmin = box_[0];
        volume_ = T{1};
        for (std::size_t i = 0; i < N; ++i) {
            if (!(box_[i] > T{0})) throw std::invalid_argument("GravityField: box sides must be positive");
            lmin = std::min(lmin, box_[i]);
            volume_ *= box_[i];
        }
        if (ewald_.eta / ewald_.alpha_L > 1.0)
            throw std::invalid_argument("GravityField: Ewald cutoff exceeds the box side (raise alpha_L or lower eta)");
        alpha_ = static_cast<T>(ewald_.alpha_L) / lmin;
        alpha2_ = alpha_ * alpha_;
        rc_ = static_cast<T>(ewald_.eta) / alpha_;
        rc2_ = rc_ * rc_;
        const T kc = T{2} * alpha_ * static_cast<T>(ewald_.eta);
        for (std::size_t i = 0; i < N; ++i)
            nmax_[i] = static_cast<int>(std::ceil(kc * box_[i] / (T{2} * std::numbers::pi_v<T>)));

        // Half space of integer vectors: first nonzero index positive.
        std::array<int, N> n{};
        for (std::size_t i = 0; i < N; ++i) n[i] = -nmax_[i];
        for (;;) {
            bool nonzero = false, positive = false;
            for (std::size_t i = 0; i < N; ++i) {
                if (n[i] != 0) { nonzero = true; positive = n[i] > 0; break; }
            }
            if (nonzero && positive) {
                KVec kv;
                T k2{};
                for (std::size_t i = 0; i < N; ++i) {
                    kv.n[i] = n[i];
                    kv.k[i] = T{2} * std::numbers::pi_v<T> * T(n[i]) / box_[i];
                    k2 += kv.k[i] * kv.k[i];
                }
                if (k2 <= kc * kc) {
                    kv.coef = reciprocal_coefficient(k2);
                    kvecs_.push_back(kv);
                }
            }
            std::size_t i = 0;
            for (; i < N; ++i) {
                if (++n[i] <= nmax_[i]) break;
                n[i] = -nmax_[i];
            }
            if (i == N) break;
        }
        structure_.assign(kvecs_.size(), std::complex<T>{});
    }

    // C(k) = pi^(N/2) / (2 Gamma(q)) * Int_{1/alpha^2}^inf u^(N/2-q) exp(-k^2 u/4 - eps^2/u) du.
    // With eps = 0 this is (4/k^2)^a Gamma(a, k^2/4alpha^2) times the prefactor, a = N/2 - q + 1.
    T reciprocal_coefficient(T k2) const {
        const double q = (law_.exponent + 1) / 2.0;
        const double u0 = 1.0 / alpha2_;
        const double sc = 4.0 / k2;
        const double e = 0.5 * N - q;   // exponent of u
        static constexpr double edges[] = {0, 1, 3, 7, 15, 31, 63};
        const auto& gl = gravity_detail::gauss_legendre_12();
        double integral = 0;
        for (int p = 0; p < 6; ++p) {
            const double a = edges[p], b = edges[p + 1];
            for (const auto& [x, w] : gl) {
                const double v = 0.5 * (b - a) * x + 0.5 * (b + a);
                const double u = u0 + sc * v;
                integral += 0.5 * (b - a) * w * std::exp(-v) * std::pow(u, e) * std::exp(-eps2_ / u);
            }
        }
        integral *= sc * std::exp(-k2 * u0 / 4.0);
        return static_cast<T>(std::pow(std::numbers::pi, 0.5 * N) / (2.0 * std::tgamma(q)) * integral);
    }

    T inv_power(T s) const {   // s^(-(p+1)/2)
        switch (law_.exponent) {
            case 1:  return T{1} / s;
            case 2:  return T{1} / (s * std::sqrt(s));
            case 3:  return T{1} / (s * s);
            default: return std::pow(s, -T(law_.exponent + 1) / T{2});
        }
    }

    T free_potential(T s) const {
        if (law_.exponent == 1) return T{0.5} * std::log(s);
        const T q = T(law_.exponent + 1) / T{2};
        return std::pow(s, T{1} - q) / (T{2} * (T{1} - q));
    }

    // Sum of potential_real over the images of one source that reach x.
    T pair_potential_real(Vec<T, N> d) const {
        if (!periodic_) return potential_real(d);
        T phi{};
        for (std::size_t i = 0; i < N; ++i) d[i] -= box_[i] * std::round(d[i] / box_[i]);
        std::array<int, N> n;
        n.fill(-1);
        for (;;) {
            Vec<T, N> dd;
            for (std::size_t i = 0; i < N; ++i) dd[i] = d[i] + T(n[i]) * box_[i];
            phi += potential_real(dd);
            std::size_t i = 0;
            for (; i < N; ++i) {
                if (++n[i] <= 1) break;
                n[i] = -1;
            }
            if (i == N) break;
        }
        return phi;
    }

    void phase_tables(const Vec<T, N>& x, Tables& tab) const {
        for (std::size_t i = 0; i < N; ++i) {
            const int m = nmax_[i];
            tab[i].resize(2 * m + 1);
            const T th = T{2} * std::numbers::pi_v<T> * x[i] / box_[i];
            const std::complex<T> base{std::cos(th), std::sin(th)};
            tab[i][m] = {T{1}, T{0}};
            for (int n = 1; n <= m; ++n) {
                tab[i][m + n] = tab[i][m + n - 1] * base;
                tab[i][m - n] = std::conj(tab[i][m + n]);
            }
        }
    }

    std::complex<T> phase(const Tables& tab, const KVec& kv) const {
        std::complex<T> p = tab[0][kv.n[0] + nmax_[0]];
        for (std::size_t i = 1; i < N; ++i) p *= tab[i][kv.n[i] + nmax_[i]];
        return p;
    }

    void build_structure_factors() {
        std::fill(structure_.begin(), structure_.end(), std::complex<T>{});
        Tables tab;
        for (const auto& s : sources_) {
            phase_tables(s.position, tab);
            for (std::size_t v = 0; v < kvecs_.size(); ++v)
                structure_[v] += s.mass * std::conj(phase(tab, kvecs_[v]));
        }
    }
};

// ── Fixed-field force for the single-body integrators ─────────
// F = coupling * body.mass * g(x): a callable of the force signature, so
// `verlet_step(body, GravityPull{&field, c}, dt)` and `composite_force(...)`
// take it as they take `PointGravity`. `Field` is a `GravityField` or a
// `GravityTree`: anything with `acceleration(x)`. `coupling` is the body's
// gravitational charge over its inertial mass.
template<typename Field>
struct GravityPull {
    using T = typename Field::ScalarType;
    static constexpr std::size_t N = Field::dimension;

    const Field* field = nullptr;
    T coupling{1};

    Vec<T, N> operator()(const PointMass<N, T>& body, T /*t*/) const {
        return Vec<T, N>{field->acceleration(body.state.position) * (coupling * body.mass)};
    }
};

} // namespace spatium::physics::mechanics
