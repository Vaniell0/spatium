#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <algorithm>
#  include <array>
#  include <cmath>
#  include <compare>
#  include <cstddef>
#  include <format>
#  include <limits>
#  include <type_traits>
#endif

SPATIUM_EXPORT namespace spatium {
inline namespace algebra {

// A limit as arithmetic, not as extrapolation.
//
// `Series<T,N>` is a truncated Laurent series in an infinitesimal eps > 0:
//
//   f = c0 eps^v + c1 eps^(v+1) + ... + c(N-1) eps^(v+N-1) + O(eps^p)
//
// with the valuation v (the exponent of the leading term, negative for a
// pole) and the precision p (the order up to which the terms are known) kept
// beside the coefficients. It is a `Scalar`, so any function written once over
// `Scalar` and calling math through ADL -- the rule of this codebase, the same
// one `Dual` relies on -- runs on it unchanged, and what comes out is the
// behaviour of the function at a point, exactly:
//
//   x = a + eps       f(x) is f's expansion at a; `limit` reads it
//   x = 1 / eps       the same at infinity: a point at infinity is a chart,
//                     and no inf ever enters an operation
//
// `Dual` is the case N = 2 with v = 0 and the derivative as the second term;
// a Series is what is left when the derivative is not the only thing wanted.
// What it answers where extrapolation guesses: sin(x)/x at 0 is 1, not a
// division of two zeros; 1/x - 1/sin(x) cancels a pole exactly and returns the
// 0 that is left; a pole comes back as a pole with its order and its sign.
//
// Honest about what it cannot hold. A Laurent series has integer exponents and
// no logarithms: exp(1/x) and log(x) at 0 are not Laurent, and the series says
// so (`Status::NonLaurent`) instead of returning a number. Fractional powers
// are met by the ramification: sqrt(x) is a series in t where x = t^2, so
// `limit` retries with x = a + t^r, r in {1, 2, 3, 4, 6}, and the order it
// reports is an exponent over r. What a cancellation consumes is recorded in
// the precision, not hidden: (sin x - x) / x^3 keeps fewer terms than sin x.
//
// Order is that of eps -> 0+: a < b when the leading term of b - a is
// positive, which is a total order and what an `if (x < y)` in f should mean
// next to a limit from the right. Equality is equality of the known terms.

namespace series_detail {

inline constexpr int kExact = 1 << 28;     // "known to every order": a polynomial, an exact zero

// p + d on orders, where kExact absorbs.
constexpr int padd(int p, int d) {
    if (p >= kExact || d >= kExact || d <= -kExact) return kExact;
    const long long r = static_cast<long long>(p) + d;
    return r >= kExact ? kExact : static_cast<int>(r);
}

enum class Status { Ok, NonLaurent, Indeterminate };

constexpr Status worse(Status a, Status b) { return a == Status::Ok ? b : a; }

}  // namespace series_detail

template<Scalar T = double, int N = 8>
class Series {
    static_assert(N >= 2, "a Series keeps at least two terms");
public:
    using coefficient_type = T;
    using Status = series_detail::Status;
    static constexpr int terms = N;
    static constexpr int kExact = series_detail::kExact;

    // The exact zero.
    Series() = default;
    Series(T v) { if (!(v == T{0})) { c_[0] = v; val_ = 0; } }
    template<class U> requires (std::is_arithmetic_v<U> && !std::is_same_v<U, T>)
    Series(U v) : Series(T(v)) {}

    // coef * eps^exponent, exact.
    static Series monomial(T coef, int exponent) {
        Series s;
        if (!(coef == T{0})) { s.c_[0] = coef; s.val_ = exponent; }
        return s;
    }
    // a + eps^r, the variable of a limit at a (r > 1 for fractional powers).
    static Series variable(T a, int ramification = 1) { return Series(a) + monomial(T{1}, ramification); }
    // a - eps^r: the same from the left.
    static Series variable_from_left(T a, int ramification = 1) { return Series(a) - monomial(T{1}, ramification); }
    // +-eps^-r: x -> +-infinity as eps -> 0.
    static Series at_infinity(int ramification = 1, bool negative = false) {
        return monomial(negative ? T{-1} : T{1}, -ramification);
    }

    Status status() const { return status_; }
    bool ok() const { return status_ == Status::Ok; }
    // No known nonzero term: either the zero, or O(eps^known_to()).
    bool is_zero() const { return val_ == prec_; }
    bool is_exact_zero() const { return val_ == kExact && prec_ == kExact; }
    // Exponent of the leading term (for a zero, the order below which it is zero).
    int valuation() const { return val_; }
    // Terms of exponent >= known_to() are not known.
    int known_to() const { return prec_; }
    // Coefficient of eps^e; 0 below the leading term and where unknown.
    T coefficient(int e) const {
        const long long i = static_cast<long long>(e) - val_;
        if (i < 0 || i >= N || e >= prec_) return T{0};
        return c_[static_cast<std::size_t>(i)];
    }
    T leading() const { return c_[0]; }

    // A single known term and nothing else, exactly.
    bool is_exact_monomial() const {
        if (!ok() || prec_ != kExact || is_zero()) return false;
        for (int i = 1; i < N; ++i) if (!(c_[static_cast<std::size_t>(i)] == T{0})) return false;
        return true;
    }
    bool is_exact_constant() const {
        return ok() && (is_exact_zero() || (is_exact_monomial() && val_ == 0));
    }

    // What the expression tends to as eps -> 0+: the constant term, 0 where the
    // leading exponent is positive, +-inf for a pole, NaN where it is not a series.
    T limit_value() const {
        using L = std::numeric_limits<T>;
        if (!ok()) return L::quiet_NaN();
        if (is_zero() || val_ > 0) return T{0};
        if (val_ == 0) return c_[0];
        return c_[0] > T{0} ? L::infinity() : -L::infinity();
    }

    // ── arithmetic ────────────────────────────────────────────

    Series operator-() const {
        Series s = *this;
        for (auto& x : s.c_) x = -x;
        return s;
    }
    friend Series operator+(const Series& a, const Series& b) { return add(a, b, false); }
    friend Series operator-(const Series& a, const Series& b) { return add(a, b, true); }
    friend Series operator*(const Series& a, const Series& b) { return mul(a, b); }
    friend Series operator/(const Series& a, const Series& b) { return div(a, b); }
    Series& operator+=(const Series& o) { return *this = *this + o; }
    Series& operator-=(const Series& o) { return *this = *this - o; }
    Series& operator*=(const Series& o) { return *this = *this * o; }
    Series& operator/=(const Series& o) { return *this = *this / o; }

    // Equal as far as both are known; ordered by the sign of the leading term
    // of the difference (eps -> 0+). A series that is not one is unordered.
    bool operator==(const Series& o) const {
        const Series d = *this - o;
        return d.ok() && d.is_zero();
    }
    std::partial_ordering operator<=>(const Series& o) const {
        const Series d = *this - o;
        if (!d.ok()) return std::partial_ordering::unordered;
        if (d.is_zero()) return std::partial_ordering::equivalent;
        return d.c_[0] < T{0} ? std::partial_ordering::less : std::partial_ordering::greater;
    }

    // ── building blocks for the functions below ───────────────

    static constexpr int kWork = 2 * N;       // terms computed before the window is cut

    // Normalise coefficients raw[i] of eps^(base+i), known below `prec`: strip
    // leading zeros, keep a window of N, and lower the precision where nonzero
    // terms had to be dropped.
    template<std::size_t L>
    static Series make(int base, const std::array<T, L>& raw, int prec, Status st = Status::Ok) {
        Series s;
        s.status_ = st;
        prec = std::min(prec, kExact);
        int first = -1;
        for (int i = 0; i < static_cast<int>(L) && base + i < prec; ++i)
            if (!(raw[static_cast<std::size_t>(i)] == T{0})) { first = i; break; }
        if (first < 0) { s.val_ = s.prec_ = prec; return s; }       // O(eps^prec), or the exact zero
        s.val_ = base + first;
        s.prec_ = prec;
        for (int i = 0; i < N; ++i) {
            const int idx = first + i;
            if (idx < static_cast<int>(L) && base + idx < prec) s.c_[static_cast<std::size_t>(i)] = raw[static_cast<std::size_t>(idx)];
        }
        for (int idx = first + N; idx < static_cast<int>(L) && base + idx < prec; ++idx)
            if (!(raw[static_cast<std::size_t>(idx)] == T{0})) { s.prec_ = std::min(s.prec_, base + first + N); break; }
        return s;
    }
    // A result that is not a series.
    static Series invalid(Status st) {
        Series s;
        s.status_ = st;
        s.val_ = s.prec_ = 0;
        s.c_[0] = std::numeric_limits<T>::quiet_NaN();
        return s;
    }

private:
    std::array<T, N> c_{};
    int val_ = kExact;
    int prec_ = kExact;
    Status status_ = Status::Ok;

    static Series add(const Series& a, const Series& b, bool negate_b) {
        const Status st = series_detail::worse(a.status_, b.status_);
        if (st != Status::Ok) return invalid(st);
        std::array<T, kWork> raw{};
        int prec = std::min(a.prec_, b.prec_);
        const int v = std::min(a.val_, b.val_);
        bool lost = false;
        const auto put = [&](const Series& s, bool neg) {
            for (int i = 0; i < N; ++i) {
                const int e = s.val_ + i;
                if (e >= prec) break;
                const T& x = s.c_[static_cast<std::size_t>(i)];
                if (x == T{0}) continue;
                const int k = e - v;
                if (k >= kWork) { lost = true; continue; }
                raw[static_cast<std::size_t>(k)] = neg ? T(raw[static_cast<std::size_t>(k)] - x) : T(raw[static_cast<std::size_t>(k)] + x);
            }
        };
        put(a, false);
        put(b, negate_b);
        if (lost) prec = std::min(prec, v + kWork);
        return make(v, raw, prec);
    }

    static Series mul(const Series& a, const Series& b) {
        const Status st = series_detail::worse(a.status_, b.status_);
        if (st != Status::Ok) return invalid(st);
        std::array<T, kWork> raw{};
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                raw[static_cast<std::size_t>(i + j)] = raw[static_cast<std::size_t>(i + j)] + a.c_[static_cast<std::size_t>(i)] * b.c_[static_cast<std::size_t>(j)];
        const int prec = std::min(series_detail::padd(b.prec_, a.val_), series_detail::padd(a.prec_, b.val_));
        return make(std::min(a.val_ + b.val_, kExact), raw, prec);
    }

    static Series div(const Series& a, const Series& b) {
        const Status st = series_detail::worse(a.status_, b.status_);
        if (st != Status::Ok) return invalid(st);
        if (b.is_zero()) return invalid(Status::Indeterminate);       // a series cannot be divided by a zero
        const int base = a.val_ - b.val_;
        std::array<T, kWork> q{};
        const T b0 = b.c_[0];
        for (int k = 0; k < kWork; ++k) {
            T acc = k < N ? a.c_[static_cast<std::size_t>(k)] : T{0};
            for (int j = 1; j <= k && j < N; ++j) acc = acc - b.c_[static_cast<std::size_t>(j)] * q[static_cast<std::size_t>(k - j)];
            q[static_cast<std::size_t>(k)] = acc / b0;
        }
        int prec = std::min(series_detail::padd(a.prec_, -b.val_), series_detail::padd(b.prec_, a.val_ - 2 * b.val_));
        if (!b.is_exact_monomial()) prec = std::min(prec, base + kWork);   // 1/(1-eps) never ends
        return make(base, q, a.is_zero() ? series_detail::padd(a.prec_, -b.val_) : prec);
    }
};

// ── functions, found by ADL like Dual's ───────────────────────

namespace series_detail {

// Coefficient of eps^k of f, k >= 0, as the work array the recurrences read.
template<Scalar T, int N>
std::array<T, Series<T, N>::kWork> work_coefficients(const Series<T, N>& f) {
    std::array<T, Series<T, N>::kWork> g{};
    for (int k = 0; k < Series<T, N>::kWork; ++k) g[static_cast<std::size_t>(k)] = f.coefficient(k);
    return g;
}

// exp, sin, cos of a series that has no pole: not a series where it does.
template<Scalar T, int N>
bool analytic_at_zero(const Series<T, N>& f) {
    return f.ok() && (f.is_exact_zero() || std::min(f.valuation(), f.known_to()) >= 0);
}

}  // namespace series_detail

template<Scalar T, int N>
Series<T, N> exp(const Series<T, N>& f) {
    using std::exp;
    if (!f.ok()) return f;
    if (!series_detail::analytic_at_zero(f)) return Series<T, N>::invalid(series_detail::Status::NonLaurent);   // exp(pole): essential
    if (f.is_exact_constant()) return Series<T, N>(exp(f.coefficient(0)));
    constexpr int M = Series<T, N>::kWork;
    const auto g = series_detail::work_coefficients(f);
    std::array<T, M> e{};
    e[0] = exp(g[0]);
    for (int n = 1; n < M; ++n) {
        T acc{0};
        for (int k = 1; k <= n; ++k) acc = acc + T(k) * g[static_cast<std::size_t>(k)] * e[static_cast<std::size_t>(n - k)];
        e[static_cast<std::size_t>(n)] = acc / T(n);
    }
    return Series<T, N>::make(0, e, std::min(f.known_to(), M));
}

template<Scalar T, int N>
Series<T, N> sin(const Series<T, N>& f) {
    using std::sin, std::cos;
    if (!f.ok()) return f;
    if (!series_detail::analytic_at_zero(f)) return Series<T, N>::invalid(series_detail::Status::NonLaurent);   // sin(1/x)
    if (f.is_exact_constant()) return Series<T, N>(sin(f.coefficient(0)));
    constexpr int M = Series<T, N>::kWork;
    const auto g = series_detail::work_coefficients(f);
    std::array<T, M> S{}, C{}, out{};
    S[0] = T{0}; C[0] = T{1};
    for (int n = 1; n < M; ++n) {
        T s{0}, c{0};
        for (int k = 1; k <= n; ++k) {
            s = s + T(k) * g[static_cast<std::size_t>(k)] * C[static_cast<std::size_t>(n - k)];
            c = c + T(k) * g[static_cast<std::size_t>(k)] * S[static_cast<std::size_t>(n - k)];
        }
        S[static_cast<std::size_t>(n)] = s / T(n);
        C[static_cast<std::size_t>(n)] = -c / T(n);
    }
    const T s0 = sin(g[0]), c0 = cos(g[0]);
    for (int n = 0; n < M; ++n)
        out[static_cast<std::size_t>(n)] = s0 * C[static_cast<std::size_t>(n)] + c0 * S[static_cast<std::size_t>(n)];
    return Series<T, N>::make(0, out, std::min(f.known_to(), M));
}

template<Scalar T, int N>
Series<T, N> cos(const Series<T, N>& f) {
    using std::sin, std::cos;
    if (!f.ok()) return f;
    if (!series_detail::analytic_at_zero(f)) return Series<T, N>::invalid(series_detail::Status::NonLaurent);
    if (f.is_exact_constant()) return Series<T, N>(cos(f.coefficient(0)));
    constexpr int M = Series<T, N>::kWork;
    const auto g = series_detail::work_coefficients(f);
    std::array<T, M> S{}, C{}, out{};
    S[0] = T{0}; C[0] = T{1};
    for (int n = 1; n < M; ++n) {
        T s{0}, c{0};
        for (int k = 1; k <= n; ++k) {
            s = s + T(k) * g[static_cast<std::size_t>(k)] * C[static_cast<std::size_t>(n - k)];
            c = c + T(k) * g[static_cast<std::size_t>(k)] * S[static_cast<std::size_t>(n - k)];
        }
        S[static_cast<std::size_t>(n)] = s / T(n);
        C[static_cast<std::size_t>(n)] = -c / T(n);
    }
    const T s0 = sin(g[0]), c0 = cos(g[0]);
    for (int n = 0; n < M; ++n)
        out[static_cast<std::size_t>(n)] = c0 * C[static_cast<std::size_t>(n)] - s0 * S[static_cast<std::size_t>(n)];
    return Series<T, N>::make(0, out, std::min(f.known_to(), M));
}

template<Scalar T, int N>
Series<T, N> tan(const Series<T, N>& f) { return sin(f) / cos(f); }

template<Scalar T, int N>
Series<T, N> sinh(const Series<T, N>& f) { return (exp(f) - exp(-f)) / Series<T, N>(T{2}); }

template<Scalar T, int N>
Series<T, N> cosh(const Series<T, N>& f) { return (exp(f) + exp(-f)) / Series<T, N>(T{2}); }

// log of a series that starts at a positive constant: log(eps) is not a series.
template<Scalar T, int N>
Series<T, N> log(const Series<T, N>& f) {
    using std::log;
    if (!f.ok()) return f;
    if (f.is_zero() || f.valuation() != 0) return Series<T, N>::invalid(series_detail::Status::NonLaurent);
    if (f.is_exact_constant()) return Series<T, N>(log(f.coefficient(0)));
    constexpr int M = Series<T, N>::kWork;
    const auto g = series_detail::work_coefficients(f);
    std::array<T, M> L{};
    L[0] = log(g[0]);
    for (int n = 1; n < M; ++n) {
        T acc{0};
        for (int k = 1; k < n; ++k) acc = acc + T(k) * L[static_cast<std::size_t>(k)] * g[static_cast<std::size_t>(n - k)];
        L[static_cast<std::size_t>(n)] = (g[static_cast<std::size_t>(n)] - acc / T(n)) / g[0];
    }
    return Series<T, N>::make(0, L, std::min(f.known_to(), M));
}

// f^n for a real n: (c0 eps^v (1 + h))^n = c0^n eps^(v n) (1 + h)^n, which is a
// Laurent series only where v n is an integer; otherwise the caller raises the
// ramification (`limit` does).
template<Scalar T, int N>
Series<T, N> pow(const Series<T, N>& f, T n) {
    using std::pow;
    if (!f.ok()) return f;
    if (n == T{0}) return Series<T, N>(T{1});
    if (f.is_zero()) {
        if (n > T{0} && f.is_exact_zero()) return f;
        return Series<T, N>::invalid(n > T{0} ? series_detail::Status::NonLaurent : series_detail::Status::Indeterminate);
    }
    const double vn = static_cast<double>(f.valuation()) * primal_double(n);
    const double m_round = std::round(vn);
    if (std::abs(vn - m_round) > 1e-9) return Series<T, N>::invalid(series_detail::Status::NonLaurent);
    const int m = static_cast<int>(m_round);
    constexpr int M = Series<T, N>::kWork;
    std::array<T, M> P{}, c{};
    for (int k = 0; k < M; ++k) c[static_cast<std::size_t>(k)] = k < N ? f.coefficient(f.valuation() + k) : T{0};
    P[0] = pow(c[0], n);
    for (int k = 1; k < M; ++k) {
        T acc{0};
        for (int j = 1; j <= k; ++j)
            acc = acc + (n * T(j) - T(k - j)) * c[static_cast<std::size_t>(j)] * P[static_cast<std::size_t>(k - j)];
        P[static_cast<std::size_t>(k)] = acc / (T(k) * c[0]);
    }
    int prec = series_detail::padd(f.known_to(), m - f.valuation());
    if (!f.is_exact_monomial()) prec = std::min(prec, m + M);
    return Series<T, N>::make(m, P, prec);
}

template<Scalar T, int N>
Series<T, N> pow(const Series<T, N>& f, int n) {
    return pow(f, T(n));
}

template<Scalar T, int N>
Series<T, N> sqrt(const Series<T, N>& f) { return pow(f, T{1} / T{2}); }

template<Scalar T, int N>
Series<T, N> abs(const Series<T, N>& f) {
    if (!f.ok() || f.is_zero()) return f;
    return f.leading() < T{0} ? -f : f;
}

using Series64 = Series<double, 8>;

// ── limits ────────────────────────────────────────────────────

enum class LimitKind {
    Finite,         // a number, possibly 0 (f decays)
    Infinite,       // a pole: +-inf, with its order
    Undetermined    // not a Laurent series in any ramification tried: log, exp(pole), sin(1/x)
};

enum class Approach { Right, Left };

// What f(x) does as x -> a (or -> infinity): f ~ coefficient * t^exponent in
// t = (x - a)^(1/ramification), i.e. like (x - a)^(exponent/ramification) at a
// and like x^(-exponent/ramification) at infinity.
template<Scalar T>
struct LimitResult {
    LimitKind kind = LimitKind::Undetermined;
    T value{};                  // the limit; +-inf for Infinite, NaN for Undetermined
    T coefficient{};            // the leading coefficient
    int exponent = 0;           // the leading exponent in t
    int ramification = 1;
    int known_to = 0;           // the series is known below this exponent in t

    bool determined() const { return kind != LimitKind::Undetermined; }
    // The order of vanishing or growth in x - a (or 1/x): negative for a pole.
    double order() const { return static_cast<double>(exponent) / ramification; }
};

namespace series_detail {

template<Scalar T, int N, typename F, typename MakeX>
LimitResult<T> limit_with(F& f, MakeX make_x) {
    const T nan = std::numeric_limits<T>::quiet_NaN();
    for (const int r : {1, 2, 3, 4, 6}) {
        const Series<T, N> y = f(make_x(r));
        if (y.status() == Status::NonLaurent) continue;       // a fractional power: a larger ramification
        LimitResult<T> out;
        out.ramification = r;
        out.known_to = y.known_to();
        if (!y.ok()) return {LimitKind::Undetermined, nan, nan, 0, r, 0};
        out.exponent = y.valuation();
        out.coefficient = y.is_zero() ? T{0} : y.leading();
        out.value = y.limit_value();
        out.kind = (!y.is_zero() && y.valuation() < 0) ? LimitKind::Infinite : LimitKind::Finite;
        return out;
    }
    return {LimitKind::Undetermined, nan, nan, 0, 1, 0};
}

}  // namespace series_detail

// The limit of f(x) as x -> a from one side. `f` is written once over a scalar:
// `[](auto x) { using std::sin; return sin(x) / x; }`.
template<int N = 8, Scalar T, typename F>
LimitResult<T> limit(F&& f, T a, Approach side = Approach::Right) {
    return series_detail::limit_with<T, N>(f, [&](int r) {
        return side == Approach::Right ? Series<T, N>::variable(a, r) : Series<T, N>::variable_from_left(a, r);
    });
}

// Both sides, and only where they agree: Finite with one value, or Infinite
// with one sign.
template<int N = 8, Scalar T, typename F>
LimitResult<T> two_sided_limit(F&& f, T a) {
    const auto r = limit<N>(f, a, Approach::Right);
    const auto l = limit<N>(f, a, Approach::Left);
    const bool agree = r.determined() && l.determined() && r.kind == l.kind && r.value == l.value;
    if (agree) return r;
    return {LimitKind::Undetermined, std::numeric_limits<T>::quiet_NaN(), T{}, 0, 1, 0};
}

// The limit as x -> +infinity (or -infinity), through the chart x = 1/eps.
template<int N = 8, Scalar T = double, typename F>
LimitResult<T> limit_at_infinity(F&& f, bool toward_negative = false) {
    return series_detail::limit_with<T, N>(f, [&](int r) {
        return Series<T, N>::at_infinity(r, toward_negative);
    });
}

} // namespace algebra
} // namespace spatium

// std::format support
template<spatium::Scalar T, int N>
struct std::formatter<spatium::Series<T, N>> {
    constexpr auto parse(auto& ctx) { return ctx.begin(); }
    auto format(const spatium::Series<T, N>& s, auto& ctx) const {
        auto out = ctx.out();
        if (!s.ok()) return std::format_to(out, "(not a Laurent series)");
        if (s.is_exact_zero()) return std::format_to(out, "0");
        bool first = true;
        for (int e = s.valuation(); e < s.valuation() + N && e < s.known_to(); ++e) {
            const T c = s.coefficient(e);
            if (c == T{0}) continue;
            out = std::format_to(out, "{}{}eps^{}", first ? "" : " + ", c, e);
            first = false;
        }
        if (s.known_to() < spatium::Series<T, N>::kExact) out = std::format_to(out, "{}O(eps^{})", first ? "" : " + ", s.known_to());
        return out;
    }
};

// numeric_limits of a Series are its coefficient type's, as a Dual's.
template<class T, int N>
struct std::numeric_limits<spatium::Series<T, N>> : std::numeric_limits<T> {
    using S = spatium::Series<T, N>;
    static S min() noexcept { return S(std::numeric_limits<T>::min()); }
    static S max() noexcept { return S(std::numeric_limits<T>::max()); }
    static S lowest() noexcept { return S(std::numeric_limits<T>::lowest()); }
    static S epsilon() noexcept { return S(std::numeric_limits<T>::epsilon()); }
    static S infinity() noexcept { return S(std::numeric_limits<T>::infinity()); }
    static S quiet_NaN() noexcept { return S(std::numeric_limits<T>::quiet_NaN()); }
};
