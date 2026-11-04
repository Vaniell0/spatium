#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <algorithm>
#  include <charconv>
#  include <cmath>
#  include <cstddef>
#  include <cstdint>
#  include <ios>
#  include <limits>
#  include <string>
#  include <system_error>
#  include <type_traits>
#endif

// `Interval<T>` -- a number known to lie between two numbers of T, whatever rounding happened on
// the way.
//
// Every operation rounds its lower end down and its upper end up, so the true value of what was
// written -- over the reals, not over T -- stays inside: a sum, a product, a square root is the
// interval of every result its arguments' values can give, widened by the rounding of T. That is
// what turns a computed number into a proved one: `1/3` in an `Interval<double>` is
// [0.33333333333333331, 0.33333333333333337], and the fifteen digits its two ends share are
// fifteen digits of the third, not of a double near it (`proved` below). The base is any
// Scalar: a double, or `Real<N>` for an enclosure N digits wide.
//
// ── How the rounding is held ─────────────────────────────────────
//
// For a built-in floating type each end moves one step outward (`std::nextafter`) after each
// operation, unless the operation was exact, which an fma or a two-sum tells. IEEE + - * / and sqrt
// are correctly rounded, half a step from the truth, so one step out holds it. exp, log, sin, cos
// and pow come from libm, which is not correctly rounded: their ends move two steps, and the claim
// rests on the libm's bound for them: glibc's on the host (under one step, its manual's table of
// errors), musl's in WebAssembly (FreeBSD's and Arm's routines, under one step as well). For a type with
// no `nextafter`, `Real<N>`, an arithmetic end moves by two of its epsilons of its own size, and a
// function's by epsilons of its argument's size (`fdown`): Boost's bound, which a test holds
// against twice the digits, since its documentation states none.
// A constant enters by its text (`Interval("0.1")`): a double that is not the decimal written is
// widened to the doubles either side of it, which hold it; a decimal base holds it exactly.
//
// ── Comparisons are three-valued ─────────────────────────────────
//
// [1, 2] < [3, 4] is true, [3, 4] < [1, 2] false, and [1, 3] < [2, 4] is neither: some values of
// the two are ordered one way and some the other. `less` and the rest say which (`Truth`). The
// operators < <= > >= == != are what generic code calls -- every Scalar algorithm in the library,
// unchanged -- and on an undecided pair they take the order of the midpoints and raise a flag, the
// way IEEE raises "inexact": `Undecided` counts them. A computation that raised none followed the
// one branch every value in its inputs follows, so its interval holds; one that raised any took a
// guessed branch somewhere, and its interval proves nothing. A square root or a logarithm of an
// interval reaching below its domain raises it too, since part of what was asked has no value.
// The library's Scalar asks for a total order; this one is total on what it decides and says
// where it did not, rather than weakening every algorithm's contract to a partial one.

SPATIUM_EXPORT namespace spatium {

// How many comparisons of intervals were undecided on this thread, ever. Read through `Undecided`.
inline thread_local std::uint64_t interval_undecided = 0;

// Whether any comparison of intervals since this was made was undecided.
struct Undecided {
    std::uint64_t start = interval_undecided;
    bool raised() const { return interval_undecided != start; }
};

enum class Truth : std::uint8_t { No, Yes, Undecided };

namespace interval_detail {

template<class T>
inline constexpr bool builtin = std::is_floating_point_v<T>;

// One end moved outward by `steps` units of T's rounding.
template<class T>
T down(const T& x, int steps = 1) {
    if constexpr (builtin<T>) {
        T y = x;
        for (int i = 0; i < steps; ++i) y = std::nextafter(y, -std::numeric_limits<T>::infinity());
        return y;
    } else {
        using std::abs;
        return x - (abs(x) * std::numeric_limits<T>::epsilon() * T(2 * steps) + std::numeric_limits<T>::min());
    }
}
template<class T>
T up(const T& x, int steps = 1) {
    if constexpr (builtin<T>) {
        T y = x;
        for (int i = 0; i < steps; ++i) y = std::nextafter(y, std::numeric_limits<T>::infinity());
        return y;
    } else {
        using std::abs;
        return x + (abs(x) * std::numeric_limits<T>::epsilon() * T(2 * steps) + std::numeric_limits<T>::min());
    }
}

// A function's end moved by its error. For a built-in type, `steps` units of the result: glibc's
// bound, which holds near a sine's zero too, its argument reduced exactly. For another base, `steps`
// epsilons of the larger of the result and `scale`: Boost's sine near its zero is off by epsilons
// of the argument, not of the result, so sin(3.141592653589793) at fifty digits has thirty-three
// right, and an exponential loses as many as its argument has digits before the point.
template<class T> T fdown(const T& r, const T& scale, int steps) {
    if constexpr (builtin<T>) return down(r, steps);
    else {
        using std::abs;
        const T m = abs(r) < scale ? scale : abs(r);
        return r - (m * std::numeric_limits<T>::epsilon() * T(2 * steps) + std::numeric_limits<T>::min());
    }
}
template<class T> T fup(const T& r, const T& scale, int steps) {
    if constexpr (builtin<T>) return up(r, steps);
    else {
        using std::abs;
        const T m = abs(r) < scale ? scale : abs(r);
        return r + (m * std::numeric_limits<T>::epsilon() * T(2 * steps) + std::numeric_limits<T>::min());
    }
}

template<class T> bool is_nan(const T& x) { return !(x == x); }
template<class T> bool is_inf(const T& x) {
    return x == std::numeric_limits<T>::infinity() || x == -std::numeric_limits<T>::infinity();
}

// The smaller and larger of two ends, a NaN skipped: inf * 0 and inf / inf at an end of an
// interval are a value no point of it reaches, and the other products bound it.
template<class T> T lower(const T& a, const T& b) { return is_nan(a) ? b : is_nan(b) ? a : (b < a ? b : a); }
template<class T> T upper(const T& a, const T& b) { return is_nan(a) ? b : is_nan(b) ? a : (a < b ? b : a); }

// An operation's end rounded outward, and left where it is when the operation was exact: a sum or a
// product of doubles whose rounding error, found exactly (Knuth's two-sum, an fma), is zero. So
// 2 + 2 is [4, 4], not the 4 and its neighbours. Far below the normal range an error can vanish by
// underflow, so there the end moves regardless. A base with no fma always moves.
template<class T> bool tiny(const T& x) {
    using std::abs;
    return abs(x) < std::numeric_limits<T>::min() * T(4503599627370496.0);   // 2^52 above the least normal
}
template<class T> T add(const T& a, const T& b, bool to_up) {
    const T s = a + b;
    if constexpr (builtin<T>) {
        const T bb = s - a;
        if ((a - (s - bb)) + (b - bb) == T(0) && !tiny(s)) return s;
    }
    return to_up ? up(s) : down(s);
}
// An end times an end, where zero times an unbounded end is zero: the interval reaches the
// infinity, no value in it is one.
template<class T> T mul(const T& a, const T& b, bool to_up) {
    if (a == T(0) || b == T(0)) return T(0);
    const T p = a * b;
    if constexpr (builtin<T>)
        if (std::fma(a, b, -p) == T(0) && !tiny(p)) return p;
    return to_up ? up(p) : down(p);
}
template<class T> T div(const T& a, const T& b, bool to_up) {
    const T q = a / b;
    if constexpr (builtin<T>)
        if (std::fma(q, b, -a) == T(0) && !tiny(q) && !tiny(a)) return q;
    return to_up ? up(q) : down(q);
}

template<class T> T pi() {
    if constexpr (builtin<T>) return static_cast<T>(3.14159265358979323846264338327950288L);
    else { using std::atan; return T(4) * atan(T(1)); }
}

} // namespace interval_detail

template<Scalar T = double>
struct Interval {
    T lo{};
    T hi{};

    constexpr Interval() = default;
    constexpr Interval(T v) : lo(v), hi(v) {}   // NOLINT: a number is the interval of itself
    template<class U> requires (std::is_integral_v<U>)
    constexpr Interval(U v) : lo(T(v)), hi(T(v)) {}   // NOLINT: exact while |v| < 2^53 in a double
    constexpr Interval(T l, T h) : lo(l), hi(h) {}

    // A decimal as written, held: exact where T holds it, else between the two numbers of T beside it.
    explicit Interval(const char* text) {
        if constexpr (interval_detail::builtin<T>) {
            T x{};
            const char* end = text;
            while (*end) ++end;
            std::from_chars(text, end, x);
            lo = hi = x;
            if (!exactly(x, text, end)) lo = interval_detail::down(x), hi = interval_detail::up(x);
        } else {
            lo = hi = T(text);
        }
    }

    static Interval empty() { return {std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN()}; }
    static Interval entire() { return {-std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity()}; }

    T mid() const { return lo / T(2) + hi / T(2); }
    T width() const { return hi - lo; }
    bool is_empty() const { return interval_detail::is_nan(lo) || interval_detail::is_nan(hi); }
    bool contains(const T& x) const { return lo <= x && x <= hi; }

    Interval operator-() const { return {-hi, -lo}; }

    friend Interval operator+(const Interval& a, const Interval& b) {
        using namespace interval_detail;
        return {add(a.lo, b.lo, false), add(a.hi, b.hi, true)};
    }
    friend Interval operator-(const Interval& a, const Interval& b) {
        using namespace interval_detail;
        return {add(a.lo, T(-b.hi), false), add(a.hi, T(-b.lo), true)};
    }
    friend Interval operator*(const Interval& a, const Interval& b) {
        using namespace interval_detail;
        const T l[4] = {mul(a.lo, b.lo, false), mul(a.lo, b.hi, false), mul(a.hi, b.lo, false), mul(a.hi, b.hi, false)};
        const T h[4] = {mul(a.lo, b.lo, true), mul(a.lo, b.hi, true), mul(a.hi, b.lo, true), mul(a.hi, b.hi, true)};
        return {lower(lower(l[0], l[1]), lower(l[2], l[3])), upper(upper(h[0], h[1]), upper(h[2], h[3]))};
    }
    // A divisor that holds zero leaves nothing to bound the quotient by: the whole line, which is
    // still true, rather than a guess. A divisor that is zero alone has no quotient at all.
    friend Interval operator/(const Interval& a, const Interval& b) {
        using namespace interval_detail;
        if (b.lo == T(0) && b.hi == T(0)) return empty();
        if (!(T(0) < b.lo || b.hi < T(0))) return entire();
        const T l[4] = {div(a.lo, b.lo, false), div(a.lo, b.hi, false), div(a.hi, b.lo, false), div(a.hi, b.hi, false)};
        const T h[4] = {div(a.lo, b.lo, true), div(a.lo, b.hi, true), div(a.hi, b.lo, true), div(a.hi, b.hi, true)};
        return {lower(lower(l[0], l[1]), lower(l[2], l[3])), upper(upper(h[0], h[1]), upper(h[2], h[3]))};
    }

    Interval& operator+=(const Interval& o) { return *this = *this + o; }
    Interval& operator-=(const Interval& o) { return *this = *this - o; }
    Interval& operator*=(const Interval& o) { return *this = *this * o; }
    Interval& operator/=(const Interval& o) { return *this = *this / o; }

    friend Truth less(const Interval& a, const Interval& b) {
        return a.hi < b.lo ? Truth::Yes : !(a.lo < b.hi) ? Truth::No : Truth::Undecided;
    }
    friend Truth less_equal(const Interval& a, const Interval& b) {
        return !(b.lo < a.hi) ? Truth::Yes : b.hi < a.lo ? Truth::No : Truth::Undecided;
    }
    friend Truth equal(const Interval& a, const Interval& b) {
        if (a.lo == a.hi && b.lo == b.hi && a.lo == b.lo) return Truth::Yes;
        return a.hi < b.lo || b.hi < a.lo ? Truth::No : Truth::Undecided;
    }

    // What generic code branches on: the decided answer, or the midpoints' order with the flag raised.
    friend bool operator<(const Interval& a, const Interval& b) { return decide(less(a, b), a.mid() < b.mid()); }
    friend bool operator>(const Interval& a, const Interval& b) { return decide(less(b, a), b.mid() < a.mid()); }
    friend bool operator<=(const Interval& a, const Interval& b) { return decide(less_equal(a, b), !(b.mid() < a.mid())); }
    friend bool operator>=(const Interval& a, const Interval& b) { return decide(less_equal(b, a), !(a.mid() < b.mid())); }
    // Undecided, two intervals are taken equal when their ends are, so a copy equals what it copies.
    friend bool operator==(const Interval& a, const Interval& b) { return decide(equal(a, b), a.lo == b.lo && a.hi == b.hi); }
    friend bool operator!=(const Interval& a, const Interval& b) { return !(a == b); }

private:
    static bool decide(Truth t, bool guess) {
        if (t == Truth::Undecided) {
            ++interval_undecided;
            return guess;
        }
        return t == Truth::Yes;
    }

    // Whether the double `x` read from the decimal [text, end) is that decimal: when its exact
    // expansion, which is finite, has the same digits as the text.
    static bool exactly(T x, const char* text, const char* end) {
        if (!std::isfinite(x)) return false;
        char exact[1200];
        const auto r = std::to_chars(exact, exact + sizeof exact, x, std::chars_format::scientific, 1100);
        if (r.ec != std::errc{}) return false;
        const auto digits = [](const char* b, const char* e, int& exponent) {
            std::string d;
            int point = -1, lead = 0;
            bool any = false;
            const char* p = b;
            if (p != e && (*p == '-' || *p == '+')) ++p;
            for (; p != e && *p != 'e' && *p != 'E'; ++p) {
                if (*p == '.') { point = static_cast<int>(d.size()) + lead; continue; }
                if (*p == '0' && !any) { ++lead; continue; }
                any = true, d += *p;
            }
            int e10 = 0;
            if (p != e) std::from_chars(p + 1 + (p[1] == '+'), e, e10);
            const int whole = point < 0 ? static_cast<int>(d.size()) + lead : point;
            exponent = e10 + whole - lead;
            while (!d.empty() && d.back() == '0') d.pop_back();
            return d;
        };
        int e1 = 0, e2 = 0;
        const std::string a = digits(text, end, e1), b = digits(exact, r.ptr, e2);
        return a == b && (a.empty() || e1 == e2);
    }
};

static_assert(Scalar<Interval<double>>);
static_assert(Scalar<Interval<float>>);

// ── Functions, as enclosures ─────────────────────────────────────

template<Scalar T> Interval<T> abs(const Interval<T>& x) {
    if (!(x.lo < T(0))) return x;
    if (!(T(0) < x.hi)) return -x;
    return {T(0), x.hi < -x.lo ? -x.lo : x.hi};
}
template<Scalar T> Interval<T> min(const Interval<T>& a, const Interval<T>& b) {
    return {b.lo < a.lo ? b.lo : a.lo, b.hi < a.hi ? b.hi : a.hi};
}
template<Scalar T> Interval<T> max(const Interval<T>& a, const Interval<T>& b) {
    return {a.lo < b.lo ? b.lo : a.lo, a.hi < b.hi ? b.hi : a.hi};
}

template<Scalar T> Interval<T> sqrt(const Interval<T>& x) {
    using std::sqrt;
    using namespace interval_detail;
    if (x.hi < T(0)) return Interval<T>::empty();
    T l = x.lo;
    if (l < T(0)) ++interval_undecided, l = T(0);
    const auto end = [](const T& v, bool to_up) {
        const T r = sqrt(v);
        if constexpr (builtin<T>)
            if (std::fma(r, r, -v) == T(0) && !tiny(v)) return r;
        return to_up ? up(r) : down(r);
    };
    const T s = end(l, false);
    return {s < T(0) ? T(0) : s, end(x.hi, true)};
}

template<Scalar T> Interval<T> exp(const Interval<T>& x) {
    using std::exp;
    using namespace interval_detail;
    using std::abs;
    const T a = T(exp(x.lo)), b = T(exp(x.hi));
    const T l = fdown(a, T(a * (abs(x.lo) + T(1))), 2);
    return {l < T(0) ? T(0) : l, fup(b, T(b * (abs(x.hi) + T(1))), 2)};
}

namespace interval_detail {
// log's ends, where an end at or below zero is the logarithm's limit there: what a power reads.
template<Scalar T> Interval<T> log_ends(const Interval<T>& x) {
    using std::log;
    return {x.lo > T(0) ? fdown(T(log(x.lo)), T(1), 2) : -std::numeric_limits<T>::infinity(), fup(T(log(x.hi)), T(1), 2)};
}
} // namespace interval_detail

template<Scalar T> Interval<T> log(const Interval<T>& x) {
    if (!(T(0) < x.hi)) return Interval<T>::empty();
    if (!(T(0) < x.lo)) ++interval_undecided;
    return interval_detail::log_ends(x);
}

namespace interval_detail {
// Whether some phase + 2 pi k may lie on [lo, hi], k an integer: decided by where the ends fall in
// periods, with a margin for the rounding of that division, so that "may" is never "no" wrongly.
template<Scalar T> bool may_hit(const T& lo, const T& hi, const T& phase) {
    using std::abs, std::floor, std::ceil;
    const T period = T(2) * pi<T>();
    const T eps = std::numeric_limits<T>::epsilon();
    const T a = (lo - phase) / period, b = (hi - phase) / period;
    const T slack = T(16) * eps * (abs(lo) + abs(hi) + T(1));
    return !(floor(b + slack) < ceil(a - slack));
}
// sin (shift 0) or cos (shift pi/2 back): a maximum at pi/2 - shift, a minimum at -pi/2 - shift.
template<Scalar T> Interval<T> wave(const Interval<T>& x, bool cosine) {
    using std::sin, std::cos, std::abs;
    const T big = T(1099511627776);   // 2^40: past it a period is too few units of T to place an end in
    if (!(x.lo == x.lo) || !(x.hi == x.hi) || is_inf(x.lo) || is_inf(x.hi) || !(x.hi - x.lo < T(6)) ||
        big < abs(x.lo) || big < abs(x.hi))
        return {T(-1), T(1)};
    const T half = pi<T>() / T(2);
    const T top = cosine ? T(0) : half, bottom = cosine ? pi<T>() : -half;
    const T a = cosine ? T(cos(x.lo)) : T(sin(x.lo)), b = cosine ? T(cos(x.hi)) : T(sin(x.hi));
    const T scale = T(1) + (abs(x.lo) < abs(x.hi) ? abs(x.hi) : abs(x.lo));
    T l = may_hit(x.lo, x.hi, bottom) ? T(-1) : fdown(a < b ? a : b, scale, 2);
    T h = may_hit(x.lo, x.hi, top) ? T(1) : fup(a < b ? b : a, scale, 2);
    return {l < T(-1) ? T(-1) : l, T(1) < h ? T(1) : h};
}
} // namespace interval_detail

template<Scalar T> Interval<T> sin(const Interval<T>& x) { return interval_detail::wave(x, false); }
template<Scalar T> Interval<T> cos(const Interval<T>& x) { return interval_detail::wave(x, true); }

// x^n for a whole n is a product with x's sign kept, an even one with |x|, so [-1, 2]^2 is [0, 4]
// and not the [-2, 4] multiplying would give; any other power is exp(y log x), defined for x >= 0.
template<Scalar T> Interval<T> pow(const Interval<T>& x, const Interval<T>& y) {
    using std::pow, std::floor, std::abs;
    using namespace interval_detail;
    if (y.lo == y.hi && floor(y.lo) == y.lo && abs(y.lo) < T(1 << 30)) {
        const T n = y.lo;
        if (n == T(0)) return {T(1), T(1)};
        if (n < T(0)) return Interval<T>{T(1)} / pow(x, Interval<T>{-n});
        const bool even = floor(n / T(2)) * T(2) == n;
        // one number to a whole power: by squaring, so a power that is exact, 10^20, stays one number
        if (x.lo == x.hi) {
            Interval<T> r{T(1)}, b = x;
            for (long long k = static_cast<long long>(n); k > 0; k >>= 1) {
                if (k & 1) r = r * b;
                if (k > 1) b = b * b;
            }
            return r;
        }
        const Interval<T> base = even ? abs(x) : x;
        const T pl = T(pow(base.lo, n)), ph = T(pow(base.hi, n));
        const T l = fdown(pl, T(abs(pl) * n), 2);
        return {even && l < T(0) ? T(0) : l, fup(ph, T(abs(ph) * n), 2)};
    }
    if (x.hi < T(0)) return Interval<T>::empty();
    if (x.lo < T(0)) ++interval_undecided;
    const Interval<T> base{x.lo < T(0) ? T(0) : x.lo, x.hi};
    if (base.hi == T(0)) return y.lo > T(0) ? Interval<T>{T(0)} : Interval<T>::empty();
    return exp(y * log_ends(base));
}

// ── What an interval proves, in decimal ──────────────────────────
//
// The most significant digits n for which both ends round to the same n-digit decimal D: then
// every value between them is within half a unit of D's last digit, which is what "correct to n
// digits" means. An interval about zero proves no digit; `bound` is then how far from zero it
// reaches. `text` is D, plain where its exponent is moderate.
struct Proved {
    int digits = 0;
    std::string text;
    double bound = 0;     // the largest |value| the interval allows
    bool exact = false;   // the interval is one number: the value is it, and `text` its shortest digits
};

namespace interval_detail {
template<class T> std::string scientific(const T& x, int digits) {
    if constexpr (builtin<T>) {
        char buf[64];
        const auto r = std::to_chars(buf, buf + sizeof buf, x, std::chars_format::scientific, digits - 1);
        return std::string(buf, r.ptr);
    } else {
        return x.str(digits - 1, std::ios_base::scientific);
    }
}
// "4.2982e+00" as "4.2982", "1.25e-03" as "0.00125"; a far exponent stays as it is.
inline std::string plain(const std::string& s) {
    const auto e = s.find_first_of("eE");
    if (e == std::string::npos) return s;
    int exponent = 0;
    std::from_chars(s.data() + e + 1 + (s[e + 1] == '+'), s.data() + s.size(), exponent);
    if (exponent < -6 || exponent > 20) return s;
    const bool negative = s[0] == '-';
    std::string d;
    for (std::size_t i = negative; i < e; ++i)
        if (s[i] != '.') d += s[i];
    std::string out;
    if (exponent < 0) out = "0." + std::string(static_cast<std::size_t>(-exponent - 1), '0') + d;
    else if (static_cast<std::size_t>(exponent) + 1 >= d.size())
        out = d + std::string(static_cast<std::size_t>(exponent) + 1 - d.size(), '0');
    else out = d.substr(0, static_cast<std::size_t>(exponent) + 1) + "." + d.substr(static_cast<std::size_t>(exponent) + 1);
    return negative ? "-" + out : out;
}
} // namespace interval_detail

template<Scalar T> Proved proved(const Interval<T>& x) {
    using namespace interval_detail;
    using std::abs;
    Proved out;
    if (x.is_empty()) return out;
    const T reach = abs(x.lo) < abs(x.hi) ? abs(x.hi) : abs(x.lo);
    out.bound = static_cast<double>(reach);
    if (is_inf(x.lo) || is_inf(x.hi) || (x.lo < T(0) && T(0) < x.hi)) return out;
    if (x.lo == x.hi) {
        if constexpr (builtin<T>) {
            char buf[64];
            const auto r = std::to_chars(buf, buf + sizeof buf, x.lo);
            out.text = plain(std::string(buf, r.ptr));
        } else {
            out.text = plain(x.lo.str(0, std::ios_base::scientific));
        }
        if (const auto e = out.text.find_first_of("eE"), point = out.text.find('.'); point != std::string::npos) {
            const std::size_t stop = e == std::string::npos ? out.text.size() : e;
            std::size_t last = stop;
            while (last > point + 1 && out.text[last - 1] == '0') --last;
            if (last == point + 1) --last;
            out.text.erase(last, stop - last);
        }
        bool lead = true;
        for (char c : out.text) {
            if (c == 'e' || c == 'E') break;
            if (c < '0' || c > '9' || (lead && c == '0')) continue;
            lead = false, ++out.digits;
        }
        out.exact = true;
        return out;
    }
    for (int n = std::numeric_limits<T>::max_digits10; n > 0; --n) {
        const std::string a = scientific(x.lo, n);
        if (a == scientific(x.hi, n)) {
            out.digits = n;
            out.text = plain(a);
            return out;
        }
    }
    return out;
}

} // namespace spatium
