#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/complex.hpp>
#  include <spatium/core/interval.hpp>
#  include <algorithm>
#  include <array>
#  include <limits>
#  include <cmath>
#  include <vector>
#endif

// The Riemann zeta function, proved. zeta(s) is a complex interval that holds the true value: the
// Euler-Maclaurin sum with its remainder bounded, every operation rounded outward. On the critical
// line, Hardy's Z(t) = e^{i theta(t)} zeta(1/2 + it) is real, so a zero of zeta there is a change
// of Z's sign, proved where Z is proved positive on one side and negative on the other. How many
// zeros lie below a height T is N(T) = theta(T)/pi + 1 + S(T); where Re zeta stays positive on the
// way from 2 + iT to 1/2 + iT, |S(T)| < 1/2 (Backlund), so N(T) is the one whole number near
// theta(T)/pi + 1. When the sign changes found between two such heights are as many as N says,
// they are all the zeros there are between them: each simple, each on the critical line.

SPATIUM_EXPORT namespace spatium {
inline namespace algebra {

namespace zeta_detail {

// B_2, B_4, ..., B_60 as fractions.
inline constexpr const char* bernoulli[30][2] = {
    {"1", "6"}, {"-1", "30"}, {"1", "42"}, {"-1", "30"}, {"5", "66"}, {"-691", "2730"}, {"7", "6"},
    {"-3617", "510"}, {"43867", "798"}, {"-174611", "330"}, {"854513", "138"}, {"-236364091", "2730"},
    {"8553103", "6"}, {"-23749461029", "870"}, {"8615841276005", "14322"}, {"-7709321041217", "510"},
    {"2577687858367", "6"}, {"-26315271553053477373", "1919190"}, {"2929993913841559", "6"},
    {"-261082718496449122051", "13530"}, {"1520097643918070802691", "1806"},
    {"-27833269579301024235023", "690"}, {"596451111593912163277961", "282"},
    {"-5609403368997817686249127547", "46410"}, {"495057205241079648212477525", "66"},
    {"-801165718135489957347924991853", "1590"}, {"29149963634884862421418123812691", "798"},
    {"-2479392929313226753685415739663229", "870"}, {"84483613348880041862046775994036021", "354"},
    {"-1215233140483755572040304994079820246041491", "56786730"}};

// B_2k, k = 1 .. 30, read from the fractions once
template<Scalar T> Interval<T> bernoulli_number(int k) {
    static const auto table = [] {
        std::array<Interval<T>, 30> b{};
        for (int i = 0; i < 30; ++i) b[i] = Interval<T>(bernoulli[i][0]) / Interval<T>(bernoulli[i][1]);
        return b;
    }();
    return table[k - 1];
}

template<Scalar T> Interval<T> pi() {
    const T p = interval_detail::pi<T>();
    return {interval_detail::down(p), interval_detail::up(p)};
}

template<Scalar T> Interval<T> square(const Interval<T>& x) { return pow(x, Interval<T>(2)); }
template<Scalar T> Interval<T> modulus(const Complex<Interval<T>>& z) { return sqrt(square(z.re) + square(z.im)); }

// arg z where Re z > 0
template<Scalar T> Interval<T> arg_right(const Complex<Interval<T>>& z) { return atan(z.im / z.re); }

template<class R> inline constexpr bool is_interval = false;
template<class T> inline constexpr bool is_interval<Interval<T>> = true;

inline double modulus(const Complex<double>& z) { return std::hypot(z.re, z.im); }
inline double low(double x) { return x; }
inline double high(double x) { return x; }
template<Scalar T> T low(const Interval<T>& x) { return x.lo; }
template<Scalar T> T high(const Interval<T>& x) { return x.hi; }

// B_2k in R: the interval, or its middle in a double
template<class R> R bernoulli_as(int k) {
    if constexpr (is_interval<R>) return bernoulli_number<decltype(R{}.lo)>(k);
    else return bernoulli_number<R>(k).mid();
}

// n^-s
template<class R> Complex<R> power_minus(long n, const Complex<R>& s) {
    using std::log, std::exp, std::cos, std::sin;
    const R l = log(R(n)), a = exp(-(s.re * l)), angle = s.im * l;
    return {a * cos(angle), -(a * sin(angle))};
}

// Im log Gamma(z) for Re z > 0, on the branch real on the positive axis: moved right until |z| is
// at least 16 by log Gamma(z) = log Gamma(z + 1) - log z, then Stirling's series to ten terms. Its
// remainder is at most the first term left out times sec^{2m+2}(arg z / 2), which is below 2^{m+1}
// where Re z > 0.
template<Scalar T> Interval<T> im_log_gamma(Complex<Interval<T>> z) {
    using I = Interval<T>;
    constexpr int m = 10;
    I shifted{0};
    while (modulus(z).lo < T(16)) shifted = shifted - arg_right(z), z.re = z.re + I(1);
    const I a = arg_right(z), l = log(square(z.re) + square(z.im)) / I(2);
    I im = (z.re - I(T(0.5))) * a + z.im * l - z.im;
    const Complex<I> w = Complex<I>{I(1)} / z, w2 = w * w;
    Complex<I> p = w;
    for (int j = 1; j <= m; ++j) {
        im = im + bernoulli_number<T>(j) / I(2 * j * (2 * j - 1)) * p.im;
        p = p * w2;
    }
    const T r = (abs(bernoulli_number<T>(m + 1)) / I((2 * m + 2) * (2 * m + 1)) * modulus(p) * I(1 << (m + 1))).hi;
    return im + I(-r, r) + shifted;
}

} // namespace zeta_detail

// zeta(s): the sum over n < N of n^-s, the integral and half the term at N, then 29 terms of
// Euler-Maclaurin. N is (|s| + 59) / pi or more, so the terms fall by four each. In intervals the
// value is proved: every operation rounded outward and the remainder after the 29 terms bounded by
// |s + 59| / (Re s + 59) times the first term left out. In doubles it is the same sum, for a
// picture. Defined for Re s > -58 and s away from 1; elsewhere every value is allowed (NaN in a
// double).
template<Scalar R>
Complex<R> zeta(const Complex<R>& s) {
    using C = Complex<R>;
    using zeta_detail::modulus, zeta_detail::low, zeta_detail::high;
    constexpr int M = 29;
    if (!(low(s.re) > -2 * M)) {
        if constexpr (zeta_detail::is_interval<R>) return {R::entire(), R::entire()};
        else return {std::numeric_limits<R>::quiet_NaN(), std::numeric_limits<R>::quiet_NaN()};
    }
    const double size = static_cast<double>(high(modulus(s)));
    const long N = std::max(10L, static_cast<long>(std::ceil((size + 2 * M + 1) / 3.14159)) + 1);
    C sum{R(0)};
    for (long n = 1; n < N; ++n) sum += zeta_detail::power_minus(N - n, s);   // smallest first
    const R n_(N);
    const C p = zeta_detail::power_minus(N, s);    // N^-s
    sum += p * n_ / (s - C{R(1)}) + p * R(0.5);
    C rising = s;                                   // s (s + 1) ... (s + 2k - 2)
    R size_rising = modulus(s);                     // its modulus, the factors' moduli multiplied: its square overflows
    R inverse = R(1) / n_;                          // N^{1 - 2k}
    const R inverse2 = inverse * inverse;
    R factorial(2);                                 // (2k)!
    for (int k = 1; k <= M; ++k) {
        sum += rising * p * (zeta_detail::bernoulli_as<R>(k) / factorial * inverse);
        rising = rising * (s + C{R(2 * k - 1)}) * (s + C{R(2 * k)});
        if constexpr (zeta_detail::is_interval<R>)
            size_rising = size_rising * modulus(s + C{R(2 * k - 1)}) * modulus(s + C{R(2 * k)});
        inverse = inverse * inverse2;
        factorial = factorial * R((2 * k + 1) * (2 * k + 2));
    }
    if constexpr (zeta_detail::is_interval<R>) {
        using std::abs;
        const R tail = size_rising * modulus(s + C{R(2 * M + 1)}) * abs(zeta_detail::bernoulli_as<R>(M + 1)) / factorial *
                       modulus(p) * inverse / (R(s.re.lo) + R(2 * M + 1));
        const auto r = tail.hi;
        sum.re = sum.re + R(-r, r);
        sum.im = sum.im + R(-r, r);
    }
    return sum;
}

// The Riemann-Siegel theta: Im log Gamma(1/4 + it/2) - (t/2) log pi.
template<Scalar T>
Interval<T> riemann_siegel_theta(const Interval<T>& t) {
    using I = Interval<T>;
    return zeta_detail::im_log_gamma(Complex<I>{I(T(0.25)), t / I(2)}) - t / I(2) * log(zeta_detail::pi<T>());
}

// Hardy's Z(t) = e^{i theta(t)} zeta(1/2 + it), real, |Z(t)| = |zeta(1/2 + it)|.
template<Scalar T>
Interval<T> hardy_z(const Interval<T>& t) {
    using I = Interval<T>;
    const I th = riemann_siegel_theta(t);
    const Complex<I> z = zeta(Complex<I>{I(T(0.5)), t});
    return cos(th) * z.re - sin(th) * z.im;
}

// A Gram point g_n, where theta(g_n) = n pi; past t = 7 theta rises, so there is one for each
// n >= -1. Gram's law holds there when (-1)^n Z(g_n) > 0: the point is good.
template<Scalar T = double>
struct GramPoint {
    long n = 0;
    Interval<T> t{};
    Truth good = Truth::Undecided;
};

template<Scalar T>
std::vector<GramPoint<T>> gram_points(T t0, T t1) {
    using I = Interval<T>;
    using std::log, std::abs;
    std::vector<GramPoint<T>> out;
    const T pi = interval_detail::pi<T>();
    const T from = t0 < T(7) ? T(7) : t0;
    const auto theta = [](T t) { return riemann_siegel_theta(I(t)).mid(); };
    // theta'(t) is near log(t / 2 pi) / 2; Newton's step on it, in doubles, then the point proved
    // between two numbers where theta is below and above n pi
    long n = static_cast<long>(std::floor(static_cast<double>(theta(from) / pi))) + 1;
    T g = from;
    for (;;) {
        const I target = I(n) * zeta_detail::pi<T>();
        for (int i = 0; i < 60; ++i) {
            const T slope = log((g < T(10) ? T(10) : g) / (T(2) * pi)) / T(2);
            const T step = (theta(g) - T(n) * pi) / slope;
            g = g - step;
            if (abs(step) < T(1e-14) * g) break;
        }
        if (t1 < g) break;
        T d = T(4) * std::numeric_limits<T>::epsilon() * g;
        GramPoint<T> p{n, I(g), Truth::Undecided};
        for (int i = 0; i < 12; ++i, d = d * T(8)) {
            if (riemann_siegel_theta(I(g - d)).hi < target.lo && target.hi < riemann_siegel_theta(I(g + d)).lo) {
                p.t = I(g - d, g + d);
                const I z = hardy_z(p.t);
                p.good = n % 2 == 0 ? less(I(0), z) : less(z, I(0));
                break;
            }
        }
        if (t0 <= g) out.push_back(p);
        ++n;
        g = g + T(2) * pi / log((g < T(10) ? T(10) : g) / (T(2) * pi));
    }
    return out;
}

// The zeros of zeta on the critical line between t0 and t1, each the ordinate t of a zero
// 1/2 + it in an interval where Z changes sign. `counted` is where their number is known: N(T) at
// both ends by Backlund's bound; `complete` when as many were found there as N says, so these are
// every zero of zeta with its ordinate there, each simple and on the line.
template<Scalar T = double>
struct CriticalZeros {
    std::vector<Interval<T>> zeros;
    Interval<T> counted = Interval<T>::empty();
    long expected = -1;   // N at the top of `counted` minus N at its bottom, or -1 not known
    bool complete = false;
    std::size_t evaluations = 0;   // of Z
};

namespace zeta_detail {

// N(T) = the number of zeros with 0 < Im < T, or -1 when Re zeta is not proved positive on
// [1/2, 2] + iT (beyond 2 it is: |zeta - 1| <= zeta(2) - 1 < 1).
template<Scalar T> long zeros_below(T height, std::size_t evaluations) {
    using I = Interval<T>;
    const Box<T, 1> segment{I(T(0.5), T(2))};
    const auto proof = prove_on(segment, [&](const Box<T, 1>& b) {
        return less(I(0), zeta(Complex<I>{b[0], I(height)}).re);
    }, Budget{.evaluations = evaluations, .smallest = 0x1p-40});
    if (proof.verdict != Truth::Yes) return -1;
    const I x = riemann_siegel_theta(I(height)) / pi<T>() + I(1);
    using std::round;
    const T k = round(x.mid());
    return x.hi - T(0.5) < k && k < x.lo + T(0.5) ? static_cast<long>(k) : -1;
}

template<Scalar T> T spacing(T t) {
    using std::log;
    const T pi = interval_detail::pi<T>();
    return T(2) * pi / log((t < T(20) ? T(20) : t) / (T(2) * pi));
}

} // namespace zeta_detail

template<Scalar T>
CriticalZeros<T> critical_zeros(T t0, T t1, std::size_t evaluations = 200000) {
    using I = Interval<T>;
    CriticalZeros<T> out;
    const std::size_t per_end = 4096;
    // the ends: below t0 and above t1, or at 10, below the first zero, when t0 is lower; of twelve
    // heights a third of a spacing apart, those where Re zeta(1/2 + iT) is largest tried first, as
    // Backlund's bound holds most easily there
    const auto end = [&](T from, T direction, T& at) -> long {
        std::vector<std::pair<T, T>> heights;
        for (int j = 0; j < 12; ++j) {
            const T h = from + direction * T(j) * zeta_detail::spacing(from) / T(3);
            if (h >= T(1)) heights.emplace_back(-zeta(Complex<double>{0.5, static_cast<double>(h)}).re, h);
        }
        std::sort(heights.begin(), heights.end());
        for (const auto& [minus_re, h] : heights) {
            if (!(minus_re < T(0))) break;
            if (const long n = zeta_detail::zeros_below(h, per_end); n >= 0) return at = h, n;
        }
        return -1;
    };
    T lo = t0, hi = t1;
    long below = -1;
    if (t0 < T(10)) {
        if (zeta_detail::zeros_below(T(10), per_end) == 0) lo = T(0), below = 0;
    } else {
        below = end(t0, T(-1), lo);
    }
    const long above = end(t1, T(1), hi);
    if (below >= 0 && above >= 0) out.expected = above - below, out.counted = I(lo, hi);
    const T start = below == 0 && lo == T(0) ? T(10) : lo;
    // Z's sign on a grid, finer until the changes are as many as expected
    const auto sign = [&](T t) {
        ++out.evaluations;
        const I z = hardy_z(I(t));
        return z.lo > T(0) ? 1 : z.hi < T(0) ? -1 : 0;
    };
    std::vector<std::pair<T, T>> brackets;
    for (T step = zeta_detail::spacing(hi) / T(8);; step = step / T(2)) {
        brackets.clear();
        T last = start;
        int s = sign(start);
        for (T t = start + step; ; t = t + step) {
            const T at = t < hi ? t : hi;
            const int z = sign(at);
            if (z != 0) {
                if (s != 0 && z != s) brackets.emplace_back(last, at);
                s = z, last = at;
            }
            if (!(at < hi)) break;
        }
        if (out.expected < 0 || static_cast<long>(brackets.size()) >= out.expected || out.evaluations > evaluations) break;
    }
    // each change narrowed by halving while Z's sign at the middle is proved
    for (auto [a, b] : brackets) {
        int sa = sign(a);
        for (int i = 0; i < 64; ++i) {
            const T m = a / T(2) + b / T(2);
            if (!(a < m && m < b)) break;
            const int sm = sign(m);
            if (sm == 0) break;
            (sm == sa ? a : b) = m;
        }
        if (t0 <= b && a <= t1) out.zeros.push_back(I(a, b));
    }
    out.complete = out.expected >= 0 && static_cast<long>(brackets.size()) == out.expected;
    return out;
}

} // namespace algebra
} // namespace spatium
