#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/complex.hpp>
#  include <spatium/algebra/dual.hpp>
#  include <spatium/algebra/quadrature.hpp>
#  include <spatium/algebra/series.hpp>
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

template<class R> inline constexpr bool is_interval = false;
template<class T> inline constexpr bool is_interval<Interval<T>> = true;

// pi and x^2 in R: an interval's are proved, the others' are a number's
template<class R> R pi_as() {
    if constexpr (is_interval<R>) return pi<decltype(R{}.lo)>();
    else return R(3.14159265358979323846);
}
template<class R> R sq(const R& x) {
    if constexpr (is_interval<R>) return square(x);
    else return x * x;
}

inline double modulus(const Complex<double>& z) { return std::hypot(z.re, z.im); }
inline double low(double x) { return x; }
inline double high(double x) { return x; }
template<Scalar T> T low(const Interval<T>& x) { return x.lo; }
template<Scalar T> T high(const Interval<T>& x) { return x.hi; }
// a dual number, to differentiate the sum in doubles: its value decides as the double does
template<class T> double low(const Dual<T>& x) { return low(x.value); }
template<class T> double high(const Dual<T>& x) { return high(x.value); }
template<class T> double modulus(const Complex<Dual<T>>& z) { return std::hypot(low(z.re), low(z.im)); }

// B_2k in R: the interval, or its middle as a number of R
template<class R> R bernoulli_as(int k) {
    if constexpr (is_interval<R>) return bernoulli_number<decltype(R{}.lo)>(k);
    else return R(bernoulli_number<double>(k).mid());
}

// n^-s
template<class R> Complex<R> power_minus(long n, const Complex<R>& s) {
    using std::log, std::exp, std::cos, std::sin;
    const R l = log(R(n)), a = exp(-(s.re * l)), angle = s.im * l;
    return {a * cos(angle), -(a * sin(angle))};
}

// log Gamma(z) for Re z > 0, on the branch real on the positive axis: moved right until |z| is at
// least 16 by log Gamma(z) = log Gamma(z + 1) - log z, then Stirling's series to ten terms. Its
// remainder is at most the first term left out times sec^{2m+2}(arg z / 2), which is below 2^{m+1}
// where Re z > 0; in intervals it is added to both parts.
template<class R> Complex<R> log_gamma(Complex<R> z) {
    using std::atan, std::log;
    using C = Complex<R>;
    constexpr int m = 10;
    const auto log_of = [](const C& w) { return C{log(sq(w.re) + sq(w.im)) / R(2), atan(w.im / w.re)}; };
    C shifted{R(0)};
    while (low(modulus(z)) < 16) shifted -= log_of(z), z.re = z.re + R(1);
    C out = (z - C{R(0.5)}) * log_of(z) - z + C{log(R(2) * pi_as<R>()) / R(2)};
    const C w = C{R(1)} / z, w2 = w * w;
    C p = w;
    for (int j = 1; j <= m; ++j) {
        out += p * (bernoulli_as<R>(j) / R(2 * j * (2 * j - 1)));
        p = p * w2;
    }
    if constexpr (is_interval<R>) {
        using std::abs;
        const auto r = (abs(bernoulli_as<R>(m + 1)) / R((2 * m + 2) * (2 * m + 1)) * modulus(p) * R(1 << (m + 1))).hi;
        out += C{R(-r, r), R(-r, r)};
    }
    return out + shifted;
}

// zeta(s): the sum over n < N of n^-s, the integral and half the term at N, then 29 terms of
// Euler-Maclaurin. N is (|s| + 59) / pi or more, so the terms fall by four each. In intervals the
// remainder after the 29 terms is bounded by |s + 59| / (Re s + 59) times the first term left out.
// Defined for Re s > -58; left of the line its terms cancel, so `zeta` takes it only for Re s >= 0.
template<class R>
Complex<R> euler_maclaurin(const Complex<R>& s) {
    using C = Complex<R>;
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
    // the k-th term is B_2k r_k, r_k = s (s + 1) ... (s + 2k - 2) N^{-s+1-2k} / (2k)!, each r from the
    // one before by the ratio of the two: the rising product and the power of N, kept apart, leave
    // the doubles past t = 1e5 (one overflows, the other underflows)
    const R n2 = n_ * n_;
    C r = s * p / (n_ * R(2));
    R size_r = modulus(s) * modulus(p) / (n_ * R(2));   // |r_k| as the factors' moduli multiplied
    for (int k = 1; k <= M; ++k) {
        sum += r * zeta_detail::bernoulli_as<R>(k);
        const R step = R((2 * k + 1) * (2 * k + 2)) * n2;
        r = r * (s + C{R(2 * k - 1)}) * (s + C{R(2 * k)}) / step;
        if constexpr (zeta_detail::is_interval<R>)
            size_r = size_r * modulus(s + C{R(2 * k - 1)}) * modulus(s + C{R(2 * k)}) / step;
    }
    if constexpr (zeta_detail::is_interval<R>) {
        using std::abs;
        // after the M terms: |r_{M+1}| |s + 2M + 1| |B_{2M+2}| / (Re s + 2M + 1)
        const R tail = size_r * modulus(s + C{R(2 * M + 1)}) * abs(zeta_detail::bernoulli_as<R>(M + 1)) /
                       (R(s.re.lo) + R(2 * M + 1));
        const auto r = tail.hi;
        sum.re = sum.re + R(-r, r);
        sum.im = sum.im + R(-r, r);
    }
    return sum;
}

// zeta(s) = chi(s) zeta(1 - s), chi(s) = 2^s pi^{s-1} sin(pi s / 2) Gamma(1 - s), for Re s < 0, where
// zeta(1 - s) is the sum's to take. sin(pi s / 2) grows as e^{pi |t| / 2} and Gamma(1 - s) falls as
// much, so the factor is moved from one to the other: sin(pi s / 2) e^{-pi |t| / 2} is two
// exponentials of at most 1, and the rest of chi one exponential of its logarithm. At s = -2n the
// sine is zero: the trivial zeros.
template<class R>
Complex<R> reflected(const Complex<R>& s) {
    using std::abs, std::exp, std::log, std::cos, std::sin;
    using C = Complex<R>;
    // pi s / 2 less k pi, k the whole number nearest s / 2, subtracted where that is exact; the sine
    // changes sign with each pi taken. So at s = -2n it is the sine at 0, not at a rounded -n pi.
    const R half = s.re / R(2);
    const double k = std::round(low(half));
    const R pi = pi_as<R>(), height = abs(s.im), h = pi * (half - R(k));
    // (e^{i pi s / 2} - e^{-i pi s / 2}) / 2i, each exponential times e^{-pi |t| / 2}
    const R up = exp(-(pi * (s.im + height)) / R(2)), down = exp(pi * (s.im - height) / R(2));
    const C d = C{cos(h) * up, sin(h) * up} - C{cos(h) * down, -(sin(h) * down)};
    const R odd = R(std::fmod(k, 2.0) == 0 ? 0.5 : -0.5);
    const C sine{d.im * odd, -(d.re * odd)};
    const C reflection{R(1) - s.re, -s.im};
    const C lg = log_gamma(reflection);
    const R l2 = log(R(2)), lp = log(pi);
    const R mag = exp(s.re * l2 + (s.re - R(1)) * lp + lg.re + pi * height / R(2)), angle = s.im * (l2 + lp) + lg.im;
    return C{mag * cos(angle), mag * sin(angle)} * sine * euler_maclaurin(reflection);
}

} // namespace zeta_detail

// zeta(s), everywhere but the pole at 1: Euler-Maclaurin for Re s >= 0, and left of that the
// functional equation, which hands the sum 1 - s. In intervals the value is proved: every operation
// rounded outward, every remainder bounded. In doubles it is the same computation, for a picture.
template<Scalar R>
Complex<R> zeta(const Complex<R>& s) {
    if (zeta_detail::high(s.re) < 0) return zeta_detail::reflected(s);
    return zeta_detail::euler_maclaurin(s);
}

// The Riemann-Siegel theta: Im log Gamma(1/4 + it/2) - (t/2) log pi.
template<Scalar T>
Interval<T> riemann_siegel_theta(const Interval<T>& t) {
    using I = Interval<T>;
    return zeta_detail::log_gamma(Complex<I>{I(T(0.25)), t / I(2)}).im - t / I(2) * log(zeta_detail::pi<T>());
}

// Hardy's Z(t) = e^{i theta(t)} zeta(1/2 + it), real, |Z(t)| = |zeta(1/2 + it)|.
template<Scalar T>
Interval<T> hardy_z(const Interval<T>& t) {
    using I = Interval<T>;
    const I th = riemann_siegel_theta(t);
    const Complex<I> z = zeta(Complex<I>{I(T(0.5)), t});
    return cos(th) * z.re - sin(th) * z.im;
}

// ── Z(t) by the Riemann-Siegel formula ────────────────────────────
//
// Z(t) = 2 sum_{n <= N} n^{-1/2} cos(theta(t) - t log n)
//        + (-1)^{N-1} (t / 2 pi)^{-1/4} sum_{k <= 4} C_k(p) (t / 2 pi)^{-k/2},
// N = floor(sqrt(t / 2 pi)) and p its fraction, the C_k from Psi(p) = cos(2 pi (p^2 - p - 1/16)) /
// cos(2 pi p) and its derivatives to the twelfth (Edwards 1974, 7.4). sqrt(t / 2 pi) terms where
// Euler-Maclaurin takes t / pi. In doubles, computed, not proved: what is left after C_4 is of order
// t^{-11/4}, and the phases t log n carry t times the rounding of a double.

// theta(t) in doubles, from the same Stirling series
inline double riemann_siegel_theta(double t) {
    return zeta_detail::log_gamma(Complex<double>{0.25, t / 2}).im - t / 2 * std::log(3.14159265358979323846);
}

namespace zeta_detail {

// x less the nearest multiple of 2 pi, in long double, where 2 pi has 64 bits: a phase t log n near
// 3e9 keeps its angle to about 3e-10, where in a double it kept 1e-6 (measured against mpmath at 30
// digits: Z off by 6e-6 at t = 2.6e8).
inline double reduced(long double x) {
    constexpr long double two_pi = 6.283185307179586476925286766559005768L;
    return static_cast<double>(x - two_pi * std::nearbyint(x / two_pi));
}

// theta(t) less a multiple of 2 pi: for t >= 200 its asymptotic series to t^{-9} in long double (what
// is left is below 1e-25 there), below that from the Stirling series in doubles
inline double theta_reduced(double t) {
    if (!(t >= 200)) return reduced(riemann_siegel_theta(t));
    constexpr long double pi = 3.141592653589793238462643383279502884L;
    const long double T = t, u = 1 / (T * T);
    const long double tail = (1.0L / 48 + u * (7.0L / 5760 + u * (31.0L / 80640 + u * (127.0L / 430080 + u * 511.0L / 1216512)))) / T;
    return reduced(T / 2 * std::log(T / (2 * pi)) - T / 2 - pi / 8 + tail);
}

} // namespace zeta_detail

namespace zeta_detail {

// sin(y) / y of a series, by its Taylor series in y^2: no division, so nothing lost where y is small
template<int N> Series<double, N> sinc(const Series<double, N>& y) {
    const Series<double, N> y2 = y * y;
    std::array<double, 16> c{};
    c[0] = 1;
    for (int k = 1; k < 16; ++k) c[k] = -c[k - 1] / ((2 * k) * (2 * k + 1));
    Series<double, N> out(c[15]);
    for (int k = 14; k >= 0; --k) out = out * y2 + Series<double, N>(c[k]);
    return out;
}

// Psi and its derivatives to the twelfth at p. Psi is symmetric about 1/2, and on [0, 1/2], with
// u = p - 1/4, it is sin(pi u (1 - 2u)) / sin(2 pi u) = (1/2 - u) sinc(pi u (1 - 2u)) / sinc(2 pi u):
// the two zeros that cancel at p = 1/4 taken out, the sinc below at least 2 / pi. Its Taylor series
// at u gives the derivatives; dividing Psi's series at 1/2 instead loses 2^n of the n-th.
inline std::array<double, 13> psi_derivatives(double p) {
    using S = Series<double, 13>;
    const bool right = p > 0.5;
    const S u = S::variable((right ? 1 - p : p) - 0.25);
    const double pi = 3.14159265358979323846;
    const S psi = (S(0.5) - u) * sinc(u * (S(1.0) - S(2.0) * u) * S(pi)) / sinc(u * S(2 * pi));
    std::array<double, 13> d{};
    double factorial = 1;
    for (int n = 0; n <= 12; ++n) {
        if (n > 0) factorial *= n;
        d[n] = psi.coefficient(n) * factorial * (right && n % 2 ? -1 : 1);
    }
    return d;
}

} // namespace zeta_detail

namespace zeta_detail {

// C_0 .. C_4 at p, from Psi's derivatives
inline std::array<double, 5> riemann_siegel_coefficients_exact(double p) {
    const double pi = 3.14159265358979323846, pi2 = pi * pi, pi4 = pi2 * pi2, pi6 = pi4 * pi2, pi8 = pi4 * pi4;
    const auto d = psi_derivatives(p);
    return {d[0], -d[3] / (96 * pi2), d[2] / (64 * pi2) + d[6] / (18432 * pi4),
            -d[1] / (64 * pi2) - d[5] / (3840 * pi4) - d[9] / (5308416 * pi6),
            d[0] / (128 * pi2) + 19 * d[4] / (24576 * pi4) + 11 * d[8] / (5898240 * pi6) + d[12] / (2038431744 * pi8)};
}

// The same from a table of 4096 steps over [0, 1], built once, by Lagrange's cubic through the four
// nearest nodes: the derivatives to the twelfth cost 8 us a point, the formula's other terms a tenth
// of that, so a moment at t = 1e8 spent half its time here. The error is of h^4, below 1e-13.
inline std::array<double, 5> riemann_siegel_coefficients(double p) {
    constexpr int steps = 4096;
    static const auto table = [] {
        std::vector<std::array<double, 5>> out(steps + 3);
        for (int i = 0; i < steps + 3; ++i) out[i] = riemann_siegel_coefficients_exact((i - 1) / double(steps));
        return out;
    }();
    const double x = p * steps;
    const int i = std::clamp(static_cast<int>(x), 0, steps - 1);
    const double u = x - i;   // between node i and i + 1, which are table rows i + 1 and i + 2
    const double w0 = -u * (u - 1) * (u - 2) / 6, w1 = (u + 1) * (u - 1) * (u - 2) / 2, w2 = -(u + 1) * u * (u - 2) / 2,
                 w3 = (u + 1) * u * (u - 1) / 6;
    std::array<double, 5> c{};
    for (int k = 0; k < 5; ++k) c[k] = w0 * table[i][k] + w1 * table[i + 1][k] + w2 * table[i + 2][k] + w3 * table[i + 3][k];
    return c;
}

// What the formula adds to 2 sum_{n <= N} n^{-1/2} cos(theta - t log n):
// (-1)^{N-1} (t / 2 pi)^{-1/4} sum_{k <= 4} C_k(p) (t / 2 pi)^{-k/2}
inline double riemann_siegel_rest(double t) {
    const double pi = 3.14159265358979323846, a = std::sqrt(t / (2 * pi));
    const long N = static_cast<long>(a);
    const auto c = riemann_siegel_coefficients(a - static_cast<double>(N));
    double rest = 0;
    for (int k = 4; k >= 0; --k) rest = rest / a + c[k];
    return (N % 2 ? 1 : -1) * rest / std::sqrt(a);
}

inline long riemann_siegel_terms(double t) { return static_cast<long>(std::sqrt(t / (2 * 3.14159265358979323846))); }

} // namespace zeta_detail

// A bound on |Z(t) - riemann_siegel_z(t)| from where the formula is cut, for t >= 200: what is left
// after C_2 is at most 0.011 t^{-7/4} (Gabcke 1979, as Gourdon 2004 quotes it), and the sum here goes
// on with the C_3 and C_4 terms, which are added to the bound as they are. The rounding of the doubles
// is not in it. Infinite below 200, where the bound is not known.
inline double riemann_siegel_truncation(double t) {
    if (!(t >= 200)) return std::numeric_limits<double>::infinity();
    const double pi = 3.14159265358979323846, a = std::sqrt(t / (2 * pi));
    const long N = static_cast<long>(a);
    const auto c = zeta_detail::riemann_siegel_coefficients(a - static_cast<double>(N));
    const double beyond = std::abs((c[3] + c[4] / a) / (a * a * a)) / std::sqrt(a);
    return 0.011 * std::pow(t, -1.75) + beyond;
}

inline double riemann_siegel_z(double t) {
    const long N = zeta_detail::riemann_siegel_terms(t);
    // theta - t log n reduced in long double, theta itself already reduced
    const long double theta = zeta_detail::theta_reduced(t), T = t;
    double sum = 0;
    for (long n = 1; n <= N; ++n)
        sum += std::cos(zeta_detail::reduced(theta - T * std::log(static_cast<long double>(n)))) / std::sqrt(static_cast<double>(n));
    return 2 * sum + zeta_detail::riemann_siegel_rest(t);
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

// ── zeta' and what it shows ───────────────────────────────────────
//
// zeta and its first two derivatives in doubles: the same sum run in a dual number of dual
// numbers along Re s, whose derivative is the complex one since zeta is analytic. Pictures'
// values, as zeta in doubles is: the remainder bound is zeta's, not its derivatives'.
struct ZetaJet {
    Complex<double> value, first, second;
};

inline ZetaJet zeta_jet(Complex<double> s) {
    using D1 = Dual<double>;
    using D2 = Dual<D1>;
    const auto z = zeta(Complex<D2>{D2(D1(s.re, 1), D1(1, 0)), D2(D1(s.im, 0), D1(0, 0))});
    return {{z.re.value.value, z.im.value.value}, {z.re.value.deriv, z.im.value.deriv}, {z.re.deriv.deriv, z.im.deriv.deriv}};
}

// Newton's method on zeta' from s: s - zeta'(s) / zeta''(s) until a step is below 1e-13 of where
// it stands. Where it settles and in how many steps is one pixel of the Newton fractal: a start in
// the core of a basin settles in a few, one near its boundary wanders.
// A run that leaves the strip |Re s| <= 12 has `escaped` (+1 to the right, -1 to the left): right of
// it zeta' is near -log 2 / 2^s, no zero, and each step moves 1 / log 2 further right.
struct NewtonRun {
    Complex<double> at;
    int steps = 0;
    bool settled = false;
    int escaped = 0;
};

inline NewtonRun zeta_prime_newton(Complex<double> s, int most = 64) {
    NewtonRun out{s};
    for (; out.steps < most; ++out.steps) {
        const ZetaJet j = zeta_jet(out.at);
        const Complex<double> step = j.first / j.second;
        if (!std::isfinite(step.re) || !std::isfinite(step.im)) return out;
        out.at = out.at - step;
        if (std::abs(out.at.re) > 12) {
            out.escaped = out.at.re > 0 ? 1 : -1;
            ++out.steps;
            return out;
        }
        if (zeta_detail::modulus(step) <= 1e-13 * std::max(1.0, zeta_detail::modulus(out.at))) {
            out.settled = std::abs(out.at.re - 1) > 1e-3 || std::abs(out.at.im) > 1e-3;   // not at the pole
            ++out.steps;
            return out;
        }
    }
    return out;
}

// The zeros of zeta' in a box, found by Newton from a grid of starts `per_unit` to a unit of the
// box's side, each once (two within 1e-9 are one): computed, not proved, and not proved to be all.
inline std::vector<Complex<double>> zeta_prime_zeros(Complex<double> lo, Complex<double> hi, int per_unit = 4) {
    std::vector<Complex<double>> out;
    const int w = std::max(2, static_cast<int>((hi.re - lo.re) * per_unit)), h = std::max(2, static_cast<int>((hi.im - lo.im) * per_unit));
    for (int j = 0; j < h; ++j)
        for (int i = 0; i < w; ++i) {
            const NewtonRun r = zeta_prime_newton({lo.re + (hi.re - lo.re) * (i + 0.5) / w, lo.im + (hi.im - lo.im) * (j + 0.5) / h});
            const auto& z = r.at;
            if (!r.settled || z.re < lo.re || z.re > hi.re || z.im < lo.im || z.im > hi.im) continue;
            if (std::ranges::none_of(out, [&](const Complex<double>& y) { return zeta_detail::modulus(y - z) < 1e-9; })) out.push_back(z);
        }
    std::ranges::sort(out, {}, [](const Complex<double>& z) { return z.im; });
    return out;
}

// The integral of |zeta(1/2 + it)|^{2k} over [t0, t1], in doubles: Gauss-Kronrod on each unit of t,
// |zeta| = |Z| by Riemann-Siegel from t = 1000 (where it meets Euler-Maclaurin to 1e-10) and by
// Euler-Maclaurin below. Each unit's tolerance is `relative` times the mean of the units before it:
// the quadrature's tolerance is absolute, Z carries a rounding of about 1e-11, so one fixed at the
// rule's own never settles, and one from the unit just before falls with it where |zeta| passes near
// 0 (2100 evaluations a unit, against 82). Computed, its error estimated, not proved.
struct ZetaMoment {
    double value = 0, error = 0;
    long evaluations = 0;
};

inline ZetaMoment zeta_moment(int k, double t0, double t1, double relative = 1e-10) {
    ZetaMoment out;
    const auto f = [k](double t) {
        const double z = t >= 1000 ? riemann_siegel_z(t) : zeta_detail::modulus(zeta(Complex<double>{0.5, t}));
        return std::pow(std::abs(z), 2 * k);
    };
    // the scale before any unit is done: one pass of the rule with no tolerance to meet
    const double first = integrate(f, Finite<double>{t0, std::min(t0 + 1, t1)}, IntegralOptions{.tolerance = 1e300, .method = Method::GaussKronrod}).value;
    long units = 0;
    for (double a = t0; a < t1; a += 1, ++units) {
        const double mean = std::max(units > 0 ? out.value / units : first, 1.0);
        const auto r = integrate(f, Finite<double>{a, std::min(a + 1, t1)}, IntegralOptions{.tolerance = relative * mean, .method = Method::GaussKronrod});
        out.value += r.value, out.error += r.error_estimate, out.evaluations += r.evaluations;
    }
    return out;
}

// M_k(T), over [0, T]
inline ZetaMoment zeta_moment(int k, double T) { return zeta_moment(k, 0, T); }

namespace zeta_detail {

// Gauss-Legendre with m points on [0, 1]: the roots of P_m by Newton's method from Tricomi's
// first guess, and the weights from P_m' there
inline void gauss_legendre(int m, std::vector<double>& x, std::vector<double>& w) {
    const double pi = 3.14159265358979323846;
    x.assign(m, 0), w.assign(m, 0);
    for (int i = 0; i < m; ++i) {
        double z = std::cos(pi * (i + 0.75) / (m + 0.5)), dp = 0;
        for (int it = 0; it < 100; ++it) {
            double p0 = 1, p1 = z;
            for (int j = 2; j <= m; ++j) {
                const double p2 = ((2 * j - 1) * z * p1 - (j - 1) * p0) / j;
                p0 = p1, p1 = p2;
            }
            dp = m * (z * p1 - p0) / (z * z - 1);
            const double step = p1 / dp;
            z -= step;
            if (std::abs(step) < 1e-16) break;
        }
        x[i] = (1 - z) / 2, w[i] = 1 / ((1 - z * z) * dp * dp);
    }
}

} // namespace zeta_detail

// The same integral by a fixed rule, for heights where Riemann-Siegel holds (t0 >= 1000): Gauss-
// Legendre with `nodes` points on each unit of t, 0 to take k log(t1 / 2 pi) + 16, which is past the
// band |Z|^{2k} oscillates in (its frequencies are below k log(t / 2 pi) a unit). Every node steps
// by exactly 1 from one unit to the next, so each phase t log n of the sum goes on by a product with
// n^{-i} instead of a cosine, and is taken afresh every 64 units. Sixteen units are taken in one pass
// over the terms, each term's phases turned sixteen times while they are in the nearest cache: one
// unit a pass read every phase from memory each time, and twelve threads waited on it (a block 3.5
// times as long as alone). The loops over a term's nodes have no sum across them, so they run in
// vector lanes. A last part of a unit is summed directly. `error` bounds what the cut of the
// Riemann-Siegel formula costs: with |Z - Z~| <= d (`riemann_siegel_truncation`), |Z^{2k} - Z~^{2k}| <=
// 2k (|Z~| + d)^{2k-1} d, summed with the rule's weights. The rule's own error and the rounding are not
// in it: those are held by the same rule at twice the points, and by phases in wider arithmetic.
inline ZetaMoment zeta_moment_fixed(int k, double t0, double t1, int nodes = 0) {
    ZetaMoment out;
    if (!(t0 >= 1000) || !(t0 < t1)) return out;
    const double pi = 3.14159265358979323846;
    const int m = nodes > 0 ? nodes : static_cast<int>(k * std::log(t1 / (2 * pi))) + 16;
    std::vector<double> x, w;
    zeta_detail::gauss_legendre(m, x, w);
    const long top = zeta_detail::riemann_siegel_terms(t1) + 1;
    std::vector<double> amp(top + 1), rc(top + 1), rs(top + 1);
    std::vector<long double> lg(top + 1);   // the seeds' phases t log n, reduced in long double
    for (long n = 1; n <= top; ++n) {
        lg[n] = std::log(static_cast<long double>(n));
        const double l = static_cast<double>(lg[n]);
        amp[n] = 1 / std::sqrt(static_cast<double>(n)), rc[n] = std::cos(l), rs[n] = std::sin(l);
    }
    constexpr int G = 16;   // units a pass takes; divides 64
    // cos and sin of t log n at the pass's first unit, node j of term n at [n * m + j]
    std::vector<double> c(static_cast<std::size_t>(top + 1) * m), s(c.size()), a(G * m), b(G * m), t(G * m);
    std::vector<long> terms(G * m);
    const long whole = static_cast<long>(std::floor(t1 - t0));
    long seeded = 0;   // the terms whose phases are those of the pass's first unit
    for (long u0 = 0; u0 < whole; u0 += G) {
        const int count = static_cast<int>(std::min<long>(G, whole - u0));
        long most = 0, fewest = top;
        for (int g = 0; g < count; ++g)
            for (int j = 0; j < m; ++j) {
                const int i = g * m + j;
                t[i] = t0 + static_cast<double>(u0 + g) + x[j], terms[i] = zeta_detail::riemann_siegel_terms(t[i]), a[i] = b[i] = 0;
                most = std::max(most, terms[i]), fewest = std::min(fewest, terms[i]);
            }
        if (u0 % 64 == 0) seeded = 0;
        for (long n = seeded + 1; n <= most; ++n)
            for (int j = 0; j < m; ++j) {
                const double phase = zeta_detail::reduced(static_cast<long double>(t[j]) * lg[n]);
                c[n * m + j] = std::cos(phase), s[n * m + j] = std::sin(phase);
            }
        seeded = std::max(seeded, most);
        for (long n = 1; n <= most; ++n) {
            double* cn = &c[n * m];
            double* sn = &s[n * m];
            const double an = amp[n], r = rc[n], q = rs[n];
            for (int g = 0; g < count; ++g) {
                double* ag = &a[g * m];
                double* bg = &b[g * m];
                if (n <= fewest)
                    for (int j = 0; j < m; ++j) ag[j] += an * cn[j], bg[j] += an * sn[j];
                else
                    for (int j = 0; j < m; ++j)
                        if (n <= terms[g * m + j]) ag[j] += an * cn[j], bg[j] += an * sn[j];
                for (int j = 0; j < m; ++j) {
                    const double cj = cn[j], sj = sn[j];
                    cn[j] = cj * r - sj * q, sn[j] = sj * r + cj * q;
                }
            }
        }
        for (int i = 0; i < count * m; ++i) {
            const double th = zeta_detail::theta_reduced(t[i]);
            const double z = 2 * (std::cos(th) * a[i] + std::sin(th) * b[i]) + zeta_detail::riemann_siegel_rest(t[i]);
            const double d = riemann_siegel_truncation(t[i]);
            out.value += w[i % m] * std::pow(std::abs(z), 2 * k);
            out.error += w[i % m] * 2 * k * std::pow(std::abs(z) + d, 2 * k - 1) * d;
        }
        out.evaluations += static_cast<long>(count) * m;
    }
    if (const double base = t0 + static_cast<double>(whole), width = t1 - base; width > 0) {
        double sum = 0, bound = 0;
        for (int j = 0; j < m; ++j) {
            const double tj = base + width * x[j], z = std::abs(riemann_siegel_z(tj)), d = riemann_siegel_truncation(tj);
            sum += w[j] * std::pow(z, 2 * k), bound += w[j] * 2 * k * std::pow(z + d, 2 * k - 1) * d;
        }
        out.value += width * sum, out.error += width * bound, out.evaluations += m;
    }
    return out;
}

} // namespace algebra
} // namespace spatium
