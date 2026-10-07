#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/calculus.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <array>
#  include <cmath>
#  include <cstddef>
#  include <cstdint>
#  include <limits>
#  include <random>
#endif

// The integral of f over a box by sampling, for the dimensions where nested quadrature
// stops being possible (the cost of the door of `integrate` grows as n^K).
//
//   monte_carlo(f, lo, hi)         random points; error ~ 1/sqrt(n)
//   quasi_monte_carlo(f, lo, hi)   a low-discrepancy sequence (Sobol by default, Halton on request),
//                                  randomly shifted; error ~ (log n)^K / n for a smooth f, and the
//                                  estimate comes with its own error
//
// Both return the `IntegralResult` of the rest of the library, with the difference that its
// `error_estimate` is a standard error -- one sigma of the estimate, not a bound: the true
// error is within it about two times in three. The status is `Failed` if f is ever not finite,
// `Suspicious` where there are too few points for a standard error to mean anything, and
// otherwise `Converged`, which for a sampling estimator says only that an estimate with an
// error was formed. `evaluations` is the number of points.
//
// For the quasi-random one the points are a low-discrepancy sequence, randomised in each of
// `shifts` independent copies; the standard error is the spread of the copies' means, which is
// what makes an error estimate possible for a deterministic sequence at all. The default is
// Sobol (Gray-code order, the Joe-Kuo direction numbers, up to 40 dimensions) with a random
// digital shift -- an XOR of every point with one random word per coordinate, which keeps the
// net structure that Sobol's quality is made of. Halton (the first K primes, a uniform shift
// mod 1 as Cranley and Patterson) is kept for comparison, up to 32 dimensions; it degrades
// quickly with dimension, its large bases being poorly distributed on few points.
//
// Reproducible: the generator is a seeded `std::mt19937_64` and uniform numbers are formed
// from its bits, not from `std::uniform_real_distribution`, whose output differs between
// standard libraries.

SPATIUM_EXPORT namespace spatium {

enum class Sequence { Sobol, Halton };

struct MonteCarloOptions {
    std::size_t samples = std::size_t{1} << 16;
    std::uint64_t seed = 1;
    int shifts = 16;            // quasi_monte_carlo: independent random shifts
    Sequence sequence = Sequence::Sobol;
};

namespace monte_carlo_detail {

inline double uniform01(std::mt19937_64& rng) {
    return static_cast<double>(rng() >> 11) * (1.0 / 9007199254740992.0);    // 53 bits
}

// the i-th Halton point's coordinate in base b (radical inverse), i >= 1
inline double radical_inverse(std::uint64_t i, std::uint32_t b) {
    double f = 1.0, r = 0.0;
    const double inv = 1.0 / static_cast<double>(b);
    while (i > 0) {
        f *= inv;
        r += f * static_cast<double>(i % b);
        i /= b;
    }
    return r;
}

inline constexpr std::array<std::uint32_t, 32> kPrimes = {
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53,
    59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131};

// Direction numbers of Joe and Kuo, "Constructing Sobol sequences with better two-dimensional
// projections", SIAM J. Sci. Comput. 30 (2008), file new-joe-kuo-6.21201, dimensions 2 to 40:
// the degree s of the primitive polynomial, its coefficients a (the middle ones, as an integer)
// and the s initial numbers m_i. Dimension 1 is the van der Corput sequence.
struct SobolRow { int s; std::uint32_t a; std::array<std::uint32_t, 8> m; };
inline constexpr std::array<SobolRow, 39> kSobol = {{
    {1, 0, {1}}, {2, 1, {1, 3}}, {3, 1, {1, 3, 1}}, {3, 2, {1, 1, 1}}, {4, 1, {1, 1, 3, 3}},
    {4, 4, {1, 3, 5, 13}}, {5, 2, {1, 1, 5, 5, 17}}, {5, 4, {1, 1, 5, 5, 5}}, {5, 7, {1, 1, 7, 11, 19}},
    {5, 11, {1, 1, 5, 1, 1}}, {5, 13, {1, 1, 1, 3, 11}}, {5, 14, {1, 3, 5, 5, 31}},
    {6, 1, {1, 3, 3, 9, 7, 49}}, {6, 13, {1, 1, 1, 15, 21, 21}}, {6, 16, {1, 3, 1, 13, 27, 49}},
    {6, 19, {1, 1, 1, 15, 7, 5}}, {6, 22, {1, 3, 1, 15, 13, 25}}, {6, 25, {1, 1, 5, 5, 19, 61}},
    {7, 1, {1, 3, 7, 11, 23, 15, 103}}, {7, 4, {1, 3, 7, 13, 13, 15, 69}}, {7, 7, {1, 1, 3, 13, 7, 35, 63}},
    {7, 8, {1, 3, 5, 9, 1, 25, 53}}, {7, 14, {1, 3, 1, 13, 9, 35, 107}}, {7, 19, {1, 3, 1, 5, 27, 61, 31}},
    {7, 21, {1, 1, 5, 11, 19, 41, 61}}, {7, 28, {1, 3, 5, 3, 3, 13, 69}}, {7, 31, {1, 1, 7, 13, 1, 19, 1}},
    {7, 32, {1, 3, 7, 5, 13, 19, 59}}, {7, 37, {1, 1, 3, 9, 25, 29, 41}}, {7, 41, {1, 3, 5, 13, 23, 1, 55}},
    {7, 42, {1, 3, 7, 3, 13, 59, 17}}, {7, 50, {1, 3, 1, 3, 5, 53, 69}}, {7, 55, {1, 1, 5, 5, 23, 33, 13}},
    {7, 56, {1, 1, 7, 7, 1, 61, 123}}, {7, 59, {1, 1, 7, 9, 13, 61, 49}}, {7, 62, {1, 3, 3, 5, 3, 55, 33}},
    {8, 14, {1, 3, 1, 15, 31, 13, 49, 245}}, {8, 21, {1, 3, 5, 15, 31, 59, 63, 97}},
    {8, 22, {1, 3, 1, 11, 11, 11, 77, 249}}}};

inline constexpr std::size_t kMaxSobolDimensions = 40;

// The 32 direction numbers V[1..32] (as 32-bit fractions) of Sobol dimension `dim` (0-based).
inline std::array<std::uint32_t, 33> sobol_directions(std::size_t dim) {
    std::array<std::uint32_t, 33> v{};
    if (dim == 0) {
        for (int i = 1; i <= 32; ++i) v[i] = std::uint32_t{1} << (32 - i);
        return v;
    }
    const SobolRow& row = kSobol[dim - 1];
    const int s = row.s;
    for (int i = 1; i <= s; ++i) v[i] = row.m[static_cast<std::size_t>(i - 1)] << (32 - i);
    for (int i = s + 1; i <= 32; ++i) {
        v[i] = v[i - s] ^ (v[i - s] >> s);
        for (int k = 1; k <= s - 1; ++k)
            if ((row.a >> (s - 1 - k)) & 1u) v[i] ^= v[i - k];
    }
    return v;
}

template<Scalar T>
IntegralResult<T> failed(long evaluations) {
    const T nan = T(std::numeric_limits<double>::quiet_NaN());
    return {nan, nan, evaluations, IntegralStatus::Failed};
}

template<std::size_t K, Scalar T>
T box_volume(const Vec<T, K>& lo, const Vec<T, K>& hi) {
    T v{1};
    for (std::size_t i = 0; i < K; ++i) v = T(v * (hi[i] - lo[i]));
    return v;
}

}  // namespace monte_carlo_detail

template<std::size_t K, Scalar T, typename F>
IntegralResult<T> monte_carlo(F&& f, const Vec<T, K>& lo, const Vec<T, K>& hi, MonteCarloOptions o = {}) {
    using std::sqrt; using std::abs;
    std::mt19937_64 rng(o.seed);
    const T volume = monte_carlo_detail::box_volume<K, T>(lo, hi);
    T mean{0}, m2{0};
    std::size_t n = 0;
    Vec<T, K> x{};
    for (std::size_t s = 0; s < o.samples; ++s) {
        for (std::size_t i = 0; i < K; ++i)
            x[i] = T(lo[i] + (hi[i] - lo[i]) * T(monte_carlo_detail::uniform01(rng)));
        const T fx = f(x);
        if (!std::isfinite(primal_double(fx))) return monte_carlo_detail::failed<T>(static_cast<long>(n) + 1);
        ++n;
        const T delta = T(fx - mean);
        mean = T(mean + delta / T(static_cast<double>(n)));
        m2 = T(m2 + delta * (fx - mean));
    }
    if (n < 2) return monte_carlo_detail::failed<T>(static_cast<long>(n));
    const T variance = T(abs(m2) / T(static_cast<double>(n - 1)));
    const T standard_error = T(volume * sqrt(T(variance / T(static_cast<double>(n)))));
    return {T(volume * mean), standard_error, static_cast<long>(n),
            n >= 30 ? IntegralStatus::Converged : IntegralStatus::Suspicious};
}

template<std::size_t K, Scalar T, typename F>
IntegralResult<T> quasi_monte_carlo(F&& f, const Vec<T, K>& lo, const Vec<T, K>& hi, MonteCarloOptions o = {}) {
    using std::sqrt; using std::abs;
    const bool sobol = o.sequence == Sequence::Sobol;
    if (sobol ? K > monte_carlo_detail::kMaxSobolDimensions : K > monte_carlo_detail::kPrimes.size())
        return monte_carlo_detail::failed<T>(0);                // more dimensions than the sequence has
    const int R = o.shifts < 2 ? 2 : o.shifts;
    const std::size_t per_shift = o.samples / static_cast<std::size_t>(R) > 0 ? o.samples / static_cast<std::size_t>(R) : 1;
    std::mt19937_64 rng(o.seed);
    const T volume = monte_carlo_detail::box_volume<K, T>(lo, hi);

    std::array<std::array<std::uint32_t, 33>, 40> dirs{};
    if (sobol) for (std::size_t i = 0; i < K; ++i) dirs[i] = monte_carlo_detail::sobol_directions(i);

    std::array<double, 40> shift{};
    std::array<std::uint32_t, 40> word{}, state{};
    T mean_of_means{0}, m2{0};
    Vec<T, K> x{};
    long evaluations = 0;
    for (int r = 0; r < R; ++r) {
        for (std::size_t i = 0; i < K; ++i) {
            shift[i] = monte_carlo_detail::uniform01(rng);
            word[i] = static_cast<std::uint32_t>(rng() >> 32);       // the digital shift
            state[i] = 0;
        }
        T sum{0};
        for (std::size_t s_idx = 0; s_idx < per_shift; ++s_idx) {
            if (sobol) {
                if (s_idx > 0) {                                      // Gray-code step: flip by the lowest zero bit of s_idx - 1
                    std::uint64_t value = s_idx - 1;
                    int c = 1;
                    while (value & 1u) { value >>= 1; ++c; }
                    for (std::size_t i = 0; i < K; ++i) state[i] ^= dirs[i][c];
                }
                for (std::size_t i = 0; i < K; ++i) {
                    const double u = static_cast<double>(state[i] ^ word[i]) * (1.0 / 4294967296.0);
                    x[i] = T(lo[i] + (hi[i] - lo[i]) * T(u));
                }
            } else {
                for (std::size_t i = 0; i < K; ++i) {
                    double u = monte_carlo_detail::radical_inverse(s_idx + 1, monte_carlo_detail::kPrimes[i]) + shift[i];
                    u -= std::floor(u);
                    x[i] = T(lo[i] + (hi[i] - lo[i]) * T(u));
                }
            }
            const T fx = f(x);
            ++evaluations;
            if (!std::isfinite(primal_double(fx))) return monte_carlo_detail::failed<T>(evaluations);
            sum = T(sum + fx);
        }
        const T m = T(sum / T(static_cast<double>(per_shift)));
        const T delta = T(m - mean_of_means);
        mean_of_means = T(mean_of_means + delta / T(static_cast<double>(r + 1)));
        m2 = T(m2 + delta * (m - mean_of_means));
    }
    const T variance_of_means = T(abs(m2) / T(static_cast<double>(R - 1)));
    const T standard_error = T(volume * sqrt(T(variance_of_means / T(static_cast<double>(R)))));
    return {T(volume * mean_of_means), standard_error, evaluations,
            per_shift >= 30 ? IntegralStatus::Converged : IntegralStatus::Suspicious};
}

} // namespace spatium
