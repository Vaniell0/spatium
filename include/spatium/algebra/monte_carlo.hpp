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
//   quasi_monte_carlo(f, lo, hi)   a Halton sequence, randomly shifted; error ~ (log n)^K / n
//                                  for a smooth f, and the estimate comes with its own error
//
// Both return the `IntegralResult` of the rest of the library, with the difference that its
// `error_estimate` is a standard error -- one sigma of the estimate, not a bound: the true
// error is within it about two times in three. The status is `Failed` if f is ever not finite,
// `Suspicious` where there are too few points for a standard error to mean anything, and
// otherwise `Converged`, which for a sampling estimator says only that an estimate with an
// error was formed. `evaluations` is the number of points.
//
// For the quasi-random one the points are the Halton points of the first K primes, shifted
// by an independent uniform vector (mod 1) in each of `shifts` copies; the standard error is
// the spread of the copies' means, which is what makes an error estimate possible for a
// deterministic low-discrepancy sequence at all (Cranley-Patterson). Up to 32 dimensions.
//
// Reproducible: the generator is a seeded `std::mt19937_64` and uniform numbers are formed
// from its bits, not from `std::uniform_real_distribution`, whose output differs between
// standard libraries.

SPATIUM_EXPORT namespace spatium {

struct MonteCarloOptions {
    std::size_t samples = std::size_t{1} << 16;
    std::uint64_t seed = 1;
    int shifts = 16;            // quasi_monte_carlo: independent random shifts
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
    static_assert(K >= 1 && K <= 32, "the Halton sequence here has 32 bases");
    const int R = o.shifts < 2 ? 2 : o.shifts;
    const std::size_t per_shift = o.samples / static_cast<std::size_t>(R) > 0 ? o.samples / static_cast<std::size_t>(R) : 1;
    std::mt19937_64 rng(o.seed);
    const T volume = monte_carlo_detail::box_volume<K, T>(lo, hi);

    std::array<double, 32> shift{};
    T mean_of_means{0}, m2{0};
    Vec<T, K> x{};
    long evaluations = 0;
    for (int r = 0; r < R; ++r) {
        for (std::size_t i = 0; i < K; ++i) shift[i] = monte_carlo_detail::uniform01(rng);
        T sum{0};
        for (std::size_t s = 1; s <= per_shift; ++s) {
            for (std::size_t i = 0; i < K; ++i) {
                double u = monte_carlo_detail::radical_inverse(s, monte_carlo_detail::kPrimes[i]) + shift[i];
                u -= std::floor(u);
                x[i] = T(lo[i] + (hi[i] - lo[i]) * T(u));
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
