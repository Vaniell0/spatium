#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/error.hpp>
#  include <spatium/discrete/finite_set.hpp>
#  include <algorithm>
#  include <cstddef>
#  include <limits>
#  include <vector>
#endif

// ── Counting + generation ─────────────────────────────────────────
//
// factorial/binomial_coefficient/permutations_count and the
// k_combinations generator below. Deliberately just counting and
// generation -- no probability distributions or anything statistical
// belongs in this file; that's a separate, larger piece the project
// has explicitly not scoped in yet.
//
// Integral-vs-Scalar genericity: n and k are always plain
// `unsigned long long` -- they are *counts* (which/how-many), not
// values living in whatever numeric type the result is expressed in.
// This mirrors algebra/functions.hpp's power(g, a, int n): the
// exponent there stays a bare `int` even though the base/result is a
// generic Group::ElementType, for the same reason. What *is*
// templated on Scalar is the result type T (default std::uint64_t):
// callers who need exact counts beyond 64 bits (e.g. factorial(30) or
// binomial_coefficient(1000, 500), both far outside uint64_t range)
// can instantiate with T = Real50/Real100 and get the same
// multiplicative algorithm running in arbitrary precision, no
// separate code path. The overflow check below only fires for
// *bounded* integral T (std::numeric_limits<T>::is_bounded &&
// is_integer) -- Real50/Real100 and plain double skip it, since
// neither has the fixed-width wraparound this guards against.

SPATIUM_EXPORT namespace spatium::discrete {

// n! -- product 1*2*...*n (0! = 1! = 1, by convention).
template<Scalar T = std::uint64_t>
Result<T> factorial(unsigned long long n) {
    T result{1};
    for (unsigned long long i = 2; i <= n; ++i) {
        T next = T(i);
        if constexpr (std::numeric_limits<T>::is_bounded && std::numeric_limits<T>::is_integer) {
            if (result > std::numeric_limits<T>::max() / next)
                return std::unexpected(Error{ErrorCode::OutOfDomain,
                    "factorial(n) overflows T"});
        }
        result = T(result * next);
    }
    return result;
}

// nCr = n! / (k! (n-k)!), via the incremental multiplicative formula
// (result_{i+1} = result_i * (n-i) / (i+1), which is exact at every
// step -- no need to compute the three separate, much larger
// factorials and divide). C(n, k) = 0 for k > n, the standard
// combinatorial convention, not an error condition.
//
// Bounded integral T needs a wider scratch accumulator than T itself:
// unlike factorial/permutations_count (plain increasing products, no
// division), this formula's *intermediate* result_i * (n-i) -- before
// the following /(i+1) brings it back down -- can temporarily run
// several times larger than the final coefficient. Checking overflow
// against T's own width at every step (the naive approach) rejects a
// real, non-hypothetical range of valid (n,k) even at the default
// T = std::uint64_t: e.g. C(63,29) = 759510004936100355 fits easily
// under 2^64-1 but trips a per-step uint64 guard -- confirmed by
// brute-force simulation against Python's math.comb (2418 spurious
// rejections for n up to 2000 with the naive check, 0 with the fix
// below checked up to n=250 including every true-overflow boundary
// in that range). Accumulating in unsigned __int128 (available on
// both compilers this project targets, GCC >= 15 and Clang >= 19)
// and only bounds-checking against T's own max at the very end fixes
// this for any T up to 64 bits; T itself may still be narrower (the
// overflow check on the final cast handles that) or may be an
// unbounded type (Real50/Real100), which skips this path entirely.
template<Scalar T = std::uint64_t>
Result<T> binomial_coefficient(unsigned long long n, unsigned long long k) {
    if (k > n) return T{0};
    k = std::min(k, n - k); // C(n,k) == C(n,n-k); fewer multiplications

    if constexpr (std::numeric_limits<T>::is_bounded && std::numeric_limits<T>::is_integer) {
#if defined(__SIZEOF_INT128__)
        using Wide = unsigned __int128;
#else
        using Wide = unsigned long long;
#endif
        Wide result = 1;
        for (unsigned long long i = 0; i < k; ++i) {
            Wide num = Wide(n - i);
            if (num != 0 && result > std::numeric_limits<Wide>::max() / num)
                return std::unexpected(Error{ErrorCode::OutOfDomain,
                    "binomial_coefficient(n,k) overflows the working precision"});
            result = result * num / Wide(i + 1);
        }
        if (result > Wide(std::numeric_limits<T>::max()))
            return std::unexpected(Error{ErrorCode::OutOfDomain,
                "binomial_coefficient(n,k) overflows T"});
        return T(result);
    } else {
        T result{1};
        for (unsigned long long i = 0; i < k; ++i)
            result = T(result * T(n - i) / T(i + 1));
        return result;
    }
}

// nPr = n! / (n-k)! -- the number of ways to arrange k items chosen
// from n distinguishable items, order mattering. Same 0-for-k>n
// convention as binomial_coefficient.
template<Scalar T = std::uint64_t>
Result<T> permutations_count(unsigned long long n, unsigned long long k) {
    if (k > n) return T{0};
    T result{1};
    for (unsigned long long i = 0; i < k; ++i) {
        T next = T(n - i);
        if constexpr (std::numeric_limits<T>::is_bounded && std::numeric_limits<T>::is_integer) {
            if (result != T{0} && result > std::numeric_limits<T>::max() / next)
                return std::unexpected(Error{ErrorCode::OutOfDomain,
                    "permutations_count(n,k) overflows T"});
        }
        result = T(result * next);
    }
    return result;
}

// All size-k subsets of `set`, as FiniteSet<T> values, generated in
// lexicographic order of the underlying (sorted, per FiniteSet's own
// invariant) element indices. Exactly binomial_coefficient(set.size(), k)
// results -- checked as a regression in tests/test_combinatorics.cpp,
// not just asserted here. k > set.size() yields an empty vector (no
// combinations exist), matching binomial_coefficient's own
// k > n == 0 convention.
template<typename T>
std::vector<FiniteSet<T>> k_combinations(const FiniteSet<T>& set, std::size_t k) {
    std::vector<FiniteSet<T>> result;
    auto n = set.size();
    if (k > n) return result;

    if (k == 0) {
        result.emplace_back(); // the empty set is the sole 0-combination
        return result;
    }

    std::vector<std::size_t> indices(k);
    for (std::size_t i = 0; i < k; ++i) indices[i] = i;

    while (true) {
        std::vector<T> subset;
        subset.reserve(k);
        for (auto idx : indices) subset.push_back(set.elements[idx]);
        result.emplace_back(std::move(subset));

        // Advance to the next combination: find the rightmost index
        // not already at its maximum position, bump it, then reset
        // everything to its right to run consecutively after it.
        auto i = static_cast<std::ptrdiff_t>(k) - 1;
        while (i >= 0 && indices[static_cast<std::size_t>(i)] == n - k + static_cast<std::size_t>(i))
            --i;
        if (i < 0) break; // last combination was the lexicographically final one
        ++indices[static_cast<std::size_t>(i)];
        for (auto j = static_cast<std::size_t>(i) + 1; j < k; ++j)
            indices[j] = indices[j - 1] + 1;
    }
    return result;
}

} // namespace spatium::discrete
