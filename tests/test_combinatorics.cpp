#include <catch2/catch_test_macros.hpp>
#include <spatium/core/precision.hpp>
#include <spatium/discrete/combinatorics.hpp>
#include <algorithm>
#include <cstdint>
#include <set>

using namespace spatium;
using namespace spatium::discrete;

// ── factorial ─────────────────────────────────────────────────────

TEST_CASE("factorial: known small values", "[combinatorics]") {
    CHECK(*factorial(0ULL) == 1ULL);
    CHECK(*factorial(1ULL) == 1ULL);
    CHECK(*factorial(5ULL) == 120ULL);
    CHECK(*factorial(10ULL) == 3628800ULL);
    // 20! = 2432902008176640000, the largest factorial that still fits
    // in an unsigned 64-bit result (21! overflows).
    CHECK(*factorial(20ULL) == 2432902008176640000ULL);
}

TEST_CASE("factorial: overflow is reported, not silently wrapped", "[combinatorics]") {
    // 8! = 40320 fits in uint16_t (max 65535); 9! = 362880 does not.
    auto ok = factorial<std::uint16_t>(8);
    REQUIRE(ok.has_value());
    CHECK(*ok == 40320);

    auto overflowed = factorial<std::uint16_t>(9);
    CHECK_FALSE(overflowed.has_value());
}

// ── binomial_coefficient ──────────────────────────────────────────

TEST_CASE("binomial_coefficient: known small values", "[combinatorics]") {
    CHECK(*binomial_coefficient(5ULL, 2ULL) == 10ULL);
    CHECK(*binomial_coefficient(10ULL, 5ULL) == 252ULL);
    CHECK(*binomial_coefficient(6ULL, 0ULL) == 1ULL);
    CHECK(*binomial_coefficient(6ULL, 6ULL) == 1ULL);
}

TEST_CASE("binomial_coefficient: k > n is zero by convention, not an error", "[combinatorics]") {
    auto r = binomial_coefficient(6ULL, 7ULL);
    REQUIRE(r.has_value());
    CHECK(*r == 0ULL);
}

TEST_CASE("binomial_coefficient: symmetry C(n,k) == C(n,n-k)", "[combinatorics]") {
    CHECK(*binomial_coefficient(12ULL, 5ULL) == *binomial_coefficient(12ULL, 7ULL));
}

TEST_CASE("binomial_coefficient: overflow is reported for bounded integral T", "[combinatorics]") {
    // C(20,10) = 184756, which does not fit in uint16_t (max 65535).
    auto overflowed = binomial_coefficient<std::uint16_t>(20, 10);
    CHECK_FALSE(overflowed.has_value());
}

TEST_CASE("binomial_coefficient: does not spuriously overflow uint64_t on values that fit",
          "[combinatorics]") {
    // Regression for a real bug found while writing these tests: a
    // naive per-step overflow check on the incremental multiplicative
    // formula (result = result * (n-i) / (i+1)) rejects this input
    // even though the true value fits easily under 2^64-1, because the
    // *intermediate* product before the /(i+1) division briefly runs
    // several times larger than the final coefficient. Verified
    // independently against Python's math.comb(63, 29).
    auto r = binomial_coefficient<std::uint64_t>(63, 29);
    REQUIRE(r.has_value());
    CHECK(*r == 759510004936100355ULL);
}

TEST_CASE("binomial_coefficient: Real50 reaches counts uint64_t cannot hold", "[combinatorics]") {
    // C(100,50) = 100891344545564193334812497256, far past uint64_t's
    // ~1.8e19 ceiling -- computed independently via Python's
    // math.comb(100, 50) as ground truth, not derived from this file.
    auto overflowed = binomial_coefficient<std::uint64_t>(100, 50);
    CHECK_FALSE(overflowed.has_value());

    auto big = binomial_coefficient<Real50>(100, 50);
    REQUIRE(big.has_value());
    Real50 expected("100891344545564193334812497256");
    Real50 diff = *big - expected;
    if (diff < Real50(0)) diff = -diff;
    CHECK(diff < Real50("1e10")); // well within Real50's working precision
}

// ── permutations_count ────────────────────────────────────────────

TEST_CASE("permutations_count: known small values", "[combinatorics]") {
    CHECK(*permutations_count(5ULL, 2ULL) == 20ULL);
    CHECK(*permutations_count(5ULL, 5ULL) == *factorial(5ULL));
    CHECK(*permutations_count(5ULL, 0ULL) == 1ULL);
}

TEST_CASE("permutations_count: k > n is zero by convention", "[combinatorics]") {
    auto r = permutations_count(5ULL, 6ULL);
    REQUIRE(r.has_value());
    CHECK(*r == 0ULL);
}

// ── k_combinations ────────────────────────────────────────────────

TEST_CASE("k_combinations: count matches binomial_coefficient exactly", "[combinatorics]") {
    FiniteSet<int> set{1, 2, 3, 4, 5, 6};
    for (std::size_t k = 0; k <= set.size(); ++k) {
        auto combos = k_combinations(set, k);
        auto expected = *binomial_coefficient(set.size(), k);
        CHECK(combos.size() == expected);
    }
}

TEST_CASE("k_combinations: no duplicate combinations are produced", "[combinatorics]") {
    FiniteSet<int> set{1, 2, 3, 4, 5};
    auto combos = k_combinations(set, 3);

    std::set<std::vector<int>> unique;
    for (auto& c : combos) unique.insert(c.elements);
    CHECK(unique.size() == combos.size());
}

TEST_CASE("k_combinations: every result is a genuine subset with the right size", "[combinatorics]") {
    FiniteSet<int> set{1, 2, 3, 4, 5};
    auto combos = k_combinations(set, 3);
    for (auto& c : combos) {
        CHECK(c.size() == 3);
        for (auto& elem : c.elements)
            CHECK(set.contains(elem));
    }
}

TEST_CASE("k_combinations: k == 0 yields exactly the empty set", "[combinatorics]") {
    FiniteSet<int> set{1, 2, 3};
    auto combos = k_combinations(set, 0);
    REQUIRE(combos.size() == 1);
    CHECK(combos[0].empty());
}

TEST_CASE("k_combinations: k > size yields no combinations", "[combinatorics]") {
    FiniteSet<int> set{1, 2, 3};
    auto combos = k_combinations(set, 4);
    CHECK(combos.empty());
}

TEST_CASE("k_combinations: small hand-checkable example", "[combinatorics]") {
    FiniteSet<int> set{1, 2, 3};
    auto combos = k_combinations(set, 2);
    REQUIRE(combos.size() == 3);

    std::vector<std::vector<int>> expected{{1, 2}, {1, 3}, {2, 3}};
    std::vector<std::vector<int>> actual;
    for (auto& c : combos) actual.push_back(c.elements);
    std::sort(actual.begin(), actual.end());
    CHECK(actual == expected);
}
