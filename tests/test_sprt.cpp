// The sequential test of who is stronger (statistics/sprt.hpp).
//
// The ratio and its bounds are held to their formulas; then games are drawn: between equal players
// the test takes the first as stronger about as rarely as alpha says, and a player stronger by the
// test's own margin it takes as stronger about as often as 1 - beta says. The case it must fail: a
// fixed 24 games scored at 0.55, the gate it replaces, takes an equal player as stronger once in
// four.

#include <spatium/statistics/sprt.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <random>

using spatium::statistics::Sprt;
using Catch::Matchers::WithinAbs;

namespace {

// How often the test takes the first player as stronger when a game is a win with chance `p`, and
// the mean number of games it plays, over `trials` runs stopped at `cap` games.
struct Run {
    double taken = 0, games = 0;
};
Run run(const Sprt& test, double p, int trials, int cap, std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::bernoulli_distribution win(p);
    Run r;
    for (int t = 0; t < trials; ++t) {
        int w = 0, l = 0, v = 0;
        while (v == 0 && w + l < cap) {
            (win(rng) ? w : l) += 1;
            v = test.verdict(w, l);
        }
        r.taken += v > 0;
        r.games += w + l;
    }
    r.taken /= trials, r.games /= trials;
    return r;
}

} // namespace

TEST_CASE("sprt: the ratio and the bounds are Wald's", "[sprt]") {
    const Sprt test{.elo = 35};
    const double p1 = Sprt::win_chance(35);
    CHECK_THAT(Sprt::win_chance(0), WithinAbs(0.5, 1e-15));
    CHECK_THAT(test.llr(7, 3), WithinAbs(7 * std::log(p1 / 0.5) + 3 * std::log((1 - p1) / 0.5), 1e-12));
    CHECK_THAT(test.upper(), WithinAbs(std::log(19.0), 1e-12));
    CHECK_THAT(test.lower(), WithinAbs(-std::log(19.0), 1e-12));
    CHECK(test.verdict(0, 0) == 0);
    CHECK(test.verdict(200, 100) == 1);
    CHECK(test.verdict(100, 200) == -1);
}

TEST_CASE("sprt: an equal player is rarely taken as stronger, a stronger one usually is", "[sprt]") {
    const Sprt test{.elo = 70};
    const Run equal = run(test, 0.5, 4000, 100000, 1);
    const Run stronger = run(test, Sprt::win_chance(70), 4000, 100000, 2);
    CHECK(equal.taken < 0.065);
    CHECK(stronger.taken > 0.935);
    // and it needs games in the hundreds at this margin, not thousands
    CHECK(equal.games < 300);
    CHECK(stronger.games < 300);

    // the case it must fail: the fixed gate takes an equal player as stronger once in four
    std::mt19937_64 rng(3);
    std::bernoulli_distribution win(0.5);
    int passed = 0;
    for (int t = 0; t < 4000; ++t) {
        int w = 0;
        for (int g = 0; g < 24; ++g)
            w += win(rng);
        passed += w / 24.0 >= 0.55;
    }
    CHECK(passed / 4000.0 > 0.24);
    CHECK(passed / 4000.0 > 4 * equal.taken);
}
