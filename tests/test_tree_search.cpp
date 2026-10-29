// The tree search of rsc/core, guided by a learned prior and value, and the games it plays
// against itself.
//
// The method is tried on a game that has an exact answer: Nim with a pile and takes of 1 to 3,
// the player who takes the last stone winning. A pile is lost for the side to move exactly when it
// is a multiple of four, so the winning move from any other pile is the one that leaves one, and
// "the search plays well" is a count of piles where it takes that. Each piece is held to a second
// route: the search to the exact answer, and what the games teach to the same exact answer on
// piles it never started from. The learned parts are the library's (a linear softmax for the
// prior, a tanh net for the value), fitted by statistics::fit_pass with AdaGrad and Adam. Each
// check gets a case it must fail: a search that tries nothing, a game whose turn is two moves
// searched as if the sides alternated, and nets that have not played a game.

#include <rsc/core/tree_search.hpp>
#include <spatium/algebra/analysis/optimize.hpp>
#include <spatium/algebra/models/network.hpp>
#include <spatium/algebra/models/softmax.hpp>
#include <spatium/statistics/fit.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <random>
#include <vector>

using namespace rsc;
using Catch::Matchers::WithinAbs;

namespace {

constexpr std::size_t row_size = 5, state_size = 5, hidden = 8;

// The learned parts the search asks for: a softmax over the move rows for the prior, a tanh net
// over the state's row for the value.
struct NimNets {
    spatium::LinearSoftmax<row_size> policy;
    spatium::TanhNet<state_size, hidden> net;
    void prior(std::span<const std::array<double, row_size>> rows, std::span<double> out) const {
        policy.distribution(rows, out);
    }
    double value(const std::array<double, state_size>& f) const { return net.value(f); }
};

// One pass over a batch: the policy toward the visits the search made, by AdaGrad, and the net
// toward how each game ended, by Adam.
struct Trainer {
    spatium::AdaGrad<row_size> policy_rule;
    spatium::Adam<spatium::TanhNet<state_size, hidden>::weight_count> net_rule;
    template <class A> void fit(NimNets& nets, std::span<const Sample<A>> batch, std::mt19937_64& rng) {
        using spatium::statistics::Fit;
        spatium::statistics::fit_pass(nets.policy.weights(), policy_rule, batch, 0.5, rng, [&](const Sample<A>& s) {
            const std::span<const std::array<double, row_size>> rows(s.rows);
            const std::span<const double> target(s.target);
            return Fit<std::array<double, row_size>>{nets.policy.cross_entropy(rows, target), nets.policy.gradient(rows, target)};
        });
        spatium::statistics::fit_pass(nets.net.weights(), net_rule, batch, 0.01, rng, [&](const Sample<A>& s) {
            return Fit<spatium::TanhNet<state_size, hidden>::Weights>{nets.net.loss(s.features, s.outcome),
                                                                       nets.net.gradient(s.features, s.outcome)};
        });
    }
};

struct Nim {
    struct State {
        int pile;
    };
    using Move = int;
    using Row = std::array<double, row_size>;
    using Features = std::array<double, state_size>;

    static std::vector<Move> moves(const State& s) {
        std::vector<Move> out;
        for (int take = 1; take <= 3 && take <= s.pile; ++take)
            out.push_back(take);
        return out;
    }
    static State after(const State& s, Move take) { return {s.pile - take}; }
    // With no stones the player before took the last one: the side to move has lost.
    static std::optional<double> outcome(const State&, bool no_moves) {
        return no_moves ? std::optional<double>(-1.0) : std::nullopt;
    }
    // The pile as a count, and its remainder mod 4 as a one-hot: a way of writing the pile down.
    static Features features(const State& s) {
        Features f{};
        f[static_cast<std::size_t>(s.pile % 4)] = 1;
        f[4] = s.pile / 30.0;
        return f;
    }
    static std::vector<Row> rows(const State& s, std::span<const Move> moves) {
        std::vector<Row> out;
        for (Move take : moves) {
            Row r{};
            r[static_cast<std::size_t>((s.pile - take) % 4)] = 1;
            r[4] = take / 3.0;
            out.push_back(r);
        }
        return out;
    }
};
static_assert(SearchGame<Nim>);
static_assert(Evaluator<NimNets, Nim>);

// In how many of the piles 5..largest that are not multiples of four the search takes the stones
// that leave a multiple of four.
double winning_share(const NimNets& nets, int simulations, std::uint64_t seed, int largest = 30) {
    std::mt19937_64 rng(seed);
    int right = 0, total = 0;
    SearchOptions options;
    options.simulations = simulations;
    for (int pile = 5; pile <= largest; ++pile) {
        if (pile % 4 == 0)
            continue;
        const auto result = search<Nim>({pile}, nets, options, rng);
        if (result.moves.empty())
            continue;
        const int take = result.moves[result.pick(0, rng)];
        right += (pile - take) % 4 == 0;
        ++total;
    }
    return static_cast<double>(right) / total;
}

// The same game with a turn of two moves: first the stones are chosen (the side keeps the move),
// then the choice is confirmed (the side passes). Played right it is Nim, so what the search finds
// must be Nim's answer; played with the sides assumed to alternate it is a different game.
struct SplitNim {
    struct State {
        int pile;
        int pending; // stones chosen and not yet taken
    };
    using Move = int; // 1..3 to choose, 0 to confirm
    using Row = std::array<double, row_size>;
    using Features = std::array<double, state_size>;

    static std::vector<Move> moves(const State& s) {
        if (s.pending > 0)
            return {0};
        std::vector<Move> out;
        for (int take = 1; take <= 3 && take <= s.pile; ++take)
            out.push_back(take);
        return out;
    }
    static State after(const State& s, Move m) {
        return m == 0 ? State{s.pile - s.pending, 0} : State{s.pile, m};
    }
    static bool keeps_side(const State& s, Move) { return s.pending == 0; }
    static std::optional<double> outcome(const State&, bool no_moves) {
        return no_moves ? std::optional<double>(-1.0) : std::nullopt;
    }
    static Features features(const State& s) {
        Features f{};
        f[static_cast<std::size_t>((s.pile - s.pending) % 4)] = 1;
        f[4] = s.pending / 3.0;
        return f;
    }
    static std::vector<Row> rows(const State& s, std::span<const Move> moves) {
        std::vector<Row> out;
        for (Move m : moves) {
            Row r{};
            r[static_cast<std::size_t>((s.pile - m) % 4)] = m > 0 ? 1 : 0;
            r[4] = m / 3.0;
            out.push_back(r);
        }
        return out;
    }
};
static_assert(SearchGame<SplitNim>);

} // namespace

TEST_CASE("search: with no learning at all, enough simulations find the winning take", "[rsc][search]") {
    NimNets blank;
    // (piles to 12: terminal results reach the root from a few moves deep)
    REQUIRE(winning_share(blank, 600, 5, 12) > 0.95);
    // the same search with no simulations tries nothing: it can only take the first of the moves
    REQUIRE(winning_share(blank, 1, 5, 12) < 0.6);
}

TEST_CASE("search: the standing of each move is what the simulations under it found", "[rsc][search]") {
    NimNets blank;
    std::mt19937_64 rng(3);
    SearchOptions options;
    options.simulations = 600;
    // Five stones: taking one leaves four, which loses for the player who must move; taking two
    // or three leaves a pile that wins for them.
    const auto r = search<Nim>({5}, blank, options, rng);
    REQUIRE(r.values.size() == 3);
    CHECK(r.values[0] > 0.8);
    CHECK(r.values[1] < 0);
    CHECK(r.values[2] < 0);
    // the mean over the visits is the standing of the root
    double sum = 0;
    int visits = 0;
    for (std::size_t i = 0; i < r.values.size(); ++i)
        sum += r.values[i] * r.visits[i], visits += r.visits[i];
    CHECK_THAT(r.searched, WithinAbs(sum / visits, 1e-12));
    CHECK(r.searched > 0.5);
    // one simulation away from nothing, the root stands where the net puts it
    options.simulations = 0;
    const auto none = search<Nim>({5}, blank, options, rng);
    CHECK(none.searched == none.value);
    // a finished game is exact
    CHECK(search<Nim>({0}, blank, options, rng).searched == -1.0);
}

TEST_CASE("search: a result forced by every reply is proved, and its line is the most visited", "[rsc][search]") {
    NimNets blank;
    std::mt19937_64 rng(5);
    SearchOptions options;
    options.simulations = 600;
    options.prove = true;
    // Five stones: take one, whatever follows take what is left: a win in three plies.
    const auto win = search<Nim>({5}, blank, options, rng);
    CHECK(win.proof == 3);
    REQUIRE(!win.line.empty());
    CHECK(win.line.front() == 1);
    CHECK(win.line.size() <= 12);
    // Four stones: whatever is taken, the other takes the rest: a loss in two plies.
    CHECK(search<Nim>({4}, blank, options, rng).proof == -2);
    // One stone is a win in one ply.
    CHECK(search<Nim>({1}, blank, options, rng).proof == 1);
    // Nothing is claimed without the option, nor with too few simulations to try every reply.
    options.prove = false;
    CHECK(search<Nim>({5}, blank, options, rng).proof == 0);
    options.prove = true;
    options.simulations = 2;
    CHECK(search<Nim>({12}, blank, options, rng).proof == 0);
}

TEST_CASE("search: playing itself, it learns to win Nim from a blank start", "[rsc][search]") {
    NimNets nets;
    nets.net = spatium::TanhNet<state_size, hidden>::random(3);
    Trainer trainer;
    // a few simulations are not enough for the blank nets to be right
    const double before = winning_share(nets, 6, 17);

    std::mt19937_64 rng(7);
    SearchOptions options;
    options.simulations = 24;
    options.noise = 0.25;
    for (int iteration = 0; iteration < 40; ++iteration) {
        std::vector<Sample<Nim>> batch;
        for (int g = 0; g < 12; ++g) {
            const int pile = 8 + static_cast<int>(rng() % 20);
            auto game = self_play<Nim>({pile}, nets, options, 3, rng);
            for (auto& s : game.samples)
                batch.push_back(std::move(s));
        }
        for (int pass = 0; pass < 2; ++pass)
            trainer.fit<Nim>(nets, std::span<const Sample<Nim>>(batch), rng);
    }
    const double after = winning_share(nets, 6, 17);
    REQUIRE(after > 0.9);
    REQUIRE(after > before + 0.2);
}

TEST_CASE("search: the policy's prior decides where the first simulations go", "[rsc][search]") {
    NimNets nets;
    // rows end with take / 3: a weight there makes taking three the policy's favourite
    nets.policy.weights() = {0, 0, 0, 0, 6};
    SearchOptions options;
    options.simulations = 3;
    std::mt19937_64 rng(1);
    const auto result = search<Nim>({20}, nets, options, rng);
    REQUIRE(result.moves == std::vector<int>{1, 2, 3});
    REQUIRE(result.visits[2] >= 2);
    REQUIRE(result.visits[2] > result.visits[0]);
    // with the opposite weight it is taking one
    nets.policy.weights() = {0, 0, 0, 0, -6};
    const auto other = search<Nim>({20}, nets, options, rng);
    REQUIRE(other.visits[0] >= 2);
}

TEST_CASE("search: the games differ by their first moves only while they explore", "[rsc][search]") {
    NimNets nets;
    SearchOptions options;
    options.simulations = 12;
    auto distinct_second_states = [&](int explore_plies) {
        std::mt19937_64 rng(4);
        std::vector<double> seen;
        for (int g = 0; g < 30; ++g) {
            const auto game = self_play<Nim>({20}, nets, options, explore_plies, rng);
            const double f = game.samples[1].features[4];
            if (std::ranges::find(seen, f) == seen.end())
                seen.push_back(f);
        }
        return seen.size();
    };
    REQUIRE(distinct_second_states(0) == 1);
    REQUIRE(distinct_second_states(4) >= 2);
}

TEST_CASE("search: the policy is fitted to what the search found, not to uniform", "[rsc][search]") {
    // the probability the bare policy (no search) gives to the winning take, over piles 5..30
    auto policy_share = [](const NimNets& nets) {
        double sum = 0;
        int count = 0;
        for (int pile = 5; pile <= 30; ++pile) {
            if (pile % 4 == 0)
                continue;
            const auto moves = Nim::moves({pile});
            const auto rows = Nim::rows({pile}, std::span<const int>(moves));
            std::vector<double> p(rows.size());
            nets.policy.distribution(std::span<const Nim::Row>(rows), std::span<double>(p));
            for (std::size_t i = 0; i < moves.size(); ++i)
                if ((pile - moves[i]) % 4 == 0)
                    sum += p[i];
            ++count;
        }
        return sum / count;
    };
    NimNets nets;
    nets.net = spatium::TanhNet<state_size, hidden>::random(1);
    Trainer trainer;
    const double before = policy_share(nets);
    std::mt19937_64 rng(7);
    SearchOptions options;
    options.simulations = 24;
    options.noise = 0.25;
    for (int iteration = 0; iteration < 40; ++iteration) {
        std::vector<Sample<Nim>> batch;
        for (int g = 0; g < 12; ++g) {
            const int pile = 8 + static_cast<int>(rng() % 20);
            for (auto& s : self_play<Nim>({pile}, nets, options, 3, rng).samples)
                batch.push_back(std::move(s));
        }
        trainer.fit<Nim>(nets, std::span<const Sample<Nim>>(batch), rng);
    }
    REQUIRE(before < 0.45);
    REQUIRE(policy_share(nets) > 0.8);
}

TEST_CASE("search: a turn of two moves keeps its sign", "[rsc][search]") {
    // A choice keeps the side, a confirmation passes it; a game with no hook always passes it.
    CHECK(!changes_side<SplitNim>({9, 0}, 2));
    CHECK(changes_side<SplitNim>({9, 2}, 0));
    CHECK(changes_side<Nim>({9}, 2));
    NimNets blank;
    std::mt19937_64 rng(31);
    SearchOptions options;
    options.simulations = 1500;
    int right = 0, total = 0;
    for (int pile = 5; pile <= 11; ++pile) {
        if (pile % 4 == 0)
            continue;
        const SearchResult<SplitNim> result = search<SplitNim>({pile, 0}, blank, options, rng);
        const int chosen = result.moves[result.pick(0, rng)];
        right += (pile - chosen) % 4 == 0;
        ++total;
    }
    // The sides assumed to alternate, the search values the confirmation as the other side's and
    // so prefers the wrong choice in about half the piles; with the hook it is Nim.
    CHECK(right == total);
}

TEST_CASE("search: the outcome of a sample follows the sides through a turn", "[rsc][search]") {
    NimNets nets;
    std::mt19937_64 rng(8);
    SearchOptions options;
    options.simulations = 40;
    const SelfPlayGame<SplitNim> game = self_play<SplitNim>({9, 0}, nets, options, 100, rng);
    REQUIRE(game.samples.size() >= 4);
    const auto& s = game.samples;
    // The last move took the last stone: its side won.
    CHECK(s.back().outcome == 1.0);
    for (std::size_t i = 0; i + 1 < s.size(); ++i) {
        const bool next_is_confirmation = s[i + 1].features[4] > 0;
        CHECK(s[i].kept_side == next_is_confirmation);
        CHECK(s[i].outcome == (s[i].kept_side ? s[i + 1].outcome : -s[i + 1].outcome));
    }
}
