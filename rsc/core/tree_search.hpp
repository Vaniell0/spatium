#pragma once

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <numeric>
#include <optional>
#include <random>
#include <span>
#include <utility>
#include <vector>

// A tree search guided by a learned prior and value, and the games it plays against itself.
//
// The prior says which moves look worth trying; the value says how a state stands for the side
// to move. A search uses both, and what it finds is better than either: the number of times it
// tried each move is a sharper choice than the prior's, and the result of the game it ends in is
// a truer value than the guess. Playing itself, the search makes those two targets; fitting the
// learned parts to them and searching the next games with what was fitted is the whole method.
//
// Nothing here knows a game or a model. A game is an adapter (`SearchGame`): its states and
// moves, the moves of a state, the state after a move, how a finished or cut game stands for the
// side to move, a row of numbers for each move and a row for the state. The sides alternate,
// unless the adapter says that a move keeps the side (`keeps_side(state, move)`): a game whose
// turn is several moves passes the move on only with the last of them. The learned parts are an
// `Evaluator`; how they are fitted to the samples `self_play` returns is the caller's.
//
// Like the rest of rsc/core, this includes nothing of the library it serves.

namespace rsc {

// ── The game ──────────────────────────────────────────────────────

template <class A>
concept SearchGame = requires(const typename A::State& s, const typename A::Move& m) {
    typename A::State;
    typename A::Move;
    typename A::Row;
    typename A::Features;
    { A::moves(s) } -> std::same_as<std::vector<typename A::Move>>;
    { A::after(s, m) } -> std::same_as<typename A::State>;
    // How the game stands for the side to move when it has no moves (`no_moves`) or has been
    // cut: +1 won, 0 drawn, -1 lost; nothing while it goes on.
    { A::outcome(s, true) } -> std::same_as<std::optional<double>>;
    {
        A::rows(s, std::span<const typename A::Move>{})
    } -> std::same_as<std::vector<typename A::Row>>;
    { A::features(s) } -> std::same_as<typename A::Features>;
};

// What the search asks of the learned parts: a prior over a state's moves, written into the
// caller's buffer (one entry per row, summing to one), and how the state stands for the side to
// move, in [-1, 1].
template <class E, class A>
concept Evaluator = SearchGame<A>
    && requires(const E& e, std::span<const typename A::Row> rows, std::span<double> out,
                const typename A::Features& f) {
        e.prior(rows, out);
        { e.value(f) } -> std::convertible_to<double>;
    };

// ── The search ────────────────────────────────────────────────────

// Whether a move passes the turn to the other side; true unless the adapter says it keeps it.
template <SearchGame A>
bool changes_side(const typename A::State& s, const typename A::Move& m) {
    if constexpr (requires {
                      { A::keeps_side(s, m) } -> std::convertible_to<bool>;
                  })
        return !A::keeps_side(s, m);
    else
        return true;
}

struct SearchOptions {
    int simulations = 100;
    double c_puct = 1.5; // how much the policy's prior counts against what was found
    double noise = 0;    // the share of Dirichlet noise mixed into the root's priors
    double noise_alpha = 0.3;
    // Carry forced results up the tree (`SearchResult::proof`). It costs a pass over the edges of
    // each node on the path of every simulation, so training leaves it off.
    bool prove = false;
};

template <SearchGame A> struct SearchResult {
    std::vector<typename A::Move> moves;
    std::vector<typename A::Row> rows;
    std::vector<int> visits;
    double value = 0; // the root's standing for the side to move, as the net alone gives it
    // What the simulations found: the mean standing under each move for the side to move at the
    // root (0 for a move never tried), and the mean over all of them. Equal to `value` when no
    // simulation was run.
    std::vector<double> values;
    double searched = 0;
    // The most visited line from the root, the best move first, as far as the simulations went
    // (at most twelve moves).
    std::vector<typename A::Move> line;
    // With `SearchOptions::prove`: a result the simulations proved by playing out every reply.
    // Above zero the side to move at the root wins by force in that many moves (plies), below zero
    // it loses in that many, zero when nothing is proved.
    int proof = 0;

    // The move chosen: the most visited at temperature 0, else drawn with probability
    // proportional to visits^(1/temperature).
    std::size_t pick(double temperature, std::mt19937_64& rng) const {
        if (temperature <= 0)
            return static_cast<std::size_t>(std::ranges::max_element(visits) - visits.begin());
        std::vector<double> w(visits.size());
        for (std::size_t i = 0; i < w.size(); ++i)
            w[i] = std::pow(static_cast<double>(visits[i]), 1 / temperature);
        return std::discrete_distribution<std::size_t>(w.begin(), w.end())(rng);
    }

    // The visits as a distribution over the moves.
    std::vector<double> distribution() const {
        double total = 0;
        for (int v : visits)
            total += v;
        std::vector<double> p(visits.size());
        for (std::size_t i = 0; i < p.size(); ++i)
            p[i] = total > 0 ? visits[i] / total : 1.0 / static_cast<double>(p.size());
        return p;
    }
};

template <SearchGame A, Evaluator<A> E>
SearchResult<A> search(
    const typename A::State& root, const E& evaluator, const SearchOptions& options,
    std::mt19937_64& rng
) {
    using Move = typename A::Move;
    struct Edge {
        double prior = 0, total = 0;
        int visits = 0, child = -1;
        bool flips = true; // the child's side to move is the other one
    };
    struct Node {
        typename A::State state;
        std::vector<Move> moves;
        std::vector<Edge> edges;
        std::optional<double> terminal;
        bool expanded = false;
        int visits = 0;
        int win = -1, lose = -1; // plies to a forced win or loss for the side to move; -1: none
    };

    std::vector<Node> nodes;
    nodes.reserve(static_cast<std::size_t>(options.simulations) + 1);
    nodes.push_back({root, {}, {}, std::nullopt, false, 0, -1, -1});
    SearchResult<A> result;
    std::vector<double> prior; // the evaluator's prior for the node being expanded

    // Fills in a node's moves and priors; returns its standing for the side to move.
    auto expand = [&](int id, bool is_root) {
        Node& node = nodes[static_cast<std::size_t>(id)];
        node.expanded = true;
        node.moves = A::moves(node.state);
        if (const auto o = A::outcome(node.state, node.moves.empty())) {
            node.terminal = *o;
            if (*o <= -1.0)
                node.lose = 0;
            else if (*o >= 1.0)
                node.win = 0;
            return *o;
        }
        const auto rows = A::rows(node.state, std::span<const Move>(node.moves));
        prior.resize(rows.size());
        evaluator.prior(std::span<const typename A::Row>(rows), std::span<double>(prior));
        node.edges.resize(node.moves.size());
        for (std::size_t i = 0; i < node.edges.size(); ++i)
            node.edges[i].prior = prior[i];
        if (is_root) {
            result.moves = node.moves;
            result.rows = rows;
            if (options.noise > 0) {
                std::gamma_distribution<double> g(options.noise_alpha, 1);
                std::vector<double> dirichlet(node.edges.size());
                double total = 0;
                for (double& d : dirichlet)
                    total += d = g(rng);
                for (std::size_t i = 0; i < node.edges.size(); ++i)
                    node.edges[i].prior = (1 - options.noise) * node.edges[i].prior +
                                          options.noise * dirichlet[i] / total;
            }
        }
        return static_cast<double>(evaluator.value(A::features(node.state)));
    };

    result.value = expand(0, true);
    if (nodes[0].terminal) {
        result.value = result.searched = *nodes[0].terminal;
        return result;
    }

    std::vector<std::pair<int, std::size_t>> path; // (node, edge taken)
    for (int sim = 0; sim < options.simulations; ++sim) {
        path.clear();
        int at = 0;
        double v;
        while (true) {
            Node& node = nodes[static_cast<std::size_t>(at)];
            if (!node.expanded) {
                v = expand(at, false);
                break;
            }
            if (node.terminal) {
                v = *node.terminal;
                break;
            }
            const double scale =
                options.c_puct * std::sqrt(static_cast<double>(node.visits) + 1.0);
            std::size_t best = 0;
            double best_score = -1e300;
            for (std::size_t i = 0; i < node.edges.size(); ++i) {
                const Edge& e = node.edges[i];
                const double q = e.visits ? e.total / e.visits : 0.0;
                const double score = q + scale * e.prior / (1.0 + e.visits);
                if (score > best_score)
                    best_score = score, best = i;
            }
            Edge& edge = node.edges[best];
            if (edge.child < 0) {
                edge.child = static_cast<int>(nodes.size());
                edge.flips = changes_side<A>(node.state, node.moves[best]);
                const auto next = A::after(node.state, node.moves[best]);
                nodes.push_back({next, {}, {}, std::nullopt, false, 0, -1, -1});
            }
            path.push_back({at, best});
            at = nodes[static_cast<std::size_t>(path.back().first)]
                     .edges[path.back().second]
                     .child;
        }
        // v is the leaf's standing for its side to move; each step up is the other side, unless
        // the move kept the side.
        for (auto it = path.rbegin(); it != path.rend(); ++it) {
            Node& node = nodes[static_cast<std::size_t>(it->first)];
            Edge& edge = node.edges[it->second];
            if (edge.flips)
                v = -v;
            edge.total += v;
            ++edge.visits;
            ++node.visits;
        }
        if (options.prove) {
            for (auto it = path.rbegin(); it != path.rend(); ++it) {
                Node& node = nodes[static_cast<std::size_t>(it->first)];
                int best_win = -1, worst = -1;
                bool all = !node.edges.empty();
                for (const Edge& e : node.edges) {
                    if (e.child < 0) {
                        all = false;
                        continue;
                    }
                    const Node& c = nodes[static_cast<std::size_t>(e.child)];
                    const int w = e.flips ? c.lose : c.win, l = e.flips ? c.win : c.lose;
                    if (w >= 0 && (best_win < 0 || w + 1 < best_win))
                        best_win = w + 1;
                    if (l >= 0)
                        worst = std::max(worst, l + 1);
                    else
                        all = false;
                }
                node.win = best_win;
                node.lose = best_win < 0 && all ? worst : -1;
            }
        }
    }
    result.visits.resize(nodes[0].edges.size());
    result.values.assign(nodes[0].edges.size(), 0.0);
    double sum = 0;
    int tried = 0;
    for (std::size_t i = 0; i < result.visits.size(); ++i) {
        const Edge& e = nodes[0].edges[i];
        result.visits[i] = e.visits;
        if (e.visits)
            result.values[i] = e.total / e.visits;
        sum += e.total;
        tried += e.visits;
    }
    result.searched = tried ? sum / tried : result.value;
    if (options.prove)
        result.proof = nodes[0].win >= 0 ? nodes[0].win : nodes[0].lose >= 0 ? -nodes[0].lose : 0;
    for (int at = 0; result.line.size() < 12;) {
        const Node& node = nodes[static_cast<std::size_t>(at)];
        std::size_t top = 0;
        for (std::size_t i = 1; i < node.edges.size(); ++i)
            if (node.edges[i].visits > node.edges[top].visits)
                top = i;
        if (node.edges.empty() || node.edges[top].visits == 0)
            break;
        result.line.push_back(node.moves[top]);
        at = node.edges[top].child;
    }
    return result;
}

// ── Playing itself ────────────────────────────────────────────────

// What one position of a game says to learn: the rows of its moves and how often the search
// tried each, the features of the state, and how the game ended for the side that was to move.
template <SearchGame A> struct Sample {
    std::vector<typename A::Row> rows;
    std::vector<double> target;
    typename A::Features features;
    double outcome = 0;
    bool kept_side = false; // the move played from this state kept the side to move
};

template <SearchGame A> struct SelfPlayGame {
    std::vector<Sample<A>> samples;
    double first_player_result = 0; // +1 if the side that moved first won
    int plies = 0;
};

// One game from `start`, every move searched. The first `explore_plies` moves are drawn in
// proportion to the visits (so games differ), the rest are the most visited.
template <SearchGame A, Evaluator<A> E>
SelfPlayGame<A> self_play(
    const typename A::State& start, const E& evaluator, const SearchOptions& options,
    int explore_plies, std::mt19937_64& rng
) {
    SelfPlayGame<A> game;
    typename A::State state = start;
    double final_standing = 0; // for the side to move at the end
    for (int ply = 0;; ++ply) {
        const auto result = search<A>(state, evaluator, options, rng);
        if (result.moves.empty()) {
            const auto moves = A::moves(state);
            final_standing = A::outcome(state, moves.empty()).value_or(0.0);
            game.plies = ply;
            break;
        }
        const double total = std::accumulate(result.visits.begin(), result.visits.end(), 0.0);
        if (total == 0) { // a search with no simulations: nothing to learn from, play on
            final_standing = 0;
            game.plies = ply;
            break;
        }
        const std::size_t move = result.pick(ply < explore_plies ? 1.0 : 0.0, rng);
        game.samples.push_back(
            {result.rows, result.distribution(), A::features(state), 0.0,
             !changes_side<A>(state, result.moves[move])}
        );
        state = A::after(state, result.moves[move]);
        // a cut or a finished game is recognised by the adapter's own outcome
        const auto next_moves = A::moves(state);
        if (const auto o = A::outcome(state, next_moves.empty())) {
            final_standing = *o;
            game.plies = ply + 1;
            break;
        }
    }
    // The last state's standing is for its side to move; the sample before it was the other
    // side's unless its move kept the side, and so on back.
    double standing = final_standing;
    for (std::size_t i = game.samples.size(); i-- > 0;) {
        if (!game.samples[i].kept_side)
            standing = -standing;
        game.samples[i].outcome = standing;
    }
    game.first_player_result = standing;
    return game;
}

} // namespace rsc
