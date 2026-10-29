// chess_zero -- the tree search of rsc/core learning chess from nothing, by playing itself.
//
//   chess_zero [--iterations N] [--games G] [--sims S] [--threads T] [--seed X]
//              [--eval-games E] [--eval-sims S] [--eval-every K] [--cap PLIES]
//
// What it knows of chess is the legal moves of a position and the vectors each piece moves by: a
// move is a row of numbers (examples/chess_rows.hpp), a position is thirteen numbers (the
// material, the pieces' reach, how far the pawns have gone, how central the pieces are) and the
// two learned parts are a softmax over the rows (algebra/models/softmax.hpp, 13 weights) and a
// tanh net from the position's numbers to its standing (algebra/models/network.hpp, 13 * 32 +
// 2 * 32 + 1 weights). It begins with none: a softmax that likes every move alike and a net with
// small random weights.
//
// One iteration: G games are played by searching every move (S simulations of the tree), the
// first twelve moves drawn in proportion to the visits so that games differ; the softmax is then
// fitted to the visit counts by AdaGrad and the net to how each game ended by Adam, one pass of
// statistics::fit_pass each, over the new positions and as many older ones drawn at random. A
// game ends where the rules end it (discrete::outcome), or is cut at 120 plies and decided by
// material if a side is ahead by three points; without that cut nearly every game between
// players this weak is a draw and the net has nothing to learn from.
//
// Every few iterations the player is measured against three fixed ones, E games each, half as
// each colour, searching with `--eval-sims` simulations and playing the most visited move: a
// random player, a greedy one (the most it can take, less what it leaves hanging, a mate above
// all), and a search two plies deep over material alone. Points are wins plus half the draws.
// The weights live in memory for the run.

#include <rsc/core/tree_search.hpp>
#include <spatium/algebra/analysis/optimize.hpp>
#include <spatium/algebra/models/network.hpp>
#include <spatium/algebra/models/softmax.hpp>
#include <spatium/discrete/position.hpp>
#include <spatium/statistics/fit.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <optional>
#include <print>
#include <random>
#include <span>
#include <string>
#include <thread>
#include <vector>

#include "chess_rows.hpp"

using namespace spatium::discrete;
using namespace chess_rows;

namespace {

constexpr std::size_t state_columns = 13, hidden_units = 32;
using Net = spatium::TanhNet<state_columns, hidden_units>;

// The two learned parts, as the search asks for them.
struct ChessNets {
    spatium::LinearSoftmax<columns> policy;
    Net net = Net::random(1);
    void prior(std::span<const Row> rows, std::span<double> out) const { policy.distribution(rows, out); }
    double value(const std::array<double, state_columns>& f) const { return net.value(f); }
    static constexpr std::size_t weight_count = columns + Net::weight_count;
};

struct Chess {
    struct State {
        Position pos;
        int ply = 0;
    };
    using Move = spatium::discrete::Move;
    using Row = chess_rows::Row;
    using Features = std::array<double, state_columns>;

    // The ply at which a game is cut and decided by material.
    static inline int cap = 120;

    static std::vector<Move> moves(const State& s) {
        const MoveList list = legal_moves(s.pos);
        return {list.begin(), list.end()};
    }
    static State after(const State& s, const Move& m) { return {s.pos.made(m), s.ply + 1}; }

    // The rules' end of a game (mate, stalemate, fifty moves, a dead position), else the cut. A
    // third repetition needs the game's history, which a search state does not carry.
    static std::optional<double> outcome(const State& s, bool) {
        switch (spatium::discrete::outcome(s.pos)) {
        case Outcome::Checkmate:
            return -1.0;
        case Outcome::Ongoing:
            break;
        default:
            return 0.0;
        }
        if (s.ply >= cap) {
            const Color me = s.pos.side_to_move();
            const double lead = material(s.pos, me) - material(s.pos, opposite(me));
            return lead >= 3 ? 1.0 : lead <= -3 ? -1.0 : 0.0;
        }
        return std::nullopt;
    }

    static std::vector<Row> rows(const State& s, std::span<const Move> moves) {
        std::vector<Row> out;
        out.reserve(moves.size());
        for (const Move& m : moves)
            out.push_back(row_of(s.pos, m));
        return out;
    }

    static Features features(const State& s) {
        const Position& p = s.pos;
        const Color me = p.side_to_move(), you = opposite(me);
        Features f{};
        constexpr std::array<PieceType, 5> kinds{
            PieceType::Pawn, PieceType::Knight, PieceType::Bishop, PieceType::Rook, PieceType::Queen
        };
        constexpr std::array<double, 5> most{8, 2, 2, 2, 1};
        for (std::size_t i = 0; i < 5; ++i)
            f[i] = (p.pieces(me, kinds[i]).count() - p.pieces(you, kinds[i]).count()) / most[i];
        f[5] = material(p, me) / 39.0;
        f[6] = material(p, you) / 39.0;
        f[7] = p.in_check() ? 1.0 : 0.0;
        f[8] = legal_moves(p).size / 40.0;
        auto advance = [&](Color c) {
            double sum = 0;
            for (Square s : p.pieces(c, PieceType::Pawn))
                sum += (c == Color::White ? s.rank() : 7 - s.rank()) - 1;
            return sum / 40.0;
        };
        auto centre = [&](Color c) {
            double sum = 0;
            for (PieceType t :
                 {PieceType::Knight, PieceType::Bishop, PieceType::Rook, PieceType::Queen})
                for (Square s : p.pieces(c, t))
                    sum += 1.0 - std::max(std::abs(s.file() - 3.5), std::abs(s.rank() - 3.5)) / 3.5;
            return sum / 8.0;
        };
        f[9] = advance(me);
        f[10] = advance(you);
        f[11] = centre(me);
        f[12] = centre(you);
        return f;
    }
};
static_assert(rsc::SearchGame<Chess>);
static_assert(rsc::Evaluator<ChessNets, Chess>);

// ── The players it is measured against ────────────────────────────

using Mover = std::function<Move(const Position&, int, std::mt19937_64&)>;

Move random_mover(const Position& p, int, std::mt19937_64& rng) {
    const MoveList moves = legal_moves(p);
    return moves.moves[rng() % static_cast<std::size_t>(moves.size)];
}

Move greedy_mover(const Position& p, int, std::mt19937_64& rng) {
    const MoveList moves = legal_moves(p);
    std::vector<Move> best;
    double top = -1e300;
    for (const Move& m : moves) {
        const Judged j = judge(p, m);
        const double v = 100 * j.mate + j.captured + j.promoted - j.hanging + 0.5 * j.check;
        if (v > top)
            top = v, best.clear();
        if (v == top)
            best.push_back(m);
    }
    return best[rng() % best.size()];
}

// Material and mates, `depth` plies, alpha-beta.
int negamax(const Position& p, int depth, int alpha, int beta) {
    const MoveList moves = legal_moves(p);
    if (!moves.size)
        return p.in_check() ? -1000 - depth : 0;
    if (depth == 0) {
        const Color me = p.side_to_move();
        return static_cast<int>(material(p, me) - material(p, opposite(me)));
    }
    int best = -100000;
    for (const Move& m : moves) {
        best = std::max(best, -negamax(p.made(m), depth - 1, -beta, -alpha));
        alpha = std::max(alpha, best);
        if (alpha >= beta)
            break;
    }
    return best;
}

Move material_mover(const Position& p, int, std::mt19937_64& rng) {
    const MoveList moves = legal_moves(p);
    std::vector<Move> best;
    int top = -100000;
    for (const Move& m : moves) {
        const int v = -negamax(p.made(m), 1, -100000, 100000);
        if (v > top)
            top = v, best.clear();
        if (v == top)
            best.push_back(m);
    }
    return best[rng() % best.size()];
}

Mover agent_mover(const ChessNets& nets, int simulations) {
    return [&nets, simulations](const Position& p, int ply, std::mt19937_64& rng) {
        rsc::SearchOptions options;
        options.simulations = simulations;
        const auto result = rsc::search<Chess>({p, ply}, nets, options, rng);
        return result.moves[result.pick(0, rng)];
    };
}

// One game, `white` against `black`. +1 if white wins, -1 if black does.
int referee(const Mover& white, const Mover& black, std::mt19937_64& rng) {
    Chess::State state{Position::initial(), 0};
    while (true) {
        const auto moves = Chess::moves(state);
        if (const auto o = Chess::outcome(state, moves.empty())) {
            const double for_white = state.pos.side_to_move() == Color::White ? *o : -*o;
            return for_white > 0 ? 1 : for_white < 0 ? -1 : 0;
        }
        const Move m = (state.ply % 2 == 0 ? white : black)(state.pos, state.ply, rng);
        state = Chess::after(state, m);
    }
}

struct Points {
    int win = 0, draw = 0, loss = 0;
    double value() const {
        const int n = win + draw + loss;
        return n ? (win + 0.5 * draw) / n : 0.0;
    }
};

Points match(const Mover& player, const Mover& opponent, int games, std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    Points p;
    for (int g = 0; g < games; ++g) {
        const bool white = g % 2 == 0;
        const int r = white ? referee(player, opponent, rng) : -referee(opponent, player, rng);
        (r > 0 ? p.win : r < 0 ? p.loss : p.draw) += 1;
    }
    return p;
}

void evaluate(const ChessNets& nets, int games, int simulations, std::uint64_t seed) {
    constexpr int eval_cap = 200;
    const int saved = Chess::cap;
    Chess::cap = eval_cap;
    const Mover agent = agent_mover(nets, simulations);
    const Points r = match(agent, random_mover, games, seed);
    const Points g = match(agent, greedy_mover, games, seed + 1);
    const Points m = match(agent, material_mover, games, seed + 2);
    std::println(
        "    points over {} games: random {:.3f} ({}-{}-{})  greedy {:.3f} ({}-{}-{})  "
        "material 2 plies {:.3f} ({}-{}-{})",
        games, r.value(), r.win, r.draw, r.loss, g.value(), g.win, g.draw, g.loss, m.value(), m.win,
        m.draw, m.loss
    );
    Chess::cap = saved;
}

// One pass over `sample`: the softmax toward the visits, the net toward the results. Returns the
// net's mean squared error before its steps.
struct Trainer {
    spatium::AdaGrad<columns> policy_rule;
    spatium::Adam<Net::weight_count> net_rule;

    double fit(ChessNets& nets, std::span<const rsc::Sample<Chess>> sample, std::mt19937_64& rng) {
        using spatium::statistics::Fit;
        spatium::statistics::fit_pass(nets.policy.weights(), policy_rule, sample, 0.3, rng,
                                      [&](const rsc::Sample<Chess>& s) {
                                          const std::span<const Row> rows(s.rows);
                                          const std::span<const double> target(s.target);
                                          return Fit<std::array<double, columns>>{
                                              nets.policy.cross_entropy(rows, target),
                                              nets.policy.gradient(rows, target)};
                                      });
        return spatium::statistics::fit_pass(nets.net.weights(), net_rule, sample, 0.003, rng,
                                             [&](const rsc::Sample<Chess>& s) {
                                                 return Fit<Net::Weights>{nets.net.loss(s.features, s.outcome),
                                                                          nets.net.gradient(s.features, s.outcome)};
                                             });
    }
};

} // namespace

int main(int argc, char** argv) {
    int iterations = 60, games = 24, simulations = 48, threads = 2, eval_games = 20,
        eval_simulations = 48, eval_every = 10, train_cap = 120;
    std::uint64_t seed = 1;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&]() -> std::string { return i + 1 < argc ? argv[++i] : ""; };
        if (a == "--iterations")
            iterations = std::atoi(next().c_str());
        else if (a == "--games")
            games = std::atoi(next().c_str());
        else if (a == "--sims")
            simulations = std::atoi(next().c_str());
        else if (a == "--threads")
            threads = std::max(1, std::atoi(next().c_str()));
        else if (a == "--eval-games")
            eval_games = std::atoi(next().c_str());
        else if (a == "--eval-sims")
            eval_simulations = std::atoi(next().c_str());
        else if (a == "--eval-every")
            eval_every = std::max(1, std::atoi(next().c_str()));
        else if (a == "--cap")
            train_cap = std::atoi(next().c_str());
        else if (a == "--seed")
            seed = std::strtoull(next().c_str(), nullptr, 10);
        else {
            std::println(stderr, "unknown argument {}", a);
            return 2;
        }
    }

    ChessNets nets;
    nets.net = Net::random(seed);
    Trainer trainer;
    std::println("{} weights: {} in the softmax, {} in the net", ChessNets::weight_count, columns,
                 Net::weight_count);
    std::println("untrained:");
    evaluate(nets, eval_games, eval_simulations, seed + 5000);

    std::vector<rsc::Sample<Chess>> replay;
    constexpr std::size_t replay_limit = 30000;
    std::mt19937_64 rng(seed);
    rsc::SearchOptions options;
    options.simulations = simulations;
    options.noise = 0.25;
    Chess::cap = train_cap;
    const auto begin = std::chrono::steady_clock::now();
    for (int iteration = 1; iteration <= iterations; ++iteration) {
        std::vector<std::vector<rsc::Sample<Chess>>> per_thread(static_cast<std::size_t>(threads));
        std::vector<std::array<int, 4>> tally(static_cast<std::size_t>(threads)); // w, d, l, plies
        {
            std::vector<std::jthread> pool;
            for (int t = 0; t < threads; ++t)
                pool.emplace_back([&, t] {
                    std::mt19937_64 local(seed * 1000003 + static_cast<std::uint64_t>(iteration) * 131 +
                                          static_cast<std::uint64_t>(t));
                    const int mine = games / threads + (t < games % threads ? 1 : 0);
                    for (int g = 0; g < mine; ++g) {
                        auto game = rsc::self_play<Chess>(Chess::State{Position::initial(), 0}, nets,
                                                          options, 12, local);
                        auto& tl = tally[static_cast<std::size_t>(t)];
                        (game.first_player_result > 0   ? tl[0]
                         : game.first_player_result < 0 ? tl[2]
                                                        : tl[1]) += 1;
                        tl[3] += game.plies;
                        for (auto& s : game.samples)
                            per_thread[static_cast<std::size_t>(t)].push_back(std::move(s));
                    }
                });
        }
        std::vector<rsc::Sample<Chess>> fresh;
        std::array<int, 4> total{};
        for (std::size_t t = 0; t < per_thread.size(); ++t) {
            for (auto& s : per_thread[t])
                fresh.push_back(std::move(s));
            for (int k = 0; k < 4; ++k)
                total[static_cast<std::size_t>(k)] += tally[t][static_cast<std::size_t>(k)];
        }
        const double error = trainer.fit(nets, std::span<const rsc::Sample<Chess>>(fresh), rng);
        // an equal number of older positions, drawn at random, so that the net does not forget
        if (!replay.empty()) {
            std::vector<rsc::Sample<Chess>> old;
            for (std::size_t i = 0; i < fresh.size() && i < replay.size(); ++i)
                old.push_back(replay[rng() % replay.size()]);
            trainer.fit(nets, std::span<const rsc::Sample<Chess>>(old), rng);
        }
        for (auto& s : fresh)
            replay.push_back(std::move(s));
        if (replay.size() > replay_limit)
            replay.erase(replay.begin(),
                         replay.begin() + static_cast<std::ptrdiff_t>(replay.size() - replay_limit));
        const double seconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
        std::println("iteration {:3}: {} games ({} white wins, {} draws, {} black wins), {:.0f} plies a "
                     "game, value error {:.3f}, {} positions kept, {:.0f} s",
                     iteration, games, total[0], total[1], total[2],
                     double(total[3]) / std::max(games, 1), error, replay.size(), seconds);
        if (iteration % eval_every == 0 || iteration == iterations)
            evaluate(nets, eval_games, eval_simulations, seed + 5000);
    }
    return 0;
}
