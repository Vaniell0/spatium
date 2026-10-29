// chess_agent -- an agent that knows how each piece moves and nothing else, learning to choose
// among the legal moves it is handed.
//
//   chess_agent [--observe G] [--play G] [--distill G] [--reinforce G] [--seed S]
//
// The agent is `LinearSoftmax` (algebra/models/softmax.hpp), trained by `AdaGrad`
// (algebra/analysis/optimize.hpp), over rows of numbers (discrete/move_features.hpp). A row has two
// blocks. The geometric one is made from the piece's vectors and the two squares of the move, on
// the same code that serves a box of any number of axes (here two): how far along each axis, how
// many steps of the piece's vector, whether it slid, how far from the centre it ended. The other is
// this program's own idea of what a move is worth, the part a game's author writes: what is
// captured, whether it checks or mates, whether the piece is left where it can be taken, what it
// promotes to. The agent starts with every weight zero, so it plays at random.
//
// Three ways to learn from lists of legal moves, each measured:
//
//   --observe G   a hand-written player (the most material won, a piece left hanging counted
//                 against it, a mate above all) plays G games against a random one, and each of
//                 its choices is shown to the agent with the list it chose from. Afterwards the
//                 agent is asked to choose in positions it has not seen and is compared with the
//                 teacher, where the teacher has something to prefer.
//   --play G      the untrained agent, the agent trained by `--observe`, and the teacher, each
//                 against a random player for G games.
//   --distill G   the agent trained by `--observe`, distilled into trees of depth 1 to 8 grown on
//                 the moves of G random games and measured on G others: the divergence of the
//                 tree's distribution from the agent's, how often their best moves agree, and the
//                 tree's description length in bits, beside the agent's own.
//   --reinforce G the agent plays a random player and, after each game, every move it chose is
//                 made likelier if the game was won and less likely if it was lost. No teacher.
//
// A game ends where the rules end it (`discrete::outcome`: mate, stalemate, the fifty-move rule,
// a dead position, a third repetition), and is otherwise cut at 200 plies and decided by material
// if one side is ahead by three points.

#include <spatium/algebra/analysis/optimize.hpp>
#include <spatium/algebra/models/softmax.hpp>
#include <spatium/algebra/models/tree.hpp>
#include <spatium/discrete/move_features.hpp>
#include <spatium/statistics/information.hpp>
#include <spatium/statistics/tree_fit.hpp>
#include <spatium/discrete/position.hpp>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <print>
#include <random>
#include <string>
#include <vector>

#include "chess_rows.hpp"

using namespace spatium::discrete;
using namespace chess_rows;

namespace {

constexpr int ply_limit = 200;
using Agent = spatium::LinearSoftmax<columns>;
using Rule = spatium::AdaGrad<columns>;

// One step of `rule` toward the move chosen from `rows`, scaled by `weight` (negative moves
// away from it).
void learn(Agent& agent, Rule& rule, const std::vector<Row>& rows, std::size_t chosen, double rate,
           double weight = 1) {
    auto g = agent.gradient(rows, chosen);
    for (double& v : g)
        v *= weight;
    rule.step(agent.weights(), g, rate);
}

// The hand-written player: material won, less the piece it leaves to be taken, a mate above all.
double teacher_value(const Position& p, const Move& m) {
    const Judged j = judge(p, m);
    return 100 * j.mate + j.captured + j.promoted - j.hanging + 0.5 * j.check;
}

std::size_t teacher_choice(const Position& p, const MoveList& moves, std::mt19937_64& rng) {
    std::vector<std::size_t> best;
    double top = -1e300;
    for (int i = 0; i < moves.size; ++i) {
        const double v = teacher_value(p, moves.moves[static_cast<std::size_t>(i)]);
        if (v > top)
            top = v, best.clear();
        if (v == top)
            best.push_back(static_cast<std::size_t>(i));
    }
    return best[rng() % best.size()];
}

enum class Result { Win, Draw, Loss };

// The result for `first`: what the rules say at the end (a mate is a loss for the side to move,
// every other end a draw), or at the cut the material, three points deciding.
Result result_for(const Position& p, Outcome end, Color first) {
    if (end == Outcome::Checkmate)
        return p.side_to_move() == first ? Result::Loss : Result::Win;
    if (end != Outcome::Ongoing)
        return Result::Draw;
    const double lead = material(p, first) - material(p, opposite(first));
    return lead >= 3 ? Result::Win : lead <= -3 ? Result::Loss : Result::Draw;
}

struct Chosen {
    std::vector<Row> rows;
    std::size_t index;
};

enum class Player { Random, Agent, Teacher };

struct Score {
    int win = 0, draw = 0, loss = 0, mates = 0;
    void add(Result r, bool mate) {
        (r == Result::Win ? win : r == Result::Draw ? draw : loss) += 1;
        mates += mate && r == Result::Win;
    }
    int games() const { return win + draw + loss; }
    double points() const { return games() ? (win + 0.5 * draw) / games() : 0.0; }
};

// One game of `first` against `second`. The agent's choices (the lists and the index) are kept
// when `record` is given. Returns the result for `first`.
Result play(
    Player first, Player second, const Agent& agent, double temperature, std::mt19937_64& rng,
    std::vector<Chosen>* record = nullptr, bool* mated = nullptr
) {
    std::vector<Position> game{Position::initial()};
    std::uniform_real_distribution<double> u(0, 1);
    Outcome end = Outcome::Ongoing;
    for (int ply = 0; ply < ply_limit; ++ply) {
        const Position& p = game.back();
        if ((end = outcome(std::span<const Position>(game))) != Outcome::Ongoing)
            break;
        const MoveList moves = legal_moves(p);
        const Player who = (ply % 2 == 0) ? first : second;
        std::size_t pick = rng() % static_cast<std::size_t>(moves.size);
        if (who == Player::Agent) {
            std::vector<Row> rows = rows_of(p, moves);
            pick = agent.sample(rows, u(rng), temperature);
            if (record)
                record->push_back({std::move(rows), pick});
        } else if (who == Player::Teacher) {
            pick = teacher_choice(p, moves, rng);
        }
        game.push_back(p.made(moves.moves[pick]));
    }
    if (end == Outcome::Ongoing)
        end = outcome(std::span<const Position>(game));
    if (mated)
        *mated = end == Outcome::Checkmate;
    // `first` moves on even plies, and is White.
    return result_for(game.back(), end, Color::White);
}

Score versus_random(
    Player who, const Agent& agent, double temperature, int games, std::uint64_t seed
) {
    std::mt19937_64 rng(seed);
    Score score;
    for (int g = 0; g < games; ++g) {
        const bool white = g % 2 == 0;
        bool mated = false;
        Result r = white ? play(who, Player::Random, agent, temperature, rng, nullptr, &mated)
                         : play(Player::Random, who, agent, temperature, rng, nullptr, &mated);
        if (!white)
            r = r == Result::Win ? Result::Loss : r == Result::Loss ? Result::Win : r;
        score.add(r, mated);
    }
    return score;
}

void print_score(const char* name, const Score& s) {
    std::println(
        "{:<22} {:4} won ({:3} by mate)  {:4} drawn  {:4} lost   points {:.3f}", name, s.win,
        s.mates, s.draw, s.loss, s.points()
    );
}

// How often the agent's greedy choice is one of the teacher's best, in the positions where the
// teacher has something to prefer (its best value is above zero) and in all of them.
struct Agreement {
    int decisive = 0, decisive_same = 0, all = 0, all_same = 0;
};

Agreement agreement(const Agent& agent, int games, std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    Agreement a;
    for (int g = 0; g < games; ++g) {
        Position p = Position::initial();
        for (int ply = 0; ply < ply_limit; ++ply) {
            const MoveList moves = legal_moves(p);
            if (!moves.size)
                break;
            double top = -1e300;
            for (const Move& m : moves)
                top = std::max(top, teacher_value(p, m));
            const std::vector<Row> rows = rows_of(p, moves);
            const bool same = teacher_value(p, moves.moves[agent.choose(rows)]) == top;
            ++a.all, a.all_same += same;
            if (top > 0)
                ++a.decisive, a.decisive_same += same;
            p = p.made(moves.moves[rng() % static_cast<std::size_t>(moves.size)]);
        }
    }
    return a;
}

int observe(int games, std::uint64_t seed, double rate, Agent& agent) {
    Rule rule;
    std::println("\nobserving the teacher against a random player");
    std::println(
        "games watched   moves shown   agrees where it matters   (positions)   agrees overall"
    );
    std::mt19937_64 rng(seed);
    int shown = 0, failures = 0;
    auto report = [&](int watched) {
        const Agreement a = agreement(agent, 20, seed + 4242);
        std::println(
            "{:13}   {:11}   {:21.1f}%   {:11}   {:13.1f}%", watched, shown,
            a.decisive ? 100.0 * a.decisive_same / a.decisive : 0.0, a.decisive,
            100.0 * a.all_same / a.all
        );
        return a.decisive ? double(a.decisive_same) / a.decisive : 0.0;
    };
    report(0);
    double last = 0;
    int next_report = 1;
    for (int g = 1; g <= games; ++g) {
        Position p = Position::initial();
        const bool teacher_white = g % 2 == 0;
        for (int ply = 0; ply < ply_limit; ++ply) {
            const MoveList moves = legal_moves(p);
            if (!moves.size)
                break;
            const bool teachers_turn = (ply % 2 == 0) == teacher_white;
            std::size_t pick = rng() % static_cast<std::size_t>(moves.size);
            if (teachers_turn) {
                pick = teacher_choice(p, moves, rng);
                learn(agent, rule, rows_of(p, moves), pick, rate);
                ++shown;
            }
            p = p.made(moves.moves[pick]);
        }
        if (g == next_report || g == games) {
            last = report(g);
            next_report *= 4;
        }
    }
    if (games > 0 && last < 0.9) {
        std::println("the agent does not agree with the teacher where it matters");
        ++failures;
    }
    return failures;
}

// The positions of `games` games of random play, each as the rows of its legal moves and the
// index of the teacher's choice there.
struct Seen {
    std::vector<Row> rows;
    std::size_t teacher;
};
std::vector<Seen> positions(int games, std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::vector<Seen> out;
    for (int g = 0; g < games; ++g) {
        Position p = Position::initial();
        for (int ply = 0; ply < ply_limit; ++ply) {
            const MoveList moves = legal_moves(p);
            if (!moves.size)
                break;
            out.push_back({rows_of(p, moves), teacher_choice(p, moves, rng)});
            p = p.made(moves.moves[rng() % static_cast<std::size_t>(moves.size)]);
        }
    }
    return out;
}

// The agent distilled into trees of growing depth, each grown on the agent's scores of the moves
// of one set of positions and measured on another: the mean KL(agent | tree) of the two
// distributions over the legal moves, how often the tree's best move is the agent's and the
// teacher's, and the tree's description length. The agent's own line is first: its weights as
// the description, nothing lost.
void distill(const Agent& agent, int games, std::uint64_t seed) {
    using Tree = spatium::ScoreTree<columns>;
    const auto train = positions(games, seed), held = positions(games, seed + 777);
    std::vector<Row> rows;
    std::vector<double> scores;
    for (const Seen& s : train)
        for (const Row& r : s.rows)
            rows.push_back(r), scores.push_back(agent.score(r));
    std::println("\ndistilling the agent into trees: {} moves of {} positions to grow on, {} positions to measure on",
                 rows.size(), train.size(), held.size());
    std::println("model            bits    KL(agent | model)   agrees with agent   with teacher");
    auto measure = [&](auto&& score_of) {
        double kl = 0;
        int with_agent = 0, with_teacher = 0;
        for (const Seen& s : held) {
            std::vector<double> mine(s.rows.size()), p(s.rows.size()), q(s.rows.size());
            for (std::size_t i = 0; i < s.rows.size(); ++i)
                mine[i] = score_of(s.rows[i]);
            agent.distribution(s.rows, p);
            spatium::softmax(std::span<const double>(mine), std::span<double>(q));
            kl += spatium::statistics::kl_divergence(std::span<const double>(p), std::span<const double>(q));
            const auto best = static_cast<std::size_t>(std::ranges::max_element(mine) - mine.begin());
            with_agent += best == agent.choose(s.rows);
            with_teacher += best == s.teacher;
        }
        const double n = double(held.size());
        return std::array{kl / n, with_agent / n, with_teacher / n};
    };
    const auto own = measure([&](const Row& r) { return agent.score(r); });
    std::println("{:<14} {:6}    {:17.4f}   {:16.1f}%   {:11.1f}%", "linear agent", columns * 32, own[0],
                 100 * own[1], 100 * own[2]);
    for (std::size_t depth : {1u, 2u, 4u, 6u, 8u}) {
        const Tree tree = spatium::statistics::fit_tree<columns>(std::span<const Row>(rows),
                                                                  std::span<const double>(scores),
                                                                  spatium::statistics::TreeOptions{.depth = depth, .min_leaf = 8});
        const auto m = measure([&](const Row& r) { return tree.score(r); });
        std::println("tree, depth {:<2} {:6.0f}    {:17.4f}   {:16.1f}%   {:11.1f}%", depth, tree.bits(), m[0],
                     100 * m[1], 100 * m[2]);
    }
}

void reinforce(int games, std::uint64_t seed) {
    constexpr double rate = 0.05;
    std::println("\nreinforcing from nothing, against a random player (no teacher)");
    std::println("games   points over the last block (sampling at temperature 1)");
    Agent agent;
    Rule rule;
    std::mt19937_64 rng(seed);
    double baseline = 0;
    const int block = std::max(games / 10, 1);
    Score recent;
    for (int g = 1; g <= games; ++g) {
        const bool white = g % 2 == 0;
        std::vector<Chosen> record;
        bool mated = false;
        Result r = white ? play(Player::Agent, Player::Random, agent, 1.0, rng, &record, &mated)
                         : play(Player::Random, Player::Agent, agent, 1.0, rng, &record, &mated);
        if (!white)
            r = r == Result::Win ? Result::Loss : r == Result::Loss ? Result::Win : r;
        recent.add(r, mated);
        const double reward = r == Result::Win ? 1.0 : r == Result::Loss ? -1.0 : 0.0;
        for (const Chosen& c : record)
            learn(agent, rule, c.rows, c.index, rate, reward - baseline);
        baseline += 0.02 * (reward - baseline);
        if (g % block == 0) {
            std::println("{:5}   {:.3f}", g, recent.points());
            recent = {};
        }
    }
    print_score("after, greedy-ish", versus_random(Player::Agent, agent, 0.2, 200, seed + 7));
    std::println("weights (geometric block, then captured, check, mate, hanging, promoted, moved):"
    );
    for (std::size_t i = 0; i < columns; ++i)
        std::print("{:8.3f}", agent.weights()[i]);
    std::println("");
}

} // namespace

int main(int argc, char** argv) {
    int observe_games = 0, play_games = 0, reinforce_games = 0, distill_games = 0;
    std::uint64_t seed = 1;
    double rate = 10;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&]() -> std::string { return i + 1 < argc ? argv[++i] : ""; };
        if (a == "--observe")
            observe_games = std::atoi(next().c_str());
        else if (a == "--play")
            play_games = std::atoi(next().c_str());
        else if (a == "--distill")
            distill_games = std::atoi(next().c_str());
        else if (a == "--reinforce")
            reinforce_games = std::atoi(next().c_str());
        else if (a == "--rate")
            rate = std::atof(next().c_str());
        else if (a == "--seed")
            seed = std::strtoull(next().c_str(), nullptr, 10);
        else {
            std::println(stderr, "unknown argument {}", a);
            return 2;
        }
    }
    if (!observe_games && !play_games && !reinforce_games && !distill_games)
        observe_games = 256, play_games = 200, reinforce_games = 2000;

    int failures = 0;
    Agent trained;
    if (observe_games)
        failures += observe(observe_games, seed, rate, trained);
    if (play_games) {
        std::println("\n{} games against a random player, half as each colour", play_games);
        print_score(
            "untrained agent", versus_random(Player::Agent, Agent{}, 0.2, play_games, seed + 1)
        );
        print_score(
            "agent that watched", versus_random(Player::Agent, trained, 0.2, play_games, seed + 1)
        );
        print_score(
            "the teacher", versus_random(Player::Teacher, Agent{}, 0.2, play_games, seed + 1)
        );
    }
    if (distill_games)
        distill(trained, distill_games, seed + 31);
    if (reinforce_games)
        reinforce(reinforce_games, seed);
    std::println("{}", failures ? "FAILED" : "checks passed");
    return failures ? 1 : 0;
}
