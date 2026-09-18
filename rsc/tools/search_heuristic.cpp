// Does a trained policy earn its keep as a search heuristic?
//
// This is the question ROADMAP's "RSC as search, not classification" entry
// proposes and nobody has asked with a number. The tempting version of it
// -- "can the model solve the task" -- is the wrong one, because an
// uninformed breadth-first search also solves the task, every time, and
// optimally. A heuristic is not there to make the answer reachable. It is
// there to make it cheaper, and the only way that claim can be false or
// true is in nodes expanded.
//
// So three searches run on identical problems:
//
//   breadth-first     uninformed, expands by depth, always optimal
//   best-first        ordered by the policy's own log-probability, so the
//                     sequences it considers likely are tried first
//   greedy rollout    the policy alone, no search at all -- the degenerate
//                     case, and the one the capacity sweep measured
//
// Reported per search: how often it found the target inside a fixed node
// budget, how many nodes that took, and how often what it found was
// actually shortest. Optimality is separated from success on purpose --
// best-first is free to return a longer chain than breadth-first, and a
// heuristic that saves nodes by accepting worse answers has not saved
// anything, it has changed the question.
//
// The cost side is already known from dispatch_cost: one policy
// evaluation is about 600 ns for a hidden-64 network, or 8 ns once
// distilled into per-domain trees. Against that, a node expansion here is
// five integer operations and a map probe. The break-even is therefore not
// obvious in either direction, which is the reason to measure rather than
// assume.
//
// Build: part of the rsc tools target.

#include "toy_search.hpp"

#include <algorithm>
#include <map>
#include <print>
#include <queue>
#include <random>
#include <vector>

namespace {

constexpr std::size_t kNodeBudget = 20000;

struct Outcome {
    bool solved = false;
    std::size_t expanded = 0;
    std::size_t length = 0;
};

// Uninformed, by depth. The reference every other row is measured against,
// and the only one guaranteed to return a shortest chain.
Outcome breadth_first(const toy::Problem& p, std::size_t max_depth) {
    struct Item { long long state; std::size_t depth; };
    std::deque<Item> queue{{p.start, 0}};
    std::map<long long, bool> seen{{p.start, true}};
    Outcome out;

    while (!queue.empty() && out.expanded < kNodeBudget) {
        const auto [s, d] = queue.front();
        queue.pop_front();
        ++out.expanded;
        if (d >= max_depth) continue;
        for (int a = 0; a < toy::kNumActions; ++a) {
            const long long n = toy::apply_action(s, a);
            if (std::llabs(n) > toy::kStateLimit || seen.contains(n)) continue;
            if (n == p.target) { out.solved = true; out.length = d + 1; return out; }
            seen[n] = true;
            queue.push_back({n, d + 1});
        }
    }
    return out;
}

// Ordered by accumulated -log probability under the policy, so the most
// plausible continuations are expanded first. Depth is still bounded, so
// this cannot "win" by searching deeper than breadth-first was allowed to.
Outcome best_first(const rsc::Dispatcher& model, const toy::Problem& p, std::size_t max_depth) {
    struct Item {
        double cost;          // accumulated -log p
        long long state;
        std::size_t depth;
        bool operator>(const Item& o) const { return cost > o.cost; }
    };
    std::priority_queue<Item, std::vector<Item>, std::greater<Item>> frontier;
    frontier.push({0.0, p.start, 0});
    std::map<long long, double> best{{p.start, 0.0}};
    Outcome out;

    while (!frontier.empty() && out.expanded < kNodeBudget) {
        const Item cur = frontier.top();
        frontier.pop();
        ++out.expanded;
        if (cur.depth >= max_depth) continue;

        auto cache = model.forward_cached(toy::features(cur.state, p.target));
        auto probs = rsc::softmax(cache.logits);

        for (int a = 0; a < toy::kNumActions; ++a) {
            const long long n = toy::apply_action(cur.state, a);
            if (std::llabs(n) > toy::kStateLimit) continue;
            const double step = -std::log(std::max(probs[static_cast<std::size_t>(a)], 1e-12));
            const double cost = cur.cost + step;
            if (n == p.target) { out.solved = true; out.length = cur.depth + 1; return out; }
            auto it = best.find(n);
            if (it != best.end() && it->second <= cost) continue;
            best[n] = cost;
            frontier.push({cost, n, cur.depth + 1});
        }
    }
    return out;
}

// The policy with no search behind it. One expansion per step, which is
// what makes it the cheapest row by a wide margin and also the one that
// gives up most often.
Outcome greedy(const rsc::Dispatcher& model, const toy::Problem& p, std::size_t max_steps,
               std::mt19937_64& rng) {
    const toy::Episode ep = toy::rollout(model, p, max_steps, rng, /*greedy=*/true);
    return Outcome{ep.solved, ep.length, ep.length};
}

}  // namespace

int main() {
    constexpr std::size_t kInputDim = 6;
    constexpr std::size_t kHidden = 128;
    constexpr std::size_t kEpisodes = 60000;
    constexpr std::size_t kBatch = 32;
    constexpr double kLearningRate = 0.05;
    constexpr std::size_t kProblems = 400;

    std::println("Is a trained policy worth its cost as a search heuristic?");
    std::println("hidden {}, {} training episodes, {} problems per row, node budget {}",
                 kHidden, kEpisodes, kProblems, kNodeBudget);
    std::println("");

    for (std::size_t walk : {3u, 4u, 5u, 6u}) {
        rsc::Dispatcher model(kInputDim, kHidden, toy::kNumActions, /*seed=*/1234);
        toy::train(model, walk, kEpisodes, kLearningRate, kBatch, /*seed=*/7);

        std::mt19937_64 rng(99);
        struct Acc { std::size_t solved = 0, optimal = 0; double nodes = 0.0; };
        Acc bfs, bf, gr;

        for (std::size_t i = 0; i < kProblems; ++i) {
            const toy::Problem p = toy::sample_problem(rng, walk);

            const Outcome a = breadth_first(p, walk);
            const Outcome b = best_first(model, p, walk);
            const Outcome c = greedy(model, p, walk + 2, rng);

            const auto tally = [&](Acc& acc, const Outcome& o) {
                acc.nodes += static_cast<double>(o.expanded);
                if (o.solved) { ++acc.solved; if (o.length <= p.optimal) ++acc.optimal; }
            };
            tally(bfs, a); tally(bf, b); tally(gr, c);
        }

        const double n = static_cast<double>(kProblems);
        std::println("walk {} -- target reached by {} random actions", walk, walk);
        std::println("  {:<16} | {:>8} | {:>9} | {:>12}",
                     "search", "solved", "shortest", "nodes / problem");
        const auto row = [&](const char* name, const Acc& acc) {
            std::println("  {:<16} | {:>7.1f}% | {:>8.1f}% | {:>12.1f}",
                         name,
                         100.0 * static_cast<double>(acc.solved) / n,
                         acc.solved ? 100.0 * static_cast<double>(acc.optimal) /
                                          static_cast<double>(acc.solved) : 0.0,
                         acc.nodes / n);
        };
        row("breadth-first", bfs);
        row("best-first", bf);
        row("greedy rollout", gr);
        std::println("");
    }

    std::println("A heuristic that expands fewer nodes while returning longer chains has");
    std::println("not saved work, it has answered a different question -- which is why");
    std::println("the shortest column sits next to the node count and not below it.");
    return 0;
}
