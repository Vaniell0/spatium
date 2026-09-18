// What does a learned dispatch actually cost, next to the branch it
// replaces?
//
// The question this answers is not "is the model accurate" but "can it sit
// where an `if` sits". Those are different questions and only the second
// decides whether a trained dispatcher can live in a hot path rather than
// being consulted once per scene.
//
// Timed on the same inputs:
//
//   if-cascade      the hand-written branch the model replaces
//   trees           at the depths distillation actually needs, which
//                   `distill_tree` measured rather than assumed
//   forward+softmax what the dispatcher does today
//   forward+argmax  the same, minus the softmax
//
// The tree depths are not a free choice. An earlier version of this file
// timed depth 4, called it 1.6x an `if`, and treated the question as
// settled -- before anything had checked what a depth-4 tree costs in
// accuracy. It costs 0.511 on the worst domain, which is not a tree
// anybody would ship. The depths worth timing are 8 and 12.
//
// The softmax row is here because it should not need to exist. argmax over
// a softmax is argmax over the logits -- the exponentials cannot reorder
// anything, since exp is monotonic. Every one of them at inference is
// bought and thrown away, and exp is expensive enough that the row is
// worth seeing rather than reasoning about.
//
// Build: part of the rsc tools target.

#include <dispatcher.hpp>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <print>
#include <random>
#include <vector>

namespace {

constexpr std::size_t kInputDim = 8;
constexpr std::size_t kNumOps   = 30;
constexpr std::size_t kReps     = 200000;

// Stand-in for the branch a learned dispatcher displaces: a handful of
// comparisons on the same features. Deliberately not trivial -- a single
// `if` would flatter the comparison.
std::size_t if_cascade(const std::vector<double>& x) {
    if (x[0] < -0.5) return (x[1] < 0.0) ? 0u : 1u;
    if (x[0] < 0.5)  return (x[2] < 0.0) ? 2u : 3u;
    if (x[3] > 1.0)  return (x[4] < 0.5) ? 4u : 5u;
    return (x[5] < 0.0) ? 6u : 7u;
}

// Trees at the depths the distillation actually needs, not at a depth
// chosen for looking good. `distill_tree` measured the loss against the
// network it imitates: depth 4 gives up 0.511 accuracy on the worst
// domain and is unusable, depth 8 gives up 0.092, depth 12 gives up
// 0.041. So depth 4 is here only as the number that was quoted before it
// was checked, and the honest rows are the deep ones.
//
// The comparisons are cheap; what makes a deep tree cost more than its
// depth suggests is that the node array stops fitting in L1 and the
// branches stop being predictable.
struct Tree {
    std::vector<std::size_t> feature;
    std::vector<double> threshold;
    std::vector<std::size_t> label;
    std::size_t depth = 0;
};

Tree make_tree(std::size_t depth, std::size_t num_ops, std::size_t dim, std::mt19937_64& rng) {
    const std::size_t internal = (std::size_t{1} << depth) - 1;
    Tree t;
    t.depth = depth;
    t.feature.resize(internal);
    t.threshold.resize(internal);
    t.label.resize(std::size_t{1} << depth);
    std::uniform_int_distribution<std::size_t> fd(0, dim - 1);
    std::normal_distribution<double> td(0.0, 1.0);
    for (std::size_t i = 0; i < internal; ++i) { t.feature[i] = fd(rng); t.threshold[i] = td(rng); }
    for (std::size_t i = 0; i < t.label.size(); ++i) t.label[i] = i % num_ops;
    return t;
}

std::size_t tree_eval(const Tree& t, const std::vector<double>& x) {
    // Heap layout: the children of node i are 2i+1 and 2i+2, so after
    // `depth` descents the index lands in the leaf band starting at
    // 2^depth - 1 and the offset into `label` is that band's base.
    std::size_t node = 0;
    for (std::size_t d = 0; d < t.depth; ++d)
        node = 2 * node + 1 + (x[t.feature[node]] > t.threshold[node] ? 1u : 0u);
    return t.label[node - ((std::size_t{1} << t.depth) - 1)];
}

template<class F>
double time_ns(F&& f, std::size_t reps) {
    const auto t0 = std::chrono::steady_clock::now();
    volatile std::size_t sink = 0;
    for (std::size_t i = 0; i < reps; ++i) sink = f(i);
    (void)sink;
    const auto t1 = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::nano>(t1 - t0).count() / static_cast<double>(reps);
}

}  // namespace

int main() {
    std::mt19937_64 rng(42);
    std::normal_distribution<double> nd(0.0, 1.0);

    // A pool of inputs rather than one, so the measurement is not of a
    // branch predictor that has memorised a single path.
    constexpr std::size_t kPool = 1024;
    std::vector<std::vector<double>> inputs(kPool, std::vector<double>(kInputDim));
    for (auto& v : inputs) for (auto& e : v) e = nd(rng);


    std::println("Learned dispatch against the branch it replaces");
    std::println("input dim {}, {} ops, {} reps, inputs drawn from a pool of {}",
                 kInputDim, kNumOps, kReps, kPool);
    std::println("");
    std::println("  {:<24} | {:>12} | {:>14}", "what", "ns / call", "vs if-cascade");

    const double t_if = time_ns([&](std::size_t i) {
        return if_cascade(inputs[i % kPool]);
    }, kReps);

    std::println("  {:<24} | {:>12.2f} | {:>13.1f}x", "if-cascade", t_if, 1.0);

    for (std::size_t depth : {4u, 8u, 10u, 12u}) {
        const Tree tree = make_tree(depth, kNumOps, kInputDim, rng);
        const double t_tree = time_ns([&](std::size_t i) {
            return tree_eval(tree, inputs[i % kPool]);
        }, kReps);
        std::println("  {:<24} | {:>12.2f} | {:>13.1f}x",
                     std::format("tree, depth {}", depth), t_tree, t_tree / t_if);
    }

    for (std::size_t hidden : {8u, 16u, 32u, 64u, 128u}) {
        rsc::Dispatcher model(kInputDim, hidden, kNumOps, /*seed=*/7);

        const double t_soft = time_ns([&](std::size_t i) {
            auto c = model.forward_cached(inputs[i % kPool]);
            auto p = rsc::softmax(c.logits);
            return static_cast<std::size_t>(std::max_element(p.begin(), p.end()) - p.begin());
        }, kReps);

        const double t_arg = time_ns([&](std::size_t i) {
            auto c = model.forward_cached(inputs[i % kPool]);
            return static_cast<std::size_t>(
                std::max_element(c.logits.begin(), c.logits.end()) - c.logits.begin());
        }, kReps);

        std::println("  {:<24} | {:>12.2f} | {:>13.1f}x",
                     std::format("hidden {}, +softmax", hidden), t_soft, t_soft / t_if);
        std::println("  {:<24} | {:>12.2f} | {:>13.1f}x",
                     std::format("hidden {}, argmax only", hidden), t_arg, t_arg / t_if);
    }

    std::println("");
    std::println("The softmax rows are pure waste at inference: exp is monotonic, so");
    std::println("argmax over the probabilities is argmax over the logits.");
    return 0;
}
