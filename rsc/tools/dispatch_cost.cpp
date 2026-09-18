// What does a learned dispatch actually cost, next to the branch it
// replaces?
//
// The question this answers is not "is the model accurate" but "can it sit
// where an `if` sits". Those are different questions and only the second
// decides whether a trained dispatcher can live in a hot path rather than
// being consulted once per scene.
//
// Four things are timed on the same inputs:
//
//   if-cascade      the hand-written branch the model replaces
//   forward+softmax what the dispatcher does today
//   forward+argmax  the same, minus the softmax
//   depth-4 tree    a distilled decision tree, the shape the literature
//                   reaches for when a policy has to answer in microseconds
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

// A depth-4 tree over the same features, standing in for a distilled
// policy. Sixteen leaves, each naming an op -- the shape the distillation
// literature produces, and the reason it is worth measuring is that its
// cost does not depend on the hidden layer it was distilled from.
std::size_t tree_depth4(const std::vector<double>& x, const std::vector<std::size_t>& leaf) {
    std::size_t i = 0;
    i = 2 * i + (x[0] > 0.0);
    i = 2 * i + (x[1] > 0.0);
    i = 2 * i + (x[2] > 0.0);
    i = 2 * i + (x[3] > 0.0);
    return leaf[i % leaf.size()];
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

    std::vector<std::size_t> leaf(16);
    for (std::size_t i = 0; i < leaf.size(); ++i) leaf[i] = i % kNumOps;

    std::println("Learned dispatch against the branch it replaces");
    std::println("input dim {}, {} ops, {} reps, inputs drawn from a pool of {}",
                 kInputDim, kNumOps, kReps, kPool);
    std::println("");
    std::println("  {:<24} | {:>12} | {:>14}", "what", "ns / call", "vs if-cascade");

    const double t_if = time_ns([&](std::size_t i) {
        return if_cascade(inputs[i % kPool]);
    }, kReps);

    const double t_tree = time_ns([&](std::size_t i) {
        return tree_depth4(inputs[i % kPool], leaf);
    }, kReps);

    std::println("  {:<24} | {:>12.2f} | {:>13.1f}x", "if-cascade", t_if, 1.0);
    std::println("  {:<24} | {:>12.2f} | {:>13.1f}x", "distilled tree, depth 4",
                 t_tree, t_tree / t_if);

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
