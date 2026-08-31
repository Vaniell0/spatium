#pragma once

// Step 5 of the RSC build plan: the policy-gradient (REINFORCE) training
// loop. Reward is 1 if the sampled op index matches the task's true op, 0
// otherwise -- no execution through the registry needed for this specific
// task shape, since op_index IS the thing under test (see task.hpp/
// dispatcher.hpp for why: same-arity Tier-1 ops need the observed output
// in the input for this to be learnable at all).
//
// REINFORCE with a running-average reward baseline: dL/dlogits =
// (reward - baseline) * cross_entropy_dlogits. Tried plain REINFORCE
// (advantage = raw reward, no baseline) first, exactly as the build plan
// says -- measured it empirically before deciding, not by assumption: eval
// accuracy rose to ~0.61 then regressed back toward baseline (~0.25) over
// further training, and higher learning rates didn't learn at all. Reward=0
// gives literally zero gradient with no baseline, so wrong guesses carry no
// corrective signal at all -- a known, real REINFORCE weakness, not a
// hypothetical one here. Subtracting a running-average reward gives wrong
// guesses a negative advantage (push away) and right guesses a positive one
// only when they beat the recent average (push toward) -- fixed the
// instability in the same experiment. baseline is caller-owned so it
// persists across multiple train() calls, not reset each time.

#include <dispatcher.hpp>
#include <task.hpp>
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

namespace rsc {

// One REINFORCE update from one already-sampled task. baseline is updated
// in place (exponential moving average of the reward) so the caller can
// keep training across multiple calls without losing it. Returns whether
// the sampled action matched the task's true op, for accuracy tracking.
// Takes a Task directly (not a generator) so chain.hpp can reuse this
// unchanged for a chain's evolving-state steps -- a chain step and an i.i.d.
// Tier-1 task are graded and trained identically, only how they're
// generated differs.
inline bool reinforce_update(Dispatcher& model, const Task& task, std::size_t max_in,
                              std::size_t max_out, double lr, std::mt19937_64& rng,
                              double& baseline, double baseline_lr = 0.01) {
    auto x = features(task, max_in, max_out);
    auto cache = model.forward_cached(x);
    auto probs = softmax(cache.logits);

    std::discrete_distribution<std::size_t> action_dist(probs.begin(), probs.end());
    std::size_t action = action_dist(rng);
    bool correct = (action == task.op_index);
    double reward = correct ? 1.0 : 0.0;

    double advantage = reward - baseline;
    baseline += baseline_lr * (reward - baseline);

    auto dlogits = cross_entropy_dlogits(cache.logits, action);
    for (auto& d : dlogits) d *= advantage;

    model.apply_gradients(model.backward(cache, dlogits), lr);
    return correct;
}

// Same REINFORCE gradient as reinforce_update(), but returns it instead of
// applying it immediately, and reports correctness via an out-param instead
// of the return value -- lets several samples get averaged into one
// lower-variance update before touching the weights at all. See
// train_batch()/chain.hpp's train_chains_batch() below for why this turned
// out to be necessary, not just nicer: single-sample REINFORCE (above)
// never learned chains' add-vs-multiply distinction at all (accuracy stuck
// at the ~0.50 chance level over 20k+ steps) -- reward=1 samples still
// carry a real but very high-variance gradient, and single-sample updates
// were too noisy to make progress on that specific, genuinely nonlinear
// decision. Averaging over a batch of ~32 samples fixed it completely
// (measured: chance-level -> ~0.95+ within a couple thousand steps).
// entropy_beta > 0 adds a policy-entropy bonus to the gradient (standard
// RL exploration regularizer): dLoss_entropy/dlogits[k] =
// beta * p_k * (log p_k + H), which pushes gradient descent toward a
// *less* peaked softmax whenever beta>0 (sign checked directly, not just
// asserted -- see rsc/README.md's base-model section for why this was
// added: a bigger hidden layer measurably hurt the hardest domain's
// accuracy in the unified base, consistent with faster policy-entropy
// collapse outrunning exploration on that domain specifically, not a
// capacity/architecture problem). Default 0.0 is a no-op, so every
// existing caller is unaffected.
inline Gradients reinforce_gradient(Dispatcher& model, const Task& task, std::size_t max_in,
                                     std::size_t max_out, std::mt19937_64& rng, double& baseline,
                                     bool& correct, double baseline_lr = 0.01,
                                     double entropy_beta = 0.0) {
    auto x = features(task, max_in, max_out);
    auto cache = model.forward_cached(x);
    auto probs = softmax(cache.logits);

    std::discrete_distribution<std::size_t> action_dist(probs.begin(), probs.end());
    std::size_t action = action_dist(rng);
    correct = (action == task.op_index);
    double reward = correct ? 1.0 : 0.0;

    double advantage = reward - baseline;
    baseline += baseline_lr * (reward - baseline);

    auto dlogits = cross_entropy_dlogits(cache.logits, action);
    for (auto& d : dlogits) d *= advantage;

    if (entropy_beta != 0.0) {
        double H = 0.0;
        for (double p : probs)
            if (p > 0.0) H -= p * std::log(p);
        for (std::size_t k = 0; k < probs.size(); ++k) {
            double lp = probs[k] > 0.0 ? std::log(probs[k]) : 0.0;
            dlogits[k] += entropy_beta * probs[k] * (lp + H);
        }
    }

    return model.backward(cache, dlogits);
}

inline void accumulate_gradients(Gradients& sum, const Gradients& g) {
    if (sum.dw1.empty()) { sum = g; return; }
    for (std::size_t i = 0; i < sum.dw1.size(); ++i) sum.dw1[i] += g.dw1[i];
    for (std::size_t i = 0; i < sum.db1.size(); ++i) sum.db1[i] += g.db1[i];
    for (std::size_t i = 0; i < sum.dw2.size(); ++i) sum.dw2[i] += g.dw2[i];
    for (std::size_t i = 0; i < sum.db2.size(); ++i) sum.db2[i] += g.db2[i];
}

inline void scale_gradients(Gradients& g, double s) {
    for (auto& v : g.dw1) v *= s;
    for (auto& v : g.db1) v *= s;
    for (auto& v : g.dw2) v *= s;
    for (auto& v : g.db2) v *= s;
}

inline bool train_step(Dispatcher& model, TaskGenerator& tasks, std::size_t max_in,
                        std::size_t max_out, double lr, std::mt19937_64& rng,
                        double& baseline, double baseline_lr = 0.01) {
    return reinforce_update(model, tasks.sample(), max_in, max_out, lr, rng, baseline, baseline_lr);
}

// Runs n_steps of train_step, returns the fraction of sampled actions that
// were correct over the run -- the accuracy metric step 6 checks rises.
// baseline is in/out so a caller can chain multiple train() calls.
inline double train(Dispatcher& model, TaskGenerator& tasks, std::size_t max_in,
                     std::size_t max_out, double lr, std::size_t n_steps, std::uint64_t seed,
                     double& baseline) {
    std::mt19937_64 rng(seed);
    std::size_t correct = 0;
    for (std::size_t i = 0; i < n_steps; ++i)
        if (train_step(model, tasks, max_in, max_out, lr, rng, baseline)) ++correct;
    return static_cast<double>(correct) / static_cast<double>(n_steps);
}

// Minibatched REINFORCE: average the gradient over batch_size independent
// samples, then one apply_gradients() call. See reinforce_gradient() above
// for why this exists.
inline double train_batch(Dispatcher& model, TaskGenerator& tasks, std::size_t max_in,
                           std::size_t max_out, double lr, std::size_t batch_size,
                           std::size_t n_updates, std::uint64_t seed, double& baseline) {
    std::mt19937_64 rng(seed);
    std::size_t correct = 0, total = 0;
    for (std::size_t u = 0; u < n_updates; ++u) {
        Gradients sum;
        for (std::size_t b = 0; b < batch_size; ++b) {
            bool ok;
            auto g = reinforce_gradient(model, tasks.sample(), max_in, max_out, rng, baseline, ok);
            accumulate_gradients(sum, g);
            if (ok) ++correct;
            ++total;
        }
        scale_gradients(sum, 1.0 / static_cast<double>(batch_size));
        model.apply_gradients(sum, lr);
    }
    return static_cast<double>(correct) / static_cast<double>(total);
}

// Argmax dispatch on one task, no weight update -- shared by evaluate()
// and chain.hpp's evaluation helper.
inline bool argmax_correct(const Dispatcher& model, const Task& task, std::size_t max_in,
                            std::size_t max_out) {
    auto x = features(task, max_in, max_out);
    auto logits = model.forward(x);
    std::size_t argmax = 0;
    for (std::size_t k = 1; k < logits.size(); ++k)
        if (logits[k] > logits[argmax]) argmax = k;
    return argmax == task.op_index;
}

// Held-out accuracy: same as the loop above, minus the weight update --
// for measuring where the model stands without training on the eval data.
inline double evaluate(const Dispatcher& model, TaskGenerator& tasks, std::size_t max_in,
                        std::size_t max_out, std::size_t n_samples) {
    std::size_t correct = 0;
    for (std::size_t i = 0; i < n_samples; ++i)
        if (argmax_correct(model, tasks.sample(), max_in, max_out)) ++correct;
    return static_cast<double>(correct) / static_cast<double>(n_samples);
}

} // namespace rsc
