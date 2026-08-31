#pragma once

// Step 3 of the RSC build plan: a minimal feedforward dispatcher. No
// recurrence (per "no learned hidden state by default" in the design doc),
// and small enough that step 4 can hand-derive its backward pass instead of
// needing a general reverse-mode autodiff engine.
//
// Input is deliberately [operands, observed output], not operands alone:
// for same-arity Tier-1 ops (add vs. multiply) there is no signal in the
// operands alone to tell them apart -- TaskGenerator assigned the correct
// op independently at random. Given the observed output too, "which op
// reproduces this input -> output pair" is a genuine, learnable pattern
// (target == a+b vs. target == a*b), not a coin flip. This is a mechanism
// sanity check -- does the forward/backward/training loop work at all --
// not a preview of Tier-2's actual dispatch signal (there, e.g. domain 1's
// dispatch reads a property of the input alone, like a discriminant's
// magnitude, with no "observed output" involved).
//
// Weights are flat row-major vectors, not a Spatium Matrix<T,R,C>: rows/
// cols here are runtime values (the registry's size, not a compile-time
// constant), and no dynamic linear-algebra type exists in Spatium yet.
// Hand-rolling one tiny dense matmul is small enough not to need one.

#include <registry.hpp>
#include <task.hpp>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

namespace rsc {

// Pads a task's operands and expected output into one fixed-size vector:
// [operands (max_in, zero-padded), expected_output (max_out, zero-padded)].
// A task whose expected_output is left empty (precision_task.hpp's domain,
// which has no op-agnostic "observed output" to compare against) just
// leaves those slots at zero -- capping each copy at max_in/max_out keeps
// that safe instead of relying on every caller getting the sizing right.
inline std::vector<double> features(const Task& task, std::size_t max_in, std::size_t max_out) {
    std::vector<double> x(max_in + max_out, 0.0);
    std::size_t in_n = std::min(task.inputs.size(), max_in);
    std::size_t out_n = std::min(task.expected_output.size(), max_out);
    for (std::size_t i = 0; i < in_n; ++i) x[i] = task.inputs[i];
    for (std::size_t i = 0; i < out_n; ++i) x[max_in + i] = task.expected_output[i];
    return x;
}

inline std::vector<double> softmax(const std::vector<double>& logits) {
    double m = *std::max_element(logits.begin(), logits.end());
    std::vector<double> out(logits.size());
    double sum = 0.0;
    for (std::size_t i = 0; i < logits.size(); ++i) {
        out[i] = std::exp(logits[i] - m);
        sum += out[i];
    }
    for (auto& v : out) v /= sum;
    return out;
}

// logits[target] singled out, rest lumped via log-sum-exp -- the standard
// categorical cross-entropy loss. Reused two ways: as the surrogate loss
// step 4's gradient check validates backward() against, and, unmodified, as
// step 5's REINFORCE loss when "target" is the sampled action and the whole
// loss gets scaled by the reward (reward=1 recovers ordinary supervised
// cross-entropy, which is exactly what a numerical gradient check needs).
inline double cross_entropy_loss(const std::vector<double>& logits, std::size_t target) {
    double m = *std::max_element(logits.begin(), logits.end());
    double sum = 0.0;
    for (double v : logits) sum += std::exp(v - m);
    return -(logits[target] - m) + std::log(sum);
}

// d(cross_entropy_loss)/d(logits) = softmax(logits) - one_hot(target).
inline std::vector<double> cross_entropy_dlogits(const std::vector<double>& logits,
                                                  std::size_t target) {
    auto d = softmax(logits);
    d[target] -= 1.0;
    return d;
}

struct ForwardCache {
    std::vector<double> x;       // input -- needed to compute dL/dw1
    std::vector<double> h;       // post-ReLU hidden -- needed for dL/dw2 and the ReLU mask
    std::vector<double> logits;
};

struct Gradients {
    std::vector<double> dw1, db1, dw2, db2;
};

// One hidden layer, ReLU, linear output -- logits over the registry's ops.
class Dispatcher {
public:
    Dispatcher(std::size_t input_dim, std::size_t hidden_dim, std::size_t num_ops,
               std::uint64_t seed = 0)
        : input_dim_(input_dim), hidden_dim_(hidden_dim), num_ops_(num_ops),
          w1_(hidden_dim * input_dim), b1_(hidden_dim, 0.0),
          w2_(num_ops * hidden_dim), b2_(num_ops, 0.0) {
        std::mt19937_64 rng(seed);
        std::normal_distribution<double> d1(0.0, std::sqrt(2.0 / static_cast<double>(input_dim)));
        for (auto& w : w1_) w = d1(rng);
        std::normal_distribution<double> d2(0.0, std::sqrt(2.0 / static_cast<double>(hidden_dim)));
        for (auto& w : w2_) w = d2(rng);
    }

    std::size_t input_dim() const { return input_dim_; }
    std::size_t hidden_dim() const { return hidden_dim_; }
    std::size_t num_ops() const { return num_ops_; }

    // Direct weight access -- needed by the numerical gradient check
    // (perturb one weight, re-run forward, compare to backward()'s
    // analytic gradient) and, later, by step 5's weight update.
    std::vector<double>& w1() { return w1_; }
    std::vector<double>& b1() { return b1_; }
    std::vector<double>& w2() { return w2_; }
    std::vector<double>& b2() { return b2_; }

    // Logits over num_ops(); pass through softmax() for a distribution, or
    // take the argmax directly for a hard dispatch choice.
    std::vector<double> forward(const std::vector<double>& x) const {
        return forward_cached(x).logits;
    }

    // Same forward pass, but keeps the intermediates backward() needs.
    ForwardCache forward_cached(const std::vector<double>& x) const {
        ForwardCache c;
        c.x = x;
        c.h.assign(hidden_dim_, 0.0);
        for (std::size_t j = 0; j < hidden_dim_; ++j) {
            double s = b1_[j];
            for (std::size_t i = 0; i < input_dim_; ++i) s += w1_[j * input_dim_ + i] * x[i];
            c.h[j] = std::max(0.0, s); // ReLU
        }
        c.logits.assign(num_ops_, 0.0);
        for (std::size_t k = 0; k < num_ops_; ++k) {
            double s = b2_[k];
            for (std::size_t j = 0; j < hidden_dim_; ++j) s += w2_[k * hidden_dim_ + j] * c.h[j];
            c.logits[k] = s;
        }
        return c;
    }

    // Closed-form backprop for this exact net -- matmul/ReLU/matmul, in
    // reverse. dlogits is dL/dlogits from whatever loss sits on top (e.g.
    // cross_entropy_dlogits() above, or REINFORCE's reward-scaled version
    // in step 5). Validated against finite differences, not asserted --
    // see tests/test_rsc_dispatcher.cpp.
    Gradients backward(const ForwardCache& cache, const std::vector<double>& dlogits) const {
        Gradients g;
        g.dw2.assign(num_ops_ * hidden_dim_, 0.0);
        g.db2.assign(num_ops_, 0.0);
        g.dw1.assign(hidden_dim_ * input_dim_, 0.0);
        g.db1.assign(hidden_dim_, 0.0);

        std::vector<double> dh(hidden_dim_, 0.0);
        for (std::size_t k = 0; k < num_ops_; ++k) {
            g.db2[k] = dlogits[k];
            for (std::size_t j = 0; j < hidden_dim_; ++j) {
                g.dw2[k * hidden_dim_ + j] = dlogits[k] * cache.h[j];
                dh[j] += dlogits[k] * w2_[k * hidden_dim_ + j];
            }
        }
        for (std::size_t j = 0; j < hidden_dim_; ++j) {
            double dz1 = (cache.h[j] > 0.0) ? dh[j] : 0.0; // ReLU gradient mask
            g.db1[j] = dz1;
            for (std::size_t i = 0; i < input_dim_; ++i) g.dw1[j * input_dim_ + i] = dz1 * cache.x[i];
        }
        return g;
    }

    // Plain SGD: theta -= lr * dtheta. Nothing here assumes it specifically
    // -- swap for Adam later if training turns out to need it, nothing else
    // changes.
    void apply_gradients(const Gradients& g, double lr) {
        for (std::size_t i = 0; i < w1_.size(); ++i) w1_[i] -= lr * g.dw1[i];
        for (std::size_t i = 0; i < b1_.size(); ++i) b1_[i] -= lr * g.db1[i];
        for (std::size_t i = 0; i < w2_.size(); ++i) w2_[i] -= lr * g.dw2[i];
        for (std::size_t i = 0; i < b2_.size(); ++i) b2_[i] -= lr * g.db2[i];
    }

private:
    std::size_t input_dim_, hidden_dim_, num_ops_;
    std::vector<double> w1_, b1_, w2_, b2_;
};

} // namespace rsc
