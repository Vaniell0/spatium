#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <array>
#  include <cassert>
#  include <cmath>
#  include <cstddef>
#  include <span>
#endif

SPATIUM_EXPORT namespace spatium::inline algebra {

// ── A softmax over a linear score ─────────────────────────────────
//
// A list of candidates, each a row of K numbers, and one weight per column: the score of a row
// is weight . row, and the chance of a candidate is exp(score / temperature) normalised over
// the list. The model is the function and its gradient, nothing more; how the weights move is
// the step rule's (algebra/analysis/optimize.hpp), and what a row is made of is the caller's
// (a move's, in discrete/move_features.hpp).
//
// Nothing here allocates. `distribution` writes into a buffer the caller gives; `sample`,
// `cross_entropy` and `gradient` walk the list twice, once for the largest score and the sum,
// once for each candidate's chance, so they need no buffer at all.
//
// Templated on the scalar like the rest of algebra/, so the gradient written here can be held
// against `Dual` carried through `cross_entropy`. The highest score is `choose`, which is what
// `Chooser` asks.

template <std::size_t K, Scalar T = double> class LinearSoftmax {
  public:
    using Row = std::array<T, K>;
    using Weights = std::array<T, K>;

    LinearSoftmax() = default;
    explicit LinearSoftmax(const Weights& w) : weights_(w) {}

    T score(const Row& x) const {
        T s{0};
        for (std::size_t i = 0; i < K; ++i)
            s = s + weights_[i] * x[i];
        return s;
    }

    // The chance of each candidate, written into `out` (one per candidate). At temperature 1 it
    // is the learned distribution; a large one is uniform, a small one the best candidate.
    void distribution(std::span<const Row> list, std::span<T> out, T temperature = T{1}) const {
        assert(out.size() == list.size());
        const Normaliser n = normaliser(list, temperature);
        for (std::size_t i = 0; i < list.size(); ++i)
            out[i] = n.chance(score(list[i]));
    }

    // The candidate a uniform draw `u` in [0, 1) lands on in the distribution.
    std::size_t sample(std::span<const Row> list, T u, T temperature = T{1}) const {
        const Normaliser n = normaliser(list, temperature);
        T acc{0};
        for (std::size_t i = 0; i < list.size(); ++i)
            if (u < (acc = acc + n.chance(score(list[i]))))
                return i;
        return list.size() - 1;
    }

    // The highest score, the first of equals.
    std::size_t choose(std::span<const Row> list) const {
        std::size_t best = 0;
        T top = score(list[0]);
        for (std::size_t i = 1; i < list.size(); ++i)
            if (const T s = score(list[i]); s > top)
                top = s, best = i;
        return best;
    }

    // -sum_i target[i] log p[i]: how surprised the model is by the target distribution.
    T cross_entropy(std::span<const Row> list, std::span<const T> target) const {
        const Normaliser n = normaliser(list, T{1});
        T loss{0};
        for (std::size_t i = 0; i < list.size(); ++i)
            if (target[i] != T{0})
                loss = loss - target[i] * n.log_chance(score(list[i]));
        return loss;
    }
    // The same for one candidate chosen: -log p[chosen].
    T cross_entropy(std::span<const Row> list, std::size_t chosen) const {
        return -normaliser(list, T{1}).log_chance(score(list[chosen]));
    }

    // The gradient of the cross-entropy in the weights: sum_i (p[i] - target[i]) row[i].
    Weights gradient(std::span<const Row> list, std::span<const T> target) const {
        const Normaliser n = normaliser(list, T{1});
        Weights g{};
        for (std::size_t i = 0; i < list.size(); ++i) {
            const T pull = n.chance(score(list[i])) - target[i];
            for (std::size_t j = 0; j < K; ++j)
                g[j] = g[j] + pull * list[i][j];
        }
        return g;
    }
    // The same for one candidate chosen.
    Weights gradient(std::span<const Row> list, std::size_t chosen) const {
        const Normaliser n = normaliser(list, T{1});
        Weights g{};
        for (std::size_t i = 0; i < list.size(); ++i) {
            const T pull = n.chance(score(list[i])) - (i == chosen ? T{1} : T{0});
            for (std::size_t j = 0; j < K; ++j)
                g[j] = g[j] + pull * list[i][j];
        }
        return g;
    }

    const Weights& weights() const { return weights_; }
    Weights& weights() { return weights_; }

  private:
    // exp(s / temperature - top) / total, with top the largest scaled score: no exponent overflows.
    struct Normaliser {
        T temperature, top, log_total;
        T log_chance(T s) const { return s / temperature - top - log_total; }
        T chance(T s) const {
            using std::exp;
            return exp(log_chance(s));
        }
    };

    Normaliser normaliser(std::span<const Row> list, T temperature) const {
        using std::exp;
        using std::log;
        T top = score(list[0]) / temperature;
        for (std::size_t i = 1; i < list.size(); ++i)
            if (const T s = score(list[i]) / temperature; s > top)
                top = s;
        T total{0};
        for (const Row& x : list)
            total = total + exp(score(x) / temperature - top);
        return {temperature, top, log(total)};
    }

    Weights weights_{};
};

static_assert(Chooser<LinearSoftmax<4>, std::span<const std::array<double, 4>>>);

} // namespace spatium::algebra
