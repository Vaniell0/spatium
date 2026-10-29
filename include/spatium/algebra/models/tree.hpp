#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <array>
#  include <cassert>
#  include <cmath>
#  include <cstddef>
#  include <cstdint>
#  include <span>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::inline algebra {

// ── A tree that scores a row ──────────────────────────────────────
//
// A binary tree over rows of K numbers: each inner node asks whether one column is at most a
// threshold and sends the row left if it is, right if not; each leaf holds a score. Used as the
// score of a softmax (`softmax` below) it is a policy, and distilled from another policy it is
// the second rung of the ladder of hypotheses -- simpler than a net, readable, and lowered to
// the field IR with nothing but a comparison and a product per node.
//
// The nodes are kept in one vector, the root first; how a tree is grown from a sample is
// statistics/tree_fit.hpp. `bits` is a description length: per inner node the column (log2 K
// bits, rounded up), the threshold and the two children's places are not stored, the shape
// costing one bit a node; per leaf its score. With `scalar_bits` for each number stored.

template <std::size_t K, Scalar T = double> class ScoreTree {
  public:
    using Row = std::array<T, K>;

    struct Node {
        std::int32_t column = -1; // -1 for a leaf
        T threshold{};
        std::int32_t left = -1, right = -1;
        T value{}; // the score of a leaf
    };

    ScoreTree() : nodes_{Node{}} {}
    explicit ScoreTree(std::vector<Node> nodes) : nodes_(std::move(nodes)) { assert(!nodes_.empty()); }

    T score(const Row& x) const {
        std::size_t at = 0;
        while (nodes_[at].column >= 0) {
            const Node& n = nodes_[at];
            at = static_cast<std::size_t>(x[static_cast<std::size_t>(n.column)] <= n.threshold ? n.left : n.right);
        }
        return nodes_[at].value;
    }

    std::size_t leaves() const {
        std::size_t n = 0;
        for (const Node& node : nodes_)
            n += node.column < 0;
        return n;
    }
    std::size_t depth() const { return depth_from(0); }

    double bits(double scalar_bits = 32) const {
        double column_bits = 0;
        for (std::size_t k = K - 1; k > 0; k >>= 1)
            ++column_bits;
        double total = 0;
        for (const Node& n : nodes_)
            total += 1 + (n.column >= 0 ? column_bits + scalar_bits : scalar_bits);
        return total;
    }

    const std::vector<Node>& nodes() const { return nodes_; }

  private:
    std::size_t depth_from(std::size_t at) const {
        const Node& n = nodes_[at];
        if (n.column < 0)
            return 0;
        const std::size_t l = depth_from(static_cast<std::size_t>(n.left));
        const std::size_t r = depth_from(static_cast<std::size_t>(n.right));
        return 1 + (l > r ? l : r);
    }

    std::vector<Node> nodes_;
};

// The softmax of scores at a temperature, into `out`: exp(s / T) normalised, the largest
// taken out first so that nothing overflows. The same distribution LinearSoftmax gives from
// its own scores.
template <Scalar T> void softmax(std::span<const T> scores, std::span<T> out, T temperature = T{1}) {
    using std::exp;
    assert(out.size() == scores.size());
    T top = scores[0] / temperature;
    for (const T& s : scores)
        if (s / temperature > top)
            top = s / temperature;
    T total{0};
    for (std::size_t i = 0; i < scores.size(); ++i)
        total = total + (out[i] = exp(scores[i] / temperature - top));
    for (T& p : out)
        p = p / total;
}

} // namespace spatium::algebra
