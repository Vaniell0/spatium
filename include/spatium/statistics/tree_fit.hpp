#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/models/tree.hpp>
#  include <algorithm>
#  include <cstddef>
#  include <cstdint>
#  include <numeric>
#  include <span>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::statistics {

// ── Growing a tree from a sample ──────────────────────────────────
//
// A regression tree fitted to (row, target) pairs by the greedy rule of CART (Breiman et al.
// 1984): at each node, over every column and every threshold between two neighbouring values,
// the split that most reduces the summed squared error, until the depth is reached, a node
// holds fewer than twice `min_leaf` rows, or no split reduces the error. A leaf's score is the
// mean target of its rows. Distilling a policy is this with the policy's scores as targets.
//
// The squared error on the sample never rises as the depth does: a deeper tree can always
// keep a shallower one's split and refine below it.

struct TreeOptions {
    std::size_t depth = 4;
    std::size_t min_leaf = 1;
};

template <std::size_t K, class T>
ScoreTree<K, T> fit_tree(std::span<const std::array<T, K>> rows, std::span<const T> targets,
                         TreeOptions options = {}) {
    using Tree = ScoreTree<K, T>;
    using Node = typename Tree::Node;
    std::vector<Node> nodes;
    std::vector<std::size_t> index(rows.size());
    std::iota(index.begin(), index.end(), std::size_t{0});

    const auto mean = [&](std::span<const std::size_t> at) {
        T s{0};
        for (std::size_t i : at)
            s = s + targets[i];
        return at.empty() ? T{0} : s / T(static_cast<double>(at.size()));
    };

    // Grows the node for rows[at...] at `depth`, returns its place.
    const auto grow = [&](auto&& self, std::span<std::size_t> at, std::size_t depth) -> std::int32_t {
        const auto place = static_cast<std::int32_t>(nodes.size());
        nodes.push_back(Node{-1, T{}, -1, -1, mean(at)});
        if (depth == options.depth || at.size() < 2 * options.min_leaf)
            return place;

        T total{0}, total_sq{0};
        for (std::size_t i : at)
            total = total + targets[i], total_sq = total_sq + targets[i] * targets[i];
        const T n = T(static_cast<double>(at.size()));
        const T whole = total_sq - total * total / n; // the node's summed squared error

        T best = whole;
        std::int32_t best_column = -1;
        T best_threshold{};
        std::vector<std::size_t> order(at.begin(), at.end());
        for (std::size_t c = 0; c < K; ++c) {
            std::ranges::sort(order, [&](std::size_t a, std::size_t b) { return rows[a][c] < rows[b][c]; });
            T left{0}, left_sq{0};
            for (std::size_t k = 0; k + 1 < order.size(); ++k) {
                const T y = targets[order[k]];
                left = left + y, left_sq = left_sq + y * y;
                const T here = rows[order[k]][c], next = rows[order[k + 1]][c];
                if (!(here < next) || k + 1 < options.min_leaf || order.size() - (k + 1) < options.min_leaf)
                    continue;
                const T nl = T(static_cast<double>(k + 1)), nr = n - nl;
                const T right = total - left, right_sq = total_sq - left_sq;
                const T error = (left_sq - left * left / nl) + (right_sq - right * right / nr);
                if (error < best)
                    best = error, best_column = static_cast<std::int32_t>(c), best_threshold = (here + next) / T{2};
            }
        }
        if (best_column < 0)
            return place;

        const auto c = static_cast<std::size_t>(best_column);
        const auto middle = std::partition(at.begin(), at.end(), [&](std::size_t i) { return rows[i][c] <= best_threshold; });
        const auto split = static_cast<std::size_t>(middle - at.begin());
        const std::int32_t l = self(self, at.subspan(0, split), depth + 1);
        const std::int32_t r = self(self, at.subspan(split), depth + 1);
        Node& node = nodes[static_cast<std::size_t>(place)];
        node.column = best_column;
        node.threshold = best_threshold;
        node.left = l;
        node.right = r;
        return place;
    };
    grow(grow, std::span<std::size_t>(index), 0);
    return Tree(std::move(nodes));
}

} // namespace spatium::statistics
