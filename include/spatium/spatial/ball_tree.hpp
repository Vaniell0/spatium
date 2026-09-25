#pragma once
// A tree of geodesic balls over a manifold: the acceleration structure for
// a space with no coordinates a box could be drawn in. Items are balls in
// the space's own metric (markers, bodies, the bound of anything); each
// node is a ball containing its children's, merged along geodesics
// (geodesic_ball.hpp). A geodesic ray walks it with the closed-form
// interval of that header, a point with the lower distance, so on the
// sphere and in hyperbolic space the walk is exact and needs no marching.
//
// The build splits by two pivots, using nothing but `distance`: the item
// farthest from the first, then the one farthest from that, and each item
// goes to the nearer pivot -- the metric-tree split (Uhlmann's generalised
// hyperplane), which works in any metric space because it never asks for a
// coordinate. A learned or surface-area split is a later choice over the
// same tree.
//
// Node layout as in bvh.hpp: a leaf holds `count > 0` items from `first`
// in the order array; an internal node has its left child next to it and
// its right child at `first`.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/spatial/geodesic_ball.hpp>
#  include <algorithm>
#  include <array>
#  include <cstdint>
#  include <limits>
#  include <numeric>
#  include <optional>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::spatial {

template<typename Space>
class GeodesicBallTree {
public:
    using T = typename Space::ScalarType;
    using PointType = typename Space::PointType;
    using Ball = GeodesicBall<Space>;

    struct Hit {
        std::size_t index;   // into the items the tree was built from
        T t;                 // where the ray enters that item
    };
    struct Nearest {
        std::size_t index;
        T distance;          // to the item's ball, 0 inside it
    };
    struct Node {
        Ball bound;
        std::uint32_t first{};
        std::uint32_t count{};
    };

    static GeodesicBallTree build(const Space& space, std::vector<Ball> items) {
        GeodesicBallTree tree;
        tree.space_ = space;
        tree.items_ = std::move(items);
        if (tree.items_.empty()) return tree;
        std::vector<std::uint32_t> idx(tree.items_.size());
        std::iota(idx.begin(), idx.end(), 0u);
        tree.nodes_.reserve(2 * tree.items_.size());
        tree.build_range(idx, 0, idx.size(), 1);
        tree.order_ = std::move(idx);
        return tree;
    }

    // The first item a unit-speed geodesic from p along v enters, within
    // [0, t_max].
    std::optional<Hit> ray_cast(const PointType& p, const PointType& v,
                                T t_max = std::numeric_limits<T>::infinity()) const {
        if (nodes_.empty()) return std::nullopt;
        std::optional<Hit> best;
        T best_t = t_max;
        return with_stack([&](auto& stack) {
            std::uint32_t sp = 0;
            stack[sp++] = 0;
            while (sp) {
                const auto ni = stack[--sp];
                const auto& n = nodes_[ni];
                if (!ray_interval(space_, p, v, n.bound, best_t)) continue;
                if (n.count > 0) {
                    for (std::uint32_t k = 0; k < n.count; ++k) {
                        const auto item = order_[n.first + k];
                        const auto in = ray_interval(space_, p, v, items_[item], best_t);
                        if (in && in->first <= best_t) {
                            best_t = in->first;
                            best = Hit{item, in->first};
                        }
                    }
                    continue;
                }
                const auto left = ni + 1, right = n.first;
                const auto lt = ray_interval(space_, p, v, nodes_[left].bound, best_t);
                const auto rt = ray_interval(space_, p, v, nodes_[right].bound, best_t);
                if (lt && rt) {   // near one on top
                    if (lt->first > rt->first) { stack[sp++] = left; stack[sp++] = right; }
                    else { stack[sp++] = right; stack[sp++] = left; }
                } else if (lt) {
                    stack[sp++] = left;
                } else if (rt) {
                    stack[sp++] = right;
                }
            }
            return best;
        });
    }

    std::optional<Nearest> nearest(const PointType& q) const {
        if (nodes_.empty()) return std::nullopt;
        std::optional<Nearest> best;
        T best_d = std::numeric_limits<T>::infinity();
        return with_stack([&](auto& stack) {
            std::uint32_t sp = 0;
            stack[sp++] = 0;
            while (sp) {
                const auto ni = stack[--sp];
                const auto& n = nodes_[ni];
                if (lower_distance(space_, n.bound, q) >= best_d) continue;
                if (n.count > 0) {
                    for (std::uint32_t k = 0; k < n.count; ++k) {
                        const auto item = order_[n.first + k];
                        const T d = lower_distance(space_, items_[item], q);
                        if (d < best_d) { best_d = d; best = Nearest{item, d}; }
                    }
                    continue;
                }
                const auto left = ni + 1, right = n.first;
                const T dl = lower_distance(space_, nodes_[left].bound, q);
                const T dr = lower_distance(space_, nodes_[right].bound, q);
                if (dl < dr) {
                    if (dr < best_d) stack[sp++] = right;
                    if (dl < best_d) stack[sp++] = left;
                } else {
                    if (dl < best_d) stack[sp++] = left;
                    if (dr < best_d) stack[sp++] = right;
                }
            }
            return best;
        });
    }

    const std::vector<Ball>& items() const { return items_; }
    const std::vector<Node>& nodes() const { return nodes_; }
    std::size_t depth() const { return depth_; }

private:
    static constexpr std::size_t kLeaf = 4;

    template<typename Walk>
    decltype(auto) with_stack(Walk&& walk) const {
        if (depth_ + 2 <= 64) {
            std::array<std::uint32_t, 64> stack;
            return walk(stack);
        }
        std::vector<std::uint32_t> stack(depth_ + 2);
        return walk(stack);
    }

    std::uint32_t build_range(std::vector<std::uint32_t>& idx, std::size_t begin, std::size_t end,
                              std::size_t depth) {
        depth_ = std::max(depth_, depth);
        const auto node = static_cast<std::uint32_t>(nodes_.size());
        nodes_.push_back({});
        Ball bound = items_[idx[begin]];
        for (std::size_t i = begin + 1; i < end; ++i) bound = merge(space_, bound, items_[idx[i]]);

        if (end - begin <= kLeaf) {
            nodes_[node] = {bound, static_cast<std::uint32_t>(begin), static_cast<std::uint32_t>(end - begin)};
            return node;
        }
        // Two pivots by distance alone, then each item to the nearer.
        auto farthest_from = [&](const PointType& from) {
            std::size_t best = begin;
            T best_d = T{-1};
            for (std::size_t i = begin; i < end; ++i) {
                const T d = space_.distance(from, items_[idx[i]].c);
                if (d > best_d) { best_d = d; best = i; }
            }
            return items_[idx[best]].c;
        };
        const PointType a = farthest_from(items_[idx[begin]].c);
        const PointType b = farthest_from(a);
        auto mid = std::partition(idx.begin() + static_cast<std::ptrdiff_t>(begin),
                                  idx.begin() + static_cast<std::ptrdiff_t>(end), [&](std::uint32_t i) {
                                      return space_.distance(items_[i].c, a) < space_.distance(items_[i].c, b);
                                  });
        auto split = static_cast<std::size_t>(mid - idx.begin());
        if (split == begin || split == end) split = begin + (end - begin) / 2;   // coincident centres

        build_range(idx, begin, split, depth + 1);
        const auto right = build_range(idx, split, end, depth + 1);
        nodes_[node] = {bound, right, 0};
        return node;
    }

    Space space_{};
    std::vector<Ball> items_;
    std::vector<Node> nodes_;
    std::vector<std::uint32_t> order_;
    std::size_t depth_ = 0;
};

} // namespace spatium::spatial
