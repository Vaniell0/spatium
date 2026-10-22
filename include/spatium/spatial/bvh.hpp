#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/error.hpp>
#  include <spatium/geometry/concepts.hpp>
#  include <spatium/geometry/box.hpp>
#  include <spatium/geometry/intersection.hpp>
#  include <spatium/geometry/ray_hit.hpp>
#  include <spatium/geometry/triangle.hpp>
#  include <spatium/spatial/bound.hpp>
#  include <algorithm>
#  include <array>
#  include <cstdint>
#  include <cmath>
#  include <limits>
#  include <numeric>
#  include <optional>
#  include <type_traits>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::spatial {

// `B` is what a node is bounded by -- see bound.hpp. The default is the
// axis-aligned box the tree has always used, and with it the tree builds
// and walks exactly as before; a Ball gives the same tree over spheres.
template<geometry::Bounded Shape,
         Bound B = geometry::Box<Shape::ambient_dimension, typename Shape::ScalarType>>
struct BVH {
    using T = typename Shape::ScalarType;
    static constexpr std::size_t N = Shape::ambient_dimension;
    using PointType = Vec<T, N>;
    using BoxType = geometry::Box<N, T>;
    using BoundType = B;

    // BVH itself requires `Bounded` (above) — the Shape must supply
    // `bounding_box()`. When the Shape ALSO satisfies `RayHittable`,
    // Hit carries the shape-specific barycentric weights and unit
    // normal harvested from `ray_hit`. The canonical such shape is
    // `Triangle<3, T>`. `BoundedQuadric` is the analytic counterpart:
    // RayHittable *and* Bounded, so a tree can hold exact surfaces
    // instead of the triangles they would otherwise be tessellated
    // into. Bare `Quadric` / `Torus` stay RayHittable but unbounded —
    // an infinite cylinder or cone has no bounding box to give — so
    // they remain direct-call shapes, as does any quadric that has
    // not been confined to a box.
    // Bounded shapes that lack a `ray_hit` overload fall through to
    // the legacy `intersect(ray, shape)` path; only `t` and `point`
    // are filled in that case.
    static constexpr bool has_ray_hit =
        geometry::RayHittable<Shape, typename Shape::ScalarType>;

    struct Hit {
        std::size_t index;
        T t;
        PointType point;
        T u{};              // barycentric: vertex 1 weight
        T v{};              // barycentric: vertex 2 weight (w = 1-u-v)
        PointType normal{}; // face normal (unit length) — 3D triangle only
    };

    struct NearestResult {
        std::size_t index;
        T distance;
        PointType point;
    };

    // ── Build ─────────────────────────────────────────────────

    static BVH build(std::vector<Shape> shapes) {
        BVH bvh;
        bvh.shapes_ = std::move(shapes);
        bvh.rebuild();
        return bvh;
    }

    // Build the tree again over the shapes it holds, from scratch. Shape
    // indices are untouched -- only the nodes change -- so this is safe at
    // any moment, and it is what the quality guard below calls.
    void rebuild() {
        nodes_.clear();
        prim_indices_.clear();
        free_nodes_.clear();
        leaf_of_.assign(shapes_.size(), npos);
        depth_ = 0;
        garbage_ = 0;
        cost_sum_ = 0;
        if (shapes_.empty()) {
            reference_cost_ = 0;
            reference_size_ = 0;
            return;
        }

        std::vector<B> boxes(shapes_.size());
        std::vector<PointType> centroids(shapes_.size());
        for (std::size_t i = 0; i < shapes_.size(); ++i) {
            boxes[i] = bound_of<B>(shapes_[i]);
            centroids[i] = center(boxes[i]);
        }

        std::vector<std::size_t> indices(shapes_.size());
        std::iota(indices.begin(), indices.end(), 0);

        nodes_.reserve(2 * shapes_.size());
        build_recursive(indices, 0, indices.size(), boxes, centroids, 1, npos);
        for (const auto& n : nodes_) cost_sum_ += contribution(n);
        reference_cost_ = cost();
        reference_size_ = shapes_.size();
    }

    // ── Update ────────────────────────────────────────────────
    //
    // A tree that follows a shape set which changes by a few at a time,
    // instead of being built again. Each operation uses only what the
    // bound already offers -- `merge`, `sah_measure`, `bound_of` -- so it
    // is the same for a box, a ball or any other Bound.
    //
    // Indices. `insert` appends and returns the new shape's index (the old
    // size). `remove(i)` keeps the shapes dense by moving the LAST shape
    // into slot i, and returns that shape's old index, so a caller holding
    // data parallel to the shapes does the same swap-remove on it
    // (`v[i] = std::move(v.back()); v.pop_back();`); the value is `i`
    // itself when i was the last shape and nothing moved. Every other
    // index is unchanged, and `rebuild()` leaves all of them alone.
    //
    // Layout. An updated tree is no longer in depth-first order: a node's
    // left child is `Node::left`, not `index + 1`, free slots may sit
    // between live nodes, and `prim_indices()` may hold entries no leaf's
    // range reaches. A walk from node 0 through `left` and `first` (the
    // right child) sees exactly the live tree. Code that copies the arrays
    // wholesale (the device copy of the triangle tree) wants a built tree,
    // which `rebuild()` gives back.

    std::size_t size() const { return shapes_.size(); }

    std::size_t insert(Shape shape) {
        const std::size_t idx = shapes_.size();
        shapes_.push_back(std::move(shape));
        leaf_of_.push_back(npos);
        const B b = bound_of<B>(shapes_[idx]);

        if (nodes_.empty()) {
            prim_indices_.clear();
            prim_indices_.push_back(idx);
            nodes_.push_back(Node{b, 0, 1, npos, npos});
            leaf_of_[idx] = 0;
            depth_ = 1;
            garbage_ = 0;
            cost_sum_ = contribution(nodes_[0]);
            reference_cost_ = cost();
            reference_size_ = 1;
            return idx;
        }

        // Down the cheaper side, where cheaper is the smaller growth of the
        // child's surface measure, the part of the SAH this shape adds.
        // Bounds are fixed on the way back up, from the children: for a
        // box merging the new shape into each ancestor would do, but for a
        // ball the merge of an ancestor with the shape need not contain the
        // merge of its child with it.
        std::uint32_t cur = 0;
        std::size_t level = 1;
        while (nodes_[cur].count == 0) {
            const std::uint32_t l = nodes_[cur].left, r = nodes_[cur].first;
            const double ml = sah_measure(nodes_[l].bounds), mr = sah_measure(nodes_[r].bounds);
            const double gl = double(sah_measure(merge(nodes_[l].bounds, b))) - ml;
            const double gr = double(sah_measure(merge(nodes_[r].bounds, b))) - mr;
            cur = (gl < gr || (gl == gr && ml <= mr)) ? l : r;
            ++level;
        }

        // A leaf's range cannot grow where it stands, so it moves to the
        // end of prim_indices_ with the new entry; the old range is left
        // behind and counted as garbage.
        retire(nodes_[cur]);
        const auto offset = static_cast<std::uint32_t>(prim_indices_.size());
        const std::uint32_t old_first = nodes_[cur].first, old_count = nodes_[cur].count;
        for (std::uint32_t i = 0; i < old_count; ++i)
            prim_indices_.push_back(prim_indices_[old_first + i]);
        prim_indices_.push_back(idx);
        garbage_ += old_count;
        nodes_[cur].first = offset;
        nodes_[cur].count = old_count + 1;
        nodes_[cur].bounds = merge(nodes_[cur].bounds, b);
        admit(nodes_[cur]);
        leaf_of_[idx] = cur;

        if (nodes_[cur].count > LEAF_THRESHOLD) split_leaf(cur, level);
        refit_up(nodes_[cur].parent);
        maintain();
        return idx;
    }

    std::size_t remove(std::size_t index) {
        const std::size_t last = shapes_.size() - 1;
        const std::uint32_t leaf_idx = leaf_of_[index];

        // Out of its leaf's range.
        {
            auto& leaf = nodes_[leaf_idx];
            retire(leaf);
            std::uint32_t pos = leaf.first;
            while (prim_indices_[pos] != index) ++pos;
            prim_indices_[pos] = prim_indices_[leaf.first + leaf.count - 1];
            ++garbage_;
            --leaf.count;
        }

        if (nodes_[leaf_idx].count > 0) {
            refit_leaf(leaf_idx);
            admit(nodes_[leaf_idx]);
            refit_up(nodes_[leaf_idx].parent);
        } else if (leaf_idx == 0) {
            nodes_.clear();
            prim_indices_.clear();
            free_nodes_.clear();
            cost_sum_ = 0;
            garbage_ = 0;
            depth_ = 0;
        } else {
            // The leaf is empty: its parent becomes its sibling.
            const std::uint32_t par = nodes_[leaf_idx].parent;
            const std::uint32_t sib = nodes_[par].left == leaf_idx ? nodes_[par].first
                                                                   : nodes_[par].left;
            retire(nodes_[par]);
            retire(nodes_[sib]);
            const std::uint32_t grand = nodes_[par].parent;
            nodes_[par] = nodes_[sib];
            nodes_[par].parent = grand;
            if (nodes_[par].count == 0) {
                nodes_[nodes_[par].left].parent = par;
                nodes_[nodes_[par].first].parent = par;
            } else {
                for (std::uint32_t i = 0; i < nodes_[par].count; ++i)
                    leaf_of_[prim_indices_[nodes_[par].first + i]] = par;
            }
            admit(nodes_[par]);
            free_nodes_.push_back(leaf_idx);
            free_nodes_.push_back(sib);
            refit_up(grand);
        }

        // The last shape takes the freed slot.
        if (index != last) {
            const std::uint32_t l = leaf_of_[last];
            std::uint32_t pos = nodes_[l].first;
            while (prim_indices_[pos] != last) ++pos;
            prim_indices_[pos] = index;
            leaf_of_[index] = l;
            shapes_[index] = std::move(shapes_[last]);
        }
        shapes_.pop_back();
        leaf_of_.pop_back();
        maintain();
        return last;
    }

    // ── Quality ───────────────────────────────────────────────
    //
    // The SAH cost of the live tree, the quantity `build` minimises:
    // (sum over internal nodes of their measure + sum over leaves of
    // count * measure) / measure of the root. The operations above keep it
    // up to date in O(depth); `measured_cost()` recomputes it from the
    // nodes. Updates pick the cheaper side but never reshape what is
    // already there, so it creeps up. Once it passes `degradation_limit()`
    // times what a fresh build had (scaled by log2 of the size, since a
    // tree over more shapes costs more) the tree is built again, with
    // every shape index kept.
    T cost() const {
        if (nodes_.empty()) return T{0};
        const double root = sah_measure(nodes_[0].bounds);
        return root > 0 ? static_cast<T>(cost_sum_ / root) : T{0};
    }

    T measured_cost() const {
        if (nodes_.empty()) return T{0};
        double sum = 0;
        walk_live([&](const Node& n) { sum += contribution(n); });
        const double root = sah_measure(nodes_[0].bounds);
        return root > 0 ? static_cast<T>(sum / root) : T{0};
    }

    double degradation_limit() const { return degradation_limit_; }
    void set_degradation_limit(double f) { degradation_limit_ = f; }
    std::size_t rebuilds() const { return rebuilds_; }

    // Whether the cost has passed the limit.
    bool degraded() const {
        if (shapes_.size() < 2 * LEAF_THRESHOLD || !(reference_cost_ > T{0})) return false;
        const double scale = (1.0 + std::log2(static_cast<double>(shapes_.size()))) /
                             (1.0 + std::log2(static_cast<double>(reference_size_)));
        return static_cast<double>(cost()) >
               degradation_limit_ * scale * static_cast<double>(reference_cost_);
    }

    // ── Queries ───────────────────────────────────────────────

    // `t_max` ends the ray: a hit beyond it is not a hit. A shadow ray
    // stops at the light, a swept query at the end of its step.
    std::optional<Hit> ray_cast(const geometry::Ray<N, T>& ray,
                                T t_max = std::numeric_limits<T>::max()) const {
        if (nodes_.empty()) return std::nullopt;

        std::optional<Hit> best;
        T best_t = t_max;

        return with_stack([&](auto& stack) {
            std::uint32_t sp = 0;
            stack[sp++] = 0;

            while (sp) {
                auto node_idx = stack[--sp];
                auto& node = nodes_[node_idx];

                if (!ray_interval(ray, node.bounds, best_t)) continue;

                if (node.count > 0) {
                    for (std::uint32_t i = 0; i < node.count; ++i) {
                        auto idx = prim_indices_[node.first + i];
                        if constexpr (has_ray_hit) {
                            // Unqualified, after pulling the library's own
                            // overloads into scope: this is the two-step that
                            // makes the extension point in ray_hit.hpp true.
                            // A qualified `geometry::ray_hit` finds only what
                            // lives in that namespace, so a user shape with
                            // its own overload beside it satisfied
                            // `RayHittable` -- which looks the name up
                            // unqualified and finds theirs by ADL -- while
                            // `BVH<TheirShape>` failed to compile. The concept
                            // said yes and the use said no.
                            using geometry::ray_hit;
                            auto h = ray_hit(ray, shapes_[idx]);
                            if (h && h->t >= T{0} && h->t <= best_t) {
                                best_t = h->t;
                                best = Hit{idx, h->t, h->point,
                                           h->u, h->v, h->normal};
                            }
                        } else {
                            auto result = geometry::intersect(ray, shapes_[idx]);
                            if (result) {
                                auto t = (result.value() - ray.origin).dot(ray.direction);
                                if (t >= T{0} && t <= best_t) {
                                    best_t = t;
                                    best = Hit{idx, t, result.value()};
                                }
                            }
                        }
                    }
                } else {
                    auto left = node.left;
                    auto right = node.first;
                    auto left_t = ray_interval(ray, nodes_[left].bounds, best_t);
                    auto right_t = ray_interval(ray, nodes_[right].bounds, best_t);

                    bool left_ok = left_t.has_value();
                    bool right_ok = right_t.has_value();

                    // Push far first, near second (near popped first)
                    if (left_ok && right_ok) {
                        if (left_t->first > right_t->first) {
                            stack[sp++] = left;
                            stack[sp++] = right;
                        } else {
                            stack[sp++] = right;
                            stack[sp++] = left;
                        }
                    } else if (left_ok) {
                        stack[sp++] = left;
                    } else if (right_ok) {
                        stack[sp++] = right;
                    }
                }
            }
            return best;
        });
    }

    // Any hit in [0, t_max]. Through the same `ray_hit` as ray_cast where
    // the shape has one -- it used to call `intersect()` alone, so a tree of
    // exact quadrics could say what a ray hits first and not whether it hits
    // anything.
    bool ray_test(const geometry::Ray<N, T>& ray, T t_max = std::numeric_limits<T>::max()) const {
        if (nodes_.empty()) return false;

        return with_stack([&](auto& stack) {
            std::uint32_t sp = 0;
            stack[sp++] = 0;

            while (sp) {
                auto node_idx = stack[--sp];
                auto& node = nodes_[node_idx];

                if (!ray_interval(ray, node.bounds, t_max)) continue;

                if (node.count > 0) {
                    for (std::uint32_t i = 0; i < node.count; ++i) {
                        auto idx = prim_indices_[node.first + i];
                        if constexpr (has_ray_hit) {
                            using geometry::ray_hit;
                            auto h = ray_hit(ray, shapes_[idx]);
                            if (h && h->t >= T{0} && h->t <= t_max) return true;
                        } else {
                            auto result = geometry::intersect(ray, shapes_[idx]);
                            if (result) {
                                auto t = (result.value() - ray.origin).dot(ray.direction);
                                if (t >= T{0} && t <= t_max) return true;
                            }
                        }
                    }
                } else {
                    stack[sp++] = node.left;
                    stack[sp++] = node.first;     // right
                }
            }
            return false;
        });
    }

    std::optional<NearestResult> nearest(const PointType& p) const {
        if (nodes_.empty()) return std::nullopt;

        std::optional<NearestResult> best;
        T best_dist = std::numeric_limits<T>::max();

        return with_stack([&](auto& stack) {
            std::uint32_t sp = 0;
            stack[sp++] = 0;

            while (sp) {
                auto node_idx = stack[--sp];
                auto& node = nodes_[node_idx];

                T box_dist = lower_distance(node.bounds, p);
                if (box_dist >= best_dist) continue;

                if (node.count > 0) {
                    for (std::uint32_t i = 0; i < node.count; ++i) {
                        auto idx = prim_indices_[node.first + i];
                        auto proj = shapes_[idx].project(p);
                        T d = (p - proj).norm();
                        if (d < best_dist) {
                            best_dist = d;
                            best = NearestResult{idx, d, proj};
                        }
                    }
                } else {
                    auto left = node.left;
                    auto right = node.first;
                    T dl = lower_distance(nodes_[left].bounds, p);
                    T dr = lower_distance(nodes_[right].bounds, p);

                    if (dl < dr) {
                        if (dr < best_dist) stack[sp++] = right;
                        if (dl < best_dist) stack[sp++] = left;
                    } else {
                        if (dl < best_dist) stack[sp++] = left;
                        if (dr < best_dist) stack[sp++] = right;
                    }
                }
            }
            return best;
        });
    }

    std::vector<std::size_t> query_box(const BoxType& query) const {
        std::vector<std::size_t> result;
        if (nodes_.empty()) return result;

        return with_stack([&](auto& stack) {
            std::uint32_t sp = 0;
            stack[sp++] = 0;

            while (sp) {
                auto node_idx = stack[--sp];
                auto& node = nodes_[node_idx];

                if (!overlaps(node.bounds, query)) continue;

                if (node.count > 0) {
                    for (std::uint32_t i = 0; i < node.count; ++i) {
                        auto idx = prim_indices_[node.first + i];
                        if (shapes_[idx].bounding_box().intersects(query))
                            result.push_back(idx);
                    }
                } else {
                    stack[sp++] = node.left;
                    stack[sp++] = node.first;
                }
            }
            return std::move(result);
        });
    }

    const std::vector<Shape>& shapes() const { return shapes_; }
    std::size_t node_count() const { return nodes_.size(); }
    std::size_t depth() const { return depth_; }

    // Node layout:
    //   leaf:     first = prim_indices_ offset, count > 0
    //   internal: first = RIGHT child index, count = 0, left = LEFT child
    //             (a built tree puts it at this_node + 1, an updated one anywhere)
    //
    // Public, with the arrays behind it, because the layout is already
    // what a device wants -- flat, index-linked, walked with a stack -- and
    // a GPU copy of the tree is these two arrays converted, not a second
    // build that could split differently.
    struct Node {
        B bounds;
        std::uint32_t first{};
        std::uint32_t count{};
        std::uint32_t left{};    // internal nodes only
        std::uint32_t parent{};  // ~0 for the root
    };
    const std::vector<Node>& nodes() const { return nodes_; }
    const std::vector<std::size_t>& prim_indices() const { return prim_indices_; }

private:
    static constexpr std::uint32_t LEAF_THRESHOLD = 4;
    static constexpr int NUM_BINS = 12;

    // The traversal stack, sized by the tree it walks. A pop pushes at most
    // two, so a walk never holds more than depth + 1 entries. It was a fixed
    // array of 64 with the comment that SAH trees stay shallow, and they do
    // not: centroids spread over many orders of magnitude make the binned
    // split peel a few off per level, and the test that builds such a tree
    // overran the array. Trees up to 64 deep -- nearly all of them -- still
    // walk on a fixed array on the call's stack, exactly as before; only a
    // deeper one allocates. Each query's walk is a generic lambda so both
    // paths are the same code, and the common one costs what it always did.
    template<typename Walk>
    decltype(auto) with_stack(Walk&& walk) const {
        if (depth_ + 2 <= 64) {
            std::array<std::uint32_t, 64> stack;
            return walk(stack);
        }
        std::vector<std::uint32_t> stack(depth_ + 2);
        return walk(stack);
    }

    std::vector<Shape> shapes_;
    std::vector<Node> nodes_;
    std::vector<std::size_t> prim_indices_;
    std::size_t depth_ = 0;

    static constexpr std::uint32_t npos = ~std::uint32_t{0};

    std::vector<std::uint32_t> leaf_of_;     // shape index -> the leaf holding it
    std::vector<std::uint32_t> free_nodes_;  // slots of nodes_ no live node uses
    std::size_t garbage_ = 0;                // entries of prim_indices_ no leaf reaches
    double cost_sum_ = 0;                    // the numerator of cost()
    T reference_cost_{};                     // cost() when last built
    std::size_t reference_size_ = 0;
    double degradation_limit_ = 1.25;
    std::size_t rebuilds_ = 0;

    double contribution(const Node& n) const {
        const double a = sah_measure(n.bounds);
        return n.count == 0 ? a : a * n.count;
    }
    void retire(const Node& n) { cost_sum_ -= contribution(n); }
    void admit(const Node& n) { cost_sum_ += contribution(n); }

    // A leaf's bound from the shapes in its range; the caller accounts for
    // the cost.
    void refit_leaf(std::uint32_t i) {
        auto& n = nodes_[i];
        B b = bound_of<B>(shapes_[prim_indices_[n.first]]);
        for (std::uint32_t k = 1; k < n.count; ++k)
            b = merge(b, bound_of<B>(shapes_[prim_indices_[n.first + k]]));
        n.bounds = b;
    }

    // From `i` to the root, each bound again from its two children.
    void refit_up(std::uint32_t i) {
        while (i != npos) {
            retire(nodes_[i]);
            nodes_[i].bounds = merge(nodes_[nodes_[i].left].bounds, nodes_[nodes_[i].first].bounds);
            admit(nodes_[i]);
            i = nodes_[i].parent;
        }
    }

    std::uint32_t alloc_node() {
        if (!free_nodes_.empty()) {
            auto i = free_nodes_.back();
            free_nodes_.pop_back();
            return i;
        }
        nodes_.push_back(Node{});
        return static_cast<std::uint32_t>(nodes_.size() - 1);
    }

    // A leaf of more than LEAF_THRESHOLD shapes becomes an internal node
    // over two: cut at the median of the centroids along the axis they
    // spread over most.
    void split_leaf(std::uint32_t p, std::size_t level) {
        const std::uint32_t first = nodes_[p].first, count = nodes_[p].count;
        std::vector<std::size_t> ids(prim_indices_.begin() + first, prim_indices_.begin() + first + count);
        std::vector<B> bs;
        std::vector<PointType> cs;
        for (auto id : ids) {
            bs.push_back(bound_of<B>(shapes_[id]));
            cs.push_back(center(bs.back()));
        }
        std::size_t axis = 0;
        T widest = T{-1};
        for (std::size_t d = 0; d < N; ++d) {
            T lo = cs[0][d], hi = cs[0][d];
            for (const auto& c : cs) { lo = std::min(lo, c[d]); hi = std::max(hi, c[d]); }
            if (hi - lo > widest) { widest = hi - lo; axis = d; }
        }
        std::vector<std::size_t> ord(count);
        std::iota(ord.begin(), ord.end(), 0);
        std::sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t b) {
            return cs[a][axis] < cs[b][axis];
        });

        retire(nodes_[p]);
        garbage_ += count;
        const std::uint32_t l = alloc_node(), r = alloc_node();
        const std::size_t half = count / 2;
        for (int side = 0; side < 2; ++side) {
            const std::uint32_t node = side ? r : l;
            const std::size_t from = side ? half : 0, to = side ? count : half;
            const auto offset = static_cast<std::uint32_t>(prim_indices_.size());
            B bound = bs[ord[from]];
            for (std::size_t k = from; k < to; ++k) {
                prim_indices_.push_back(ids[ord[k]]);
                leaf_of_[ids[ord[k]]] = node;
                if (k > from) bound = merge(bound, bs[ord[k]]);
            }
            nodes_[node] = Node{bound, offset, static_cast<std::uint32_t>(to - from), npos, p};
            admit(nodes_[node]);
        }
        nodes_[p].count = 0;
        nodes_[p].left = l;
        nodes_[p].first = r;
        nodes_[p].bounds = merge(nodes_[l].bounds, nodes_[r].bounds);
        admit(nodes_[p]);
        depth_ = std::max(depth_, level + 1);
    }

    template<typename F>
    void walk_live(F&& f) const {
        if (nodes_.empty()) return;
        std::vector<std::uint32_t> stack{0};
        while (!stack.empty()) {
            const auto i = stack.back();
            stack.pop_back();
            f(nodes_[i]);
            if (nodes_[i].count == 0) {
                stack.push_back(nodes_[i].left);
                stack.push_back(nodes_[i].first);
            }
        }
    }

    void compact_prims() {
        std::vector<std::size_t> packed;
        packed.reserve(shapes_.size());
        std::vector<std::uint32_t> stack{0};
        while (!stack.empty()) {
            const auto i = stack.back();
            stack.pop_back();
            auto& n = nodes_[i];
            if (n.count == 0) {
                stack.push_back(n.left);
                stack.push_back(n.first);
            } else {
                const auto offset = static_cast<std::uint32_t>(packed.size());
                for (std::uint32_t k = 0; k < n.count; ++k) packed.push_back(prim_indices_[n.first + k]);
                n.first = offset;
            }
        }
        prim_indices_ = std::move(packed);
        garbage_ = 0;
    }

    // After an update: drop the stale entries once they outnumber the live
    // ones, and build again if the cost has degraded.
    void maintain() {
        if (nodes_.empty()) return;
        if (garbage_ > 64 && garbage_ > prim_indices_.size() / 2) compact_prims();
        if (degraded()) {
            rebuild();
            ++rebuilds_;
        }
    }


    std::uint32_t build_recursive(std::vector<std::size_t>& indices,
                                   std::size_t begin, std::size_t end,
                                   const std::vector<B>& boxes,
                                   const std::vector<PointType>& centroids,
                                   std::size_t depth, std::uint32_t parent) {
        auto node_idx = static_cast<std::uint32_t>(nodes_.size());
        nodes_.push_back({});
        nodes_[node_idx].parent = parent;
        depth_ = std::max(depth_, depth);

        B bounds = boxes[indices[begin]];
        for (std::size_t i = begin + 1; i < end; ++i)
            bounds = merge(bounds, boxes[indices[i]]);

        std::size_t count = end - begin;

        if (count <= LEAF_THRESHOLD) {
            make_leaf(node_idx, bounds, indices, begin, end);
            return node_idx;
        }

        // Centroid bounds
        BoxType cbounds{centroids[indices[begin]], centroids[indices[begin]]};
        for (std::size_t i = begin + 1; i < end; ++i) {
            auto& c = centroids[indices[i]];
            for (std::size_t d = 0; d < N; ++d) {
                cbounds.min_corner[d] = std::min(cbounds.min_corner[d], c[d]);
                cbounds.max_corner[d] = std::max(cbounds.max_corner[d], c[d]);
            }
        }

        // Binned SAH
        T best_cost = std::numeric_limits<T>::max();
        std::size_t best_axis = 0;
        std::size_t best_split_bin = 0;
        T parent_area = sah_measure(bounds);

        struct Bin { B bounds{}; std::uint32_t n{}; bool init{}; };

        for (std::size_t axis = 0; axis < N; ++axis) {
            T ext = cbounds.max_corner[axis] - cbounds.min_corner[axis];
            if (ext < epsilon<T>()) continue;

            std::array<Bin, NUM_BINS> bins{};
            T inv = T{1} / ext;
            T lo = cbounds.min_corner[axis];

            for (std::size_t i = begin; i < end; ++i) {
                auto idx = indices[i];
                int b = std::clamp(int((centroids[idx][axis] - lo) * inv * T(NUM_BINS)),
                                   0, NUM_BINS - 1);
                if (!bins[b].init) { bins[b].bounds = boxes[idx]; bins[b].init = true; }
                else bins[b].bounds = merge(bins[b].bounds, boxes[idx]);
                bins[b].n++;
            }

            // Left sweep
            std::array<T, NUM_BINS - 1> la{};
            std::array<std::uint32_t, NUM_BINS - 1> lc{};
            {
                B run{}; bool ri = false; std::uint32_t cn = 0;
                for (int i = 0; i < NUM_BINS - 1; ++i) {
                    if (bins[i].init) { run = ri ? merge(run, bins[i].bounds) : bins[i].bounds; ri = true; }
                    cn += bins[i].n;
                    lc[i] = cn;
                    la[i] = ri ? sah_measure(run) : T{0};
                }
            }

            // Right sweep + cost
            {
                B run{}; bool ri = false; std::uint32_t cn = 0;
                for (int i = NUM_BINS - 1; i >= 1; --i) {
                    if (bins[i].init) { run = ri ? merge(run, bins[i].bounds) : bins[i].bounds; ri = true; }
                    cn += bins[i].n;
                    T ra = ri ? sah_measure(run) : T{0};
                    T cost = T{1} + (lc[i - 1] * la[i - 1] + cn * ra) / parent_area;
                    if (cost < best_cost) {
                        best_cost = cost;
                        best_axis = axis;
                        best_split_bin = static_cast<std::size_t>(i);
                    }
                }
            }
        }

        // SAH says leaf cheaper? Make leaf.
        if (best_cost >= static_cast<T>(count)) {
            make_leaf(node_idx, bounds, indices, begin, end);
            return node_idx;
        }

        // Partition
        T smin = cbounds.min_corner[best_axis];
        T sext = cbounds.max_corner[best_axis] - smin;
        T inv = T{1} / sext;

        auto mid = std::partition(indices.begin() + begin, indices.begin() + end,
            [&](std::size_t idx) {
                int b = std::clamp(int((centroids[idx][best_axis] - smin) * inv * T(NUM_BINS)),
                                   0, NUM_BINS - 1);
                return static_cast<std::size_t>(b) < best_split_bin;
            });

        auto mid_pos = static_cast<std::size_t>(mid - indices.begin());
        if (mid_pos == begin || mid_pos == end)
            mid_pos = begin + count / 2;

        // Left child built immediately after current: index = node_idx + 1
        build_recursive(indices, begin, mid_pos, boxes, centroids, depth + 1, node_idx);
        auto right = build_recursive(indices, mid_pos, end, boxes, centroids, depth + 1, node_idx);

        // Internal: first = right child index, count = 0
        nodes_[node_idx] = {bounds, right, 0, node_idx + 1, parent};

        return node_idx;
    }

    void make_leaf(std::uint32_t node_idx, const B& bounds,
                   const std::vector<std::size_t>& indices,
                   std::size_t begin, std::size_t end) {
        auto offset = static_cast<std::uint32_t>(prim_indices_.size());
        for (std::size_t i = begin; i < end; ++i) {
            prim_indices_.push_back(indices[i]);
            leaf_of_[indices[i]] = node_idx;
        }
        nodes_[node_idx] = {bounds, offset, static_cast<std::uint32_t>(end - begin), npos,
                            nodes_[node_idx].parent};
    }
};

} // namespace spatium::spatial
