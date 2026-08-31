#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/error.hpp>
#  include <spatium/geometry/concepts.hpp>
#  include <spatium/geometry/box.hpp>
#  include <spatium/geometry/intersection.hpp>
#  include <spatium/geometry/ray_hit.hpp>
#  include <spatium/geometry/triangle.hpp>
#  include <algorithm>
#  include <array>
#  include <cstdint>
#  include <limits>
#  include <numeric>
#  include <optional>
#  include <type_traits>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::spatial {

template<geometry::Bounded Shape>
struct BVH {
    using T = typename Shape::ScalarType;
    static constexpr std::size_t N = Shape::ambient_dimension;
    using PointType = Vec<T, N>;
    using BoxType = geometry::Box<N, T>;

    // BVH itself requires `Bounded` (above) — the Shape must supply
    // `bounding_box()`. When the Shape ALSO satisfies `RayHittable`,
    // Hit carries the shape-specific barycentric weights and unit
    // normal harvested from `ray_hit`. The canonical such shape is
    // `Triangle<3, T>`. `Quadric` / `Torus` are RayHittable but not
    // currently Bounded, so they cannot enter a BVH yet — they
    // shine via direct `ray_quadric` / `ray_torus` calls instead.
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
        if (bvh.shapes_.empty()) return bvh;

        std::vector<BoxType> boxes(bvh.shapes_.size());
        std::vector<PointType> centroids(bvh.shapes_.size());
        for (std::size_t i = 0; i < bvh.shapes_.size(); ++i) {
            boxes[i] = bvh.shapes_[i].bounding_box();
            centroids[i] = boxes[i].centroid();
        }

        std::vector<std::size_t> indices(bvh.shapes_.size());
        std::iota(indices.begin(), indices.end(), 0);

        bvh.nodes_.reserve(2 * bvh.shapes_.size());
        bvh.build_recursive(indices, 0, indices.size(), boxes, centroids);

        return bvh;
    }

    // ── Queries ───────────────────────────────────────────────

    std::optional<Hit> ray_cast(const geometry::Ray<N, T>& ray) const {
        if (nodes_.empty()) return std::nullopt;

        std::optional<Hit> best;
        T best_t = std::numeric_limits<T>::max();

        std::array<std::uint32_t, STACK_DEPTH> stack;
        std::uint32_t sp = 0;
        stack[sp++] = 0;

        while (sp) {
            auto node_idx = stack[--sp];
            auto& node = nodes_[node_idx];

            auto box_hit = geometry::intersect_parameters(ray, node.bounds);
            if (!box_hit || box_hit->first > best_t) continue;

            if (node.count > 0) {
                for (std::uint32_t i = 0; i < node.count; ++i) {
                    auto idx = prim_indices_[node.first + i];
                    if constexpr (has_ray_hit) {
                        auto h = geometry::ray_hit(ray, shapes_[idx]);
                        if (h && h->t >= T{0} && h->t < best_t) {
                            best_t = h->t;
                            best = Hit{idx, h->t, h->point,
                                       h->u, h->v, h->normal};
                        }
                    } else {
                        auto result = geometry::intersect(ray, shapes_[idx]);
                        if (result) {
                            auto t = (result.value() - ray.origin).dot(ray.direction);
                            if (t >= T{0} && t < best_t) {
                                best_t = t;
                                best = Hit{idx, t, result.value()};
                            }
                        }
                    }
                }
            } else {
                auto left = node_idx + 1;
                auto right = node.first;
                auto left_t = geometry::intersect_parameters(ray, nodes_[left].bounds);
                auto right_t = geometry::intersect_parameters(ray, nodes_[right].bounds);

                bool left_ok = left_t && left_t->first <= best_t;
                bool right_ok = right_t && right_t->first <= best_t;

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
    }

    bool ray_test(const geometry::Ray<N, T>& ray) const {
        if (nodes_.empty()) return false;

        std::array<std::uint32_t, STACK_DEPTH> stack;
        std::uint32_t sp = 0;
        stack[sp++] = 0;

        while (sp) {
            auto node_idx = stack[--sp];
            auto& node = nodes_[node_idx];

            auto box_hit = geometry::intersect_parameters(ray, node.bounds);
            if (!box_hit) continue;

            if (node.count > 0) {
                for (std::uint32_t i = 0; i < node.count; ++i) {
                    auto idx = prim_indices_[node.first + i];
                    if (geometry::intersect(ray, shapes_[idx])) return true;
                }
            } else {
                stack[sp++] = node_idx + 1;   // left
                stack[sp++] = node.first;     // right
            }
        }
        return false;
    }

    std::optional<NearestResult> nearest(const PointType& p) const {
        if (nodes_.empty()) return std::nullopt;

        std::optional<NearestResult> best;
        T best_dist = std::numeric_limits<T>::max();

        std::array<std::uint32_t, STACK_DEPTH> stack;
        std::uint32_t sp = 0;
        stack[sp++] = 0;

        while (sp) {
            auto node_idx = stack[--sp];
            auto& node = nodes_[node_idx];

            T box_dist = node.bounds.distance(p);
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
                auto left = node_idx + 1;
                auto right = node.first;
                T dl = nodes_[left].bounds.distance(p);
                T dr = nodes_[right].bounds.distance(p);

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
    }

    std::vector<std::size_t> query_box(const BoxType& query) const {
        std::vector<std::size_t> result;
        if (nodes_.empty()) return result;

        std::array<std::uint32_t, STACK_DEPTH> stack;
        std::uint32_t sp = 0;
        stack[sp++] = 0;

        while (sp) {
            auto node_idx = stack[--sp];
            auto& node = nodes_[node_idx];

            if (!node.bounds.intersects(query)) continue;

            if (node.count > 0) {
                for (std::uint32_t i = 0; i < node.count; ++i) {
                    auto idx = prim_indices_[node.first + i];
                    if (shapes_[idx].bounding_box().intersects(query))
                        result.push_back(idx);
                }
            } else {
                stack[sp++] = node_idx + 1;
                stack[sp++] = node.first;
            }
        }
        return result;
    }

    const std::vector<Shape>& shapes() const { return shapes_; }
    std::size_t node_count() const { return nodes_.size(); }

private:
    static constexpr std::uint32_t LEAF_THRESHOLD = 4;
    static constexpr int NUM_BINS = 12;
    // Traversal stack size: 64 supports up to 2^64 leaves in a balanced BVH.
    // Near-first push is at most 2 per level → depth ≤ STACK_DEPTH/2 worst-case
    // for pathological unbalanced trees; 64 is safe for SAH-built trees up to 10M+ prims.
    static constexpr std::size_t STACK_DEPTH = 64;

    // Node layout:
    //   leaf:     first = prim_indices_ offset, count > 0
    //   internal: first = RIGHT child index, count = 0
    //             LEFT child is always this_node + 1 (built immediately after)
    struct Node {
        BoxType bounds;
        std::uint32_t first{};
        std::uint32_t count{};
    };

    std::vector<Shape> shapes_;
    std::vector<Node> nodes_;
    std::vector<std::size_t> prim_indices_;

    std::uint32_t build_recursive(std::vector<std::size_t>& indices,
                                   std::size_t begin, std::size_t end,
                                   const std::vector<BoxType>& boxes,
                                   const std::vector<PointType>& centroids) {
        auto node_idx = static_cast<std::uint32_t>(nodes_.size());
        nodes_.push_back({});

        BoxType bounds = boxes[indices[begin]];
        for (std::size_t i = begin + 1; i < end; ++i)
            bounds = bounds.union_with(boxes[indices[i]]);

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
        T parent_area = bounds.surface_measure();

        struct Bin { BoxType bounds{}; std::uint32_t n{}; bool init{}; };

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
                else bins[b].bounds = bins[b].bounds.union_with(boxes[idx]);
                bins[b].n++;
            }

            // Left sweep
            std::array<T, NUM_BINS - 1> la{};
            std::array<std::uint32_t, NUM_BINS - 1> lc{};
            {
                BoxType run{}; bool ri = false; std::uint32_t cn = 0;
                for (int i = 0; i < NUM_BINS - 1; ++i) {
                    if (bins[i].init) { run = ri ? run.union_with(bins[i].bounds) : bins[i].bounds; ri = true; }
                    cn += bins[i].n;
                    lc[i] = cn;
                    la[i] = ri ? run.surface_measure() : T{0};
                }
            }

            // Right sweep + cost
            {
                BoxType run{}; bool ri = false; std::uint32_t cn = 0;
                for (int i = NUM_BINS - 1; i >= 1; --i) {
                    if (bins[i].init) { run = ri ? run.union_with(bins[i].bounds) : bins[i].bounds; ri = true; }
                    cn += bins[i].n;
                    T ra = ri ? run.surface_measure() : T{0};
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

        // Left child built immediately after current → index = node_idx + 1
        build_recursive(indices, begin, mid_pos, boxes, centroids);
        auto right = build_recursive(indices, mid_pos, end, boxes, centroids);

        // Internal: first = right child index, count = 0
        nodes_[node_idx] = {bounds, right, 0};

        return node_idx;
    }

    void make_leaf(std::uint32_t node_idx, const BoxType& bounds,
                   const std::vector<std::size_t>& indices,
                   std::size_t begin, std::size_t end) {
        auto offset = static_cast<std::uint32_t>(prim_indices_.size());
        for (std::size_t i = begin; i < end; ++i)
            prim_indices_.push_back(indices[i]);
        nodes_[node_idx] = {bounds, offset, static_cast<std::uint32_t>(end - begin)};
    }
};

} // namespace spatium::spatial
