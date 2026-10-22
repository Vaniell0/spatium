#pragma once

// Barnes-Hut evaluation of a `GravityField` on a `spatial::BVH`.
//
// The tree is the library's BVH over the sources, each a degenerate box; the
// monopole of every node (total mass, centre of mass, size) lives in a side
// array indexed by node, computed bottom-up after a build. The BVH knows
// nothing of it, and a user of the BVH who wants no gravity pays nothing.
//
// A node is replaced by its monopole when the point is outside its box and
// size / distance < theta, `size` the longest side of the node's box and
// `distance` the distance to its centre of mass; otherwise it is opened. Free
// space: the error is that of the dropped quadrupole, growing as theta^2.
//
// Periodic: only the short-range part of the Ewald split needs a tree. The
// kernel `GravityField::kernel_real` is zero beyond the cutoff (about the
// box side), so the walk is made once per periodic image of the tree whose
// box comes within the cutoff of the point, with the same opening rule and
// the displaced centres of mass. The long-range part is the field's
// reciprocal sum over its structure factor, which is exact and already
// independent of the number of sources per point; the tree adds only the
// opening error of the short-range part, the same as in free space.
//
// The tree reads the field's sources when built. When the field's sources
// change, call `rebuild()`.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/geometry/box.hpp>
#  include <spatium/physics/mechanics/gravity_field.hpp>
#  include <spatium/spatial/bvh.hpp>
#  include <algorithm>
#  include <cstddef>
#  include <cstdint>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::physics::mechanics {

template<std::size_t N, Scalar T = double>
class GravityTree {
public:
    using ScalarType = T;
    static constexpr std::size_t dimension = N;
    using Field = GravityField<N, T>;
    using PointBox = geometry::Box<N, T>;
    using Tree = spatial::BVH<PointBox>;

    // Opening parameter used when none is given. Measured on a periodic
    // random system (docs/ROADMAP.md): about 1% mean error of the force on a
    // torus, 0.25% in free space; the table there gives the error at other theta.
    static constexpr T default_theta = T{0.4};

    // Monopole of one node, stored beside the BVH's node array.
    struct Moment {
        T mass{};
        Vec<T, N> centre{};   // centre of mass
        T size{};             // longest side of the node's box
    };

    explicit GravityTree(const Field& field, T theta = default_theta)
        : field_(&field), theta_(theta) { rebuild(); }

    T theta() const { return theta_; }
    void set_theta(T t) { theta_ = t; }
    const Tree& bvh() const { return tree_; }
    const std::vector<Moment>& moments() const { return moments_; }

    // Build the BVH again over the field's current sources and recompute the moments.
    void rebuild() {
        const auto& src = field_->sources();
        std::vector<PointBox> boxes;
        boxes.reserve(src.size());
        for (const auto& s : src) boxes.push_back(PointBox{s.position, s.position});
        tree_ = Tree::build(std::move(boxes));
        refresh();
    }

    // Recompute the moments from the tree's current nodes: after a build,
    // and after any insert/remove/refit on the BVH. Post-order over the
    // links, so it holds for an updated tree whose nodes are not in
    // depth-first order.
    void refresh() {
        const auto& nodes = tree_.nodes();
        moments_.assign(nodes.size(), Moment{});
        if (nodes.empty()) return;
        std::vector<std::uint32_t> order;
        order.reserve(nodes.size());
        std::vector<std::uint32_t> stack{0};
        while (!stack.empty()) {
            const auto i = stack.back();
            stack.pop_back();
            order.push_back(i);
            if (nodes[i].count == 0) {
                stack.push_back(nodes[i].left);
                stack.push_back(nodes[i].first);
            }
        }
        const auto& src = field_->sources();
        for (auto it = order.rbegin(); it != order.rend(); ++it) {
            const auto& n = nodes[*it];
            auto& m = moments_[*it];
            Vec<T, N> weighted{};
            T mass{};
            if (n.count > 0) {
                for (std::uint32_t k = 0; k < n.count; ++k) {
                    const auto& s = src[tree_.prim_indices()[n.first + k]];
                    mass += s.mass;
                    weighted = Vec<T, N>{weighted + s.position * s.mass};
                }
            } else {
                const auto& l = moments_[n.left];
                const auto& r = moments_[n.first];
                mass = l.mass + r.mass;
                weighted = Vec<T, N>{l.centre * l.mass + r.centre * r.mass};
            }
            m.mass = mass;
            m.centre = mass != T{0} ? Vec<T, N>{weighted / mass} : center_of(n.bounds);
            const auto ext = n.bounds.extents();
            for (std::size_t d = 0; d < N; ++d) m.size = std::max(m.size, ext[d]);
        }
    }

    // Acceleration at x (clamped as the field's law says). `skip` is the
    // index of the source that is the body asking, if it is one.
    Vec<T, N> acceleration(const Vec<T, N>& x, std::size_t skip = Field::npos) const {
        return acceleration(x, theta_, skip);
    }
    Vec<T, N> acceleration(const Vec<T, N>& x, T theta, std::size_t skip) const {
        const Field& f = *field_;
        const Vec<T, N> xw = f.wrap(x);
        Vec<T, N> sum{};
        if (!tree_.nodes().empty()) {
            if (!f.periodic()) {
                walk(xw, Vec<T, N>{}, theta, skip, sum);
            } else {
                std::array<int, N> n;
                n.fill(-1);
                for (;;) {
                    Vec<T, N> sh;
                    for (std::size_t i = 0; i < N; ++i) sh[i] = T(n[i]) * f.box()[i];
                    if (tree_.nodes()[0].bounds.distance(Vec<T, N>{xw - sh}) < f.cutoff())
                        walk(xw, sh, theta, skip, sum);
                    std::size_t i = 0;
                    for (; i < N; ++i) {
                        if (++n[i] <= 1) break;
                        n[i] = -1;
                    }
                    if (i == N) break;
                }
            }
        }
        sum = Vec<T, N>{sum + f.reciprocal_field(xw)};
        return f.clamped(Vec<T, N>{sum * (-f.law().G)});
    }

private:
    const Field* field_;
    T theta_;
    Tree tree_;
    std::vector<Moment> moments_;

    static Vec<T, N> center_of(const PointBox& b) {
        return Vec<T, N>{(b.min_corner + b.max_corner) * T{0.5}};
    }

    // One walk of the tree displaced by `sh`, adding m * kernel_real.
    void walk(const Vec<T, N>& x, const Vec<T, N>& sh, T theta, std::size_t skip,
              Vec<T, N>& sum) const {
        const Field& f = *field_;
        const auto& nodes = tree_.nodes();
        const auto& prims = tree_.prim_indices();
        const auto& src = f.sources();
        const bool shifted = sh.dot(sh) > T{0};
        const Vec<T, N> xs{x - sh};   // the point, in the tree's own frame
        const T cut = f.cutoff();
        const T cut2 = cut * cut;

        thread_local std::vector<std::uint32_t> stack;
        stack.clear();
        stack.push_back(0);
        while (!stack.empty()) {
            const auto i = stack.back();
            stack.pop_back();
            const auto& n = nodes[i];
            const auto& m = moments_[i];
            if (m.mass == T{0} && n.count == 0) continue;
            const T lower = n.bounds.distance(xs);
            if (lower >= cut) continue;   // beyond the cutoff (infinite in free space)
            if (n.count == 0) {
                if (lower > T{0}) {
                    const Vec<T, N> d{xs - m.centre};
                    const T r2 = d.dot(d);
                    if (m.size * m.size < theta * theta * r2 && r2 < cut2 * T{4}) {
                        sum = Vec<T, N>{sum + f.kernel_real(d) * m.mass};
                        continue;
                    }
                }
                stack.push_back(n.left);
                stack.push_back(n.first);
            } else {
                for (std::uint32_t k = 0; k < n.count; ++k) {
                    const std::size_t j = prims[n.first + k];
                    if (j == skip && !shifted) continue;
                    sum = Vec<T, N>{sum + f.kernel_real(Vec<T, N>{xs - src[j].position}) * src[j].mass};
                }
            }
        }
    }
};

} // namespace spatium::physics::mechanics
