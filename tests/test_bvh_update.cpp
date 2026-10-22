#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/geometry/ray_hit.hpp>
#include <spatium/geometry/ray_surface.hpp>
#include <spatium/spatial/bvh.hpp>
#include <spatium/spatial/geodesic_ball.hpp>
#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

using namespace spatium;
using namespace spatium::geometry;
using namespace spatium::spatial;
using Catch::Matchers::WithinAbs;

// ── Bounds other than the box ───────────────────────────────────
//
// A geodesic ball of the Euclidean space, wrapped so its merge and measure
// need no space argument, and a box that counts the merges asked of it.
// Both reach the tree only through the Bound concept.

namespace {

using Space3 = Euclidean<3, double>;
const Space3 space3{};

struct GeoBound {
    using ScalarType = double;
    static constexpr std::size_t ambient_dimension = 3;
    GeodesicBall<Space3> g;
};

GeoBound merge(const GeoBound& a, const GeoBound& b) { return {spatial::merge(space3, a.g, b.g)}; }
double lower_distance(const GeoBound& a, const Vec<double, 3>& p) { return spatial::lower_distance(space3, a.g, p); }
double sah_measure(const GeoBound& a) { return spatial::sah_measure(space3, a.g); }
std::optional<std::pair<double, double>> ray_interval(const Ray<3, double>& ray, const GeoBound& a, double tmax) {
    return spatial::ray_interval(ray, Ball<3, double>{a.g.c, a.g.r}, tmax);
}
Vec<double, 3> center(const GeoBound& a) { return a.g.c; }
bool overlaps(const GeoBound& a, const Box<3, double>& box) { return box.distance(a.g.c) <= a.g.r; }

struct CountingBound {
    using ScalarType = double;
    static constexpr std::size_t ambient_dimension = 3;
    Box<3, double> box;
    static inline std::size_t merges = 0;
};

CountingBound merge(const CountingBound& a, const CountingBound& b) {
    ++CountingBound::merges;
    return {a.box.union_with(b.box)};
}
double lower_distance(const CountingBound& a, const Vec<double, 3>& p) { return a.box.distance(p); }
double sah_measure(const CountingBound& a) { return a.box.surface_measure(); }
std::optional<std::pair<double, double>> ray_interval(const Ray<3, double>& ray, const CountingBound& a, double tmax) {
    return geometry::ray_interval(ray, a.box, tmax);
}
Vec<double, 3> center(const CountingBound& a) { return a.box.centroid(); }
bool overlaps(const CountingBound& a, const Box<3, double>& box) { return a.box.intersects(box); }

// Whether `outer` holds `inner`; the tree's pruning relies on it.
bool contains(const Box<3, double>& outer, const Box<3, double>& inner) {
    for (std::size_t d = 0; d < 3; ++d)
        if (inner.min_corner[d] < outer.min_corner[d] || inner.max_corner[d] > outer.max_corner[d]) return false;
    return true;
}
bool contains(const Ball<3, double>& outer, const Ball<3, double>& inner) {
    return Vec<double, 3>{inner.c - outer.c}.norm() + inner.r <= outer.r;
}
bool contains(const GeoBound& outer, const GeoBound& inner) {
    return Vec<double, 3>{inner.g.c - outer.g.c}.norm() + inner.g.r <= outer.g.r;
}
bool contains(const CountingBound& outer, const CountingBound& inner) { return contains(outer.box, inner.box); }

}  // namespace

static_assert(Bound<GeoBound>);
static_assert(Bound<CountingBound>);

namespace {

using Quad = BoundedQuadric<double>;

// The shape lives in this namespace so that `bound_of<B>` for the two
// bounds above is found by ADL, as a user's own shape and bound would be.
struct Inst : Instanced<Quad> {};

template<typename B, typename Shape>
    requires std::same_as<B, GeoBound>
B bound_of(const Shape& s) {
    const auto box = s.bounding_box();
    return B{{box.centroid(), box.extents().norm() / 2 * (1 + 1e-12)}};
}
template<typename B, typename Shape>
    requires std::same_as<B, CountingBound>
B bound_of(const Shape& s) {
    return B{s.bounding_box()};
}

const Quad unit_ball = Quad::sphere(1.0);

struct Rng {
    std::mt19937_64 g{20261022};
    double u() { return std::uniform_real_distribution<double>(0.0, 1.0)(g); }
    std::size_t below(std::size_t n) { return std::uniform_int_distribution<std::size_t>(0, n - 1)(g); }
};

Inst random_ball(Rng& r) {
    Inst i;
    i.shape = &unit_ball;
    i.translation = {r.u() * 40 - 20, r.u() * 40 - 20, r.u() * 40 - 20};
    i.scale = 0.2 + r.u() * 0.8;
    return i;
}

Ray<3, double> random_ray(Rng& r) {
    const Vec<double, 3> o{r.u() * 60 - 30, r.u() * 60 - 30, r.u() * 60 - 30};
    const Vec<double, 3> target{r.u() * 40 - 20, r.u() * 40 - 20, r.u() * 40 - 20};
    return *Ray<3, double>::from(o, Vec<double, 3>{target - o});
}

// Every shape in exactly one leaf, ranges inside prim_indices, parent
// links agreeing with child links, every bound holding its children (and a
// leaf's bound holding its shapes). Returns the number of leaves reached.
template<typename Tree>
void check_invariants(const Tree& t) {
    using B = typename Tree::BoundType;
    using spatial::bound_of;
    const auto& nodes = t.nodes();
    const auto& prims = t.prim_indices();
    if (t.size() == 0) {
        CHECK(nodes.empty());
        return;
    }
    std::vector<int> seen(t.size(), 0);
    std::vector<std::uint32_t> stack{0}, leaves;
    REQUIRE(nodes[0].parent == ~std::uint32_t{0});
    while (!stack.empty()) {
        const auto i = stack.back();
        stack.pop_back();
        const auto& n = nodes[i];
        if (n.count > 0) {
            REQUIRE(std::size_t(n.first) + n.count <= prims.size());
            for (std::uint32_t k = 0; k < n.count; ++k) {
                REQUIRE(prims[n.first + k] < t.size());
                ++seen[prims[n.first + k]];
            }
            leaves.push_back(i);
        } else {
            for (auto c : {n.left, n.first}) {
                REQUIRE(c < nodes.size());
                REQUIRE(nodes[c].parent == i);
                stack.push_back(c);
            }
        }
    }
    for (std::size_t i = 0; i < seen.size(); ++i) REQUIRE(seen[i] == 1);

    // What pruning relies on: every node on the way up from a shape holds
    // that shape's own bound. (A parent need not contain its child's
    // bound as such -- a built ball tree does not either -- only what is
    // beneath it.)
    for (auto l : leaves)
        for (std::uint32_t k = 0; k < nodes[l].count; ++k) {
            const auto sb = bound_of<B>(t.shapes()[prims[nodes[l].first + k]]);
            for (auto a = l; a != ~std::uint32_t{0}; a = nodes[a].parent) REQUIRE(contains(nodes[a].bounds, sb));
        }
}

// The tree after updates must answer as a fresh build over the same
// shapes does, and as brute force does.
template<typename Tree>
void check_against_rebuild(const Tree& t, Rng& r, int rays = 150) {
    const auto fresh = Tree::build(t.shapes());
    const auto sz = t.size();
    for (int k = 0; k < rays; ++k) {
        const auto ray = random_ray(r);
        const auto a = t.ray_cast(ray);
        const auto b = fresh.ray_cast(ray);
        REQUIRE(a.has_value() == b.has_value());
        if (a) {
            REQUIRE(a->index == b->index);
            REQUIRE(a->t == b->t);
        }
        REQUIRE(t.ray_test(ray) == fresh.ray_test(ray));

        std::optional<double> best;
        for (std::size_t i = 0; i < sz; ++i) {
            const auto h = ray_hit(ray, t.shapes()[i]);
            if (h && h->t >= 0 && (!best || h->t < *best)) best = h->t;
        }
        REQUIRE(best.has_value() == a.has_value());
        if (a) REQUIRE(*best == a->t);
    }
}

template<typename B>
void churn(std::size_t start, int ops) {
    Rng r;
    std::vector<Inst> shapes;
    for (std::size_t i = 0; i < start; ++i) shapes.push_back(random_ball(r));
    auto tree = BVH<Inst, B>::build(shapes);
    check_invariants(tree);
    for (int op = 0; op < ops; ++op) {
        if (tree.size() > 1 && r.below(2) == 0) {
            const auto i = r.below(tree.size());
            const auto last = tree.size() - 1;
            const auto moved = tree.remove(i);
            CHECK(moved == last);
            shapes[i] = shapes.back();   // the documented swap-remove
            shapes.pop_back();
        } else {
            const auto s = random_ball(r);
            CHECK(tree.insert(s) == shapes.size());
            shapes.push_back(s);
        }
        if (op % 40 == 39) {
            REQUIRE(tree.shapes().size() == shapes.size());
            for (std::size_t i = 0; i < shapes.size(); ++i)
                REQUIRE(tree.shapes()[i].translation[0] == shapes[i].translation[0]);
            check_invariants(tree);
            check_against_rebuild(tree, r);
            CHECK_THAT(tree.cost(), WithinAbs(tree.measured_cost(), 1e-6 * tree.measured_cost()));
        }
    }
}

}  // namespace

TEST_CASE("BVH update over boxes follows a churning shape set", "[bvh][update]") {
    churn<Box<3, double>>(400, 400);
}

TEST_CASE("BVH update over balls follows a churning shape set", "[bvh][update]") {
    churn<Ball<3, double>>(400, 400);
}

TEST_CASE("BVH update over geodesic balls follows a churning shape set", "[bvh][update]") {
    churn<GeoBound>(300, 300);
}

TEST_CASE("BVH update calls merge O(depth) times per operation", "[bvh][update]") {
    Rng r;
    std::vector<Inst> shapes;
    for (int i = 0; i < 4000; ++i) shapes.push_back(random_ball(r));
    auto tree = BVH<Inst, CountingBound>::build(shapes);
    tree.set_degradation_limit(1e9);   // no rebuild inside the measurement

    CountingBound::merges = 0;
    tree.insert(random_ball(r));
    const auto for_insert = CountingBound::merges;
    CountingBound::merges = 0;
    tree.remove(r.below(tree.size()));
    const auto for_remove = CountingBound::merges;
    INFO("depth " << tree.depth() << ", merges: insert " << for_insert << ", remove " << for_remove);
    CHECK(for_insert <= 4 * tree.depth() + 8);
    CHECK(for_remove <= 4 * tree.depth() + 8);
    CHECK(for_insert < 4000 / 20);
}

TEST_CASE("BVH update on empty and tiny trees", "[bvh][update]") {
    Rng r;
    using Tree = BVH<Inst>;

    SECTION("insert into an empty tree") {
        auto tree = Tree::build({});
        CHECK(tree.size() == 0);
        const auto s = random_ball(r);
        CHECK(tree.insert(s) == 0);
        CHECK(tree.size() == 1);
        check_invariants(tree);
        const Ray<3, double> ray = *Ray<3, double>::from(
            Vec<double, 3>{s.translation + Vec<double, 3>{0, 0, 10}}, Vec<double, 3>{0, 0, -1});
        REQUIRE(tree.ray_cast(ray));
        CHECK(tree.ray_cast(ray)->index == 0);
    }

    SECTION("remove the last shape, then insert again") {
        auto tree = Tree::build({random_ball(r)});
        CHECK(tree.remove(0) == 0);
        CHECK(tree.size() == 0);
        CHECK(tree.node_count() == 0);
        CHECK_FALSE(tree.ray_cast(random_ray(r)));
        CHECK_FALSE(tree.ray_test(random_ray(r)));
        CHECK(tree.cost() == 0.0);
        tree.insert(random_ball(r));
        tree.insert(random_ball(r));
        check_invariants(tree);
        CHECK(tree.size() == 2);
    }

    SECTION("emptying a larger tree one shape at a time") {
        std::vector<Inst> shapes;
        for (int i = 0; i < 60; ++i) shapes.push_back(random_ball(r));
        auto tree = Tree::build(shapes);
        while (tree.size() > 0) {
            tree.remove(r.below(tree.size()));
            check_invariants(tree);
        }
        CHECK(tree.node_count() == 0);
    }

    SECTION("growing from empty one shape at a time") {
        auto tree = Tree::build({});
        for (int i = 0; i < 300; ++i) tree.insert(random_ball(r));
        check_invariants(tree);
        check_against_rebuild(tree, r);
    }

    SECTION("many identical shapes (equal centroids) split without looping") {
        Tree tree = Tree::build({});
        Inst same = random_ball(r);
        for (int i = 0; i < 40; ++i) tree.insert(same);
        check_invariants(tree);
        const Ray<3, double> ray = *Ray<3, double>::from(
            Vec<double, 3>{same.translation + Vec<double, 3>{0, 0, 10}}, Vec<double, 3>{0, 0, -1});
        CHECK(tree.ray_cast(ray));
    }
}

TEST_CASE("BVH quality guard rebuilds a degraded tree and keeps indices", "[bvh][update]") {
    // Swap every shape for one from a flattened distribution, a few at a
    // time. The first tree knows nothing of where they land, so its cost
    // climbs above what a build over the same shapes reaches.
    auto run = [](double limit, std::vector<Inst>& shapes) {
        Rng r;
        shapes.clear();
        for (int i = 0; i < 1000; ++i) shapes.push_back(random_ball(r));
        auto tree = BVH<Inst>::build(shapes);
        tree.set_degradation_limit(limit);
        for (int i = 0; i < 2000; ++i) {
            const auto k = r.below(tree.size());
            tree.remove(k);
            shapes[k] = shapes.back();
            shapes.pop_back();
            auto s = random_ball(r);
            s.translation[2] *= 0.05;
            tree.insert(s);
            shapes.push_back(s);
        }
        return tree;
    };

    std::vector<Inst> shapes;
    const auto lax = run(1e9, shapes);
    CHECK(lax.rebuilds() == 0);
    const auto fresh = BVH<Inst>::build(lax.shapes());
    INFO("updated " << lax.cost() << ", fresh " << fresh.cost());
    CHECK(lax.cost() > fresh.cost());

    const auto guarded = run(1.0, shapes);
    CHECK(guarded.rebuilds() >= 1);
    REQUIRE(guarded.shapes().size() == shapes.size());
    for (std::size_t i = 0; i < shapes.size(); ++i)   // indices survive every rebuild
        REQUIRE(guarded.shapes()[i].translation[0] == shapes[i].translation[0]);
    CHECK(guarded.cost() < 1.1 * BVH<Inst>::build(guarded.shapes()).cost());
    CHECK_THAT(guarded.cost(), WithinAbs(guarded.measured_cost(), 1e-6 * guarded.measured_cost()));
    Rng r;
    check_invariants(guarded);
    check_against_rebuild(guarded, r);
}
