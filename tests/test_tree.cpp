// A tree that scores rows, grown by CART, and a policy distilled into it.
//
// Each claim is held against a second route: the tree grown on a step function against the
// function, cell by cell; the softmax of scores against LinearSoftmax's own distribution; a
// tree of one leaf against the closed form of the divergence from a uniform choice,
// log n - H(p); the squared error against depth, which may only fall. Each check gets a case it
// must fail: a tree one split short of the function, and the uniform choice of a single leaf.

#include <spatium/algebra/models/softmax.hpp>
#include <spatium/algebra/models/tree.hpp>
#include <spatium/statistics/information.hpp>
#include <spatium/statistics/tree_fit.hpp>

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <random>
#include <vector>

using namespace spatium;
using Catch::Matchers::WithinAbs;

namespace {

constexpr std::size_t K = 4;
using Row = std::array<double, K>;

std::vector<Row> random_rows(std::mt19937& rng, std::size_t n) {
    std::uniform_real_distribution<double> u(0, 1);
    std::vector<Row> rows(n);
    for (Row& r : rows)
        for (double& v : r)
            v = u(rng);
    return rows;
}

double squared_error(const ScoreTree<K>& tree, const std::vector<Row>& rows, const std::vector<double>& y) {
    double e = 0;
    for (std::size_t i = 0; i < rows.size(); ++i)
        e += std::pow(tree.score(rows[i]) - y[i], 2);
    return e;
}

} // namespace

TEST_CASE("tree: a hand-built tree sends each row to its leaf", "[tree]") {
    using Node = ScoreTree<K>::Node;
    // column 0 at most 0.5: then column 2 at most 0.25; leaves 1, 2, 3
    const ScoreTree<K> tree({Node{0, 0.5, 1, 4, 0}, Node{2, 0.25, 2, 3, 0}, Node{-1, 0, -1, -1, 1},
                             Node{-1, 0, -1, -1, 2}, Node{-1, 0, -1, -1, 3}});
    CHECK(tree.score({0.1, 0, 0.1, 0}) == 1);
    CHECK(tree.score({0.5, 0, 0.9, 0}) == 2); // at the threshold goes left
    CHECK(tree.score({0.7, 0, 0.0, 0}) == 3);
    CHECK(tree.leaves() == 3);
    CHECK(tree.depth() == 2);
    // two inner nodes of 1 + 2 + 32 bits (shape, column of four, threshold), three leaves of 1 + 32
    CHECK(tree.bits() == 2 * 35 + 3 * 33);
    CHECK(ScoreTree<K>{}.depth() == 0);
}

TEST_CASE("tree: CART recovers a step function, and one split short it cannot", "[tree][fit]") {
    std::mt19937 rng(5);
    const auto rows = random_rows(rng, 600);
    std::vector<double> y(rows.size());
    for (std::size_t i = 0; i < rows.size(); ++i)
        y[i] = rows[i][0] <= 0.3 ? -1.0 : (rows[i][3] <= 0.6 ? 0.5 : 2.0);
    const std::span<const Row> sr(rows);
    const std::span<const double> sy(y);
    const auto two = statistics::fit_tree<K>(sr, sy, {.depth = 2});
    CHECK(squared_error(two, rows, y) < 1e-24);
    CHECK(two.leaves() == 3); // the third split would not reduce anything
    for (const Row& r : random_rows(rng, 200)) {
        if (std::abs(r[0] - 0.3) < 0.01 || std::abs(r[3] - 0.6) < 0.01)
            continue; // the learned thresholds sit between two sample values near the true ones
        const double truth = r[0] <= 0.3 ? -1.0 : (r[3] <= 0.6 ? 0.5 : 2.0);
        REQUIRE(two.score(r) == truth);
    }
    const auto one = statistics::fit_tree<K>(sr, sy, {.depth = 1});
    CHECK(squared_error(one, rows, y) > 10);
}

TEST_CASE("tree: the squared error on the sample never rises with depth", "[tree][fit]") {
    std::mt19937 rng(8);
    const auto rows = random_rows(rng, 400);
    std::vector<double> y(rows.size());
    std::normal_distribution<double> noise(0, 0.3);
    for (std::size_t i = 0; i < rows.size(); ++i)
        y[i] = std::sin(6 * rows[i][1]) + rows[i][2] * rows[i][0] + noise(rng);
    double last = 1e300;
    for (std::size_t d = 0; d <= 8; ++d) {
        const auto t = statistics::fit_tree<K>(std::span<const Row>(rows), std::span<const double>(y), {.depth = d});
        const double e = squared_error(t, rows, y);
        CHECK(e <= last + 1e-9);
        CHECK(t.depth() <= d);
        last = e;
    }
}

TEST_CASE("tree: the softmax of scores is LinearSoftmax's distribution", "[tree][softmax]") {
    std::mt19937 rng(3);
    const LinearSoftmax<K> model({2, -1, 0.5, 3});
    const auto list = random_rows(rng, 9);
    for (double T : {1.0, 0.3, 4.0}) {
        std::vector<double> scores(list.size()), a(list.size()), b(list.size());
        for (std::size_t i = 0; i < list.size(); ++i)
            scores[i] = model.score(list[i]);
        softmax(std::span<const double>(scores), std::span<double>(a), T);
        model.distribution(list, b, T);
        for (std::size_t i = 0; i < a.size(); ++i)
            CHECK_THAT(a[i], WithinAbs(b[i], 1e-15));
    }
}

TEST_CASE("tree: a policy distilled into trees of growing depth", "[tree][distill]") {
    std::mt19937 rng(11);
    const LinearSoftmax<K> policy({3, -2, 1, 4});
    // the training sample: every row of 300 lists of eight, with the policy's score as target
    std::vector<std::vector<Row>> lists(300);
    std::vector<Row> rows;
    std::vector<double> targets;
    for (auto& l : lists) {
        l = random_rows(rng, 8);
        for (const Row& r : l)
            rows.push_back(r), targets.push_back(policy.score(r));
    }
    std::vector<std::vector<Row>> held(200);
    for (auto& l : held)
        l = random_rows(rng, 8);

    // mean KL(policy | tree) over the held-out lists, and how often the two choose alike
    auto measure = [&](const ScoreTree<K>& tree) {
        double kl = 0;
        int same = 0;
        for (const auto& l : held) {
            std::vector<double> p(l.size()), s(l.size()), q(l.size());
            policy.distribution(l, p);
            for (std::size_t i = 0; i < l.size(); ++i)
                s[i] = tree.score(l[i]);
            softmax(std::span<const double>(s), std::span<double>(q));
            kl += statistics::kl_divergence(std::span<const double>(p), std::span<const double>(q));
            std::size_t best = 0;
            for (std::size_t i = 1; i < l.size(); ++i)
                if (s[i] > s[best])
                    best = i;
            same += best == policy.choose(l);
        }
        return std::pair{kl / double(held.size()), double(same) / double(held.size())};
    };

    // one leaf: every candidate scores alike, the choice is uniform, KL = log n - H(p)
    const auto leaf = statistics::fit_tree<K>(std::span<const Row>(rows), std::span<const double>(targets), {.depth = 0});
    double closed = 0;
    for (const auto& l : held) {
        std::vector<double> p(l.size());
        policy.distribution(l, p);
        closed += std::log(double(l.size())) - statistics::entropy(std::span<const double>(p));
    }
    const auto [kl0, same0] = measure(leaf);
    CHECK_THAT(kl0, WithinAbs(closed / double(held.size()), 1e-12));
    CHECK(same0 < 0.3); // ties go to the first candidate: about chance

    double last_kl = kl0;
    for (std::size_t d : {2u, 4u, 6u, 8u}) {
        const auto tree = statistics::fit_tree<K>(std::span<const Row>(rows), std::span<const double>(targets),
                                                  {.depth = d, .min_leaf = 4});
        const auto [kl, same] = measure(tree);
        INFO("depth " << d << ": KL " << kl << ", agreement " << same << ", " << tree.bits() << " bits");
        CHECK(kl < last_kl);
        last_kl = kl;
        if (d == 8)
            CHECK(same > 0.7);
    }
}
