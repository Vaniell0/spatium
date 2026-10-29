// A softmax over a linear score, and the AdaGrad step that trains it.
//
// Each claim is held against a second route: the gradient written in softmax.hpp against Dual
// carried through the cross-entropy, column by column, and against a finite difference; the
// distribution against its own sum, temperatures and the walk of `sample`; the learner against
// a model whose weights are known, so that "it learned" is a measured agreement on lists it
// never saw. AdaGrad's first step is the rate times the sign of the gradient and no step
// changes when every gradient is scaled. Each check gets a case it must fail: a learner whose
// rate is zero, a step taken the wrong way, and plain gradient descent, whose steps do scale.

#include <spatium/algebra/analysis/optimize.hpp>
#include <spatium/algebra/dual.hpp>
#include <spatium/algebra/models/softmax.hpp>

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <random>
#include <vector>

using namespace spatium;
using Catch::Matchers::WithinAbs;

namespace {

constexpr std::size_t columns = 6;
using Model = LinearSoftmax<columns>;
using Row = Model::Row;

std::vector<Row> random_list(std::mt19937& rng, std::size_t size) {
    std::uniform_real_distribution<double> u(0, 1);
    std::vector<Row> list(size);
    for (Row& r : list)
        for (double& v : r)
            v = u(rng);
    return list;
}

std::vector<double> distribution_of(const Model& m, const std::vector<Row>& list, double temperature = 1) {
    std::vector<double> p(list.size());
    m.distribution(list, p, temperature);
    return p;
}

// One AdaGrad step toward `chosen`, scaled by `weight` (negative moves away); the cross-entropy
// before it.
double learn(Model& m, AdaGrad<columns>& rule, std::span<const Row> list, std::size_t chosen,
             double rate, double weight = 1) {
    const double surprise = m.cross_entropy(list, chosen);
    auto g = m.gradient(list, chosen);
    for (double& v : g)
        v *= weight;
    rule.step(m.weights(), g, rate);
    return surprise;
}

} // namespace

TEST_CASE("softmax: the gradient agrees with Dual and with a finite difference", "[softmax]") {
    std::mt19937 rng(3);
    std::normal_distribution<double> n(0, 1);
    for (int round = 0; round < 20; ++round) {
        Row w;
        for (double& v : w)
            v = n(rng);
        const auto list = random_list(rng, 7);
        std::vector<double> target(7, 0.0);
        if (round % 2) {
            target[rng() % 7] = 1; // one move chosen
        } else {
            double total = 0;
            for (double& t : target)
                total += t = std::uniform_real_distribution<double>(0, 1)(rng);
            for (double& t : target)
                t /= total; // a soft target
        }
        const Model model(w);
        const Row g = model.gradient(list, target);

        // Dual: the same model over Dual<double>, the derivative seeded at one weight
        using D = Dual<double>;
        std::vector<LinearSoftmax<columns, D>::Row> dual_list(list.size());
        std::vector<D> dual_target(list.size());
        for (std::size_t i = 0; i < list.size(); ++i) {
            for (std::size_t j = 0; j < columns; ++j)
                dual_list[i][j] = D(list[i][j]);
            dual_target[i] = D(target[i]);
        }
        for (std::size_t j = 0; j < columns; ++j) {
            LinearSoftmax<columns, D>::Weights dw;
            for (std::size_t k = 0; k < columns; ++k)
                dw[k] = k == j ? D::variable(w[k]) : D(w[k]);
            const D loss = LinearSoftmax<columns, D>(dw).cross_entropy(dual_list, dual_target);
            REQUIRE_THAT(g[j], WithinAbs(loss.deriv, 1e-12));
            CHECK_THAT(loss.value, WithinAbs(model.cross_entropy(list, target), 1e-12));

            Row up = w, down = w;
            up[j] += 1e-6, down[j] -= 1e-6;
            const double numeric =
                (Model(up).cross_entropy(list, target) - Model(down).cross_entropy(list, target)) / 2e-6;
            REQUIRE_THAT(g[j], WithinAbs(numeric, 1e-6));
        }
        if (round % 2) {
            const std::size_t chosen = static_cast<std::size_t>(std::ranges::find(target, 1.0) - target.begin());
            const Row one = model.gradient(list, chosen);
            for (std::size_t j = 0; j < columns; ++j)
                REQUIRE_THAT(one[j], WithinAbs(g[j], 1e-14));
            REQUIRE_THAT(model.cross_entropy(list, chosen), WithinAbs(model.cross_entropy(list, target), 1e-14));
        }
    }
}

TEST_CASE("softmax: a distribution, a sample and a choice", "[softmax]") {
    std::mt19937 rng(5);
    const Model model(Row{2, -1, 0.5, 3, 0, 1});
    const auto list = random_list(rng, 12);
    const auto p = distribution_of(model, list);
    double total = 0;
    for (double v : p) {
        total += v;
        REQUIRE(v > 0);
    }
    REQUIRE_THAT(total, WithinAbs(1.0, 1e-12));
    // low temperature is the best candidate, high is uniform
    REQUIRE(distribution_of(model, list, 1e-3)[model.choose(list)] > 0.999);
    for (double v : distribution_of(model, list, 1e6))
        REQUIRE_THAT(v, WithinAbs(1.0 / 12, 1e-4));
    // the sample walks the cumulative distribution
    REQUIRE(model.sample(list, 0.0) == 0);
    REQUIRE(model.sample(list, 0.999999999) == 11);
    double acc = 0;
    for (std::size_t i = 0; i < 12; ++i) {
        REQUIRE(model.sample(list, acc + p[i] / 2) == i);
        acc += p[i];
    }
    // scores far beyond exp's range: the largest is taken out first
    const Model huge(Row{800, 0, 0, 0, 0, 0});
    const auto q = distribution_of(huge, list);
    for (double v : q)
        REQUIRE(std::isfinite(v));
}

TEST_CASE("AdaGrad: the first step is the rate, and no step depends on the gradient's scale", "[softmax][adagrad]") {
    const std::array<double, 3> g{0.5, -2e-3, 7};
    AdaGrad<3> rule;
    std::array<double, 3> w{0, 0, 0};
    rule.step(w, g, 0.1);
    CHECK_THAT(w[0], WithinAbs(-0.1, 1e-15));
    CHECK_THAT(w[1], WithinAbs(0.1, 1e-15));
    CHECK_THAT(w[2], WithinAbs(-0.1, 1e-15));

    // the same sequence of gradients, all scaled by 1000: the same steps
    std::mt19937 rng(11);
    std::normal_distribution<double> n(0, 1);
    AdaGrad<3> a, b;
    std::array<double, 3> wa{}, wb{}, plain_a{}, plain_b{};
    for (int i = 0; i < 50; ++i) {
        std::array<double, 3> step{n(rng), n(rng), n(rng)}, scaled = step;
        for (double& v : scaled)
            v *= 1000;
        a.step(wa, step, 0.05);
        b.step(wb, scaled, 0.05);
        for (std::size_t j = 0; j < 3; ++j) {
            plain_a[j] -= 0.05 * step[j];   // plain gradient descent, for the case that must fail
            plain_b[j] -= 0.05 * scaled[j];
        }
    }
    for (std::size_t j = 0; j < 3; ++j) {
        CHECK_THAT(wa[j], WithinAbs(wb[j], 1e-12));
        CHECK(std::abs(plain_a[j] - plain_b[j]) > 1);
    }
    // a parameter whose gradient has always been zero does not move
    AdaGrad<2> idle;
    std::array<double, 2> still{1, 2};
    idle.step(still, {0, 0}, 0.5);
    CHECK(still == std::array<double, 2>{1, 2});
}

TEST_CASE("softmax: it learns the choices of a model it has only watched", "[softmax]") {
    const Model teacher(Row{3, -2, 1, 0, 4, -1});

    auto agreement = [&](const Model& student, std::mt19937& rng) {
        int same = 0;
        constexpr int lists = 400;
        for (int i = 0; i < lists; ++i) {
            const auto list = random_list(rng, 10);
            same += student.choose(list) == teacher.choose(list);
        }
        return static_cast<double>(same) / lists;
    };
    auto train = [&](double rate) {
        Model student;
        AdaGrad<columns> rule;
        std::mt19937 rng(21);
        for (int i = 0; i < 3000; ++i) {
            const auto list = random_list(rng, 10);
            learn(student, rule, list, teacher.choose(list), rate);
        }
        std::mt19937 held_out(99);
        return agreement(student, held_out);
    };

    REQUIRE(train(0.5) > 0.95);
    // chance with ten candidates is 10%; a learner that never moves stays near it
    std::mt19937 rng(99);
    REQUIRE(agreement(Model{}, rng) < 0.3);
    REQUIRE(train(0.0) < 0.3);

    SECTION("the surprise falls as it learns") {
        Model student;
        AdaGrad<columns> rule;
        std::mt19937 stream(8);
        double early = 0, late = 0;
        for (int i = 0; i < 2000; ++i) {
            const auto list = random_list(stream, 10);
            const double s = learn(student, rule, list, teacher.choose(list), 0.5);
            if (i < 100)
                early += s;
            if (i >= 1900)
                late += s;
        }
        REQUIRE(late < 0.5 * early);
    }
    SECTION("a negative weight moves away from the move") {
        std::mt19937 stream(4);
        const auto list = random_list(stream, 10);
        const std::size_t move = 3;
        const double before = distribution_of(Model{}, list)[move];
        Model push, pull;
        AdaGrad<columns> rp, rq;
        learn(push, rp, list, move, 0.5, -1.0);
        learn(pull, rq, list, move, 0.5, 1.0);
        REQUIRE(distribution_of(push, list)[move] < before);
        REQUIRE(distribution_of(pull, list)[move] > before);
    }
}
