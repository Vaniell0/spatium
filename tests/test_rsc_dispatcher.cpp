#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <dispatcher.hpp>
#include <tier1_ops.hpp>
#include <cmath>
#include <numeric>

using Catch::Matchers::WithinAbs;

TEST_CASE("features() builds a fixed-size, correctly padded input vector", "[rsc]") {
    auto reg = rsc::build_tier1_registry();
    rsc::TaskGenerator gen(reg, 3);
    std::size_t max_in = reg.max_in_size();
    std::size_t max_out = reg.max_out_size();
    REQUIRE(max_in == 7);  // lerp3
    REQUIRE(max_out == 3); // lerp3

    for (int i = 0; i < 20; ++i) {
        auto task = gen.sample();
        auto x = rsc::features(task, max_in, max_out);
        REQUIRE(x.size() == max_in + max_out);

        for (std::size_t k = 0; k < task.inputs.size(); ++k) CHECK(x[k] == task.inputs[k]);
        for (std::size_t k = task.inputs.size(); k < max_in; ++k) CHECK(x[k] == 0.0);

        for (std::size_t k = 0; k < task.expected_output.size(); ++k)
            CHECK(x[max_in + k] == task.expected_output[k]);
        for (std::size_t k = task.expected_output.size(); k < max_out; ++k)
            CHECK(x[max_in + k] == 0.0);
    }
}

TEST_CASE("softmax produces a valid probability distribution", "[rsc]") {
    auto probs = rsc::softmax({1.0, 2.0, 0.5, -1.0});
    double sum = 0.0;
    for (auto p : probs) {
        CHECK(p >= 0.0);
        CHECK(p <= 1.0);
        sum += p;
    }
    CHECK_THAT(sum, WithinAbs(1.0, 1e-12));
}

TEST_CASE("Dispatcher forward pass produces one finite logit per registered op", "[rsc]") {
    auto reg = rsc::build_tier1_registry();
    rsc::Dispatcher model(reg.max_in_size() + reg.max_out_size(), 16, reg.size(), /*seed=*/0);

    rsc::TaskGenerator gen(reg, 5);
    for (int i = 0; i < 20; ++i) {
        auto task = gen.sample();
        auto x = rsc::features(task, reg.max_in_size(), reg.max_out_size());
        auto logits = model.forward(x);
        REQUIRE(logits.size() == reg.size());
        for (double v : logits) CHECK(std::isfinite(v));

        auto probs = rsc::softmax(logits);
        double sum = std::accumulate(probs.begin(), probs.end(), 0.0);
        CHECK_THAT(sum, WithinAbs(1.0, 1e-9));
    }
}

TEST_CASE("backward() matches numerical gradients (finite differences)", "[rsc]") {
    // The actual validation for the hand-derived backward pass: perturb
    // every weight by +-eps, measure the resulting change in loss, and
    // check it agrees with backward()'s closed-form gradient. This is the
    // standard way to trust a hand-rolled backward pass before it ever
    // touches a training loop.
    struct Case { std::uint64_t seed; std::vector<double> x; std::size_t target; };
    std::vector<Case> cases{
        {17, {0.3, -1.2, 2.0, 0.7, -0.5}, 1},
        {3,  {-2.0, 1.0, 0.0, 4.5, -3.3}, 0},
        {99, {1.0, 1.0, 1.0, 1.0, 1.0}, 2},
    };

    for (const auto& c : cases) {
        rsc::Dispatcher model(5, 4, 3, c.seed);

        auto cache = model.forward_cached(c.x);
        auto dlogits = rsc::cross_entropy_dlogits(cache.logits, c.target);
        auto grads = model.backward(cache, dlogits);

        auto loss_at = [&]() { return rsc::cross_entropy_loss(model.forward(c.x), c.target); };

        auto check_param = [&](std::vector<double>& param, const std::vector<double>& analytic) {
            const double eps = 1e-6;
            REQUIRE(param.size() == analytic.size());
            for (std::size_t i = 0; i < param.size(); ++i) {
                double orig = param[i];
                param[i] = orig + eps;
                double lp = loss_at();
                param[i] = orig - eps;
                double lm = loss_at();
                param[i] = orig;
                double numeric = (lp - lm) / (2 * eps);
                CHECK_THAT(analytic[i], WithinAbs(numeric, 1e-4));
            }
        };

        check_param(model.w1(), grads.dw1);
        check_param(model.b1(), grads.db1);
        check_param(model.w2(), grads.dw2);
        check_param(model.b2(), grads.db2);
    }
}

TEST_CASE("cross_entropy_dlogits is softmax minus one-hot", "[rsc]") {
    std::vector<double> logits{1.0, 2.0, 0.5};
    auto probs = rsc::softmax(logits);
    auto d = rsc::cross_entropy_dlogits(logits, 1);
    for (std::size_t k = 0; k < logits.size(); ++k) {
        double expected = probs[k] - (k == 1 ? 1.0 : 0.0);
        CHECK_THAT(d[k], WithinAbs(expected, 1e-12));
    }
}

TEST_CASE("Dispatcher is deterministic given the same seed", "[rsc]") {
    auto reg = rsc::build_tier1_registry();
    std::size_t dim = reg.max_in_size() + reg.max_out_size();
    rsc::Dispatcher a(dim, 16, reg.size(), 99);
    rsc::Dispatcher b(dim, 16, reg.size(), 99);

    std::vector<double> x(dim, 1.0);
    CHECK(a.forward(x) == b.forward(x));
}
