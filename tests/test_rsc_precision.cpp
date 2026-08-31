#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <precision_task.hpp>
#include <train.hpp>

using Catch::Matchers::WithinAbs;

TEST_CASE("precision registry has exactly the two solve_cubic precision ops", "[rsc]") {
    auto reg = rsc::build_precision_registry();
    REQUIRE(reg.size() == 2);
    CHECK(reg[reg.index_of("solve_cubic_f64")].signature().in_size == 4);
    CHECK(reg[reg.index_of("solve_cubic_real50")].signature().in_size == 4);
}

TEST_CASE("PrecisionTaskGenerator produces a real mix of easy and delicate cubics",
          "[rsc]") {
    auto reg = rsc::build_precision_registry();
    auto gen = rsc::build_precision_task_generator(reg, 1);

    std::size_t real50_index = reg.index_of("solve_cubic_real50");
    int need_real50 = 0;
    const int n = 2000;
    for (int i = 0; i < n; ++i) {
        auto task = gen.sample();
        REQUIRE(task.inputs.size() == 4);
        REQUIRE(task.expected_output.empty()); // no observed-output trick needed here
        if (task.op_index == real50_index) ++need_real50;
    }
    // Roughly balanced by construction (50/50 easy/delicate regime) -- not
    // exact since the regime alone doesn't guarantee crossing the tolerance,
    // but should be well within a sane range, not near-all-one-class.
    double frac = static_cast<double>(need_real50) / n;
    CHECK(frac > 0.3);
    CHECK(frac < 0.7);
}

TEST_CASE("a well-separated cubic needs only double precision", "[rsc]") {
    auto reg = rsc::build_precision_registry();
    // (x+5)(x)(x-5) = x^3 - 25x, roots far apart
    std::vector<double> in{1.0, 0.0, -25.0, 0.0};
    std::vector<double> f64out(6), r50out(6);
    reg[reg.index_of("solve_cubic_f64")](in, f64out);
    reg[reg.index_of("solve_cubic_real50")](in, r50out);
    for (std::size_t i = 0; i < 6; ++i)
        CHECK_THAT(f64out[i], WithinAbs(r50out[i], 1e-9));
}

TEST_CASE("training raises precision-dispatch accuracy on held-out cubics", "[rsc]") {
    // Real finding, not assumed: the first version of this generator sampled
    // the root cluster's center (r0) from a wide [-10,10] range, and the
    // dispatcher never learned past chance (~0.50) over 1000+ batched
    // updates -- absolute root position dominated the raw coefficients'
    // scale and swamped the actual separation signal the label depends on.
    // Narrowing r0 to [-1,1] (see precision_task.hpp) removed that confound;
    // the same raw coefficients as features, same minibatched REINFORCE
    // already validated on chains, reaches ~0.99 below.
    auto reg = rsc::build_precision_registry();
    rsc::Dispatcher model(4, 16, 2, /*seed=*/0);

    auto pre_eval = rsc::build_precision_task_generator(reg, 12345);
    auto eval_acc = [&](rsc::PrecisionTaskGenerator& gen, int n) {
        int correct = 0;
        for (int i = 0; i < n; ++i)
            if (rsc::argmax_correct(model, gen.sample(), 4, 0)) ++correct;
        return static_cast<double>(correct) / n;
    };
    double pre = eval_acc(pre_eval, 2000);

    auto train_gen = rsc::build_precision_task_generator(reg, 1);
    std::mt19937_64 rng(42);
    double baseline = 0.5;
    for (int u = 0; u < 1000; ++u) {
        rsc::Gradients sum;
        for (int b = 0; b < 32; ++b) {
            auto task = train_gen.sample();
            bool ok;
            auto g = rsc::reinforce_gradient(model, task, 4, 0, rng, baseline, ok);
            rsc::accumulate_gradients(sum, g);
        }
        rsc::scale_gradients(sum, 1.0 / 32.0);
        model.apply_gradients(sum, 0.1);
    }

    auto post_eval = rsc::build_precision_task_generator(reg, 777);
    double post = eval_acc(post_eval, 2000);

    CHECK(pre < 0.65);
    CHECK(post > 0.9); // measured ~0.99; leave real margin
}
