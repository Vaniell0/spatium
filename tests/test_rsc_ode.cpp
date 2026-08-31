#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <ode_task.hpp>
#include <train.hpp>

using Catch::Matchers::WithinAbs;

TEST_CASE("ode registry has exactly ode_euler and ode_rk4", "[rsc]") {
    auto reg = rsc::build_ode_registry();
    REQUIRE(reg.size() == 2);
    CHECK(reg[reg.index_of("ode_euler")].signature().in_size == 4);
    CHECK(reg[reg.index_of("ode_rk4")].signature().in_size == 4);
}

TEST_CASE("OdeTaskGenerator produces a real mix of euler-wins and rk4-wins", "[rsc]") {
    auto reg = rsc::build_ode_registry();
    auto gen = rsc::build_ode_task_generator(reg, 1);

    std::size_t rk4_index = reg.index_of("ode_rk4");
    int need_rk4 = 0;
    const int n = 2000;
    for (int i = 0; i < n; ++i) {
        auto task = gen.sample();
        REQUIRE(task.inputs.size() == 3);
        REQUIRE(task.expected_output.empty());
        if (task.op_index == rk4_index) ++need_rk4;
    }
    double frac = static_cast<double>(need_rk4) / n;
    CHECK(frac > 0.15);
    CHECK(frac < 0.85);
}

TEST_CASE("training raises ode-dispatch accuracy on held-out problems", "[rsc]") {
    auto reg = rsc::build_ode_registry();
    rsc::Dispatcher model(3, 16, 2, /*seed=*/0);

    auto eval_acc = [&](rsc::OdeTaskGenerator& gen, int n) {
        int correct = 0;
        for (int i = 0; i < n; ++i)
            if (rsc::argmax_correct(model, gen.sample(), 3, 0)) ++correct;
        return static_cast<double>(correct) / n;
    };

    auto pre_eval = rsc::build_ode_task_generator(reg, 12345);
    double pre = eval_acc(pre_eval, 2000);

    auto train_gen = rsc::build_ode_task_generator(reg, 1);
    std::mt19937_64 rng(42);
    double baseline = 0.5;
    for (int u = 0; u < 1500; ++u) {
        rsc::Gradients sum;
        for (int b = 0; b < 32; ++b) {
            auto task = train_gen.sample();
            bool ok;
            auto g = rsc::reinforce_gradient(model, task, 3, 0, rng, baseline, ok);
            rsc::accumulate_gradients(sum, g);
        }
        rsc::scale_gradients(sum, 1.0 / 32.0);
        model.apply_gradients(sum, 0.1);
    }

    auto post_eval = rsc::build_ode_task_generator(reg, 777);
    double post = eval_acc(post_eval, 2000);

    CHECK(post > pre);
    CHECK(post > 0.8);
}
