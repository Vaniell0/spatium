#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <mesh_task.hpp>
#include <train.hpp>

using Catch::Matchers::WithinAbs;

TEST_CASE("mesh registry has exactly mesh_uniform and mesh_adapted", "[rsc]") {
    auto reg = rsc::build_mesh_registry();
    REQUIRE(reg.size() == 2);
    CHECK(reg[reg.index_of("mesh_uniform")].signature().in_size == 6);
    CHECK(reg[reg.index_of("mesh_adapted")].signature().in_size == 6);
}

TEST_CASE("MeshTaskGenerator produces a real mix of uniform-wins and adapted-wins", "[rsc]") {
    auto gen = rsc::build_mesh_task_generator(1);

    const int n = 2000;
    int need_adapted = 0;
    for (int i = 0; i < n; ++i) {
        auto task = gen.sample();
        REQUIRE(task.inputs.size() == 1);
        REQUIRE(task.expected_output.empty());
        if (task.op_index == 1) ++need_adapted;
    }
    double frac = static_cast<double>(need_adapted) / n;
    // Stratified 50/50 sampling between buckets -- close to balanced by
    // construction, not by measurement, so a loose band is enough here.
    CHECK(frac > 0.3);
    CHECK(frac < 0.7);
}

TEST_CASE("training raises mesh-dispatch accuracy on held-out problems", "[rsc]") {
    rsc::Dispatcher model(1, 16, 2, /*seed=*/0);

    auto eval_acc = [&](rsc::MeshTaskGenerator& gen, int n) {
        int correct = 0;
        for (int i = 0; i < n; ++i)
            if (rsc::argmax_correct(model, gen.sample(), 1, 0)) ++correct;
        return static_cast<double>(correct) / n;
    };

    auto pre_eval = rsc::build_mesh_task_generator(12345);
    double pre = eval_acc(pre_eval, 2000);

    auto train_gen = rsc::build_mesh_task_generator(1);
    std::mt19937_64 rng(42);
    double baseline = 0.5;
    for (int u = 0; u < 1500; ++u) {
        rsc::Gradients sum;
        for (int b = 0; b < 32; ++b) {
            auto task = train_gen.sample();
            bool ok;
            auto g = rsc::reinforce_gradient(model, task, 1, 0, rng, baseline, ok);
            rsc::accumulate_gradients(sum, g);
        }
        rsc::scale_gradients(sum, 1.0 / 32.0);
        model.apply_gradients(sum, 0.1);
    }

    auto post_eval = rsc::build_mesh_task_generator(777);
    double post = eval_acc(post_eval, 2000);

    // Measured, not aspirational: a standalone probe run to 4000 updates
    // plateaus at ~0.77-0.79 (checked every 500 updates, never crosses
    // 0.8) -- real, substantial signal over the ~0.26 pre-training
    // baseline, just not a near-ceiling domain like precision's ~0.99.
    CHECK(post > pre);
    CHECK(post > 0.7);
}
