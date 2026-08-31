#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <linear_task.hpp>
#include <train.hpp>
#include <cmath>

using Catch::Matchers::WithinAbs;

TEST_CASE("linear registry has exactly linalg_jacobi and linalg_direct", "[rsc]") {
    auto reg = rsc::build_linear_registry();
    REQUIRE(reg.size() == 2);
    CHECK(reg[reg.index_of("linalg_jacobi")].signature().in_size == 20);
    CHECK(reg[reg.index_of("linalg_direct")].signature().in_size == 20);
    CHECK(reg[reg.index_of("linalg_jacobi")].signature().out_size == 4);
}

TEST_CASE("linalg_jacobi converges for a strongly diagonally dominant system", "[rsc]") {
    auto reg = rsc::build_linear_registry();
    spatium::Matrix<double, 4, 4> A;
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 4; ++j) A(i, j) = (i == j) ? 6.0 : 1.0;  // ratio = 2.0
    spatium::Vec<double, 4> b{1.0, 2.0, -1.0, 0.5};

    std::vector<double> in(20);
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 4; ++j) in[i * 4 + j] = A(i, j);
    for (std::size_t i = 0; i < 4; ++i) in[16 + i] = b[i];

    std::vector<double> jac_out(4), dir_out(4);
    reg[reg.index_of("linalg_jacobi")](in, jac_out);
    reg[reg.index_of("linalg_direct")](in, dir_out);
    for (std::size_t i = 0; i < 4; ++i) CHECK_THAT(jac_out[i], WithinAbs(dir_out[i], 1e-6));
}

TEST_CASE("linalg_jacobi measurably fails within its budget once the matrix isn't diagonally "
          "dominant enough, while linalg_direct stays exact",
          "[rsc]") {
    auto reg = rsc::build_linear_registry();
    spatium::Matrix<double, 4, 4> A;
    // ratio = 6.0 / (3*3.0) = 0.667 < 1: outside Jacobi's guaranteed-
    // convergence regime. Checked directly (not assumed): this specific
    // matrix diverges to max_err ~139 within the 20-iteration budget, not
    // a borderline case.
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 4; ++j) A(i, j) = (i == j) ? 6.0 : 3.0;
    spatium::Vec<double, 4> b{1.0, 2.0, -1.0, 0.5};

    std::vector<double> in(20);
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 4; ++j) in[i * 4 + j] = A(i, j);
    for (std::size_t i = 0; i < 4; ++i) in[16 + i] = b[i];

    std::vector<double> jac_out(4), dir_out(4);
    reg[reg.index_of("linalg_jacobi")](in, jac_out);
    reg[reg.index_of("linalg_direct")](in, dir_out);
    double max_err = 0.0;
    for (std::size_t i = 0; i < 4; ++i) max_err = std::max(max_err, std::abs(jac_out[i] - dir_out[i]));
    CHECK(max_err > 0.05);
}

TEST_CASE("LinearTaskGenerator produces a real mix of jacobi-wins and direct-wins", "[rsc]") {
    auto reg = rsc::build_linear_registry();
    auto gen = rsc::build_linear_task_generator(reg, 1);

    std::size_t direct_index = reg.index_of("linalg_direct");
    int need_direct = 0;
    const int n = 4000;
    for (int i = 0; i < n; ++i) {
        auto task = gen.sample();
        REQUIRE(task.inputs.size() == 1);
        REQUIRE(task.expected_output.empty());
        if (task.op_index == direct_index) ++need_direct;
    }
    double frac = static_cast<double>(need_direct) / n;
    // Measured, not guessed (standalone probe, 2026-08-26): ~0.36 overall
    // under this exact sampling + tolerance combination.
    CHECK(frac > 0.2);
    CHECK(frac < 0.55);
}

TEST_CASE("training raises linear-solve-dispatch accuracy on held-out problems", "[rsc]") {
    auto reg = rsc::build_linear_registry();
    rsc::Dispatcher model(1, 16, 2, /*seed=*/0);

    auto eval_acc = [&](rsc::LinearTaskGenerator& gen, int n) {
        int correct = 0;
        for (int i = 0; i < n; ++i)
            if (rsc::argmax_correct(model, gen.sample(), 1, 0)) ++correct;
        return static_cast<double>(correct) / n;
    };

    auto pre_eval = rsc::build_linear_task_generator(reg, 12345);
    double pre = eval_acc(pre_eval, 2000);

    auto train_gen = rsc::build_linear_task_generator(reg, 1);
    std::mt19937_64 rng(42);
    double baseline = 0.5;
    for (int u = 0; u < 1000; ++u) {
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

    auto post_eval = rsc::build_linear_task_generator(reg, 777);
    double post = eval_acc(post_eval, 2000);

    CHECK(post > pre);
    CHECK(post > 0.8);
}
