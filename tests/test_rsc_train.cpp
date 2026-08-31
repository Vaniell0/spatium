#include <catch2/catch_test_macros.hpp>
#include <train.hpp>
#include <tier1_ops.hpp>

// Step 6 of the RSC build plan: does the training loop actually work, not
// just run without crashing. Everything here is seeded, so this is a
// deterministic, reproducible check, not a flaky statistical one.
//
// REINFORCE without a baseline was tried first (matching the plan literally)
// and measured, not assumed: eval accuracy rose to ~0.61 then regressed back
// toward the ~0.20-0.25 starting point over further training, and learning
// rates above 0.01 never learned at all -- reward=0 gives zero gradient with
// no baseline, so wrong guesses carry no corrective signal. A running-average
// reward baseline (see train.hpp) fixed both the regression and let lr=0.01
// train stably; this is what's tested below.

TEST_CASE("training measurably raises dispatch accuracy on held-out tasks", "[rsc]") {
    auto reg = rsc::build_tier1_registry();
    std::size_t max_in = reg.max_in_size();
    std::size_t max_out = reg.max_out_size();

    rsc::Dispatcher model(max_in + max_out, 16, reg.size(), /*seed=*/0);

    rsc::TaskGenerator pre_eval(reg, 12345);
    double pre = rsc::evaluate(model, pre_eval, max_in, max_out, 5000);

    rsc::TaskGenerator train_tasks(reg, 1);
    double reward_baseline = 1.0 / static_cast<double>(reg.size()); // chance level
    rsc::train(model, train_tasks, max_in, max_out, /*lr=*/0.01, /*n_steps=*/50000,
               /*seed=*/42, reward_baseline);

    // Different seed from pre_eval -- this is held-out, not the training data.
    rsc::TaskGenerator post_eval(reg, 777);
    double post = rsc::evaluate(model, post_eval, max_in, max_out, 5000);

    CHECK(pre < 0.35);           // untrained: near chance (1/4) for 4 ops
    CHECK(post > pre + 0.2);     // measured improvement was ~0.36; leave margin
    CHECK(post > 0.45);          // measured value was ~0.56; leave margin
}
