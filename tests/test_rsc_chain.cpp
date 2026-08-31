#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <chain.hpp>
#include <tier1_ops.hpp>
#include <sstream>

using Catch::Matchers::WithinAbs;

TEST_CASE("Registry::chainable_ops() is exactly add and multiply for Tier-1", "[rsc]") {
    auto reg = rsc::build_tier1_registry();
    auto chainable = reg.chainable_ops();
    REQUIRE(chainable.size() == 2);
    CHECK(reg[chainable[0]].signature().name == "add");
    CHECK(reg[chainable[1]].signature().name == "multiply");
}

TEST_CASE("ChainTaskGenerator produces an internally consistent, state-carrying sequence",
          "[rsc]") {
    auto reg = rsc::build_tier1_registry();
    rsc::ChainTaskGenerator gen(reg, 5);

    for (int i = 0; i < 20; ++i) {
        auto chain = gen.sample(4);
        REQUIRE(chain.steps.size() == 4);

        double state = chain.initial_state;
        for (const auto& step : chain.steps) {
            CHECK(step.inputs[0] == state); // accumulator carried forward
            // ground truth actually matches executing the op on its inputs
            std::vector<double> out(1);
            reg[step.op_index](step.inputs, out);
            CHECK_THAT(out[0], WithinAbs(step.expected_output[0], 1e-12));
            state = step.expected_output[0];
        }
        CHECK_THAT(chain.final_result(), WithinAbs(state, 1e-12));
    }
}

TEST_CASE("describe_chain reports every step readably", "[rsc]") {
    auto reg = rsc::build_tier1_registry();
    rsc::ChainTaskGenerator gen(reg, 1);
    auto chain = gen.sample(2);

    auto s = rsc::describe_chain(reg, chain);
    std::ostringstream expected_prefix;
    expected_prefix << chain.initial_state;
    CHECK(s.starts_with(expected_prefix.str()));
    for (const auto& step : chain.steps)
        CHECK(s.find(reg[step.op_index].signature().name + "(") != std::string::npos);
}

TEST_CASE("single-sample REINFORCE fails to learn chains' add-vs-multiply distinction",
          "[rsc]") {
    // Documents the actual finding, doesn't just assert the fix works in
    // isolation: with no dot3/lerp3 to inflate the number, chains restricted
    // to add/multiply are exactly the hard 50/50 case Tier-1 could hide
    // inside its overall accuracy. Single-sample updates measurably do not
    // solve it -- this is why train_chains_batch() exists.
    auto reg = rsc::build_tier1_registry();
    std::size_t max_in = reg.max_in_size(), max_out = reg.max_out_size();
    rsc::Dispatcher model(max_in + max_out, 16, reg.size(), 0);

    rsc::ChainTaskGenerator train_chains(reg, 1);
    std::mt19937_64 rng(42);
    double baseline = 0.5;
    for (int i = 0; i < 20000; ++i) {
        auto chain = train_chains.sample(3);
        rsc::train_chain(model, chain, max_in, max_out, 0.01, rng, baseline);
    }

    rsc::ChainTaskGenerator eval_chains(reg, 999);
    double acc = rsc::evaluate_chains(model, eval_chains, max_in, max_out, 500, 3);
    CHECK(acc < 0.6); // stayed near chance (measured ~0.50) -- not learning
}

TEST_CASE("minibatched REINFORCE learns chains' add-vs-multiply distinction", "[rsc]") {
    auto reg = rsc::build_tier1_registry();
    std::size_t max_in = reg.max_in_size(), max_out = reg.max_out_size();
    rsc::Dispatcher model(max_in + max_out, 16, reg.size(), 0);

    rsc::ChainTaskGenerator pre_eval(reg, 12345);
    double pre = rsc::evaluate_chains(model, pre_eval, max_in, max_out, 500, 3);

    rsc::ChainTaskGenerator train_chains(reg, 1);
    double baseline = 0.5; // chance level for 2 chainable ops
    rsc::train_chains_batch(model, train_chains, max_in, max_out, /*lr=*/0.1,
                             /*chain_length=*/3, /*batch_size=*/32, /*n_updates=*/1200,
                             /*seed=*/42, baseline);

    rsc::ChainTaskGenerator post_eval(reg, 777);
    double post = rsc::evaluate_chains(model, post_eval, max_in, max_out, 500, 3);

    CHECK(pre < 0.6);
    CHECK(post > 0.85); // measured ~0.97; leave real margin
}
