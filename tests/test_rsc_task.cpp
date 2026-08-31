#include <catch2/catch_test_macros.hpp>
#include <task.hpp>
#include <tier1_ops.hpp>

TEST_CASE("TaskGenerator samples are internally consistent", "[rsc]") {
    auto reg = rsc::build_tier1_registry();
    rsc::TaskGenerator gen(reg, /*seed=*/42);

    for (int i = 0; i < 50; ++i) {
        auto task = gen.sample();
        const auto& sig = reg[task.op_index].signature();
        CHECK(task.inputs.size() == sig.in_size);
        CHECK(task.expected_output.size() == sig.out_size);

        // The ground truth must grade as correct against itself.
        auto g = rsc::grade(task, task.op_index, task.expected_output);
        CHECK(g.passed());
    }
}

TEST_CASE("grade() fails on the wrong op index", "[rsc]") {
    auto reg = rsc::build_tier1_registry();
    rsc::TaskGenerator gen(reg, 1);
    auto task = gen.sample();

    std::size_t wrong_index = (task.op_index + 1) % reg.size();
    auto g = rsc::grade(task, wrong_index, task.expected_output);
    CHECK_FALSE(g.op_correct);
    CHECK_FALSE(g.passed());
}

TEST_CASE("grade() is a hard tolerance gate, not a soft score", "[rsc]") {
    auto reg = rsc::build_tier1_registry();
    rsc::TaskGenerator gen(reg, 2);
    auto task = gen.sample();

    auto nudged = task.expected_output;
    for (auto& v : nudged) v += 1e-10; // within default tol (1e-9)
    CHECK(rsc::grade(task, task.op_index, nudged).passed());

    auto broken = task.expected_output;
    for (auto& v : broken) v += 1.0; // wildly outside tolerance
    CHECK_FALSE(rsc::grade(task, task.op_index, broken).passed());
}

TEST_CASE("TaskGenerator is deterministic given the same seed", "[rsc]") {
    auto reg = rsc::build_tier1_registry();
    rsc::TaskGenerator a(reg, 7);
    rsc::TaskGenerator b(reg, 7);

    for (int i = 0; i < 10; ++i) {
        auto ta = a.sample();
        auto tb = b.sample();
        CHECK(ta.op_index == tb.op_index);
        CHECK(ta.inputs == tb.inputs);
        CHECK(ta.expected_output == tb.expected_output);
    }
}
