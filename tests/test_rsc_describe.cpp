#include <catch2/catch_test_macros.hpp>
#include <describe.hpp>
#include <task.hpp>
#include <tier1_ops.hpp>

TEST_CASE("describe() reads the op's own names, not positional indices", "[rsc]") {
    auto reg = rsc::build_tier1_registry();
    auto i = reg.index_of("add");

    std::array<double, 2> in{2.0, 3.0};
    double out[1];
    reg[i](in, out);

    auto s = rsc::describe(reg[i].signature(), in, out);
    CHECK(s == "add(a=2, b=3) -> sum=5");
}

TEST_CASE("describe() reports a vector op's named components", "[rsc]") {
    auto reg = rsc::build_tier1_registry();
    auto i = reg.index_of("lerp3");

    std::array<double, 7> in{0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.5};
    double out[3];
    reg[i](in, out);

    auto s = rsc::describe(reg[i].signature(), in, out);
    CHECK(s == "lerp3(a.x=0, a.y=0, a.z=0, b.x=10, b.y=0, b.z=0, t=0.5) -> x=5, y=0, z=0");
}

TEST_CASE("describe() falls back to positional labels when names are missing", "[rsc]") {
    rsc::OpSignature sig{.name = "mystery", .tier = rsc::Tier::General,
                          .in_size = 2, .out_size = 1};
    std::array<double, 2> in{1.0, 2.0};
    double out[1]{3.0};
    auto s = rsc::describe(sig, in, out);
    CHECK(s == "mystery(in0=1, in1=2) -> out0=3");
}

TEST_CASE("describe() reflects every registered Tier-1 op's actual behavior", "[rsc]") {
    // Same generator used to validate the training mechanism -- describe()
    // must never disagree with what the registry itself computed.
    auto reg = rsc::build_tier1_registry();
    rsc::TaskGenerator gen(reg, 11);

    for (int i = 0; i < 10; ++i) {
        auto task = gen.sample();
        auto s = rsc::describe(reg[task.op_index].signature(), task.inputs, task.expected_output);
        CHECK(s.starts_with(reg[task.op_index].signature().name + "("));
    }
}
