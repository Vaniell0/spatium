#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <registry.hpp>
#include <tier1_ops.hpp>
#include <array>

using Catch::Matchers::WithinAbs;

TEST_CASE("Registry dispatches ops by index, the model's own output shape", "[rsc]") {
    auto reg = rsc::build_tier1_registry();
    REQUIRE(reg.size() == 4);

    auto i = reg.index_of("add");
    std::array<double, 2> in{2.0, 3.0};
    double out[1];
    reg[i](in, out);
    CHECK(out[0] == 5.0);
}

TEST_CASE("Registry ops match the Spatium functions they wrap", "[rsc]") {
    auto reg = rsc::build_tier1_registry();

    {
        std::array<double, 2> in{4.0, 5.0};
        double out[1];
        reg[reg.index_of("multiply")](in, out);
        CHECK(out[0] == 20.0);
    }
    {
        // dot3((1,0,0),(1,0,0)) = 1
        std::array<double, 6> in{1.0, 0.0, 0.0, 1.0, 0.0, 0.0};
        double out[1];
        reg[reg.index_of("dot3")](in, out);
        CHECK(out[0] == 1.0);
    }
    {
        // lerp3((0,0,0),(10,0,0), t=0.5) = (5,0,0)
        std::array<double, 7> in{0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.5};
        double out[3];
        reg[reg.index_of("lerp3")](in, out);
        CHECK_THAT(out[0], WithinAbs(5.0, 1e-12));
        CHECK_THAT(out[1], WithinAbs(0.0, 1e-12));
        CHECK_THAT(out[2], WithinAbs(0.0, 1e-12));
    }
}

TEST_CASE("Registry rejects a bad-arity dispatch instead of reading past the span", "[rsc]") {
    auto reg = rsc::build_tier1_registry();
    auto add_index = reg.index_of("add");

    // "add" expects 2 inputs, not 3 -- must throw, not silently read in[2].
    std::array<double, 3> bad_in{1.0, 2.0, 3.0};
    double out[1];
    CHECK_THROWS_AS(reg[add_index](bad_in, out), std::invalid_argument);

    // wrong output size too
    std::array<double, 2> good_in{1.0, 2.0};
    double bad_out[2];
    CHECK_THROWS_AS(reg[add_index](good_in, bad_out), std::invalid_argument);
}

TEST_CASE("Registry index is stable and IS the dispatch-head correspondence", "[rsc]") {
    auto reg = rsc::build_tier1_registry();
    // Every op's declared index_of() must match its position -- this is
    // the exact value a trained classification head would need to produce.
    for (std::size_t i = 0; i < reg.size(); ++i)
        CHECK(reg.index_of(reg[i].signature().name) == i);
}
