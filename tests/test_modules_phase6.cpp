// Phase 6 verification — spatium.physics via `import`.
// viewer/ stays out of modules (Vulkan/GLFW/ImGui C-linkage conflicts).
#include <catch2/catch_test_macros.hpp>
#include <cstring>

import spatium.core;
import spatium.algebra;
import spatium.physics;

using spatium::physics::atomic::element;
using spatium::physics::atomic::ELEMENT_COUNT;
using spatium::physics::atomic::subshell_letter;

TEST_CASE("module spatium.physics: element lookup by Z", "[modules][phase6][physics]") {
    const auto& h = element(1);
    REQUIRE(std::strcmp(h.symbol, "H") == 0);
    REQUIRE(h.Z == 1);
    REQUIRE(h.total_electrons() == 1);

    const auto& fe = element(26);
    REQUIRE(std::strcmp(fe.symbol, "Fe") == 0);
    REQUIRE(fe.total_electrons() == 26);
}

TEST_CASE("module spatium.physics: element lookup by symbol", "[modules][phase6][physics]") {
    const auto& c = element(std::string_view{"C"});
    REQUIRE(c.Z == 6);
    REQUIRE(c.total_electrons() == 6);
}

TEST_CASE("module spatium.physics: subshell letter + 118 elements", "[modules][phase6][physics]") {
    REQUIRE(ELEMENT_COUNT == 118);
    REQUIRE(subshell_letter(0) == 's');
    REQUIRE(subshell_letter(1) == 'p');
    REQUIRE(subshell_letter(2) == 'd');
    REQUIRE(subshell_letter(3) == 'f');
    const auto& og = element(118);
    REQUIRE(std::strcmp(og.symbol, "Og") == 0);
}
