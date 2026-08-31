#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <custom_layer.hpp>
#include <cstdio>
#include <unistd.h>

using Catch::Matchers::WithinAbs;

namespace {
std::string temp_custom_layer_path() {
    return "/tmp/rsc_custom_layer_test_" + std::to_string(::getpid()) + ".txt";
}
} // namespace

TEST_CASE("save_custom_layer + load_custom_layer round-trips a keyed calibration cache", "[rsc]") {
    // Synthetic entry -- no real calibration scenario wired to this yet
    // (contact physics is blocked, see rsc/README.md), this is testing
    // the generic mechanism only.
    rsc::CustomLayer layer;
    layer.base_commit_sha = "deadbeef";
    layer.params["cloth_on_sphere_r1.0"] = {0.001, 0.5, 12.0};
    layer.params["cloth_on_torus_r0.35"] = {0.0005};

    auto path = temp_custom_layer_path();
    rsc::save_custom_layer(path, layer);
    auto loaded = rsc::load_custom_layer(path);
    std::remove(path.c_str());

    CHECK(loaded.base_commit_sha == "deadbeef");
    REQUIRE(loaded.params.size() == 2);
    REQUIRE(loaded.params.count("cloth_on_sphere_r1.0") == 1);
    auto& v = loaded.params["cloth_on_sphere_r1.0"];
    REQUIRE(v.size() == 3);
    CHECK_THAT(v[0], WithinAbs(0.001, 1e-15));
    CHECK_THAT(v[1], WithinAbs(0.5, 1e-15));
    CHECK_THAT(v[2], WithinAbs(12.0, 1e-15));
    REQUIRE(loaded.params.count("cloth_on_torus_r0.35") == 1);
    CHECK(loaded.params["cloth_on_torus_r0.35"].size() == 1);
}

TEST_CASE("custom_layer_matches_base catches a base that's moved on", "[rsc]") {
    rsc::CustomLayer layer;
    layer.base_commit_sha = "abc123";
    CHECK(rsc::custom_layer_matches_base(layer, "abc123"));
    CHECK_FALSE(rsc::custom_layer_matches_base(layer, "def456"));
}
