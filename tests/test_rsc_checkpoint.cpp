#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <checkpoint.hpp>
#include <cstdio>
#include <unistd.h>

using Catch::Matchers::WithinAbs;

namespace {
// Unique-enough temp path per test run -- no fixture framework needed
// for a single-file round trip.
std::string temp_checkpoint_path() {
    return "/tmp/rsc_checkpoint_test_" + std::to_string(::getpid()) + ".txt";
}
} // namespace

TEST_CASE("snapshot_base_registry covers all 10 ops at their documented offsets", "[rsc]") {
    auto snap = rsc::snapshot_base_registry();
    REQUIRE(snap.size() == rsc::kBaseNumOps);

    // Tier1: indices 0-3, Precision: 4-5, Rootfind: 6-7, Ode: 8-9 --
    // matches kBaseOpOffset directly, checked here so a future change to
    // either drifts loudly instead of silently.
    std::vector<std::string> names_by_index(10);
    for (auto& e : snap) names_by_index[e.global_index] = e.name;

    CHECK(names_by_index[0] == "add");
    CHECK(names_by_index[1] == "multiply");
    CHECK(names_by_index[2] == "dot3");
    CHECK(names_by_index[3] == "lerp3");
    CHECK(names_by_index[4] == "solve_cubic_f64");
    CHECK(names_by_index[5] == "solve_cubic_real50");
    CHECK(names_by_index[6] == "newton_root");
    CHECK(names_by_index[7] == "bisection_root");
    CHECK(names_by_index[8] == "ode_euler");
    CHECK(names_by_index[9] == "ode_rk4");
}

TEST_CASE("save_checkpoint + load_checkpoint round-trips weights exactly", "[rsc]") {
    rsc::Dispatcher model(rsc::kBaseFeatureDim, 8, rsc::kBaseNumOps, /*seed=*/3);
    auto snap = rsc::snapshot_base_registry();
    auto path = temp_checkpoint_path();

    rsc::save_checkpoint(path, model, "deadbeef", snap);
    auto cp = rsc::load_checkpoint(path);
    std::remove(path.c_str());

    REQUIRE(cp.w1.size() == model.w1().size());
    REQUIRE(cp.b1.size() == model.b1().size());
    REQUIRE(cp.w2.size() == model.w2().size());
    REQUIRE(cp.b2.size() == model.b2().size());
    for (std::size_t i = 0; i < cp.w1.size(); ++i) CHECK_THAT(cp.w1[i], WithinAbs(model.w1()[i], 1e-15));
    for (std::size_t i = 0; i < cp.b1.size(); ++i) CHECK_THAT(cp.b1[i], WithinAbs(model.b1()[i], 1e-15));
    for (std::size_t i = 0; i < cp.w2.size(); ++i) CHECK_THAT(cp.w2[i], WithinAbs(model.w2()[i], 1e-15));
    for (std::size_t i = 0; i < cp.b2.size(); ++i) CHECK_THAT(cp.b2[i], WithinAbs(model.b2()[i], 1e-15));
}

TEST_CASE("checkpoint_to_dispatcher reconstructs a behaviorally identical model", "[rsc]") {
    rsc::Dispatcher model(rsc::kBaseFeatureDim, 8, rsc::kBaseNumOps, /*seed=*/7);
    auto snap = rsc::snapshot_base_registry();
    auto path = temp_checkpoint_path();

    rsc::save_checkpoint(path, model, "deadbeef", snap);
    auto cp = rsc::load_checkpoint(path);
    std::remove(path.c_str());
    auto restored = rsc::checkpoint_to_dispatcher(cp);

    std::vector<double> x(rsc::kBaseFeatureDim, 0.0);
    x[0] = 1.0;
    x[5] = 2.5;
    auto original_logits = model.forward(x);
    auto restored_logits = restored.forward(x);
    REQUIRE(original_logits.size() == restored_logits.size());
    for (std::size_t i = 0; i < original_logits.size(); ++i)
        CHECK_THAT(restored_logits[i], WithinAbs(original_logits[i], 1e-12));
}

TEST_CASE("validate_checkpoint catches commit and registry mismatches, "
          "not just confirms the happy path",
          "[rsc]") {
    rsc::Dispatcher model(rsc::kBaseFeatureDim, 8, rsc::kBaseNumOps, /*seed=*/1);
    auto snap = rsc::snapshot_base_registry();
    auto path = temp_checkpoint_path();

    rsc::save_checkpoint(path, model, "abc123", snap);
    auto cp = rsc::load_checkpoint(path);
    std::remove(path.c_str());

    CHECK(rsc::validate_checkpoint(cp, "abc123", snap) == rsc::ValidationResult::Match);
    CHECK(rsc::validate_checkpoint(cp, "different-sha", snap) ==
          rsc::ValidationResult::CommitMismatch);

    // Simulate an index reshuffle: swap two ops' global_index. Same set
    // of names, same arities, different assignment -- exactly the case
    // rsc/README.md's Reproducibility section calls out by name ("index
    // reshuffles must be caught too, not just behavioral changes").
    auto reshuffled = snap;
    std::swap(reshuffled[0].global_index, reshuffled[1].global_index);
    CHECK(rsc::validate_checkpoint(cp, "abc123", reshuffled) ==
          rsc::ValidationResult::RegistryMismatch);
}
