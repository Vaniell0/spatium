// train_base — produces the actual base checkpoint, not just the
// mechanism to save/load one. Config (hidden=64, entropy_beta=0.2,
// 8000 batched updates, per-domain REINFORCE baselines) is the one
// measured in tests/test_rsc_base.cpp to raise every domain to
// ~0.81-0.99 accuracy -- this tool just runs it for real and writes the
// result to disk, pinned to the exact commit + registry layout it was
// trained against (checkpoint.hpp).
//
// Usage: train_base [output_path]  (defaults to rsc_base.checkpoint)

#include <base_task.hpp>
#include <checkpoint.hpp>
#include <train.hpp>
#include <array>
#include <cstdio>
#include <print>

int main(int argc, char** argv) {
    std::string path = argc > 1 ? argv[1] : "rsc_base.checkpoint";

    rsc::Dispatcher model(rsc::kBaseFeatureDim, 64, rsc::kBaseNumOps, /*seed=*/0);
    rsc::BaseTaskGenerator train_gen(1);
    std::mt19937_64 rng(42);
    std::array<double, 4> baselines{0.5, 0.5, 0.5, 0.5};

    std::println("training base checkpoint: hidden=64, entropy_beta=0.2, 8000 updates x 32 batch");
    for (int u = 0; u < 8000; ++u) {
        rsc::Gradients sum;
        for (int b = 0; b < 32; ++b) {
            auto task = train_gen.sample();
            auto domain = rsc::domain_of(task);
            bool ok;
            auto g = rsc::reinforce_gradient(model, task, rsc::kBaseFeatureDim, 0, rng,
                                              baselines[static_cast<std::size_t>(domain)], ok,
                                              0.01, 0.2);
            rsc::accumulate_gradients(sum, g);
        }
        rsc::scale_gradients(sum, 1.0 / 32.0);
        model.apply_gradients(sum, 0.1);
    }

    auto snap = rsc::snapshot_base_registry();
    rsc::save_checkpoint(path, model, rsc::kSpatiumCommitSha, snap);
    std::println("saved to {} (commit {})", path, rsc::kSpatiumCommitSha);

    // Load it straight back and validate, as a self-check that the
    // round trip this tool depends on actually works before reporting
    // success.
    auto cp = rsc::load_checkpoint(path);
    auto result = rsc::validate_checkpoint(cp, rsc::kSpatiumCommitSha, snap);
    if (result != rsc::ValidationResult::Match) {
        std::println("ERROR: checkpoint failed self-validation immediately after saving");
        return 1;
    }
    std::println("validated: commit + registry match on reload");
    return 0;
}
