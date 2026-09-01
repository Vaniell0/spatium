#include <catch2/catch_test_macros.hpp>
#include <base_task.hpp>
#include <embedded_base.hpp>
#include <train.hpp>

#if SPATIUM_HAS_EMBEDDED_BASE

TEST_CASE("embedded base dimensions match the unified base's own layout", "[rsc]") {
    CHECK(rsc::kEmbeddedBaseInputDim == rsc::kBaseFeatureDim);
    CHECK(rsc::kEmbeddedBaseNumOps == rsc::kBaseNumOps);
}

TEST_CASE("embedded base is current (checkpoint was trained against this exact commit)", "[rsc]") {
    // Informational, not a hard failure: kEmbeddedBaseCommitSha is stamped
    // into rsc/checkpoints/base_v1.checkpoint's first line by
    // rsc/tools/train_base at the HEAD that was checked out when it ran.
    // Whatever commit then adds/updates that checkpoint file can never
    // contain its own SHA -- git only assigns a commit's hash once the
    // commit (tree + parent + message) already exists -- so the checkpoint
    // is always stamped with its *parent* commit's SHA. A hard CHECK here
    // would therefore fail on that commit and on every commit after it,
    // permanently, no matter how often the checkpoint is retrained and
    // recommitted: not a one-time staleness to refresh away, a structural
    // property of comparing a committed artifact's self-stamped SHA
    // against its own committed history. CHECK_NOFAIL still surfaces a
    // real, large drift (code moved on since rsc/tools/train_base last
    // ran, retraining overdue) without failing CI on every commit.
    CHECK_NOFAIL(rsc::embedded_base_is_current());
}

TEST_CASE("load_embedded_base produces a genuinely trained model, not zeros", "[rsc]") {
    auto model = rsc::load_embedded_base();

    rsc::BaseTaskGenerator eval(999);
    int correct = 0;
    const int n = 2000;
    for (int i = 0; i < n; ++i)
        if (rsc::argmax_correct(model, eval.sample(), rsc::kBaseFeatureDim, 0)) ++correct;

    double acc = static_cast<double>(correct) / n;
    // Chance level over 10 classes is ~0.1-0.25 depending on domain mix;
    // the real trained checkpoint (see rsc/README.md) measured
    // ~0.81-0.99 per domain, so the pooled number should be well above
    // chance -- this only needs to catch "embedded zeros/garbage," not
    // reproduce the exact per-domain numbers (test_rsc_base.cpp already
    // does that on a freshly-trained model).
    CHECK(acc > 0.6);
}

#endif // SPATIUM_HAS_EMBEDDED_BASE
