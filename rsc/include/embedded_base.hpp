#pragma once

// Base weights compiled into the binary -- no runtime file I/O to get a
// working base dispatcher. The generated data header
// (embedded_base_data.hpp, from rsc/tools/embed_checkpoint.py, wired in
// at CMake configure time) is either real data (SPATIUM_HAS_EMBEDDED_BASE
// 1) or an empty stub (0) if no checkpoint existed yet when cmake last
// ran -- see the top-level CMakeLists.txt comment for why. The custom
// (per-deployment) layer stays real file I/O on purpose -- it has to
// update in the field without recompiling; see custom_layer.hpp.

#include <build_info.hpp>
#include <dispatcher.hpp>
#include <embedded_base_data.hpp>
#include <string_view>

namespace rsc {

#if SPATIUM_HAS_EMBEDDED_BASE

inline Dispatcher load_embedded_base() {
    Dispatcher model(kEmbeddedBaseInputDim, kEmbeddedBaseHiddenDim, kEmbeddedBaseNumOps);
    model.w1().assign(kEmbeddedBaseW1.begin(), kEmbeddedBaseW1.end());
    model.b1().assign(kEmbeddedBaseB1.begin(), kEmbeddedBaseB1.end());
    model.w2().assign(kEmbeddedBaseW2.begin(), kEmbeddedBaseW2.end());
    model.b2().assign(kEmbeddedBaseB2.begin(), kEmbeddedBaseB2.end());
    return model;
}

// The embedded weights were baked in from a checkpoint pinned to
// kEmbeddedBaseCommitSha, at whatever earlier configure produced
// rsc/checkpoints/base_v1.checkpoint. kSpatiumCommitSha (build_info.hpp)
// is captured fresh on *every* configure from the current HEAD. If they
// differ, the checked-out code has moved past what the embedded base was
// actually trained against -- flag it, don't silently keep serving a
// stale base.
inline bool embedded_base_is_current() {
    return std::string_view(kEmbeddedBaseCommitSha) == std::string_view(kSpatiumCommitSha);
}

#endif // SPATIUM_HAS_EMBEDDED_BASE

} // namespace rsc
