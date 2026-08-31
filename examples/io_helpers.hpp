// Output guard for example demos.
//
// Examples write SVG/PNG artefacts into the working directory. To avoid
// silently clobbering hand-curated outputs the helper here gates every write
// behind an existence check: if the target file is already present and the
// caller did not pass --force, the write is skipped with a stderr warning.
//
// Usage:
//   namespace fs = std::filesystem;
//   if (spatium::examples::confirm_overwrite("out.svg", force))
//       svg.save("out.svg");

#pragma once

#include <filesystem>
#if __has_include(<print>)
#  include <print>
#  define SPATIUM_HAS_STD_PRINT 1
#else
// The CUDA build's target toolchain (g++ 13 on the render VM) doesn't
// ship <print> yet -- fall back to cstdio rather than requiring every
// checkout to hand-patch this header after rsync.
#  include <cstdio>
#  define SPATIUM_HAS_STD_PRINT 0
#endif

namespace spatium::examples {

// Returns true iff the caller should proceed with the write.
//   - true  : file does not exist, or `force` was requested
//   - false : file exists and `force` is false; emits a stderr warning
inline bool confirm_overwrite(const std::filesystem::path& path, bool force) {
    std::error_code ec;
    if (force || !std::filesystem::exists(path, ec) || ec) return true;
#if SPATIUM_HAS_STD_PRINT
    std::println(stderr,
                 "warning: {} exists, skipping (pass --force to overwrite)",
                 path.string());
#else
    std::fprintf(stderr, "warning: %s exists, skipping (pass --force to overwrite)\n",
                 path.string().c_str());
#endif
    return false;
}

}  // namespace spatium::examples
