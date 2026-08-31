#pragma once

// Thin PNG-writing wrapper.
//
// Every raytracer example called `stbi_write_png()` directly at its own
// call site (8 call sites across blackhole_demo.cpp, wormhole_demo.cpp,
// collatz_demo.cpp, wave_ca_demo.cpp, primitives_demo.cpp,
// parametric_analytical_demo.cpp) -- found duplicated during the
// 2026-08-26 architecture audit. The pattern was always the same
// (stride = width * channels) and error handling was inconsistent
// (some call sites checked the return value, some didn't). One place
// to get both right.
//
// Callers still `#define STB_IMAGE_WRITE_IMPLEMENTATION` in exactly one
// translation unit -- but should NOT also separately `#include
// <spatium/vendor/stb_image_write.h>` there anymore. This header is now
// that one include; doing both recompiles stb_image_write.h's
// implementation section twice in the same TU (it has no "already
// implemented" guard of its own, only the usual header guard around its
// declarations) and fails with redefinition errors.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/vendor/stb_image_write.h>
#  include <cstdint>
#  include <filesystem>
#  include <vector>
#  if __has_include(<print>)
#    include <print>
#    define SPATIUM_HAS_STD_PRINT 1
#  else
// See examples/io_helpers.hpp's matching guard: the CUDA render VM's
// g++ 13 doesn't ship <print> yet.
#    include <cstdio>
#    define SPATIUM_HAS_STD_PRINT 0
#  endif
#endif

SPATIUM_EXPORT namespace spatium::render {

inline bool write_png(const std::filesystem::path& path, int width, int height, int channels,
                       const std::uint8_t* data) {
    if (!stbi_write_png(path.string().c_str(), width, height, channels, data, width * channels)) {
#if SPATIUM_HAS_STD_PRINT
        std::print(stderr, "failed to write {}\n", path.string());
#else
        std::fprintf(stderr, "failed to write %s\n", path.string().c_str());
#endif
        return false;
    }
    return true;
}

inline bool write_png_rgb(const std::filesystem::path& path, int width, int height,
                           const std::vector<std::uint8_t>& rgb) {
    return write_png(path, width, height, 3, rgb.data());
}

inline bool write_png_rgba(const std::filesystem::path& path, int width, int height,
                            const std::vector<std::uint8_t>& rgba) {
    return write_png(path, width, height, 4, rgba.data());
}

}  // namespace spatium::render
