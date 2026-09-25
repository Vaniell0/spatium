#pragma once
// An RGB image a renderer samples by (u, v) -- the "sampling call" half of
// the ROADMAP's "a texture is a field plus a sampling call". The field
// half, a Material or DSL field that reads one, is the next step and is
// not here: this is the image and the lookup, with the conventions a UV
// map needs pinned down once.
//
// Conventions, because each one is a place two callers otherwise disagree:
//   - (u, v) in [0, 1] covers the image once; u runs along the width, v
//     along the height, v = 0 at the image's first row;
//   - bilinear between texel centres, so a texture sampled at its own
//     texel grid returns its texels exactly;
//   - outside [0, 1] the image repeats or clamps, as asked.
//
// Loading goes through stb_image, vendored beside stb_image_write. As with
// that one, a caller defines STB_IMAGE_IMPLEMENTATION in exactly one
// translation unit before including this header.
//
// Not part of the `spatium.render` module: it includes a vendored C
// header, like `render/write_image.hpp`.
#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/error.hpp>
#  include <spatium/vendor/stb_image.h>
#  include <algorithm>
#  include <cmath>
#  include <cstdint>
#  include <filesystem>
#  include <string>
#  include <vector>
#endif

namespace spatium::render {

enum class Wrap : std::uint8_t { Repeat, Clamp };

struct Texture {
    int width = 0, height = 0;
    std::vector<std::uint8_t> rgb;   // width * height * 3, row-major, first row first

    // A texture from a generator `f(x, y) -> Vec<double, 3>` in 0..255, one
    // call per texel -- how a procedural pattern becomes an image that can
    // be written out, looked at and painted over.
    template<typename F>
    static Texture generate(int w, int h, F&& f) {
        Texture t{w, h, std::vector<std::uint8_t>(static_cast<std::size_t>(w) * h * 3)};
        for (int y = 0; y < h; ++y)
            for (int x = 0; x < w; ++x) {
                const Vec<double, 3> c = f(x, y);
                auto* px = &t.rgb[3 * (static_cast<std::size_t>(y) * w + x)];
                for (int k = 0; k < 3; ++k)
                    px[k] = static_cast<std::uint8_t>(std::clamp(c[static_cast<std::size_t>(k)], 0.0, 255.0));
            }
        return t;
    }

    static Result<Texture> load(const std::filesystem::path& path) {
        int w = 0, h = 0, channels = 0;
        std::uint8_t* data = stbi_load(path.string().c_str(), &w, &h, &channels, 3);
        if (!data)
            return std::unexpected(Error(ErrorCode::InvalidArgument,
                                         "texture: cannot read " + path.string() + ": " +
                                             (stbi_failure_reason() ? stbi_failure_reason() : "")));
        Texture t{w, h, std::vector<std::uint8_t>(data, data + static_cast<std::size_t>(w) * h * 3)};
        stbi_image_free(data);
        return t;
    }

    Vec<double, 3> texel(int x, int y) const {
        const auto* px = &rgb[3 * (static_cast<std::size_t>(y) * width + x)];
        return {static_cast<double>(px[0]), static_cast<double>(px[1]), static_cast<double>(px[2])};
    }

    // Colour at (u, v), 0..255 per channel.
    Vec<double, 3> sample(double u, double v, Wrap wrap = Wrap::Repeat) const {
        if (width <= 0 || height <= 0) return {};
        // Texel centres sit at (i + 0.5) / width; x below is in texel units
        // measured from the first centre.
        double x = u * width - 0.5, y = v * height - 0.5;
        const double fx = x - std::floor(x), fy = y - std::floor(y);
        const int x0 = static_cast<int>(std::floor(x)), y0 = static_cast<int>(std::floor(y));
        auto at = [&](int i, int n) {
            if (wrap == Wrap::Repeat) return ((i % n) + n) % n;
            return std::clamp(i, 0, n - 1);
        };
        const int xa = at(x0, width), xb = at(x0 + 1, width);
        const int ya = at(y0, height), yb = at(y0 + 1, height);
        const Vec<double, 3> top{texel(xa, ya) * (1.0 - fx) + texel(xb, ya) * fx};
        const Vec<double, 3> bottom{texel(xa, yb) * (1.0 - fx) + texel(xb, yb) * fx};
        return Vec<double, 3>{top * (1.0 - fy) + bottom * fy};
    }
};

}  // namespace spatium::render
