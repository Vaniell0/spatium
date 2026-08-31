#pragma once

// Shared colour-wheel palette for atomic visualisations.
//
// Bit-reversal index hashing: indices 0, 1, 2, 3, … land on
// 0, ½, ¼, ¾, ⅛, ⅝, … of the hue wheel, so adjacent slots are
// maximally contrasting rather than rainbow-adjacent. Used by both
// the quantum point-cloud `AtomModel` (per-subshell colour) and the
// classical `BohrModel` (per-shell colour). Lives in its own header
// so neither file silently depends on the other being included
// first.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <algorithm>
#  include <cmath>
#endif

SPATIUM_EXPORT namespace spatium::physics::atomic {

inline Vec4f orbital_palette(int idx, int total, float alpha = 0.28f) {
    int n = std::max(total, 1);
    unsigned br = 0;
    unsigned v = static_cast<unsigned>(((idx % n) + n) % n);
    unsigned bits = 0;
    for (unsigned t = static_cast<unsigned>(n - 1); t > 0; t >>= 1) ++bits;
    for (unsigned b = 0; b < bits; ++b) {
        br = (br << 1) | (v & 1u);
        v >>= 1;
    }
    float h = std::fmod(static_cast<float>(br) / static_cast<float>(1u << bits) + 0.08f, 1.0f);
    constexpr float s = 1.0f;
    constexpr float val = 1.0f;
    float h6 = h * 6.0f;
    float c = val * s;
    float x = c * (1.0f - std::abs(std::fmod(h6, 2.0f) - 1.0f));
    float m = val - c;
    float r = 0, g = 0, b = 0;
    int hi = static_cast<int>(std::floor(h6)) % 6;
    if (hi < 0) hi += 6;
    switch (hi) {
        case 0: r = c; g = x; b = 0; break;
        case 1: r = x; g = c; b = 0; break;
        case 2: r = 0; g = c; b = x; break;
        case 3: r = 0; g = x; b = c; break;
        case 4: r = x; g = 0; b = c; break;
        default: r = c; g = 0; b = x; break;
    }
    return {r + m, g + m, b + m, alpha};
}

} // namespace spatium::physics::atomic
