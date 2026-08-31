#pragma once

// Blackbody temperature -> sRGB approximation (Kelvin -> [0,255] RGB).
// A piecewise polynomial fit to the Planckian locus (the widely-used
// graphics approximation popularized by Tanner Helland's blog post
// deriving it from CIE blackbody data), not a spectral CIE integration
// -- good enough for tinting thermal emitters like an accretion disk or
// a star, not colorimetrically exact. Valid over roughly [1000K,
// 40000K]; clamped outside that range. No such utility existed anywhere
// in the tree before this (confirmed by a full-repo search during the
// 2026-08-26 GR-rewrite research pass) -- the only prior color-science
// code was examples/io_helpers.hpp's example-local HSV-to-RGB, unrelated
// to physical temperature.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <algorithm>
#  include <cmath>
#endif

SPATIUM_EXPORT namespace spatium::render {

inline Vec<double, 3> blackbody_to_rgb255(double kelvin) {
    using std::clamp, std::pow, std::log;
    double temp = clamp(kelvin, 1000.0, 40000.0) / 100.0;

    double red = (temp <= 66.0)
        ? 255.0
        : clamp(329.698727446 * pow(temp - 60.0, -0.1332047592), 0.0, 255.0);

    double green = (temp <= 66.0)
        ? clamp(99.4708025861 * log(temp) - 161.1195681661, 0.0, 255.0)
        : clamp(288.1221695283 * pow(temp - 60.0, -0.0755148492), 0.0, 255.0);

    double blue;
    if (temp >= 66.0) blue = 255.0;
    else if (temp <= 19.0) blue = 0.0;
    else blue = clamp(138.5177312231 * log(temp - 10.0) - 305.0447927307, 0.0, 255.0);

    return {red, green, blue};
}

} // namespace spatium::render
