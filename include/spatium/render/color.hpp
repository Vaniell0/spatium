#pragma once

// HSV -> sRGB, the non-physical color-picker model (hue/saturation/value),
// as opposed to spectral.hpp's physically-based blackbody fit -- the two
// solve different problems (pick N visually-distinct colors vs. tint a
// thermal emitter) and shouldn't be conflated into one file.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <cmath>
#endif

SPATIUM_EXPORT namespace spatium::render {

// h in [0,1), s/v in [0,1] -> RGB in [0,255].
inline Vec<double, 3> hsv_to_rgb255(double h, double s, double v) {
    double c = v * s;
    double hp = h * 6.0;
    double x = c * (1.0 - std::abs(std::fmod(hp, 2.0) - 1.0));
    Vec<double, 3> rgb{0.0, 0.0, 0.0};
    if      (hp < 1.0) rgb = {c, x, 0.0};
    else if (hp < 2.0) rgb = {x, c, 0.0};
    else if (hp < 3.0) rgb = {0.0, c, x};
    else if (hp < 4.0) rgb = {0.0, x, c};
    else if (hp < 5.0) rgb = {x, 0.0, c};
    else               rgb = {c, 0.0, x};
    double m = v - c;
    return Vec<double, 3>{(rgb[0] + m) * 255.0, (rgb[1] + m) * 255.0, (rgb[2] + m) * 255.0};
}

} // namespace spatium::render
