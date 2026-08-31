#pragma once

// The first few Tier-1 (library-independent) ops from the RSC design doc's
// registry: plain arithmetic plus the vector free functions Spatium already
// has in algebra/functions.hpp. Small on purpose — this exists to prove the
// Registry API end to end, not to cover Tier 1 exhaustively yet.

#include <registry.hpp>
#include <spatium/algebra/functions.hpp>
#include <spatium/algebra/vector.hpp>

namespace rsc {

inline Registry build_tier1_registry() {
    Registry reg;

    reg.add({.name = "add", .tier = Tier::General, .in_size = 2, .out_size = 1,
             .input_names = {"a", "b"}, .output_names = {"sum"}},
            [](std::span<const double> in, std::span<double> out) {
                out[0] = in[0] + in[1];
            });

    reg.add({.name = "multiply", .tier = Tier::General, .in_size = 2, .out_size = 1,
             .input_names = {"a", "b"}, .output_names = {"product"}},
            [](std::span<const double> in, std::span<double> out) {
                out[0] = in[0] * in[1];
            });

    // 3-vectors flattened into spans: a = in[0..2], b = in[3..5].
    reg.add({.name = "dot3", .tier = Tier::General, .in_size = 6, .out_size = 1,
             .input_names = {"a.x", "a.y", "a.z", "b.x", "b.y", "b.z"},
             .output_names = {"dot"}},
            [](std::span<const double> in, std::span<double> out) {
                spatium::Vec<double, 3> a{in[0], in[1], in[2]};
                spatium::Vec<double, 3> b{in[3], in[4], in[5]};
                out[0] = spatium::dot(a, b);
            });

    // a = in[0..2], b = in[3..5], t = in[6].
    reg.add({.name = "lerp3", .tier = Tier::General, .in_size = 7, .out_size = 3,
             .input_names = {"a.x", "a.y", "a.z", "b.x", "b.y", "b.z", "t"},
             .output_names = {"x", "y", "z"}},
            [](std::span<const double> in, std::span<double> out) {
                spatium::Vec<double, 3> a{in[0], in[1], in[2]};
                spatium::Vec<double, 3> b{in[3], in[4], in[5]};
                double t = in[6];
                auto r = spatium::lerp(a, b, t);
                out[0] = r[0];
                out[1] = r[1];
                out[2] = r[2];
            });

    return reg;
}

} // namespace rsc
