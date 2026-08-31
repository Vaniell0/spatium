#pragma once

// Human-readable formatting for op calls -- "add(a=2, b=3) -> sum=5" instead
// of raw indices, per the original design doc's "показывает шаги решения"
// goal. Purely descriptive: reads OpSignature's names if present, falls
// back to positional labels otherwise. This is the layer a future DSL
// front-end produces/consumes; features() (the fixed-size numeric encoding
// actually fed to the dispatcher) is untouched by anything here.

#include <registry.hpp>
#include <span>
#include <sstream>
#include <string>

namespace rsc {

namespace detail {

inline std::string slot_name(const std::vector<std::string>& names, std::size_t i,
                              const char* fallback_prefix) {
    if (i < names.size()) return names[i];
    return std::string(fallback_prefix) + std::to_string(i);
}

} // namespace detail

inline std::string describe(const OpSignature& sig, std::span<const double> in,
                             std::span<const double> out) {
    std::ostringstream os;
    os << sig.name << "(";
    for (std::size_t i = 0; i < in.size(); ++i) {
        if (i) os << ", ";
        os << detail::slot_name(sig.input_names, i, "in") << "=" << in[i];
    }
    os << ") -> ";
    for (std::size_t i = 0; i < out.size(); ++i) {
        if (i) os << ", ";
        os << detail::slot_name(sig.output_names, i, "out") << "=" << out[i];
    }
    return os.str();
}

} // namespace rsc
