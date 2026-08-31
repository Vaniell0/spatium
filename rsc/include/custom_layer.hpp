#pragma once

// Per-deployment custom layer: cached calibration coefficients keyed by
// scene/deployment class, on top of the frozen embedded base
// (embedded_base.hpp). No dispatcher weight changes here -- this is
// purely the "calibration-search-with-cache" mechanism from
// rsc/README.md's domain pipeline (e.g. XPBD contact params fit via
// calculus.hpp's minimize() for one scene, cached so it's a lookup, not
// a re-solve, next time). Real file I/O, deliberately unlike the
// embedded base: this has to update in the field without recompiling.
//
// No concrete calibration scenario is wired to this yet -- contact
// physics (the original motivating case) is blocked pending an
// ipc-toolkit backend swap (see the domain pipeline section). This is
// the generic mechanism, tested against a synthetic entry, honestly not
// yet exercised by a real one.

#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace rsc {

struct CustomLayer {
    // Which base this calibration was computed against -- a base change
    // (new commit, retrained weights) invalidates every entry here; see
    // custom_layer_matches_base().
    std::string base_commit_sha;
    std::unordered_map<std::string, std::vector<double>> params;
};

inline bool custom_layer_matches_base(const CustomLayer& layer, const std::string& current_base_sha) {
    return layer.base_commit_sha == current_base_sha;
}

inline void save_custom_layer(const std::string& path, const CustomLayer& layer) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("save_custom_layer: cannot open '" + path + "' for writing");
    f << std::setprecision(std::numeric_limits<double>::max_digits10);

    f << layer.base_commit_sha << '\n';
    f << layer.params.size() << '\n';
    for (auto& [key, values] : layer.params) {
        f << key << ' ' << values.size();
        for (double v : values) f << ' ' << v;
        f << '\n';
    }
}

inline CustomLayer load_custom_layer(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("load_custom_layer: cannot open '" + path + "' for reading");

    CustomLayer layer;
    std::getline(f, layer.base_commit_sha);

    std::size_t n_entries;
    f >> n_entries;
    for (std::size_t i = 0; i < n_entries; ++i) {
        std::string key;
        std::size_t n_values;
        f >> key >> n_values;
        std::vector<double> values(n_values);
        for (auto& v : values) f >> v;
        layer.params[key] = std::move(values);
    }
    return layer;
}

} // namespace rsc
