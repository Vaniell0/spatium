#pragma once

// Reproducibility, built for real (rsc/README.md's Reproducibility
// section designed this earlier; base_task.hpp finally produced a base
// checkpoint worth pinning). Ground truth here IS Spatium's own code
// (every domain's reward compares against a real Spatium computation),
// so a Spatium change alters both the training environment and the
// definition of "correct" at once -- loading a checkpoint recorded
// against a different commit or a reshuffled registry must never
// silently proceed.

#include <build_info.hpp>
#include <dispatcher.hpp>
#include <ode_ops.hpp>
#include <precision_ops.hpp>
#include <registry.hpp>
#include <rootfind_ops.hpp>
#include <tier1_ops.hpp>
#include <base_task.hpp>
#include <fstream>
#include <iomanip>
#include <ios>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace rsc {

struct RegistryOpEntry {
    std::string name;
    std::size_t global_index;
    std::size_t in_size;
    std::size_t out_size;

    bool operator==(const RegistryOpEntry&) const = default;
};

using RegistrySnapshot = std::vector<RegistryOpEntry>;

// Snapshot of every op across all four base domains at their global
// index offsets (kBaseOpOffset, base_task.hpp) -- catches index
// reshuffles (an op moved to a different index) as well as arity
// changes, not just "does an op with this name still exist somewhere."
inline RegistrySnapshot snapshot_base_registry() {
    RegistrySnapshot snap;
    auto add_domain = [&](const Registry& reg, std::size_t offset) {
        for (std::size_t i = 0; i < reg.size(); ++i) {
            const auto& sig = reg[i].signature();
            snap.push_back({sig.name, offset + i, sig.in_size, sig.out_size});
        }
    };
    add_domain(build_tier1_registry(), kBaseOpOffset[0]);
    add_domain(build_precision_registry(), kBaseOpOffset[1]);
    add_domain(build_rootfind_registry(), kBaseOpOffset[2]);
    add_domain(build_ode_registry(), kBaseOpOffset[3]);
    return snap;
}

struct BaseCheckpoint {
    std::string commit_sha;
    RegistrySnapshot registry;
    std::size_t input_dim, hidden_dim, num_ops;
    std::vector<double> w1, b1, w2, b2;
};

enum class ValidationResult { Match, CommitMismatch, RegistryMismatch };

inline ValidationResult validate_checkpoint(const BaseCheckpoint& cp,
                                             const std::string& current_sha,
                                             const RegistrySnapshot& current_registry) {
    if (cp.commit_sha != current_sha) return ValidationResult::CommitMismatch;
    if (cp.registry != current_registry) return ValidationResult::RegistryMismatch;
    return ValidationResult::Match;
}

// Plain-text format, one field group per line -- no JSON dependency
// anywhere in this codebase, and a base checkpoint is saved rarely
// (once per base training run), so human-readable/diffable matters more
// than compactness.
inline void save_checkpoint(const std::string& path, Dispatcher& model,
                             const std::string& commit_sha, const RegistrySnapshot& registry) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("save_checkpoint: cannot open '" + path + "' for writing");
    // max_digits10 -- exact round-trip for double, not the default 6
    // significant digits operator<< would otherwise use.
    f << std::setprecision(std::numeric_limits<double>::max_digits10);

    f << commit_sha << '\n';
    f << registry.size() << '\n';
    for (auto& e : registry)
        f << e.name << ' ' << e.global_index << ' ' << e.in_size << ' ' << e.out_size << '\n';

    f << model.input_dim() << ' ' << model.hidden_dim() << ' ' << model.num_ops() << '\n';
    auto write_vec = [&](const std::vector<double>& v) {
        f << v.size();
        for (double x : v) f << ' ' << x;
        f << '\n';
    };
    write_vec(model.w1());
    write_vec(model.b1());
    write_vec(model.w2());
    write_vec(model.b2());
}

inline BaseCheckpoint load_checkpoint(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("load_checkpoint: cannot open '" + path + "' for reading");

    BaseCheckpoint cp;
    std::getline(f, cp.commit_sha);

    std::size_t reg_size;
    f >> reg_size;
    cp.registry.resize(reg_size);
    for (auto& e : cp.registry) f >> e.name >> e.global_index >> e.in_size >> e.out_size;

    f >> cp.input_dim >> cp.hidden_dim >> cp.num_ops;
    auto read_vec = [&](std::vector<double>& v) {
        std::size_t n;
        f >> n;
        v.resize(n);
        for (auto& x : v) f >> x;
    };
    read_vec(cp.w1);
    read_vec(cp.b1);
    read_vec(cp.w2);
    read_vec(cp.b2);
    return cp;
}

// Reconstruct a usable Dispatcher from a checkpoint -- callers should
// check validate_checkpoint() first; this doesn't re-check, it just
// wires the loaded weights into a live model.
inline Dispatcher checkpoint_to_dispatcher(const BaseCheckpoint& cp) {
    Dispatcher model(cp.input_dim, cp.hidden_dim, cp.num_ops);
    model.w1() = cp.w1;
    model.b1() = cp.b1;
    model.w2() = cp.w2;
    model.b2() = cp.b2;
    return model;
}

} // namespace rsc
