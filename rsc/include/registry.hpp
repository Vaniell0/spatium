#pragma once

// The registry: the one source of truth for both the dispatch model's head
// size and the training task generator (see spatium/rsc/README.md). Every
// registered Op gets a stable index — that index IS the model's dispatch
// choice, no parser or name lookup on the hot path.
//
// v1 keeps the actual computation type-erased (std::function over spans of
// double) rather than templating Registry on a tuple of per-op signatures —
// simpler, and dispatch isn't the hot path yet. Revisit only if that cost
// actually shows up.

#include <algorithm>
#include <cstddef>
#include <functional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace rsc {

enum class Tier { General, Library };

struct OpSignature {
    std::string name;
    Tier tier;
    std::size_t in_size;
    std::size_t out_size;
    // Per-slot names, for the human/DSL-facing side only (describe(), below)
    // -- not required, not read by anything on the dispatch hot path. Empty
    // is fine; existing aggregate-init call sites that only set the four
    // fields above keep compiling unchanged.
    std::vector<std::string> input_names{};
    std::vector<std::string> output_names{};
};

class Op {
public:
    using Fn = std::function<void(std::span<const double> in, std::span<double> out)>;

    Op(OpSignature sig, Fn fn) : signature_(std::move(sig)), fn_(std::move(fn)) {}

    const OpSignature& signature() const { return signature_; }

    // Arity is checked here, not left to the caller or the wrapped Fn —
    // this is what makes a bad dispatch fail loudly and immediately instead
    // of reading past the end of a span.
    void operator()(std::span<const double> in, std::span<double> out) const {
        if (in.size() != signature_.in_size)
            throw std::invalid_argument(
                "rsc::Op '" + signature_.name + "': expected " +
                std::to_string(signature_.in_size) + " inputs, got " +
                std::to_string(in.size()));
        if (out.size() != signature_.out_size)
            throw std::invalid_argument(
                "rsc::Op '" + signature_.name + "': expected " +
                std::to_string(signature_.out_size) + " outputs, got " +
                std::to_string(out.size()));
        fn_(in, out);
    }

private:
    OpSignature signature_;
    Fn fn_;
};

class Registry {
public:
    // Returns the op's dispatch-head index — the value a trained model's
    // classification head would output to select this op.
    std::size_t add(OpSignature sig, Op::Fn fn) {
        ops_.emplace_back(std::move(sig), std::move(fn));
        return ops_.size() - 1;
    }

    std::size_t size() const { return ops_.size(); }
    const Op& operator[](std::size_t index) const { return ops_.at(index); }

    // The largest declared in/out arity across all registered ops -- lets a
    // fixed-size dispatcher pad any op's operands/output into one uniform
    // input vector.
    std::size_t max_in_size() const {
        std::size_t m = 0;
        for (const auto& op : ops_) m = std::max(m, op.signature().in_size);
        return m;
    }

    std::size_t max_out_size() const {
        std::size_t m = 0;
        for (const auto& op : ops_) m = std::max(m, op.signature().out_size);
        return m;
    }

    // Ops whose output can feed back as one of the next call's inputs
    // (in_size==2, out_size==1 -- a binary op producing exactly the kind of
    // value it also consumes). Derived from arity alone, not tagged by
    // hand: for the current Tier-1 registry this is exactly add/multiply;
    // dot3/lerp3 don't have this "state in, state out" shape and are
    // excluded automatically. What chains (multi-step composition) dispatch
    // among.
    std::vector<std::size_t> chainable_ops() const {
        std::vector<std::size_t> result;
        for (std::size_t i = 0; i < ops_.size(); ++i) {
            const auto& sig = ops_[i].signature();
            if (sig.in_size == 2 && sig.out_size == 1) result.push_back(i);
        }
        return result;
    }

    // Name lookup — for tests/debugging/tooling, not the dispatch hot path
    // (the model only ever produces an index, never a name).
    std::size_t index_of(std::string_view name) const {
        for (std::size_t i = 0; i < ops_.size(); ++i)
            if (ops_[i].signature().name == name) return i;
        throw std::invalid_argument("rsc::Registry: no op named '" + std::string(name) + "'");
    }

private:
    std::vector<Op> ops_;
};

} // namespace rsc
