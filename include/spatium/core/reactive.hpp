#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <cstddef>
#  include <cstdint>
#  include <optional>
#  include <type_traits>
#  include <utility>
#  include <vector>
#endif

// `Cell<T>` and `Memo<F>` — a pull-based reactive layer.
//
// A `Cell` is an input: `set` it, `get` it. A `Memo` is a value derived
// from other cells and memos by a function. Nothing is pushed anywhere:
// setting a cell only bumps its version, and a memo looks at the versions
// of what it read the next time somebody reads *it*. A memo nobody reads
// never runs.
//
// Dependencies are not declared, they are observed. While a memo's function
// runs, a thread-local pointer names the memo, and every `get` made by the
// function registers itself there. The list is rebuilt on each run, so a
// branch that stops reading a cell stops depending on it.
//
// ── Early cutoff ─────────────────────────────────────────────────
//
// A memo bumps its own version only when the recomputed value differs
// (by `==`) from the cached one. A chain a -> b -> c therefore stops at b
// when b came out the same, and c does not rerun.
//
// ── What it does not do ──────────────────────────────────────────
//
// Nodes hold each other's addresses, so they are neither copyable nor
// movable. Single-threaded per graph: the observer pointer is thread-local,
// the versions are not atomic. A cycle (a memo reading itself) is a
// programming error and is not detected. No virtual calls: a dependency is
// a node address plus a function pointer that yields its fresh version.

SPATIUM_EXPORT namespace spatium {

namespace reactive {

// What a running memo remembers about one thing it read.
struct Dep {
    const void* node;
    std::uint64_t seen;
    std::uint64_t (*version)(const void*);  // refreshes the node, then reads it
};

// The memo whose function is running on this thread, as its dependency list.
inline thread_local std::vector<Dep>* current_ = nullptr;

inline void record(const void* node, std::uint64_t seen,
                   std::uint64_t (*version)(const void*)) {
    if (!current_) return;
    for (const Dep& d : *current_)
        if (d.node == node) return;
    current_->push_back({node, seen, version});
}

}  // namespace reactive

template<typename T>
class Cell {
public:
    explicit Cell(T v) : value_(std::move(v)) {}
    Cell(const Cell&) = delete;
    Cell& operator=(const Cell&) = delete;

    // Bumps the version only on an actual change.
    void set(const T& v) {
        if (v == value_) return;
        value_ = v;
        ++version_;
    }

    const T& get() const {
        reactive::record(this, version_, &version_of);
        return value_;
    }

    std::uint64_t version() const noexcept { return version_; }

private:
    static std::uint64_t version_of(const void* p) {
        return static_cast<const Cell*>(p)->version_;
    }

    T value_;
    std::uint64_t version_ = 1;
};

template<typename F>
class Memo {
public:
    using value_type = std::decay_t<std::invoke_result_t<F&>>;

    explicit Memo(F f) : f_(std::move(f)) {}
    Memo(const Memo&) = delete;
    Memo& operator=(const Memo&) = delete;

    const value_type& get() const {
        refresh();
        reactive::record(this, version_, &version_of);
        return *value_;
    }

    // How many times the function has run; for tests and instrumentation.
    std::size_t runs() const noexcept { return runs_; }
    std::uint64_t version() const noexcept { return version_; }

private:
    static std::uint64_t version_of(const void* p) {
        const auto* m = static_cast<const Memo*>(p);
        m->refresh();
        return m->version_;
    }

    bool stale() const {
        if (!value_) return true;
        for (const reactive::Dep& d : deps_)
            if (d.version(d.node) != d.seen) return true;
        return false;
    }

    void refresh() const {
        if (!stale()) return;
        std::vector<reactive::Dep> fresh;
        auto* outer = std::exchange(reactive::current_, &fresh);
        struct Restore {
            std::vector<reactive::Dep>* outer;
            ~Restore() { reactive::current_ = outer; }
        } restore{outer};
        value_type next = f_();
        ++runs_;
        deps_ = std::move(fresh);
        // Early cutoff: an equal result leaves the version, and so every
        // dependent, alone.
        if (!value_ || !(*value_ == next)) {
            value_ = std::move(next);
            ++version_;
        }
    }

    mutable F f_;
    mutable std::optional<value_type> value_;
    mutable std::vector<reactive::Dep> deps_;
    mutable std::uint64_t version_ = 0;
    mutable std::size_t runs_ = 0;
};

template<typename F>
Memo<F> computed(F f) { return Memo<F>(std::move(f)); }

}  // namespace spatium
