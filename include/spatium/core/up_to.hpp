#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/error.hpp>
#  include <array>
#  include <cassert>
#  include <cstddef>
#  include <initializer_list>
#  include <utility>
#endif

// `UpTo<T, N>` — at most N values, with a count, no allocation.
//
// The name states the thing that changed. Two places in this library
// returned a fixed-size `std::array` for an answer whose size is only
// *bounded*: a quadratic has at most two roots and a degenerate one has
// one, a ray meets a quadric in at most two points and usually fewer.
// Saying "exactly N" where the truth is "at most N" is not a cosmetic
// inaccuracy — it is what forced `solve_quadratic` to divide by a
// vanishing leading coefficient rather than report a shorter answer, and
// it is why `ray_quadric` reached for a `std::vector` and a heap
// allocation on every call to say the same thing.
//
// One container rather than two, deliberately: the roots case wants
// `Complex<T>` and the hits case wants `RayHit<T>`, but the contract is
// identical and writing it twice is how two nearly-identical types drift.
//
// ── The access contract, which is the part worth reading ─────────
//
// `operator[]` is **unchecked under NDEBUG**: no branch in a hot path, so
// indexing at or past `size()` is undefined there. This has to be said
// loudly, because a reader arriving from `std::array` expects a fixed
// size in which every index is always valid, and that is precisely what
// stops being true. In Debug it carries an assert.
//
// `at(i)` returns `Result<T>` and **never throws**. Checked access
// signals; unchecked access does not. Different contracts, so different
// return types — and exceptions appear nowhere else in this library's
// numerics, so they must not appear here.
//
// Range-`for` takes its bounds from the count, so every existing
// `for (auto& r : roots)` call site keeps working unchanged and silently
// starts iterating the right number of elements.
//
// ── Poisoning ────────────────────────────────────────────────────
//
// In Debug, the slots between `size()` and `capacity()` are overwritten
// via an ADL `debug_poison(T&)` when the element type provides one —
// `Complex<T>` does, with NaN. Zero instructions under NDEBUG, and under
// the Debug CI job it turns reading past the count from undefined into
// visibly wrong. It is a debugging aid and never part of the contract,
// the same category as `_GLIBCXX_ASSERTIONS`.

SPATIUM_EXPORT namespace spatium {

// The fallback: a type with nothing meaningful to poison is left alone.
// Overload `debug_poison` in your own type's namespace to opt in — a more
// specialized overload wins over this one, which is the whole mechanism.
constexpr void debug_poison(auto&) noexcept {}

template<typename T, std::size_t N>
class UpTo {
public:
    using value_type = T;

    constexpr UpTo() { poison_from(0); }

    // Built from a braced list so a solver can still `return {a, b}` and
    // read the way it did when the answer was always full length.
    constexpr UpTo(std::initializer_list<T> init) {
        assert(init.size() <= N && "UpTo: more values than capacity");
        for (const auto& v : init) data_[size_++] = v;
        poison_from(size_);
    }

    constexpr std::size_t size() const noexcept { return size_; }
    static constexpr std::size_t capacity() noexcept { return N; }
    constexpr bool empty() const noexcept { return size_ == 0; }

    constexpr void push_back(const T& v) {
        assert(size_ < N && "UpTo: push_back past capacity");
        data_[size_++] = v;
    }

    constexpr void clear() noexcept {
        size_ = 0;
        poison_from(0);
    }

    // Unchecked under NDEBUG. See the header comment.
    constexpr T& operator[](std::size_t i) {
        assert(i < size_ && "UpTo: index past the count, not past the capacity");
        return data_[i];
    }
    constexpr const T& operator[](std::size_t i) const {
        assert(i < size_ && "UpTo: index past the count, not past the capacity");
        return data_[i];
    }

    // Checked, and never throws.
    constexpr Result<T> at(std::size_t i) const {
        if (i >= size_)
            return std::unexpected(Error{ErrorCode::OutOfDomain,
                                         "UpTo: index past the count"});
        return data_[i];
    }

    constexpr T* begin() noexcept { return data_.data(); }
    constexpr T* end() noexcept { return data_.data() + size_; }
    constexpr const T* begin() const noexcept { return data_.data(); }
    constexpr const T* end() const noexcept { return data_.data() + size_; }

    constexpr T& front() { return (*this)[0]; }
    constexpr const T& front() const { return (*this)[0]; }
    constexpr T& back() { return (*this)[size_ - 1]; }
    constexpr const T& back() const { return (*this)[size_ - 1]; }

private:
    constexpr void poison_from([[maybe_unused]] std::size_t from) {
#ifndef NDEBUG
        for (std::size_t i = from; i < N; ++i) debug_poison(data_[i]);
#endif
    }

    std::array<T, N> data_{};
    std::size_t size_ = 0;
};

} // namespace spatium
