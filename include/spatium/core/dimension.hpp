#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <cstddef>
#endif

SPATIUM_EXPORT namespace spatium {

// Sentinels of the old convention, where a dimension was a std::size_t.
inline constexpr std::size_t kDynamic = static_cast<std::size_t>(-1);    // known at run time
inline constexpr std::size_t kInfinite = static_cast<std::size_t>(-2);   // no finite basis

// The dimension of a space as a value that says which kind it is, not a number that
// has to be read against sentinels:
//
//   Finite(N)  N coordinates, stored in a Vec or a matrix: the spaces that have a basis
//   Dynamic    a number fixed only at run time
//   Infinite   no finite basis -- a space of functions
//
// Only a space that stores a basis needs a `Finite` one (`Vec`, a matrix, a Euclidean
// space); an algorithm written against exp, log and a metric never reads the dimension at
// all, which is what lets a space of functions be a space. A `std::size_t` still converts,
// both ways -- `Dimension(3)` is Finite(3), `kDynamic` is Dynamic, `kInfinite` is Infinite,
// and a Dimension reads back as the number (or the sentinel) -- so every `static constexpr
// std::size_t dimension = N;` written before this existed is still a dimension, unchanged.
// Sums follow what they mean: an infinite factor makes a product infinite, a dynamic one
// makes it dynamic, finite ones add.
struct Dimension {
    enum class Kind { Finite, Dynamic, Infinite };

    constexpr Dimension() = default;
    // N, kDynamic and kInfinite are told apart here, once.
    constexpr Dimension(std::size_t n) : n_(n) {}           // NOLINT: the conversion is the point

    static constexpr Dimension finite(std::size_t n) { return Dimension(n); }
    static constexpr Dimension dynamic() { return Dimension(kDynamic); }
    static constexpr Dimension infinite() { return Dimension(kInfinite); }

    constexpr Kind kind() const { return n_ == kDynamic ? Kind::Dynamic : n_ == kInfinite ? Kind::Infinite : Kind::Finite; }
    constexpr bool is_finite() const { return kind() == Kind::Finite; }
    constexpr bool is_dynamic() const { return kind() == Kind::Dynamic; }
    constexpr bool is_infinite() const { return kind() == Kind::Infinite; }

    // The count of a finite dimension; the sentinel otherwise (so old comparisons with
    // kDynamic still read as they did).
    constexpr operator std::size_t() const { return n_; }

    friend constexpr Dimension operator+(Dimension a, Dimension b) {
        if (a.is_infinite() || b.is_infinite()) return infinite();
        if (a.is_dynamic() || b.is_dynamic()) return dynamic();
        return finite(a.n_ + b.n_);
    }

private:
    std::size_t n_ = 0;
};

} // namespace spatium
