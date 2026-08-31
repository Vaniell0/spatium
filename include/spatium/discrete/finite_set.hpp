#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/geometry/boolean.hpp>
#  include <spatium/geometry/polygon.hpp>
#  include <spatium/geometry/triangle.hpp>
#  include <spatium/geometry/circle.hpp>
#  include <algorithm>
#  include <functional>
#  include <initializer_list>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::discrete {

// ── FiniteSet<T> ──────────────────────────────────────────────
// A finite set with standard set operations.
// Elements are stored sorted and unique.

template<typename T>
struct FiniteSet {
    std::vector<T> elements;

    FiniteSet() = default;
    FiniteSet(std::initializer_list<T> init) : elements(init) { normalize(); }
    explicit FiniteSet(std::vector<T> v) : elements(std::move(v)) { normalize(); }

    std::size_t size() const { return elements.size(); }
    bool empty() const { return elements.empty(); }

    bool contains(const T& x) const {
        return std::binary_search(elements.begin(), elements.end(), x);
    }

    // ∪ Union — `+` matches the project-wide convention
    // (`|` is reserved for intersect/pipe, see CLAUDE.md).
    FiniteSet union_with(const FiniteSet& other) const {
        std::vector<T> result;
        std::set_union(elements.begin(), elements.end(),
                       other.elements.begin(), other.elements.end(),
                       std::back_inserter(result));
        return FiniteSet(std::move(result));
    }
    FiniteSet operator+(const FiniteSet& other) const {
        return union_with(other);
    }

    // ∩ Intersection
    FiniteSet operator&(const FiniteSet& other) const {
        std::vector<T> result;
        std::set_intersection(elements.begin(), elements.end(),
                              other.elements.begin(), other.elements.end(),
                              std::back_inserter(result));
        return FiniteSet(std::move(result));
    }

    // ∖ Difference
    FiniteSet operator-(const FiniteSet& other) const {
        std::vector<T> result;
        std::set_difference(elements.begin(), elements.end(),
                            other.elements.begin(), other.elements.end(),
                            std::back_inserter(result));
        return FiniteSet(std::move(result));
    }

    // △ Symmetric difference
    FiniteSet operator^(const FiniteSet& other) const {
        std::vector<T> result;
        std::set_symmetric_difference(elements.begin(), elements.end(),
                                       other.elements.begin(), other.elements.end(),
                                       std::back_inserter(result));
        return FiniteSet(std::move(result));
    }

    // ⊆ Subset
    bool operator<=(const FiniteSet& other) const {
        return std::includes(other.elements.begin(), other.elements.end(),
                             elements.begin(), elements.end());
    }

    // ⊂ Proper subset
    bool operator<(const FiniteSet& other) const {
        return *this <= other && size() != other.size();
    }

    bool operator==(const FiniteSet& other) const = default;

    // Insert
    FiniteSet& insert(const T& x) {
        auto it = std::lower_bound(elements.begin(), elements.end(), x);
        if (it == elements.end() || *it != x)
            elements.insert(it, x);
        return *this;
    }

    // Remove
    FiniteSet& erase(const T& x) {
        auto it = std::lower_bound(elements.begin(), elements.end(), x);
        if (it != elements.end() && *it == x)
            elements.erase(it);
        return *this;
    }

    // Power set (limited to n <= 20 to prevent UB and combinatorial explosion)
    FiniteSet<FiniteSet<T>> power_set() const {
        FiniteSet<FiniteSet<T>> result;
        auto n = elements.size();
        if (n > 20)
            return result; // too large, return empty
        auto total = 1ULL << n;
        for (std::size_t mask = 0; mask < total; ++mask) {
            std::vector<T> subset;
            for (std::size_t i = 0; i < n; ++i)
                if (mask & (1ULL << i))
                    subset.push_back(elements[i]);
            result.insert(FiniteSet<T>(std::move(subset)));
        }
        return result;
    }

    // Cartesian product A × B → set of pairs
    template<typename U>
    FiniteSet<std::pair<T, U>> cartesian(const FiniteSet<U>& other) const {
        FiniteSet<std::pair<T, U>> result;
        for (const auto& a : elements)
            for (const auto& b : other.elements)
                result.insert({a, b});
        return result;
    }

    // Map: apply f to each element
    template<typename F>
    auto map(F&& f) const -> FiniteSet<decltype(f(std::declval<T>()))> {
        using U = decltype(f(std::declval<T>()));
        std::vector<U> result;
        result.reserve(elements.size());
        for (const auto& x : elements)
            result.push_back(f(x));
        return FiniteSet<U>(std::move(result));
    }

    // Filter
    FiniteSet filter(std::function<bool(const T&)> pred) const {
        std::vector<T> result;
        for (const auto& x : elements)
            if (pred(x))
                result.push_back(x);
        return FiniteSet(std::move(result));
    }

    // Comparison for FiniteSet<FiniteSet<T>> (lexicographic on elements)
    auto operator<=>(const FiniteSet& other) const {
        return elements <=> other.elements;
    }

private:
    void normalize() {
        std::sort(elements.begin(), elements.end());
        elements.erase(std::unique(elements.begin(), elements.end()), elements.end());
    }
};

// ── GeometricSet: FiniteSet with visual shape ─────────────────
// Links a FiniteSet to a geometric shape for visualization.

template<typename T, std::size_t N = 2, Scalar S = double>
struct GeometricSet {
    FiniteSet<T> set;
    geometry::Polygon<N, S> boundary;  // visual boundary for Venn diagrams etc.

    // Geometric intersection region with another GeometricSet
    auto visual_intersection(const GeometricSet& other) const {
        return geometry::intersection_region(boundary, other.boundary);
    }
};

} // namespace spatium::discrete
