#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <array>
#  include <cmath>
#  include <cstdint>
#  include <random>
#  include <utility>
#endif

SPATIUM_EXPORT namespace spatium::algebra {

// Perlin gradient noise -- Ken Perlin's 2002 "improved noise" reference
// algorithm, seeded rather than built from the original fixed table, so
// different seeds give independent, reproducible fields (what "each
// particle gets its own turbulence" needs, without any shared mutable
// state). Output lands in [-1, 1]; scale/bias at the call site if you
// want [0, 1]. A reusable library primitive rather than the ad-hoc
// per-particle RNG velocity every hand-rolled particle demo tends to
// reach for.
class PerlinNoise {
public:
    explicit PerlinNoise(std::uint32_t seed = 0) {
        for (int i = 0; i < 256; ++i) perm_[static_cast<std::size_t>(i)] = static_cast<std::uint8_t>(i);
        std::mt19937 rng(seed);
        for (int i = 255; i > 0; --i) {
            std::uniform_int_distribution<int> d(0, i);
            std::swap(perm_[static_cast<std::size_t>(i)], perm_[static_cast<std::size_t>(d(rng))]);
        }
        for (int i = 0; i < 256; ++i)
            perm_[static_cast<std::size_t>(256 + i)] = perm_[static_cast<std::size_t>(i)];
    }

    template<Scalar T>
    T operator()(T x, T y, T z) const {
        using std::floor;
        int X = static_cast<int>(floor(x)) & 255;
        int Y = static_cast<int>(floor(y)) & 255;
        int Z = static_cast<int>(floor(z)) & 255;
        x -= floor(x);
        y -= floor(y);
        z -= floor(z);
        T u = fade(x), v = fade(y), w = fade(z);

        int A = at(X) + Y, AA = at(A) + Z, AB = at(A + 1) + Z;
        int B = at(X + 1) + Y, BA = at(B) + Z, BB = at(B + 1) + Z;

        return lerp(w,
            lerp(v, lerp(u, grad(at(AA), x, y, z), grad(at(BA), x - T{1}, y, z)),
                    lerp(u, grad(at(AB), x, y - T{1}, z), grad(at(BB), x - T{1}, y - T{1}, z))),
            lerp(v, lerp(u, grad(at(AA + 1), x, y, z - T{1}), grad(at(BA + 1), x - T{1}, y, z - T{1})),
                    lerp(u, grad(at(AB + 1), x, y - T{1}, z - T{1}), grad(at(BB + 1), x - T{1}, y - T{1}, z - T{1}))));
    }

    template<Scalar T>
    T operator()(const Vec<T, 3>& p) const { return (*this)(p[0], p[1], p[2]); }

private:
    std::array<std::uint8_t, 512> perm_{};

    int at(int i) const { return perm_[static_cast<std::size_t>(i) & 511]; }

    template<Scalar T> static T fade(T t) { return t * t * t * (t * (t * T{6} - T{15}) + T{10}); }
    template<Scalar T> static T lerp(T t, T a, T b) { return a + t * (b - a); }
    template<Scalar T> static T grad(int hash, T x, T y, T z) {
        int h = hash & 15;
        T u = h < 8 ? x : y;
        T v = h < 4 ? y : (h == 12 || h == 14 ? x : z);
        return ((h & 1) == 0 ? u : -u) + ((h & 2) == 0 ? v : -v);
    }
};

} // namespace spatium::algebra
