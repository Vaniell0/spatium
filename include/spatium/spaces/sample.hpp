#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/spaces/parametric.hpp>
#  include <algorithm>
#  include <cstdint>
#  include <optional>
#  include <random>
#  include <utility>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium {

template<Scalar T = double>
struct SurfaceSample {
    T u, v;
    Vec<T, 3> position;
    Vec<T, 3> normal;
};

// sample_surface_uniform: place `count` points evenly by *surface area*
// (not by parameter (u,v)) across an analytic ParametricSurface, via
// rejection sampling weighted by the first fundamental form's area
// element sqrt(EG-F^2) -- the "sprinkles" primitive, fully analytic: no
// mesh, no geodesic-Voronoi graph. Sampling (u,v) uniformly instead would
// bunch points wherever the parametrization compresses ambient space
// (e.g. a thin-ring torus's inner rim); this doesn't.
//
// The samples are PREFIX-STABLE: they come out of one random sequence in
// order, so asking for more gives the same points followed by new ones, and
// only a new seed gives different ones. `SurfaceSampler` is that sequence
// as an object whose state can be copied and resumed, for a caller that
// wants the points after the first k without drawing the first k again.
template<Scalar T = double>
class SurfaceSampler {
public:
    // Where the sequence is: the generator and how many draws it has made.
    struct State {
        std::mt19937 rng;
        std::size_t draws = 0;
    };

    SurfaceSampler(const ParametricSurface<T>& surface, std::uint32_t seed)
        : surface_(&surface), rng_(seed) {
        auto [u_min, u_max, v_min, v_max] = surface.domain();
        // Coarse grid to bound the rejection envelope.
        constexpr int grid = 24;
        for (int j = 0; j <= grid; ++j) {
            T v = v_min + (v_max - v_min) * static_cast<T>(j) / static_cast<T>(grid);
            for (int i = 0; i <= grid; ++i) {
                T u = u_min + (u_max - u_min) * static_cast<T>(i) / static_cast<T>(grid);
                max_area_ = std::max(max_area_, surface.area_element(u, v));
            }
        }
        if (max_area_ <= T{0}) max_area_ = T{1};
        du_ = std::uniform_real_distribution<T>(u_min, u_max);
        dv_ = std::uniform_real_distribution<T>(v_min, v_max);
    }

    SurfaceSampler(const ParametricSurface<T>& surface, std::uint32_t seed, State state)
        : SurfaceSampler(surface, seed) {
        rng_ = std::move(state.rng);
        draws_ = state.draws;
    }

    State state() const { return State{rng_, draws_}; }
    std::size_t draws() const { return draws_; }

    // The parameters of the next accepted point, or nothing once `max_draws`
    // draws have been made in all (the surface accepts too rarely).
    std::optional<std::pair<T, T>> next_uv(std::size_t max_draws) {
        while (draws_ < max_draws) {
            ++draws_;
            T u = du_(rng_), v = dv_(rng_);
            if (unit_(rng_) * max_area_ <= surface_->area_element(u, v)) return std::pair<T, T>{u, v};
        }
        return std::nullopt;
    }

    SurfaceSample<T> at(T u, T v) const {
        return {u, v, surface_->evaluate(u, v), surface_->normal_at(u, v)};
    }

    std::optional<SurfaceSample<T>> next(std::size_t max_draws) {
        const auto uv = next_uv(max_draws);
        if (!uv) return std::nullopt;
        return at(uv->first, uv->second);
    }

private:
    const ParametricSurface<T>* surface_;
    std::mt19937 rng_;
    std::uniform_real_distribution<T> du_, dv_, unit_{T{0}, T{1}};
    T max_area_{0};
    std::size_t draws_ = 0;
};

// The draws `sample_surface_uniform` allows itself for `count` points
// before it gives up on a surface that accepts too rarely.
inline std::size_t sample_draw_limit(std::size_t count) { return count * 200 + 10000; }

template<Scalar T = double>
std::vector<SurfaceSample<T>> sample_surface_uniform(
    const ParametricSurface<T>& surface, std::size_t count, std::uint32_t seed = 42)
{
    SurfaceSampler<T> sampler(surface, seed);
    std::vector<SurfaceSample<T>> out;
    out.reserve(count);
    const std::size_t limit = sample_draw_limit(count);
    while (out.size() < count) {
        auto s = sampler.next(limit);
        if (!s) break;
        out.push_back(*s);
    }
    return out;
}

} // namespace spatium
