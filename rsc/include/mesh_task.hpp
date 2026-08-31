#pragma once

// Domain 5 (mesh-strategy dispatch), configured over comparison_task.hpp,
// fifth instance of that shape. See mesh_ops.hpp for the two test
// families (Torus by aspect ratio, Cone by apex proximity), the two
// candidate strategies (uniform vs. anisotropy-adapted UV step), and the
// measured reason the ground-truth metric is ambient aspect ratio, not
// face-normal deviation.
//
// Output = double (the resulting aspect ratio, always >= 1); reference =
// a constant 1.0 ("perfectly isotropic" quad), not something computed
// per-problem -- comparison_task.hpp's ReferenceFn is any
// function(Problem)->Output, and a domain where the ideal answer doesn't
// depend on the problem is just as valid an instance as
// precision_task.hpp's per-problem Real50 answer.
//
// Same precompute-once-per-level architecture as geodesic_task.hpp, for
// the same reason: computing a candidate's aspect ratio is cheap here
// (a handful of surface evaluations), but there's no reason to recompute
// it on every sample() call when the full (level, quad) population is
// finite and small enough to enumerate up front. "Level" is one concrete
// surface (a family + a shape parameter); "quad" is one cell of that
// surface's coarse UV grid -- the population is exactly enumerable, so
// feature_mean is computed exactly over it, not estimated via a warmup
// sample the way ode_task.hpp's continuous rate/dt distribution needs.
//
// Feature = log(anisotropy) - mean, same centering fix as
// geodesic_task.hpp's log(vertex_count) and rootfind_task.hpp's
// |f'(x0)|: anisotropy is always >= 1, so its raw log is always >= 0,
// which starves half of a ReLU net's hidden units at initialization the
// same way an uncentered always-nonnegative feature did there -- applied
// up front this time, not rediscovered a third time.
//
// tolerance = 0.2, not the domain-wide default 1e-3 -- measured, not
// guessed: adapted_aspect's own residual error (curvature/nonlinearity
// the local-speed linearization doesn't capture) sits around 0.0-0.09
// across every level tried, while uniform_aspect's error spans from
// near-0 (fat torus/cone patches) up to ~20 (thin-ring torus), so a
// tight tolerance made "uniform wins" essentially empty (one bucket
// collapsing to zero segfaults ComparisonTaskGenerator's stratified
// sampler below) -- 0.2 is the smallest value found, by direct probing,
// that keeps both buckets genuinely populated at every level.
//
// Ground truth uses the same stratified-bucket sampling geodesic_task.hpp
// needed: most quads need the adapted strategy at the shape parameters
// chosen here (uniform only wins outright on the least-anisotropic
// levels), so uniform sampling over quads would starve the minority
// "uniform is already fine" class the same way uniform sampling over
// levels starved Dijkstra's one winning level there.

#include <comparison_task.hpp>
#include <mesh_ops.hpp>
#include <cmath>
#include <memory>
#include <vector>

namespace rsc {

struct MeshProblem {
    std::size_t level_index;
    std::size_t quad_index;
};

using MeshOutput = double;
using MeshTaskGenerator = ComparisonTaskGenerator<MeshProblem, MeshOutput>;

inline MeshTaskGenerator build_mesh_task_generator(std::uint64_t seed = 0, double tolerance = 0.2) {
    using namespace mesh_detail;

    struct Level {
        std::vector<double> anisotropy, uniform_ar, adapted_ar;
    };
    auto levels = std::make_shared<std::vector<Level>>();

    struct LevelSpec {
        MeshFamily family;
        double shape_param;
    };
    // Torus: aspect ratio R/r sweeping fat (near-isometric) to thin-ring
    // (severe, ~measured up to ~21x in test_parametric.cpp). Cone: two
    // heights, each producing its own apex-proximity degeneracy
    // regardless of scale.
    std::vector<LevelSpec> specs{
        {MeshFamily::Torus, 1.2}, {MeshFamily::Torus, 2.0},  {MeshFamily::Torus, 4.0},
        {MeshFamily::Torus, 8.0}, {MeshFamily::Torus, 16.0}, {MeshFamily::Torus, 30.0},
        {MeshFamily::Cone, 1.5},  {MeshFamily::Cone, 4.0},
    };

    constexpr int kNu = 12, kNv = 8;

    for (auto& spec : specs) {
        auto surf = make_surface(spec.family, spec.shape_param);
        auto dom = surf.domain();
        double du = (dom.u_max - dom.u_min) / kNu;
        double dv = (dom.v_max - dom.v_min) / kNv;

        Level lv;
        for (int j = 0; j < kNv; ++j) {
            for (int i = 0; i < kNu; ++i) {
                QuadProbe q{dom.u_min + i * du, dom.v_min + j * dv, dom.u_min + (i + 1) * du,
                            dom.v_min + (j + 1) * dv};
                lv.anisotropy.push_back(
                    surf.parametrization_anisotropy(0.5 * (q.u0 + q.u1), 0.5 * (q.v0 + q.v1)));
                lv.uniform_ar.push_back(uniform_aspect(surf, q));
                lv.adapted_ar.push_back(adapted_aspect(surf, q));
            }
        }
        levels->push_back(std::move(lv));
    }

    double feature_mean = 0.0;
    std::size_t total = 0;
    for (auto& lv : *levels)
        for (double a : lv.anisotropy) {
            feature_mean += std::log(a);
            ++total;
        }
    feature_mean /= static_cast<double>(total);

    std::vector<Candidate<MeshProblem, MeshOutput>> candidates{
        {"mesh_uniform", [levels](const MeshProblem& p) {
             return (*levels)[p.level_index].uniform_ar[p.quad_index];
         }},
        {"mesh_adapted", [levels](const MeshProblem& p) {
             return (*levels)[p.level_index].adapted_ar[p.quad_index];
         }},
    };

    // Stratify by ground truth once here, at construction -- same
    // "closest wins" comparison the generic template uses at sample
    // time, just precomputed -- then sample uniformly between the two
    // buckets instead of uniformly across (level, quad) pairs.
    struct ProblemRef {
        std::size_t level_index, quad_index;
    };
    auto uniform_bucket = std::make_shared<std::vector<ProblemRef>>();
    auto adapted_bucket = std::make_shared<std::vector<ProblemRef>>();
    for (std::size_t lvl = 0; lvl < levels->size(); ++lvl) {
        auto& lv = (*levels)[lvl];
        for (std::size_t q = 0; q < lv.anisotropy.size(); ++q) {
            double uniform_err = std::abs(lv.uniform_ar[q] - 1.0);
            double adapted_err = std::abs(lv.adapted_ar[q] - 1.0);
            bool uniform_wins = uniform_err <= adapted_err + tolerance;
            (uniform_wins ? uniform_bucket : adapted_bucket)->push_back({lvl, q});
        }
    }

    auto sample_problem = [uniform_bucket, adapted_bucket](std::mt19937_64& rng) -> MeshProblem {
        std::bernoulli_distribution pick_uniform_bucket(0.5);
        auto& bucket = pick_uniform_bucket(rng) ? *uniform_bucket : *adapted_bucket;
        std::uniform_int_distribution<std::size_t> idx(0, bucket.size() - 1);
        auto ref = bucket[idx(rng)];
        return {ref.level_index, ref.quad_index};
    };

    auto extract_features = [levels, feature_mean](const MeshProblem& p) -> std::vector<double> {
        double a = (*levels)[p.level_index].anisotropy[p.quad_index];
        return {std::log(a) - feature_mean};
    };

    auto reference = [](const MeshProblem&) -> MeshOutput { return 1.0; };
    auto distance = [](double a, double b) { return std::abs(a - b); };

    return MeshTaskGenerator(std::move(candidates), sample_problem, extract_features, reference,
                              distance, tolerance, seed);
}

} // namespace rsc
