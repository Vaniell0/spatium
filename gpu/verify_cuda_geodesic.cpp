// Host-side cross-check: CUDA geodesic_kernel.cu's batch RK4 stepper
// (closed-form Christoffel symbols) vs. geodesic.hpp's generic Dual<T>-
// exact christoffel() + rk4_step(), on IDENTICAL initial rays and the
// SAME fixed dlambda/step count on both sides (the CPU side is called
// through geodesic_step() directly, not trace_ray()'s adaptive-step
// wrapper, to isolate exactly what this phase needs to prove: does the
// closed-form kernel math agree with the generic engine?). Tolerance
// comes from Phase 1's own measured fp32-vs-fp64 numbers (worst case
// ~3.5e-4 relative right at the critical curve, ~1e-6 away from it) --
// not an arbitrary stricter bound, since this kernel runs in fp32/fp64
// exactly like that probe did.
//
// NOT YET RUN. Needs an actual NVIDIA GPU (this is where Phase 2's
// local-machine gap gets closed -- see geodesic_kernel.cu's header
// comment) plus SPATIUM_CUDA=ON so gpu/CMakeLists.txt's target exists.

#include "geodesic_kernel.h"

#include <spatium/algebra/vector.hpp>
#include <spatium/physics/relativity/geodesic.hpp>
#include <spatium/physics/relativity/kerr.hpp>
#include <spatium/physics/relativity/schwarzschild.hpp>

#include <cmath>
#include <cstdio>
#include <numbers>
#include <vector>

using namespace spatium;
using namespace spatium::physics::relativity;

namespace {

// Same construction as tests/test_relativity.cpp's new_geodesic_delta_phi()
// / gpu/../fp32_precision_probe.cpp's trace(): a photon aimed inward from
// coordinate radius D with impact parameter b (E=1, L=b), confined to
// theta=pi/2 initially (Kerr's own dynamics will move it off that plane
// exactly like the real renderer's rays do -- only the INITIAL condition
// assumes equatorial launch, matching the camera setups both engines
// actually use).
Vec<double, 8> initial_state(double M, double D, double b) {
    double f0 = 1.0 - 2.0 * M / D;
    double E = 1.0, L = b;
    double ur0 = std::sqrt(std::max(0.0, E * E - f0 * L * L / (D * D)));
    return Vec<double, 8>{0.0, D, std::numbers::pi / 2.0, 0.0, E / f0, -ur0, 0.0, L / (D * D)};
}

struct Mismatch {
    int ray;
    double max_abs_diff;
};

template<typename Metric>
std::vector<Mismatch> compare(const Metric& metric, const std::vector<Vec<double, 8>>& initial,
                               double dlambda, int n_steps,
                               void (*gpu_trace)(const double*, double*, double, int, int)) {
    int n_rays = static_cast<int>(initial.size());
    std::vector<double> state_in(n_rays * 8), state_out_gpu(n_rays * 8);
    for (int i = 0; i < n_rays; ++i)
        for (int k = 0; k < 8; ++k) state_in[i * 8 + k] = initial[i][k];

    gpu_trace(state_in.data(), state_out_gpu.data(), dlambda, n_steps, n_rays);

    std::vector<Mismatch> mismatches;
    for (int i = 0; i < n_rays; ++i) {
        Vec<double, 8> cpu_state = initial[i];
        for (int step = 0; step < n_steps; ++step) cpu_state = geodesic_step(metric, cpu_state, dlambda);

        double max_abs = 0.0, scale = 0.0;
        for (int k = 0; k < 8; ++k) {
            double got = state_out_gpu[i * 8 + k];
            double want = cpu_state[k];
            max_abs = std::max(max_abs, std::abs(got - want));
            scale = std::max(scale, std::abs(want));
        }
        double rel = max_abs / std::max(1e-9, scale);
        std::printf("ray %2d: max_abs_diff=%.3e  rel=%.3e%s\n", i, max_abs, rel,
                    rel > 3.5e-4 ? "  <-- OUT OF TOLERANCE" : "");
        if (rel > 3.5e-4) mismatches.push_back({i, max_abs});
    }
    return mismatches;
}

}  // namespace

int main() {
    double M = 1.0, D = 30.0, dlambda = 0.01;
    int n_steps = 2000;  // fixed step, unlike trace_ray()'s adaptive scaling -- see header comment

    double b_crit = schwarzschild_critical_impact_parameter(M);
    std::vector<Vec<double, 8>> schw_initial;
    for (double ratio : {0.5, 0.9, 0.99, 0.999, 0.9999, 1.0001, 1.001, 1.01, 1.1, 2.0, 5.0})
        schw_initial.push_back(initial_state(M, D, ratio * b_crit));

    std::printf("=== Schwarzschild: CPU (Dual<T> exact) vs CUDA (closed-form) ===\n");
    SchwarzschildMetric<double> schw{M};
    auto schw_mismatches =
        compare(schw, schw_initial, dlambda, n_steps,
                [](const double* in, double* out, double dl, int ns, int nr) {
                    spatium_gpu_trace_schwarzschild(1.0, in, out, dl, ns, nr);
                });

    double a = 0.9;
    std::vector<Vec<double, 8>> kerr_initial;
    for (double b : {2.0, 3.0, 3.5, 3.6, 3.65, 3.7, 3.8, 4.0, 5.0, 8.0})
        kerr_initial.push_back(initial_state(M, D, b));

    std::printf("\n=== Kerr (a=%.2f): CPU (Dual<T> exact) vs CUDA (closed-form) ===\n", a);
    KerrMetric<double> kerr{M, a};
    auto kerr_mismatches =
        compare(kerr, kerr_initial, dlambda, n_steps,
                [](const double* in, double* out, double dl, int ns, int nr) {
                    spatium_gpu_trace_kerr(1.0, 0.9, in, out, dl, ns, nr);
                });

    std::printf("\nSchwarzschild mismatches: %zu/%zu, Kerr mismatches: %zu/%zu\n",
                schw_mismatches.size(), schw_initial.size(), kerr_mismatches.size(), kerr_initial.size());
    return (schw_mismatches.empty() && kerr_mismatches.empty()) ? 0 : 1;
}
