// One-thread-per-ray RK4 geodesic integration for Schwarzschild/Kerr,
// using christoffel_closed_form.hpp's straight-line expressions instead
// of the CPU engine's generic Dual<T> + general 4x4 invert() machinery
// (see that header's own comment for why: SIMT execution punishes the
// indirection/branching genericity costs, and this kernel will only
// ever run these two published metrics). Not yet ported here: the
// volumetric disk emission-absorption integration (trace_ray() in
// examples/blackhole_gr_demo.cpp) -- Phase 3, once this kernel's own
// geodesic math is confirmed to match the CPU engine bit-for-bit
// (gpu/verify_cuda_geodesic.cpp) on the actual T4.
//
// NOT YET COMPILED OR RUN ANYWHERE. Written on a machine with no NVIDIA
// GPU (nvcc itself isn't installed here -- see flake.nix's `cuda`
// devShell, deliberately not fetched locally: the toolkit build is
// several GB just to syntax-check code whose real test (numerically
// matching the CPU engine) needs an actual GPU anyway). First real
// compile + gpu/verify_cuda_geodesic.cpp's cross-check both happen
// together on the Selectel T4.
//
// Fixed dlambda per batch (no adaptive r/D_cam scaling like the CPU
// path's trace_ray() uses) -- deliberate, for this phase: an apples-to-
// apples comparison against a CPU loop using geodesic_step() with the
// same fixed dlambda isolates exactly what Phase 2 needs to prove (does
// the closed-form fp32 kernel math match the generic engine?) without
// also re-deriving the adaptive step heuristic on the device. Adaptive
// stepping is a Phase 3 concern, alongside the rest of the pixel
// pipeline.

#include "christoffel_closed_form.hpp"
#include "geodesic_kernel.h"

namespace spatium::gpu {

// Vec8 layout matches geodesic.hpp's Vec<T,8>: [0..3] = (t,r,theta,phi),
// [4..7] = d/dlambda of same.
template<typename T, void (*ChristoffelFn)(T, T, T, T, T[4][4][4])>
__device__ void geodesic_rhs_2param(T p0, T p1, const T state[8], T ds[8]) {
    T Gamma[4][4][4];
    ChristoffelFn(p0, p1, state[1], state[2], Gamma);
    for (int mu = 0; mu < 4; ++mu) ds[mu] = state[4 + mu];
    for (int lam = 0; lam < 4; ++lam) {
        T acc = T{0};
        for (int mu = 0; mu < 4; ++mu)
            for (int nu = 0; nu < 4; ++nu) acc += Gamma[lam][mu][nu] * state[4 + mu] * state[4 + nu];
        ds[4 + lam] = -acc;
    }
}

// Schwarzschild's closed-form function only takes one physical
// parameter (mass); adapt it to the two-parameter signature above so
// both metrics share one RK4 stepper template.
template<typename T>
__device__ void schwarzschild_christoffel_2param(T mass, T /*unused*/, T r, T theta,
                                                   T Gamma[4][4][4]) {
    schwarzschild_christoffel_closed_form(mass, r, theta, Gamma);
}

template<typename T, void (*ChristoffelFn)(T, T, T, T, T[4][4][4])>
__device__ void rk4_geodesic_step(T p0, T p1, T state[8], T dlambda) {
    T k1[8], k2[8], k3[8], k4[8], tmp[8];

    geodesic_rhs_2param<T, ChristoffelFn>(p0, p1, state, k1);
    for (int i = 0; i < 8; ++i) tmp[i] = state[i] + T{0.5} * dlambda * k1[i];

    geodesic_rhs_2param<T, ChristoffelFn>(p0, p1, tmp, k2);
    for (int i = 0; i < 8; ++i) tmp[i] = state[i] + T{0.5} * dlambda * k2[i];

    geodesic_rhs_2param<T, ChristoffelFn>(p0, p1, tmp, k3);
    for (int i = 0; i < 8; ++i) tmp[i] = state[i] + dlambda * k3[i];

    geodesic_rhs_2param<T, ChristoffelFn>(p0, p1, tmp, k4);

    for (int i = 0; i < 8; ++i)
        state[i] += (dlambda / T{6}) * (k1[i] + T{2} * k2[i] + T{2} * k3[i] + k4[i]);
}

template<typename T, void (*ChristoffelFn)(T, T, T, T, T[4][4][4])>
__global__ void geodesic_batch_kernel(T p0, T p1, const T* state_in, T* state_out, T dlambda,
                                       int n_steps, int n_rays) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_rays) return;

    T state[8];
    for (int k = 0; k < 8; ++k) state[k] = state_in[i * 8 + k];

    for (int step = 0; step < n_steps; ++step) rk4_geodesic_step<T, ChristoffelFn>(p0, p1, state, dlambda);

    for (int k = 0; k < 8; ++k) state_out[i * 8 + k] = state[k];
}

}  // namespace spatium::gpu

// extern "C" batch entry points -- flat arrays only across this
// boundary (no structs/STL), per the plan's own note that the
// verification host program should be able to link against this
// without pulling in the CPU engine's types.
extern "C" {

void spatium_gpu_trace_schwarzschild(double mass, const double* state_in, double* state_out,
                                      double dlambda, int n_steps, int n_rays) {
    double *d_in, *d_out;
    cudaMalloc(&d_in, sizeof(double) * 8 * n_rays);
    cudaMalloc(&d_out, sizeof(double) * 8 * n_rays);
    cudaMemcpy(d_in, state_in, sizeof(double) * 8 * n_rays, cudaMemcpyHostToDevice);

    int threads = 128;
    int blocks = (n_rays + threads - 1) / threads;
    spatium::gpu::geodesic_batch_kernel<double, spatium::gpu::schwarzschild_christoffel_2param<double>>
        <<<blocks, threads>>>(mass, 0.0, d_in, d_out, dlambda, n_steps, n_rays);

    cudaMemcpy(state_out, d_out, sizeof(double) * 8 * n_rays, cudaMemcpyDeviceToHost);
    cudaFree(d_in);
    cudaFree(d_out);
}

void spatium_gpu_trace_kerr(double mass, double spin, const double* state_in, double* state_out,
                             double dlambda, int n_steps, int n_rays) {
    double *d_in, *d_out;
    cudaMalloc(&d_in, sizeof(double) * 8 * n_rays);
    cudaMalloc(&d_out, sizeof(double) * 8 * n_rays);
    cudaMemcpy(d_in, state_in, sizeof(double) * 8 * n_rays, cudaMemcpyHostToDevice);

    int threads = 128;
    int blocks = (n_rays + threads - 1) / threads;
    spatium::gpu::geodesic_batch_kernel<double, spatium::gpu::kerr_christoffel_closed_form<double>>
        <<<blocks, threads>>>(mass, spin, d_in, d_out, dlambda, n_steps, n_rays);

    cudaMemcpy(state_out, d_out, sizeof(double) * 8 * n_rays, cudaMemcpyDeviceToHost);
    cudaFree(d_in);
    cudaFree(d_out);
}

}  // extern "C"
