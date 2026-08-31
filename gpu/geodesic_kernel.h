#pragma once

// extern "C" declarations for geodesic_kernel.cu's batch entry points --
// a plain C++ host program (gpu/verify_cuda_geodesic.cpp, compiled with
// the regular host compiler, not nvcc) includes this instead of the
// .cu file directly. Flat double arrays only, no structs/STL: state_in/
// state_out are n_rays*8 doubles each, laid out [t,r,theta,phi,
// dt,dr,dtheta,dphi] per ray, matching geodesic.hpp's Vec<T,8>
// convention.

#ifdef __cplusplus
extern "C" {
#endif

void spatium_gpu_trace_schwarzschild(double mass, const double* state_in, double* state_out,
                                      double dlambda, int n_steps, int n_rays);

void spatium_gpu_trace_kerr(double mass, double spin, const double* state_in, double* state_out,
                             double dlambda, int n_steps, int n_rays);

#ifdef __cplusplus
}
#endif
