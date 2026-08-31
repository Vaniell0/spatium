#pragma once

// extern "C" declaration for schwarzschild_render_kernel.cu's batch
// entry point -- Phase 3 of the GPU port (volumetric disk integration,
// on top of Phase 2's already-verified geodesic math). Included from
// plain C++ host code (gpu/verify_cuda_render.cpp), not from nvcc-only
// translation units.
//
// Deliberately stops short of the actual starfield lookup
// (sample_sky_color() against examples/io_helpers.hpp's procedural Sky)
// -- that's a large, CPU-only, unrelated-to-physics structure (10000
// stars + spiral/cloud procedural fields + a bucket grid) with nothing
// to do with the GR/disk math this phase is validating. Both engines
// return the PRE-sky state (accumulated emission, remaining
// transmittance, and either the escape direction or a captured flag);
// the caller composites with whatever sky it likes, identically on
// both sides, outside this comparison.

#include "disk_physics.hpp"

struct RenderParams {
    double mass, D_cam;
    double r_hat_x, r_hat_y, r_hat_z;
    double disk_time;
    double r_horizon, r_capture, r_escape;
    double r_disk_inner, r_disk_outer;
    double disk_aspect_ratio, disk_density_power;
    double disk_t0, disk_t_ref, doppler_damping;
    double absorption, emission_scale, transmittance_cutoff;
    double base_dlambda, dlambda_floor_frac;
    int max_steps;
};

struct KerrRenderParams {
    double mass, spin;
    double D_cam;
    double cam_r, cam_theta, cam_phi;
    double e_r_x, e_r_y, e_r_z;
    double e_theta_x, e_theta_y, e_theta_z;
    double e_phi_x, e_phi_y, e_phi_z;
    double ut_cam, sqrt_grr_cam, sqrt_gthth_cam, sqrt_gphph_cam;
    double g_cam_00, g_cam_03, g_cam_33;
    int prograde;
    double disk_time;
    double r_horizon, r_capture, r_escape;
    double r_disk_inner, r_disk_outer;
    double disk_aspect_ratio, disk_density_power;
    double disk_t0, disk_t_ref, doppler_damping;
    double absorption, emission_scale, transmittance_cutoff;
    double base_dlambda, dlambda_floor_frac;
    int max_steps;
};

// out layout per ray (8 doubles): accumulated.xyz, transmittance,
// exit_dir.xyz, captured_flag (1.0 = captured/absorbed, 0.0 = escaped
// or ran out of steps -- exit_dir is meaningful only when 0.0).
#ifdef __cplusplus
extern "C" {
#endif

void spatium_gpu_render_schwarzschild(RenderParams params,
                                       const spatium::gpu::DiscParticleParams* particles,
                                       int n_particles, const double* dir_world_in, double* out,
                                       int n_rays);

void spatium_gpu_render_kerr(KerrRenderParams params, const spatium::gpu::DiscParticleParams* particles,
                              int n_particles, const double* dir_world_in, double* out, int n_rays);

#ifdef __cplusplus
}
#endif
