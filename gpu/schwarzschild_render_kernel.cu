// Phase 3: the full Schwarzschild per-ray pipeline (geodesic + volumetric
// disk emission-absorption), one thread per ray -- a straight port of
// blackhole_gr_demo.cpp's Schwarzschild fast path (the tetrad-corrected
// initial state, trace_ray()'s emission-absorption loop) onto
// christoffel_closed_form.hpp's closed-form geodesic math (already
// verified against the CPU engine to ~1e-16 relative in Phase 2).
//
// Stops short of the sky lookup -- see render_kernel.h's comment.
// Cross-checked against the CPU pipeline (down to that same point) in
// gpu/verify_cuda_render.cpp.

#include "christoffel_closed_form.hpp"
#include "disk_physics.hpp"
#include "render_kernel.h"

namespace spatium::gpu {

__device__ void schw_geodesic_rhs(double mass, const double state[8], double ds[8]) {
    double Gamma[4][4][4];
    schwarzschild_christoffel_closed_form(mass, state[1], state[2], Gamma);
    for (int mu = 0; mu < 4; ++mu) ds[mu] = state[4 + mu];
    for (int lam = 0; lam < 4; ++lam) {
        double acc = 0.0;
        for (int mu = 0; mu < 4; ++mu)
            for (int nu = 0; nu < 4; ++nu) acc += Gamma[lam][mu][nu] * state[4 + mu] * state[4 + nu];
        ds[4 + lam] = -acc;
    }
}

__device__ void schw_rk4_step(double mass, double state[8], double dl) {
    double k1[8], k2[8], k3[8], k4[8], tmp[8];
    schw_geodesic_rhs(mass, state, k1);
    for (int i = 0; i < 8; ++i) tmp[i] = state[i] + 0.5 * dl * k1[i];
    schw_geodesic_rhs(mass, tmp, k2);
    for (int i = 0; i < 8; ++i) tmp[i] = state[i] + 0.5 * dl * k2[i];
    schw_geodesic_rhs(mass, tmp, k3);
    for (int i = 0; i < 8; ++i) tmp[i] = state[i] + dl * k3[i];
    schw_geodesic_rhs(mass, tmp, k4);
    for (int i = 0; i < 8; ++i) state[i] += (dl / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
}

__device__ double schw_disk_redshift(double mass, double r, double b) {
    double omega = sqrt(mass / (r * r * r));
    double ut = 1.0 / sqrt(1.0 - 3.0 * mass / r);
    return ut * (1.0 - b * omega);
}

__global__ void schwarzschild_render_batch(RenderParams p, const DiscParticleParams* particles,
                                            int n_particles, const double* dir_world_in, double* out,
                                            int n_rays) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_rays) return;

    double dwx = dir_world_in[i * 3 + 0], dwy = dir_world_in[i * 3 + 1], dwz = dir_world_in[i * 3 + 2];
    double rhx = p.r_hat_x, rhy = p.r_hat_y, rhz = p.r_hat_z;
    double f_cam = 1.0 - 2.0 * p.mass / p.D_cam;

    double dir_dot_r = dwx * rhx + dwy * rhy + dwz * rhz;
    double trx = dwx - rhx * dir_dot_r, try_ = dwy - rhy * dir_dot_r, trz = dwz - rhz * dir_dot_r;
    double t_norm = sqrt(trx * trx + try_ * try_ + trz * trz);

    double accum[3] = {0.0, 0.0, 0.0};
    double transmittance = 1.0;

    if (t_norm < 1e-9) {
        // Radial ray: either straight into the hole (no picture) or
        // straight out (nothing to accumulate) -- matches the CPU
        // early-return exactly, just without the sky sample itself.
        out[i * 8 + 0] = 0.0; out[i * 8 + 1] = 0.0; out[i * 8 + 2] = 0.0;
        out[i * 8 + 3] = 1.0;
        out[i * 8 + 4] = dwx; out[i * 8 + 5] = dwy; out[i * 8 + 6] = dwz;
        out[i * 8 + 7] = (dir_dot_r < 0.0) ? 1.0 : 0.0;
        return;
    }

    double thx = trx / t_norm, thy = try_ / t_norm, thz = trz / t_norm;
    double ut = 1.0 / sqrt(f_cam);
    double ur = dir_dot_r * sqrt(f_cam);
    double uphi = t_norm / p.D_cam;
    double b_photon = (p.D_cam * p.D_cam * uphi) / (f_cam * ut);

    // n_ray = r_hat x t_hat; line_dir = n_ray x (0,0,1)
    double nrx = rhy * thz - rhz * thy, nry = rhz * thx - rhx * thz, nrz = rhx * thy - rhy * thx;
    double ldx = nry * 1.0 - nrz * 0.0, ldy = nrz * 0.0 - nrx * 1.0, ldz = 0.0;
    double line_norm = sqrt(ldx * ldx + ldy * ldy + ldz * ldz);
    bool disk_plane_ok = line_norm > 1e-9;
    double phi_disk = 0.0;
    if (disk_plane_ok) {
        ldx /= line_norm; ldy /= line_norm; ldz /= line_norm;
        double dot_t = ldx * thx + ldy * thy + ldz * thz;
        double dot_r = ldx * rhx + ldy * rhy + ldz * rhz;
        phi_disk = atan2(dot_t, dot_r);
    }

    double state[8] = {0.0, p.D_cam, 1.5707963267948966, 0.0, ut, ur, 0.0, uphi};
    double r_min = state[1];

    for (int step = 0; step < p.max_steps; ++step) {
        double r = state[1], theta = state[2], phi = state[3];
        double dl = p.base_dlambda * fmax(p.dlambda_floor_frac, r / p.D_cam);
        schw_rk4_step(p.mass, state, dl);
        double r_new = state[1];
        if (!isfinite(r_new)) {
            out[i * 8 + 0] = accum[0]; out[i * 8 + 1] = accum[1]; out[i * 8 + 2] = accum[2];
            out[i * 8 + 3] = transmittance;
            out[i * 8 + 4] = 0.0; out[i * 8 + 5] = 0.0; out[i * 8 + 6] = 0.0;
            out[i * 8 + 7] = 1.0;
            return;
        }
        if (r_new < r_min) r_min = r_new;
        double theta_new = state[2], phi_new = state[3];

        double r_mid = 0.5 * (r + r_new), theta_mid = 0.5 * (theta + theta_new),
               phi_mid = 0.5 * (phi + phi_new);
        double height = disk_plane_ok ? r_mid * sin(phi_mid - phi_disk) * line_norm : 0.0;
        double envelope = disc_density_envelope(r_mid, height, p.r_horizon, p.r_disk_inner,
                                                 p.r_disk_outer, p.disk_aspect_ratio,
                                                 p.disk_density_power);
        if (envelope > 1e-6) {
            double texture = disc_density(r_mid, phi_mid, p.disk_time, particles, n_particles,
                                           p.r_disk_outer);
            double density = envelope * texture;

            double z = schw_disk_redshift(p.mass, r_mid, b_photon);
            double z_display = 1.0 + p.doppler_damping * (z - 1.0);
            double t_emit = disk_temperature_profile(r_mid, p.r_disk_inner, p.disk_t0);
            double t_obs = t_emit / fmax(z_display, 1e-3);
            double intensity = gpu_clamp(pow(t_obs / p.disk_t_ref, 1.2), 0.0, 1.0);
            Vec3d bb = blackbody_to_rgb255(t_obs);
            Vec3d emit_color =
                boost_saturation({bb.x * intensity, bb.y * intensity, bb.z * intensity}, 1.15);

            double optical_depth = density * p.absorption * dl;
            double emit_amount = density * p.emission_scale * dl;
            accum[0] += emit_color.x * (transmittance * emit_amount);
            accum[1] += emit_color.y * (transmittance * emit_amount);
            accum[2] += emit_color.z * (transmittance * emit_amount);
            transmittance *= exp(-optical_depth);
        }

        if (r_new <= p.r_capture || transmittance < p.transmittance_cutoff) {
            out[i * 8 + 0] = accum[0]; out[i * 8 + 1] = accum[1]; out[i * 8 + 2] = accum[2];
            out[i * 8 + 3] = transmittance;
            out[i * 8 + 4] = 0.0; out[i * 8 + 5] = 0.0; out[i * 8 + 6] = 0.0;
            out[i * 8 + 7] = 1.0;
            return;
        }

        if (r_new >= p.r_escape && state[5] > 0.0) {
            double r_e = state[1], phi_e = state[3];
            double f_esc = 1.0 - 2.0 * p.mass / r_e;
            double nr = state[5] / sqrt(f_esc);
            double nphi = r_e * state[7];
            double mag = sqrt(nr * nr + nphi * nphi);
            nr /= mag; nphi /= mag;
            double cph = cos(phi_e), sph = sin(phi_e);
            double rhx_e = rhx * cph + thx * sph, rhy_e = rhy * cph + thy * sph,
                   rhz_e = rhz * cph + thz * sph;
            double thx_e = thx * cph - rhx * sph, thy_e = thy * cph - rhy * sph,
                   thz_e = thz * cph - rhz * sph;
            out[i * 8 + 0] = accum[0]; out[i * 8 + 1] = accum[1]; out[i * 8 + 2] = accum[2];
            out[i * 8 + 3] = transmittance;
            out[i * 8 + 4] = rhx_e * nr + thx_e * nphi;
            out[i * 8 + 5] = rhy_e * nr + thy_e * nphi;
            out[i * 8 + 6] = rhz_e * nr + thz_e * nphi;
            out[i * 8 + 7] = 0.0;
            return;
        }
    }
    // Ran out of steps without resolving: fold into captured rather
    // than faking an escape along the original camera-ray direction --
    // see kerr_render_kernel.cu's matching comment for why (a visible
    // blocky discontinuity right at the resolved/unresolved boundary).
    out[i * 8 + 0] = accum[0]; out[i * 8 + 1] = accum[1]; out[i * 8 + 2] = accum[2];
    out[i * 8 + 3] = transmittance;
    out[i * 8 + 4] = 0.0; out[i * 8 + 5] = 0.0; out[i * 8 + 6] = 0.0;
    out[i * 8 + 7] = 1.0;
}

}  // namespace spatium::gpu

extern "C" void spatium_gpu_render_schwarzschild(RenderParams params,
                                                  const spatium::gpu::DiscParticleParams* particles,
                                                  int n_particles, const double* dir_world_in,
                                                  double* out, int n_rays) {
    spatium::gpu::DiscParticleParams* d_particles;
    double *d_dirs, *d_out;
    cudaMalloc(&d_particles, sizeof(spatium::gpu::DiscParticleParams) * n_particles);
    cudaMalloc(&d_dirs, sizeof(double) * 3 * n_rays);
    cudaMalloc(&d_out, sizeof(double) * 8 * n_rays);
    cudaMemcpy(d_particles, particles, sizeof(spatium::gpu::DiscParticleParams) * n_particles,
               cudaMemcpyHostToDevice);
    cudaMemcpy(d_dirs, dir_world_in, sizeof(double) * 3 * n_rays, cudaMemcpyHostToDevice);

    int threads = 128;
    int blocks = (n_rays + threads - 1) / threads;
    spatium::gpu::schwarzschild_render_batch<<<blocks, threads>>>(params, d_particles, n_particles,
                                                                    d_dirs, d_out, n_rays);

    cudaMemcpy(out, d_out, sizeof(double) * 8 * n_rays, cudaMemcpyDeviceToHost);
    cudaFree(d_particles);
    cudaFree(d_dirs);
    cudaFree(d_out);
}
