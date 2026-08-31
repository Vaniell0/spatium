// Phase 3, Kerr path: the same full per-ray pipeline as
// schwarzschild_render_kernel.cu, but with the oblate-spheroidal
// camera tetrad and genuine 3D (t,r,theta,phi) dynamics the Kerr path
// needs (theta is a real dynamical variable, not fixed at pi/2 -- see
// blackhole_gr_demo.cpp's KERR PATH header note). Reuses
// christoffel_closed_form.hpp's kerr_christoffel_closed_form (Phase 2,
// already verified to ~1e-16 relative against the CPU engine) and
// disk_physics.hpp's disc_density/disc_density_envelope/etc. unchanged
// -- only the geodesic setup and exit-tetrad differ from the
// Schwarzschild kernel, matching how the CPU source itself factors it
// (both paths delegate to one shared trace_ray(), differing only in
// initial/exit tetrad construction and the disk-height test).

#include "christoffel_closed_form.hpp"
#include "disk_physics.hpp"
#include "kerr_geometry.hpp"
#include "render_kernel.h"

namespace spatium::gpu {

__device__ void kerr_geodesic_rhs(double mass, double spin, const double state[8], double ds[8]) {
    double Gamma[4][4][4];
    kerr_christoffel_closed_form(mass, spin, state[1], state[2], Gamma);
    for (int mu = 0; mu < 4; ++mu) ds[mu] = state[4 + mu];
    for (int lam = 0; lam < 4; ++lam) {
        double acc = 0.0;
        for (int mu = 0; mu < 4; ++mu)
            for (int nu = 0; nu < 4; ++nu) acc += Gamma[lam][mu][nu] * state[4 + mu] * state[4 + nu];
        ds[4 + lam] = -acc;
    }
}

__device__ void kerr_rk4_step(double mass, double spin, double state[8], double dl) {
    double k1[8], k2[8], k3[8], k4[8], tmp[8];
    kerr_geodesic_rhs(mass, spin, state, k1);
    for (int i = 0; i < 8; ++i) tmp[i] = state[i] + 0.5 * dl * k1[i];
    kerr_geodesic_rhs(mass, spin, tmp, k2);
    for (int i = 0; i < 8; ++i) tmp[i] = state[i] + 0.5 * dl * k2[i];
    kerr_geodesic_rhs(mass, spin, tmp, k3);
    for (int i = 0; i < 8; ++i) tmp[i] = state[i] + dl * k3[i];
    kerr_geodesic_rhs(mass, spin, tmp, k4);
    for (int i = 0; i < 8; ++i) state[i] += (dl / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
}

__global__ void kerr_render_batch(KerrRenderParams p, const DiscParticleParams* particles,
                                   int n_particles, const double* dir_world_in, double* out,
                                   int n_rays) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_rays) return;

    double dwx = dir_world_in[i * 3 + 0], dwy = dir_world_in[i * 3 + 1], dwz = dir_world_in[i * 3 + 2];

    double n_r = dwx * p.e_r_x + dwy * p.e_r_y + dwz * p.e_r_z;
    double n_th = dwx * p.e_theta_x + dwy * p.e_theta_y + dwz * p.e_theta_z;
    double n_ph = dwx * p.e_phi_x + dwy * p.e_phi_y + dwz * p.e_phi_z;

    // up=(0,0,1) is always the spin axis direction, and the oblate-
    // spheroidal e_phi tetrad vector has an identically zero z-
    // component by construction -- so up.e_phi == 0 always, which
    // means there is ALWAYS exactly one screen column (same sx for
    // every sy) where n_ph (and hence uphi) is exactly 0 for every ray
    // in it. Rays launched exactly on that line sit on a genuine
    // separatrix (zero initial azimuthal rate) that neighboring rays
    // approach but never hit -- found, this session, to render as a
    // visible column of blocky capture/escape flips (reported as
    // looking "censored"). Nudged here, on n_ph rather than uphi:
    // n_ph is a dot product of two unit vectors, so its natural scale
    // is always O(1) regardless of camera distance -- uphi's own scale
    // shrinks with distance (uphi = n_ph/sqrt(g_phph), and sqrt(g_phph)
    // grows with r), so a fixed epsilon on uphi itself was miscalibrated
    // for wide shots (clamped the vast majority of rays, not just the
    // degenerate column -- a real regression caught before this render
    // was ever committed to the long run).
    constexpr double kDegeneracyEpsilon = 1e-4;
    if (fabs(n_ph) < kDegeneracyEpsilon) n_ph = kDegeneracyEpsilon;

    double ur = n_r / p.sqrt_grr_cam;
    double utheta = n_th / p.sqrt_gthth_cam;
    double uphi = n_ph / p.sqrt_gphph_cam;

    double state[8] = {0.0, p.cam_r, p.cam_theta, p.cam_phi, p.ut_cam, ur, utheta, uphi};

    double L = p.g_cam_03 * p.ut_cam + p.g_cam_33 * uphi;
    double E = -(p.g_cam_00 * p.ut_cam + p.g_cam_03 * uphi);
    double b_photon = L / E;

    double accum[3] = {0.0, 0.0, 0.0};
    double transmittance = 1.0;
    bool prograde = p.prograde != 0;

    for (int step = 0; step < p.max_steps; ++step) {
        double r = state[1], theta = state[2], phi = state[3];
        double dl = p.base_dlambda * fmax(p.dlambda_floor_frac, r / p.D_cam);
        kerr_rk4_step(p.mass, p.spin, state, dl);
        double r_new = state[1];
        if (!isfinite(r_new)) {
            out[i * 8 + 0] = accum[0]; out[i * 8 + 1] = accum[1]; out[i * 8 + 2] = accum[2];
            out[i * 8 + 3] = transmittance;
            out[i * 8 + 4] = 0.0; out[i * 8 + 5] = 0.0; out[i * 8 + 6] = 0.0;
            out[i * 8 + 7] = 1.0;
            return;
        }
        double theta_new = state[2], phi_new = state[3];

        double r_mid = 0.5 * (r + r_new), theta_mid = 0.5 * (theta + theta_new),
               phi_mid = 0.5 * (phi + phi_new);
        double height = r_mid * cos(theta_mid);
        double envelope = disc_density_envelope(r_mid, height, p.r_horizon, p.r_disk_inner,
                                                 p.r_disk_outer, p.disk_aspect_ratio,
                                                 p.disk_density_power);
        if (envelope > 1e-6) {
            double texture = disc_density(r_mid, phi_mid, p.disk_time, particles, n_particles,
                                           p.r_disk_outer);
            double density = envelope * texture;

            double z = kerr_disk_redshift_factor(p.mass, p.spin, r_mid, b_photon, prograde);
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
            double r_e = state[1], theta_e = state[2], phi_e = state[3];
            double g_rr, g_thth, g_phph;
            kerr_metric_diag(p.mass, p.spin, r_e, theta_e, &g_rr, &g_thth, &g_phph);
            double nr_e = state[5] * sqrt(g_rr);
            double nth_e = state[6] * sqrt(g_thth);
            double nph_e = state[7] * sqrt(g_phph);
            double mag = sqrt(nr_e * nr_e + nth_e * nth_e + nph_e * nph_e);
            nr_e /= mag; nth_e /= mag; nph_e /= mag;
            SpheroidalTetradD tet_e = spheroidal_tetrad(r_e, theta_e, phi_e, p.spin);
            out[i * 8 + 0] = accum[0]; out[i * 8 + 1] = accum[1]; out[i * 8 + 2] = accum[2];
            out[i * 8 + 3] = transmittance;
            out[i * 8 + 4] = tet_e.e_r.x * nr_e + tet_e.e_theta.x * nth_e + tet_e.e_phi.x * nph_e;
            out[i * 8 + 5] = tet_e.e_r.y * nr_e + tet_e.e_theta.y * nth_e + tet_e.e_phi.y * nph_e;
            out[i * 8 + 6] = tet_e.e_r.z * nr_e + tet_e.e_theta.z * nth_e + tet_e.e_phi.z * nph_e;
            out[i * 8 + 7] = 0.0;
            return;
        }
    }
    // Ran out of steps without resolving: a ray still this close to the
    // photon sphere after max_steps is, for practical purposes,
    // indistinguishable from captured (it's either about to plunge, or
    // would need many more orbits to escape than the step budget
    // allows either way). Previously this returned "escaped along the
    // ORIGINAL camera-ray direction" -- a fake sky sample that has
    // nothing to do with where the ray actually got to, which produced
    // a visible hard-edged, blocky discontinuity right at the
    // resolved/unresolved boundary (reported as looking like a
    // "censored" strip of squares). Folding it into captured blends
    // smoothly into the shadow instead.
    out[i * 8 + 0] = accum[0]; out[i * 8 + 1] = accum[1]; out[i * 8 + 2] = accum[2];
    out[i * 8 + 3] = transmittance;
    out[i * 8 + 4] = 0.0; out[i * 8 + 5] = 0.0; out[i * 8 + 6] = 0.0;
    out[i * 8 + 7] = 1.0;
}

}  // namespace spatium::gpu

extern "C" void spatium_gpu_render_kerr(KerrRenderParams params,
                                         const spatium::gpu::DiscParticleParams* particles,
                                         int n_particles, const double* dir_world_in, double* out,
                                         int n_rays) {
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
    spatium::gpu::kerr_render_batch<<<blocks, threads>>>(params, d_particles, n_particles, d_dirs,
                                                           d_out, n_rays);

    cudaMemcpy(out, d_out, sizeof(double) * 8 * n_rays, cudaMemcpyDeviceToHost);
    cudaFree(d_particles);
    cudaFree(d_dirs);
    cudaFree(d_out);
}
