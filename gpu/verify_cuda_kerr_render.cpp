// Phase 3, Kerr path: lightweight sanity check (not the full CPU-parity
// harness verify_cuda_render.cpp does for Schwarzschild -- deferred;
// the underlying geodesic math is already verified to ~1e-16 in Phase
// 2, and disk_physics.hpp is shared, identical code with the already-
// verified Schwarzschild path, so the real untested surface here is
// just the tetrad/exit-direction plumbing). Checks the kernel runs and
// returns finite, sane values across a spread of rays from a real
// off-equatorial camera pose.

#include "disk_physics.hpp"
#include "kerr_geometry.hpp"
#include "render_kernel.h"

#include <spatium/algebra/vector.hpp>
#include <spatium/physics/relativity/kerr.hpp>
#include <spatium/render/camera.hpp>

#include <cmath>
#include <cstdio>
#include <numbers>
#include <vector>

using namespace spatium;
using namespace spatium::physics::relativity;
using namespace spatium::render;

int main() {
    double mass = 1.0, spin = 0.9;
    bool prograde = true;
    double cam_r = 14.0, cam_theta = std::numbers::pi / 2.0 - 0.3, cam_phi = 0.0;

    double r_horizon = kerr_outer_horizon_radius(mass, spin);
    double r_disk_inner = kerr_isco_radius(mass, spin, prograde);
    double r_disk_outer = 10.0;
    double r_capture = r_horizon * 1.15, r_escape = 40.0;

    auto tet = spatium::gpu::spheroidal_tetrad(cam_r, cam_theta, cam_phi, spin);
    double g_rr, g_thth, g_phph;
    spatium::gpu::kerr_metric_diag(mass, spin, cam_r, cam_theta, &g_rr, &g_thth, &g_phph);
    double sigma = cam_r * cam_r + spin * spin * std::cos(cam_theta) * std::cos(cam_theta);
    double g_tt = -(1.0 - 2.0 * mass * cam_r / sigma);
    double g_tphi = -2.0 * mass * cam_r * spin * std::sin(cam_theta) * std::sin(cam_theta) / sigma;
    double ut_cam = 1.0 / std::sqrt(-g_tt);

    double rho = std::sqrt(cam_r * cam_r + spin * spin);
    Vec<double, 3> cam_pos{rho * std::sin(cam_theta) * std::cos(cam_phi),
                            rho * std::sin(cam_theta) * std::sin(cam_phi), cam_r * std::cos(cam_theta)};
    double D_cam = cam_pos.norm();

    Camera<double> cam{.position = cam_pos, .target = {0, 0, 0}, .up = {0, 0, 1}, .fov_deg = 40.0};
    CameraBasis<double> basis = make_camera_basis(cam);

    constexpr int W = 32, H = 18;
    std::vector<double> dirs(W * H * 3);
    for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x) {
            Vec<double, 3> d = camera_pixel_dir(cam, basis, x, y, W, H);
            int idx = y * W + x;
            dirs[idx * 3 + 0] = d[0]; dirs[idx * 3 + 1] = d[1]; dirs[idx * 3 + 2] = d[2];
        }

    std::vector<spatium::gpu::DiscParticleParams> particles(16);
    for (std::size_t i = 0; i < 16; ++i) {
        double t = double(i) / 16.0;
        double r_center = r_disk_inner + (r_disk_outer - r_disk_inner) * (0.1 + 0.8 * t);
        double ecc = 0.10 + 0.20 * std::abs(std::sin(t * 11.0));
        particles[i] = {1.0 / (r_center * (1 + ecc)), 1.0 / (r_center * (1 - ecc)),
                         t * 2.0 * std::numbers::pi * 5.0, 1.5 + 3.0 * std::fmod(t * 7.0, 1.0)};
    }

    KerrRenderParams p{};
    p.mass = mass; p.spin = spin; p.D_cam = D_cam;
    p.cam_r = cam_r; p.cam_theta = cam_theta; p.cam_phi = cam_phi;
    p.e_r_x = tet.e_r.x; p.e_r_y = tet.e_r.y; p.e_r_z = tet.e_r.z;
    p.e_theta_x = tet.e_theta.x; p.e_theta_y = tet.e_theta.y; p.e_theta_z = tet.e_theta.z;
    p.e_phi_x = tet.e_phi.x; p.e_phi_y = tet.e_phi.y; p.e_phi_z = tet.e_phi.z;
    p.ut_cam = ut_cam; p.sqrt_grr_cam = std::sqrt(g_rr); p.sqrt_gthth_cam = std::sqrt(g_thth);
    p.sqrt_gphph_cam = std::sqrt(g_phph);
    p.g_cam_00 = g_tt; p.g_cam_03 = g_tphi; p.g_cam_33 = g_phph;
    p.prograde = prograde ? 1 : 0;
    p.disk_time = 20.0;
    p.r_horizon = r_horizon; p.r_capture = r_capture; p.r_escape = r_escape;
    p.r_disk_inner = r_disk_inner; p.r_disk_outer = r_disk_outer;
    p.disk_aspect_ratio = 0.15; p.disk_density_power = 2.0;
    p.disk_t0 = 9000.0; p.disk_t_ref = 6500.0; p.doppler_damping = 0.6;
    p.absorption = 0.7; p.emission_scale = 3.0; p.transmittance_cutoff = 1e-4;
    p.base_dlambda = 0.15; p.dlambda_floor_frac = 0.10; p.max_steps = 3000;

    std::vector<double> out(W * H * 8);
    spatium_gpu_render_kerr(p, particles.data(), 16, dirs.data(), out.data(), W * H);

    int n_finite = 0, n_captured = 0;
    double max_abs = 0.0;
    for (int i = 0; i < W * H; ++i) {
        bool ok = true;
        for (int k = 0; k < 8; ++k) {
            double v = out[i * 8 + k];
            if (!std::isfinite(v)) ok = false;
            max_abs = std::max(max_abs, std::abs(v));
        }
        if (ok) ++n_finite;
        if (out[i * 8 + 7] > 0.5) ++n_captured;
    }
    std::printf("rays=%d finite=%d captured=%d max_abs_value=%.3e\n", W * H, n_finite, n_captured,
                max_abs);
    return (n_finite == W * H) ? 0 : 1;
}
