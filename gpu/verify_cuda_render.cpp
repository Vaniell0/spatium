// Phase 3 cross-check: schwarzschild_render_kernel.cu's full per-ray
// pipeline (geodesic + volumetric disk emission-absorption) vs. a CPU
// reference built from the SAME real library calls
// blackhole_gr_demo.cpp's trace_ray()/Schwarzschild path uses
// (geodesic_step() with the Dual<T>-exact Christoffel symbols,
// disk_redshift_factor(), blackbody_to_rgb255()) -- trace_ray() itself
// lives in the example file, not a shared header, so the orchestration
// loop below is a fresh, careful transcription of it (matching
// blackhole_gr_demo.cpp exactly as of this session), not a shared
// implementation; the physics primitives it calls ARE the real,
// already-tested library functions, not reimplemented.
//
// Deliberately stops at the same point render_kernel.h's comment
// describes: pre-sky-lookup (accumulated, transmittance, exit/capture).
// A representative frame (a plain fixed camera pose, not the full
// choreography -- that arc was already visually validated in Phase 0;
// this phase validates the PER-RAY MATH, which doesn't depend on which
// camera pose produced the ray) at 160x90, one ray per pixel center (no
// AA -- supersampling would average away the exact per-ray comparison
// this test needs).

#include "disk_physics.hpp"
#include "render_kernel.h"

#include <spatium/algebra/vector.hpp>
#include <spatium/physics/relativity/accretion_disk.hpp>
#include <spatium/physics/relativity/geodesic.hpp>
#include <spatium/physics/relativity/schwarzschild.hpp>
#include <spatium/render/camera.hpp>
#include <spatium/render/spectral.hpp>

#include <cmath>
#include <cstdio>
#include <numbers>
#include <vector>

using namespace spatium;
using namespace spatium::physics::relativity;
using namespace spatium::render;

namespace {

constexpr double kMass = 1.0;
constexpr double kRHorizon = 2.0, kRCapture = 2.3, kREscape = 40.0;
constexpr double kRDiskInner = 6.0, kRDiskOuter = 10.0;
constexpr double kDiskAspectRatio = 0.15, kDiskDensityPower = 2.0;
constexpr double kDiskT0 = 9000.0, kDiskTRef = 6500.0, kDopplerDamping = 0.6;
constexpr double kAbsorption = 0.7, kEmissionScale = 3.0, kTransmittanceCutoff = 1e-4;
constexpr double kBaseDlambda = 0.15, kDlambdaFloorFrac = 0.10;
constexpr int kMaxSteps = 3000;
constexpr double kDCam = 14.0;
constexpr double kDiskTime = 20.0;  // an arbitrary mid-sequence-ish value

std::vector<spatium::gpu::DiscParticleParams> make_disc_particles() {
    constexpr std::size_t N = 16;
    std::vector<spatium::gpu::DiscParticleParams> particles(N);
    for (std::size_t i = 0; i < N; ++i) {
        double t = double(i) / double(N);
        double r_center = kRDiskInner + (kRDiskOuter - kRDiskInner) * (0.1 + 0.8 * t);
        double ecc = 0.10 + 0.20 * std::abs(std::sin(t * 11.0));
        double r_min = r_center * (1.0 - ecc);
        double r_max = r_center * (1.0 + ecc);
        particles[i] = {1.0 / r_max, 1.0 / r_min, t * 2.0 * std::numbers::pi * 5.0,
                         1.5 + 3.0 * std::fmod(t * 7.0, 1.0)};
    }
    return particles;
}

double smoothstep01_cpu(double t) {
    t = std::clamp(t, 0.0, 1.0);
    return t * t * (3.0 - 2.0 * t);
}

double disc_density_envelope_cpu(double r, double height) {
    double scale_height = kDiskAspectRatio * r;
    double vertical = std::exp(-(height * height) / (2.0 * scale_height * scale_height));
    double radial = std::pow(r / kRDiskInner, -kDiskDensityPower);
    double inner_taper = std::clamp((r - kRHorizon) / (kRDiskInner - kRHorizon), 0.0, 1.0);
    double outer_taper = 1.0 - smoothstep01_cpu((r - kRDiskOuter) / (0.5 * kRDiskOuter));
    return vertical * radial * inner_taper * outer_taper;
}

double streak_noise_cpu(double x, double y) {
    double s = std::sin(x * 4.0 + y * 2.3) + 0.6 * std::sin(x * 7.3 - y * 5.1 + 1.3) +
               0.4 * std::sin(x * 13.0 + y * 9.0 - 0.7);
    return std::clamp(0.5 + 0.25 * s, 0.0, 1.0);
}

double disc_density_cpu(double r, double phi, double p_t,
                         const std::vector<spatium::gpu::DiscParticleParams>& particles) {
    double density = 0.0;
    for (const auto& p : particles) {
        double u_avg = (p.u1 + p.u2) * 0.5;
        double dphi_dt = u_avg * std::sqrt(0.5 * u_avg);
        double phi_p = dphi_dt * p_t + p.phi0;
        double a = std::fmod(phi - phi_p, 2.0 * std::numbers::pi);
        if (a < 0.0) a += 2.0 * std::numbers::pi;
        double s = std::sin(p.dtheta_dphi * (a + phi_p));
        double r_p = 1.0 / (p.u1 + (p.u2 - p.u1) * s * s);
        double dx = (a - std::numbers::pi) / std::numbers::pi;
        double dy = (r_p - r) * 0.5;
        double dist = std::sqrt(dx * dx + dy * dy);
        double falloff_t = std::clamp(1.0 - dist, 0.0, 1.0);
        double falloff = falloff_t * falloff_t * (3.0 - 2.0 * falloff_t);
        double noise = streak_noise_cpu(dx * (r / kRDiskOuter) * 7.0, dy * 4.0);
        density += falloff * noise;
    }
    return std::max(density, 0.05);
}

double disk_temperature_profile_cpu(double r) {
    double taper = std::max(0.0, 1.0 - std::sqrt(kRDiskInner / r));
    return kDiskT0 * std::pow(r / kRDiskInner, -0.75) * std::pow(taper, 0.25);
}

Vec<double, 3> boost_saturation_cpu(Vec<double, 3> c, double factor) {
    double gray = (c[0] + c[1] + c[2]) / 3.0;
    Vec<double, 3> b{gray + (c[0] - gray) * factor, gray + (c[1] - gray) * factor,
                      gray + (c[2] - gray) * factor};
    return {std::clamp(b[0], 0.0, 255.0), std::clamp(b[1], 0.0, 255.0), std::clamp(b[2], 0.0, 255.0)};
}

struct RayResult {
    Vec<double, 3> accum;
    double transmittance;
    Vec<double, 3> exit_dir;
    bool captured;
};

// Transcription of trace_ray() + the Schwarzschild ray_color lambda
// setup from blackhole_gr_demo.cpp (this session's version) -- see file
// header comment.
RayResult cpu_reference(const SchwarzschildMetric<double>& metric, const Vec<double, 3>& dir_world,
                         const Vec<double, 3>& r_hat, double f_cam,
                         const std::vector<spatium::gpu::DiscParticleParams>& particles) {
    double dir_dot_r = dir_world.dot(r_hat);
    Vec<double, 3> t_raw{dir_world - r_hat * dir_dot_r};
    double t_norm = t_raw.norm();
    if (t_norm < 1e-9)
        return {{0, 0, 0}, 1.0, dir_world, dir_dot_r < 0.0};

    Vec<double, 3> t_hat{t_raw / t_norm};
    double ut = 1.0 / std::sqrt(f_cam);
    double ur = dir_dot_r * std::sqrt(f_cam);
    double uphi = t_norm / kDCam;
    double b_photon = (kDCam * kDCam * uphi) / (f_cam * ut);

    Vec<double, 3> n_disk{0.0, 0.0, 1.0};
    Vec<double, 3> n_ray{r_hat.cross(t_hat)};
    Vec<double, 3> line_dir{n_ray.cross(n_disk)};
    double line_norm = line_dir.norm();
    bool disk_plane_ok = line_norm > 1e-9;
    double phi_disk = 0.0;
    if (disk_plane_ok) {
        line_dir = Vec<double, 3>{line_dir / line_norm};
        phi_disk = std::atan2(line_dir.dot(t_hat), line_dir.dot(r_hat));
    }

    Vec<double, 8> state{0.0, kDCam, std::numbers::pi / 2.0, 0.0, ut, ur, 0.0, uphi};
    Vec<double, 3> accumulated{0.0, 0.0, 0.0};
    double transmittance = 1.0;

    for (int step = 0; step < kMaxSteps; ++step) {
        double r = state[1], theta = state[2], phi = state[3];
        double dl = kBaseDlambda * std::max(kDlambdaFloorFrac, r / kDCam);
        state = geodesic_step(metric, state, dl);
        double r_new = state[1];
        if (!std::isfinite(r_new)) return {accumulated, transmittance, {0, 0, 0}, true};
        double theta_new = state[2], phi_new = state[3];

        double r_mid = 0.5 * (r + r_new), theta_mid = 0.5 * (theta + theta_new),
               phi_mid = 0.5 * (phi + phi_new);
        (void)theta_mid;
        double height = disk_plane_ok ? r_mid * std::sin(phi_mid - phi_disk) * line_norm : 0.0;
        double envelope = disc_density_envelope_cpu(r_mid, height);
        if (envelope > 1e-6) {
            double texture = disc_density_cpu(r_mid, phi_mid, kDiskTime, particles);
            double density = envelope * texture;

            double z = disk_redshift_factor(kMass, r_mid, b_photon);
            double z_display = 1.0 + kDopplerDamping * (z - 1.0);
            double t_emit = disk_temperature_profile_cpu(r_mid);
            double t_obs = t_emit / std::max(z_display, 1e-3);
            double intensity = std::clamp(std::pow(t_obs / kDiskTRef, 1.2), 0.0, 1.0);
            Vec<double, 3> emit_color = boost_saturation_cpu(
                Vec<double, 3>{blackbody_to_rgb255(t_obs) * intensity}, 1.15);

            double optical_depth = density * kAbsorption * dl;
            double emit_amount = density * kEmissionScale * dl;
            accumulated = Vec<double, 3>{accumulated + emit_color * (transmittance * emit_amount)};
            transmittance *= std::exp(-optical_depth);
        }

        if (r_new <= kRCapture || transmittance < kTransmittanceCutoff)
            return {accumulated, transmittance, {0, 0, 0}, true};

        if (r_new >= kREscape && state[5] > 0.0) {
            double f_esc = 1.0 - 2.0 * kMass / r_new;
            double nr = state[5] / std::sqrt(f_esc);
            double nphi = r_new * state[7];
            double mag = std::sqrt(nr * nr + nphi * nphi);
            nr /= mag; nphi /= mag;
            double cph = std::cos(phi_new), sph = std::sin(phi_new);
            Vec<double, 3> r_hat_exit{r_hat * cph + t_hat * sph};
            Vec<double, 3> t_hat_exit{t_hat * cph - r_hat * sph};
            Vec<double, 3> exit_dir{r_hat_exit * nr + t_hat_exit * nphi};
            return {accumulated, transmittance, exit_dir, false};
        }
    }
    return {accumulated, transmittance, dir_world, false};
}

}  // namespace

int main() {
    Camera<double> cam{.position = {kDCam, 0.0, 0.0}, .target = {0.0, 0.0, 0.0},
                        .up = {0.0, 0.0, 1.0}, .fov_deg = 40.0};
    CameraBasis<double> basis = make_camera_basis(cam);
    Vec<double, 3> r_hat{cam.position / kDCam};
    double f_cam = 1.0 - 2.0 * kMass / kDCam;

    constexpr int W = 160, H = 90;
    std::vector<spatium::gpu::DiscParticleParams> particles = make_disc_particles();

    std::vector<double> dir_world_flat(W * H * 3);
    std::vector<Vec<double, 3>> dirs(W * H);
    for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x) {
            Vec<double, 3> d = camera_pixel_dir(cam, basis, x, y, W, H);
            int idx = y * W + x;
            dirs[idx] = d;
            dir_world_flat[idx * 3 + 0] = d[0];
            dir_world_flat[idx * 3 + 1] = d[1];
            dir_world_flat[idx * 3 + 2] = d[2];
        }

    RenderParams params{};
    params.mass = kMass;
    params.D_cam = kDCam;
    params.r_hat_x = r_hat[0]; params.r_hat_y = r_hat[1]; params.r_hat_z = r_hat[2];
    params.disk_time = kDiskTime;
    params.r_horizon = kRHorizon; params.r_capture = kRCapture; params.r_escape = kREscape;
    params.r_disk_inner = kRDiskInner; params.r_disk_outer = kRDiskOuter;
    params.disk_aspect_ratio = kDiskAspectRatio; params.disk_density_power = kDiskDensityPower;
    params.disk_t0 = kDiskT0; params.disk_t_ref = kDiskTRef; params.doppler_damping = kDopplerDamping;
    params.absorption = kAbsorption; params.emission_scale = kEmissionScale;
    params.transmittance_cutoff = kTransmittanceCutoff;
    params.base_dlambda = kBaseDlambda; params.dlambda_floor_frac = kDlambdaFloorFrac;
    params.max_steps = kMaxSteps;

    std::vector<double> gpu_out(W * H * 8);
    spatium_gpu_render_schwarzschild(params, particles.data(), static_cast<int>(particles.size()),
                                      dir_world_flat.data(), gpu_out.data(), W * H);

    SchwarzschildMetric<double> metric{kMass};
    int n_mismatch = 0, n_captured_cpu = 0, n_captured_gpu = 0;
    double max_abs_diff = 0.0, max_rel_diff = 0.0;
    for (int idx = 0; idx < W * H; ++idx) {
        RayResult ref = cpu_reference(metric, dirs[idx], r_hat, f_cam, particles);
        if (ref.captured) ++n_captured_cpu;
        bool gpu_captured = gpu_out[idx * 8 + 7] > 0.5;
        if (gpu_captured) ++n_captured_gpu;

        if (ref.captured != gpu_captured) {
            ++n_mismatch;
            continue;
        }
        double got[4] = {gpu_out[idx * 8 + 0], gpu_out[idx * 8 + 1], gpu_out[idx * 8 + 2],
                          gpu_out[idx * 8 + 3]};
        double want[4] = {ref.accum[0], ref.accum[1], ref.accum[2], ref.transmittance};
        double abs_diff = 0.0, scale = 1.0;
        for (int k = 0; k < 4; ++k) {
            abs_diff = std::max(abs_diff, std::abs(got[k] - want[k]));
            scale = std::max(scale, std::abs(want[k]));
        }
        if (!ref.captured) {
            for (int k = 0; k < 3; ++k) {
                double g = gpu_out[idx * 8 + 4 + k];
                double w = (k == 0 ? ref.exit_dir[0] : k == 1 ? ref.exit_dir[1] : ref.exit_dir[2]);
                abs_diff = std::max(abs_diff, std::abs(g - w));
            }
        }
        double rel = abs_diff / scale;
        max_abs_diff = std::max(max_abs_diff, abs_diff);
        max_rel_diff = std::max(max_rel_diff, rel);
        if (rel > 1e-6) {
            ++n_mismatch;
            if (n_mismatch <= 5)
                std::printf("mismatch at pixel %d (%d,%d): abs=%.3e rel=%.3e\n", idx, idx % W, idx / W,
                            abs_diff, rel);
        }
    }

    std::printf("W=%d H=%d rays=%d\n", W, H, W * H);
    std::printf("captured: cpu=%d gpu=%d\n", n_captured_cpu, n_captured_gpu);
    std::printf("max_abs_diff=%.3e max_rel_diff=%.3e\n", max_abs_diff, max_rel_diff);
    std::printf("mismatches: %d/%d\n", n_mismatch, W * H);
    return n_mismatch == 0 ? 0 : 1;
}
