#pragma once

// Device-portable transcription of blackhole_gr_demo.cpp's disk physics
// (disc_density/disc_density_envelope/disk_temperature_profile/
// blackbody_to_rgb255/boost_saturation/streak_noise) -- Phase 3 of the
// GPU port. Plain `double` throughout, matching the demo file itself
// (it isn't templated either): this is a concrete port of one specific,
// already-tuned application, not new general library material.
//
// Every function here is a straight transcription with no behavior
// changes -- correctness rests on being IDENTICAL to the CPU source,
// not on independent derivation (unlike christoffel_closed_form.hpp,
// there's no symbolic ground truth to derive this from; it's a tuned
// visual model). Cross-checked component-by-component against the CPU
// versions in gpu/verify_cuda_render.cpp.

#ifndef SPATIUM_BUILDING_MODULE
#  include <cmath>
#endif

#if defined(__CUDACC__)
#  define SPATIUM_HOST_DEVICE __host__ __device__
#else
#  define SPATIUM_HOST_DEVICE
#endif

namespace spatium::gpu {

struct DiscParticleParams {
    double u1, u2, phi0, dtheta_dphi;
};

struct Vec3d {
    double x, y, z;
};

SPATIUM_HOST_DEVICE inline double gpu_clamp(double v, double lo, double hi) {
    return v < lo ? lo : (v > hi ? hi : v);
}

SPATIUM_HOST_DEVICE inline double smoothstep01(double t) {
    t = gpu_clamp(t, 0.0, 1.0);
    return t * t * (3.0 - 2.0 * t);
}

SPATIUM_HOST_DEVICE inline double streak_noise(double x, double y) {
    using std::sin;
    double s = sin(x * 4.0 + y * 2.3) + 0.6 * sin(x * 7.3 - y * 5.1 + 1.3) +
               0.4 * sin(x * 13.0 + y * 9.0 - 0.7);
    return gpu_clamp(0.5 + 0.25 * s, 0.0, 1.0);
}

SPATIUM_HOST_DEVICE inline double disc_density(double r, double phi, double p_t,
                                                const DiscParticleParams* particles, int n_particles,
                                                double r_disk_outer) {
    using std::fmod;
    using std::sin;
    using std::sqrt;
    constexpr double pi = 3.14159265358979323846;
    double density = 0.0;
    for (int i = 0; i < n_particles; ++i) {
        const DiscParticleParams& p = particles[i];
        double u_avg = (p.u1 + p.u2) * 0.5;
        double dphi_dt = u_avg * sqrt(0.5 * u_avg);
        double phi_p = dphi_dt * p_t + p.phi0;
        double a = fmod(phi - phi_p, 2.0 * pi);
        if (a < 0.0) a += 2.0 * pi;
        double s = sin(p.dtheta_dphi * (a + phi_p));
        double r_p = 1.0 / (p.u1 + (p.u2 - p.u1) * s * s);
        double dx = (a - pi) / pi;
        double dy = (r_p - r) * 0.5;
        double dist = sqrt(dx * dx + dy * dy);
        double falloff_t = gpu_clamp(1.0 - dist, 0.0, 1.0);
        double falloff = falloff_t * falloff_t * (3.0 - 2.0 * falloff_t);
        double noise = streak_noise(dx * (r / r_disk_outer) * 7.0, dy * 4.0);
        density += falloff * noise;
    }
    double floor = 0.05;
    return density > floor ? density : floor;
}

SPATIUM_HOST_DEVICE inline double disc_density_envelope(double r, double height, double r_horizon,
                                                          double r_disk_inner, double r_disk_outer,
                                                          double disk_aspect_ratio,
                                                          double disk_density_power) {
    using std::exp;
    using std::pow;
    double scale_height = disk_aspect_ratio * r;
    double vertical = exp(-(height * height) / (2.0 * scale_height * scale_height));
    double radial = pow(r / r_disk_inner, -disk_density_power);
    double inner_taper = gpu_clamp((r - r_horizon) / (r_disk_inner - r_horizon), 0.0, 1.0);
    double outer_taper = 1.0 - smoothstep01((r - r_disk_outer) / (0.5 * r_disk_outer));
    return vertical * radial * inner_taper * outer_taper;
}

SPATIUM_HOST_DEVICE inline double disk_temperature_profile(double r, double r_disk_inner,
                                                             double disk_t0) {
    using std::pow;
    using std::sqrt;
    double taper = 1.0 - sqrt(r_disk_inner / r);
    if (taper < 0.0) taper = 0.0;
    return disk_t0 * pow(r / r_disk_inner, -0.75) * pow(taper, 0.25);
}

SPATIUM_HOST_DEVICE inline Vec3d blackbody_to_rgb255(double kelvin) {
    using std::log;
    using std::pow;
    double temp = gpu_clamp(kelvin, 1000.0, 40000.0) / 100.0;

    double red =
        (temp <= 66.0) ? 255.0 : gpu_clamp(329.698727446 * pow(temp - 60.0, -0.1332047592), 0.0, 255.0);

    double green = (temp <= 66.0)
                       ? gpu_clamp(99.4708025861 * log(temp) - 161.1195681661, 0.0, 255.0)
                       : gpu_clamp(288.1221695283 * pow(temp - 60.0, -0.0755148492), 0.0, 255.0);

    double blue;
    if (temp >= 66.0) blue = 255.0;
    else if (temp <= 19.0) blue = 0.0;
    else blue = gpu_clamp(138.5177312231 * log(temp - 10.0) - 305.0447927307, 0.0, 255.0);

    return {red, green, blue};
}

SPATIUM_HOST_DEVICE inline Vec3d boost_saturation(Vec3d c, double factor) {
    double gray = (c.x + c.y + c.z) / 3.0;
    Vec3d boosted{gray + (c.x - gray) * factor, gray + (c.y - gray) * factor,
                  gray + (c.z - gray) * factor};
    return {gpu_clamp(boosted.x, 0.0, 255.0), gpu_clamp(boosted.y, 0.0, 255.0),
            gpu_clamp(boosted.z, 0.0, 255.0)};
}

}  // namespace spatium::gpu
