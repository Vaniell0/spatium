// Rigorous-GR black hole ray-tracer (Schwarzschild or Kerr), with a real
// (tilted, lensed) accretion disk and camera choreography. Every
// perceptually-tunable constant below is CLI-overridable -- see
// print_usage() -- so visual iteration doesn't require recompiling.
//
// Built alongside examples/blackhole_demo.cpp, not replacing it.
// blackhole_demo.cpp integrates the plane-confined photon equation
// u(phi) directly (exact for Schwarzschild given a correct impact
// parameter, but derives that parameter from a flat-space approximation,
// b = D*sin(alpha), and reconstructs the escape direction by rotating
// the radial vector rather than using the photon's actual outgoing
// velocity). This file replaces both of those specific approximations
// with the exact tetrad-corrected values, and integrates the full
// 4-coordinate geodesic equation via physics/relativity/geodesic.hpp
// (Dual<T>-exact Christoffel symbols, RK4 through algebra/ode.hpp) --
// the metric-agnostic engine built this session and cross-validated in
// tests/test_relativity.cpp against blackhole_demo.cpp's own u(phi)
// integration (agreement to ~2e-5 relative).
//
// WHY A SEPARATE FILE: this renderer is meaningfully heavier per ray
// (exact Christoffel symbols recomputed at every RK4 stage, via a
// general 4x4 metric inverse) than blackhole_demo.cpp's closed-form
// ODE. Profiling during this build found the metric inverse was the
// dominant per-step cost -- fixed once, generally, by adding
// algebra/linear_solve.hpp's invert() (Gauss-Jordan, one elimination
// pass) instead of four independent solve_direct() calls -- and that a
// FIXED integration step is wasteful (curvature near the photon sphere
// needs a small step; the same step size is needlessly fine everywhere
// else). The step size below scales with r (shrinking near the BH,
// growing far away), verified against blackhole_demo.cpp's own u(phi)
// solution to stay within ~2e-4 relative agreement while cutting typical
// per-ray step counts by an order of magnitude.
//
// SCHWARZSCHILD FAST PATH (--spin 0, the default): a camera at coordinate
// radius D is modeled as a local static observer with tetrad
// e_r_hat = sqrt(f)*d/dr (f = 1-2M/D). A photon whose LOCAL direction has
// radial/tangential components (n_r, n_t) (the flat dot-product
// decomposition of the camera ray against the radial direction -- exact
// by construction, that IS the local tetrad) has 4-velocity u^t=1/sqrt(f),
// u^r=n_r*sqrt(f), u^phi=n_t/D. This lets the whole ray stay confined to
// its own orbital plane (any plane through the origin is "equatorial" by
// Schwarzschild's full spherical symmetry), so theta is fixed at pi/2 and
// the disk-plane crossing is found via the ray's own orbital-plane normal
// against the disk normal -- see height_of in ray_color_schwarzschild().
//
// KERR PATH (--spin != 0): Kerr has only axisymmetry about the spin axis,
// not full spherical symmetry, so a general camera ray is NOT confined to
// any fixed plane -- theta becomes a real, independent dynamical
// variable, and the whole state must be genuinely 3D. The camera's flat
// choreography position (already a 3-vector) is placed into
// Boyer-Lindquist coordinates via the EXACT oblate-spheroidal embedding
// (X,Y,Z) = (sqrt(r^2+a^2)sin(theta)cos(phi), ..., r*cos(theta)) -- not
// an approximation, that IS the coordinate system's defining relation --
// then the camera ray direction is decomposed against the metric's own
// orthonormal spatial tetrad (e_r_hat, e_theta_hat, e_phi_hat), which is
// exactly diagonal because g_r_theta = g_r_phi = g_theta_phi = 0
// identically in Boyer-Lindquist form (only g_t_phi is off-diagonal).
// This tetrad construction ignores frame-dragging at the camera itself
// (a true static observer, not a ZAMO) -- valid wherever g_tt<0, and the
// resulting error is of order a*M/r^3, negligible at this render's camera
// distances (>=14M) for any spin up to extremal. The disk, ISCO and
// photon-orbit radii DO use the exact Kerr (Bardeen-Press-Teukolsky 1972)
// formulas throughout, since that's where the real spin-dependent
// physics -- the tighter, narrower photon ring the user asked about --
// actually lives. ray_color_schwarzschild() and ray_color_kerr() both
// delegate their shared emission-absorption integration to trace_ray(),
// so the two paths can't silently drift apart on the disk physics; they
// differ only in how the initial/exit tetrads and disk-height test are
// built, which is where Schwarzschild's extra symmetry actually pays off.
//
// ACCRETION DISK / HALO: NOT a single hit-test against a thin plane, and
// NOT a separately painted "glow" near the photon sphere -- both of
// those were tried and both read as artificial (direct comparison
// against a reference Interstellar frame, plus a pointer to the actual
// technique real GRRT (general-relativistic radiative-transfer) codes
// use -- Porth et al. 2017's BHAC paper, arxiv.org/abs/1611.09720, and
// the structured-torus GRRT literature -- made the fix concrete). The
// disk is a genuine 3D volumetric density field, integrated via
// emission-absorption radiative transfer AT EVERY STEP of the geodesic,
// not just once: disc_density_envelope() gives the standard isothermal-
// hydrostatic-equilibrium profile (Gaussian falloff in height above the
// midplane, scale height growing linearly with r -- a flared disk --
// times a radial power law, no hard cutoff), and disc_density() layers
// real streak/turbulence structure on top, adapted from Eric Bruneton's
// open-source real-time black hole renderer
// (github.com/ebruneton/black_hole_shader, model.glsl's
// DefaultDiscColor()): many overlapping elliptical, precessing "clumps"
// of material summed together, not a smooth gradient with noise
// multiplied on top. disk_redshift_factor()/kerr_disk_redshift_factor()
// (gravitational + orbital Doppler, the standard thin-disk result also
// used by Luminet 1979 and the "Interstellar" visualization paper, James
// et al. 2015) still provide the exact GR frequency shift at every sample
// point, damped before it reaches the display color (DOPPLER_DAMPING)
// the same way Double Negative documented dialing down Gargantua's
// physically-exact brightness asymmetry for legibility -- though damped,
// not erased: a fully-flattened asymmetry read as "Doppler effect is
// absent" in review, so this stays well short of 1.0 without going back
// to the undamped, "confusing" full value either.
//
// THE MULTI-ORBIT PHOTON RING: the bright, thin ring exactly at the
// shadow's edge comes from photons that wind several times near the
// photon sphere before finally crossing the disk plane -- each winding
// is a higher-order lensed image of the disk, at the disk's own real
// radii, not material actually sitting at the photon sphere. In a review
// pass this ring read as "weakly expressed, blending into the general
// mass": the root cause was transmittance decaying too fast (a real
// consequence of ABSORPTION and the transmittance floor, not a missing
// feature) -- by the time a ray reached its second or third disk-plane
// crossing near the photon sphere, transmittance had already collapsed
// below the cutoff, extinguishing exactly the crossings that produce the
// ring before they could contribute. Lowering ABSORPTION and the
// transmittance floor lets those later windings actually survive to add
// their own (highly redshifted/lensed) light.
//
// CAMERA CHOREOGRAPHY: an establishing distance (~1/3 frame width, not
// forced close), a flyby to a near-peak distance (the ~4/5-width
// dominant framing happens HERE, as a moment, not as the baseline), and
// a continuous "weak focus" throughout -- the look direction blends the
// live direction-to-BH with one frozen early in the sequence, so the
// black hole's screen position is mostly a consequence of the camera's
// own arcing motion (natural parallax) rather than the camera hunting to
// re-center it. The elevation angle stays fairly oblique (ELEV_PEAK_DEG,
// default 30 -- lowered from an earlier, more overhead 42, since a
// steeper/more-overhead view of the disk visibly suppresses both the
// Doppler brightness asymmetry and the ring's lensed wrap, which is
// exactly the "massive"/"Interstellar" look being aimed for). The flyby's
// closest point still clears R_DISK_OUTER with real margin (a camera
// positioned at or inside the disk's own radial extent registers a
// spurious near-face-on disk hit within the first few integration steps
// for almost every ray, hiding the shadow entirely -- found by rendering
// a mid-choreography frame during tuning), and stays well outside the
// photon sphere rather than attempting the plan's "optional" photon-
// sphere flythrough -- a camera that close would put most of the frame
// inside the expensive near-critical regime for a marginal shot.
//
// Output: blackhole_gr_frames/frame_%04d.png -- assemble separately:
//   ffmpeg -framerate 30 -i blackhole_gr_frames/frame_%04d.png -pix_fmt yuv420p blackhole_gr.mp4

// Only `#define` here -- render/write_image.hpp is the one place that
// includes <spatium/vendor/stb_image_write.h> in this file.
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "io_helpers.hpp"

#include <spatium/algebra/vector.hpp>
#include <spatium/physics/relativity/accretion_disk.hpp>
#include <spatium/physics/relativity/geodesic.hpp>
#include <spatium/physics/relativity/kerr.hpp>
#include <spatium/physics/relativity/schwarzschild.hpp>
#include <spatium/render/camera.hpp>
#include <spatium/render/parallel_for_rows.hpp>
#include <spatium/render/sky.hpp>
#include <spatium/render/spectral.hpp>
#include <spatium/render/supersample.hpp>
#include <spatium/render/write_image.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <numbers>
#include <print>
#include <string_view>
#include <vector>

using spatium::Matrix;
using spatium::Vec;
using spatium::physics::relativity::disk_redshift_factor;
using spatium::physics::relativity::geodesic_step;
using spatium::physics::relativity::keplerian_omega;
using spatium::physics::relativity::kerr_disk_redshift_factor;
using spatium::physics::relativity::kerr_isco_radius;
using spatium::physics::relativity::kerr_outer_horizon_radius;
using spatium::physics::relativity::kerr_photon_orbit_radius;
using spatium::physics::relativity::killing_angular_momentum;
using spatium::physics::relativity::killing_energy;
using spatium::physics::relativity::KerrMetric;
using spatium::physics::relativity::schwarzschild_isco_radius;
using spatium::physics::relativity::schwarzschild_photon_sphere_radius;
using spatium::physics::relativity::SchwarzschildMetric;
using spatium::render::make_starfield;
using spatium::render::sample_sky_color;
using spatium::render::Sky;
using spatium::render::blackbody_to_rgb255;
using spatium::render::Camera;
using spatium::render::make_camera_basis;
using spatium::render::camera_ray_dir;
using spatium::render::parallel_for_rows;
using spatium::render::supersample_pixel_hdr;
using spatium::render::write_png_rgb;

namespace {

// ---------------------------------------------------------------------
// Tunable render constants. All plain (non-constexpr) so main() can
// override them from argv before anything below reads them --
// recompute_derived() must run exactly once, after CLI parsing and
// before the frame loop, to keep the M-dependent radii and the "in units
// of M" camera/disk distances consistent with a possibly-overridden mass.
// ---------------------------------------------------------------------

int W = 1280;
int H = 720;
// rays/pixel = AA_LEVEL^2. A bump to 3 (9 rays/pixel) was tried after a
// review complaint about missing anti-aliasing, but combined with the
// lower ABSORPTION/TRANSMITTANCE_CUTOFF below (each ray now runs longer
// before an early exit) it pushed a single 480x270 test frame past 6
// minutes on 12 cores -- for a >=750-frame final render that's weeks,
// not hours. Back to 2 (supersample.hpp's own kDefaultAA) as the
// default; --aa 3+ is still there for an occasional high-quality still.
int AA_LEVEL = 2;

double M_BH = 1.0;
double SPIN = 0.0;      // Kerr a; 0 selects the Schwarzschild fast path exactly
bool PROGRADE = true;   // disk corotation direction relative to the spin axis

// g_rr diverges at the horizon (Delta=0 for Kerr, 1-2M/r=0 for
// Schwarzschild), and so does the Christoffel tensor built from it -- a
// razor-thin capture margin let the adaptive step overshoot into that
// numerically unstable region before the per-step capture check could
// fire, corrupting the state to NaN (found via a probe at a closer camera
// distance). 1.15x gives the integrator real room to detect capture
// before things blow up -- visually indistinguishable, since nothing
// this deep is ever visible regardless of the exact cutoff. All three
// are recomputed by recompute_derived(), not read before it runs.
double R_HORIZON = 2.0;
double R_CAPTURE = 2.3;
double R_ESCAPE = 40.0;

// Step size at coordinate radius r is BASE_DLAMBDA * max(FLOOR_FRAC,
// r/D_cam) -- shrinks near the BH (strong curvature), grows far away
// (gentle curvature, safe to take large RK4 steps). Verified via a
// standalone probe against blackhole_demo.cpp's u(phi) solution before
// picking these constants: ~1e-4 to 2e-4 relative agreement across
// b = 5.3M..20M, typical per-ray step counts in the 600-1600 range (an
// order of magnitude below a naive fixed-step budget).
constexpr double BASE_DLAMBDA = 0.15;
constexpr double DLAMBDA_FLOOR_FRAC = 0.10;
constexpr long MAX_GEODESIC_STEPS = 3000;

double R_PHOTON_SPHERE = 3.0;  // recomputed: Schwarzschild 3M, or Kerr's spin/prograde-dependent value
double R_DISK_INNER = 6.0;     // recomputed: the ISCO, same spin dependence
// "In units of M" until recompute_derived() scales it in place -- see
// that function's own comment.
// Widened from an earlier 10 -- reported as reading small/lacking a
// sense of scale; a broader disk against the shadow communicates mass
// more directly than a size-accurate but visually modest ring does.
double R_DISK_OUTER = 16.0;

// Capped well below a real thin disk's actual peak (which reaches tens
// of thousands of K and reads as blue-white) -- the Interstellar
// reference frame shows no blue anywhere, just warm white fading to
// gold/brown, so the blackbody curve is kept on its warm side throughout.
double DISK_T0 = 9000.0;     // Kelvin scale just outside the ISCO
double DISK_T_REF = 6500.0;  // brightness-normalization reference

// disk_redshift_factor()/kerr_disk_redshift_factor()'s exact GR frequency
// shift is damped before it reaches display color -- see file header.
double DOPPLER_DAMPING = 0.6;
// Optical-thickness knobs for the volumetric emission-absorption
// integration -- see file header's MULTI-ORBIT PHOTON RING note for why
// ABSORPTION and the transmittance floor both came down from their
// original, ring-extinguishing values.
double ABSORPTION = 0.7;
double EMISSION_SCALE = 3.0;
double TRANSMITTANCE_CUTOFF = 1e-4;

// A true blackbody colorimetrically desaturates toward white above
// ~6600K -- physically correct, but it read as a flat grey/white sheet
// against the starfield rather than a disk with visible color. This is
// tone-mapping for display, same category of choice as the intensity
// exponent above, not a physics question: push each disk pixel's RGB
// away from its own greyscale value before display so the hot inner
// disk keeps a visible blue-white cast instead of crushing to neutral
// white, and the cooler outer disk's orange reads more vividly.
Vec<double, 3> boost_saturation(Vec<double, 3> c, double factor) {
    double gray = (c[0] + c[1] + c[2]) / 3.0;
    Vec<double, 3> boosted{gray + (c[0] - gray) * factor, gray + (c[1] - gray) * factor,
                            gray + (c[2] - gray) * factor};
    return Vec<double, 3>{std::clamp(boosted[0], 0.0, 255.0), std::clamp(boosted[1], 0.0, 255.0),
                           std::clamp(boosted[2], 0.0, 255.0)};
}

// Baseline framing stays at the establishing-shot distance (~1/3 frame
// width) -- the dominant ~4/5-width framing an earlier cut tried to
// force as a permanent camera position instead comes from the
// choreography itself: the camera flies past close to the black hole
// at its nearest approach, exactly the way an actual flythrough would.
double FOV_DEG = 40.0;

// D_PEAK (closest approach) must clear R_DISK_OUTER with real margin: a
// camera positioned at or inside the disk's own radial extent registers
// a disk hit within the first few integration steps for almost every
// ray, painting a faint near-face-on smear over the whole frame and
// hiding the shadow entirely -- found by rendering a mid-choreography
// frame where the camera distance had been left at 10M against a 12M
// disk. All three, plus R_DISK_OUTER above, are "in units of M" until
// recompute_derived() scales them in place.
double D_APPROACH_START = 45.0;
double D_PEAK = 14.0;
double D_DEPART_END = 30.0;
double ELEV_START_DEG = 12.0;
// Lowered from an earlier 42 -- see file header CAMERA CHOREOGRAPHY note.
double ELEV_PEAK_DEG = 30.0;
// A modest sweep, not a near-full orbit: this choreography is a single
// flyby pass along an arc, not an orbit around a stationary subject --
// too wide a sweep would fight the fixed-ish look direction below.
double AZIMUTH_TOTAL_DEG = 80.0;

// "Weak focus": the camera's look direction is a blend of the live
// direction to the black hole and a direction frozen early in the
// sequence, weighted toward the frozen one. This is deliberately NOT a
// hard switch from full tracking to no tracking -- with real (if
// gentle) tracking blended in throughout, the black hole's screen
// position is mostly a consequence of the camera's own arcing motion
// (natural parallax as it passes close and then continues on), the way
// it would look through a spacecraft window on an actual flyby: small
// and off-center early, swelling to dominate the frame at closest
// approach, then sliding across and out of frame as the camera moves
// past -- not a subject the camera hunts to keep centered.
double FOCUS_BLEND = 0.22;
double FOCUS_FREEZE_T = 0.12;

int N_FRAMES = 150;
// Raised from an earlier 0.6 (tuned for "half an orbit over the full
// sequence") -- reported as barely-perceptible rotation, since the
// camera is also far away (small disk, little time elapsed) for most
// of that half-orbit budget. ~2 (up from ~0.5) full orbits for the
// fastest (innermost) particles gives real, visible motion right
// through the close-approach section where the disk is actually big
// enough on screen to show it.
double DISK_TIME_PER_FRAME = 2.0;

// Raised from an earlier 0.15 -- a thicker disk reads as more
// voluminous/massive, alongside the widened R_DISK_OUTER above.
double DISK_ASPECT_RATIO = 0.22;
// Lowered from an earlier 2.5 -- a shallower falloff keeps visible
// material further out, addressing a review complaint that the ring
// should be wider/more spread rather than a tight, sharply-cut band.
double DISK_DENSITY_POWER = 2.0;

// Called exactly once in main(), after CLI parsing, before anything
// below reads these values. D_APPROACH_START/D_PEAK/D_DEPART_END and
// R_DISK_OUTER are CLI-set as dimensionless multiples of M (matching
// this file's existing "* M_BH" convention from before these became
// runtime-configurable) and get scaled here to actual radii exactly
// once -- calling this twice would double-scale them.
void recompute_derived() {
    if (SPIN != 0.0) {
        R_HORIZON = kerr_outer_horizon_radius(M_BH, SPIN);
        R_PHOTON_SPHERE = kerr_photon_orbit_radius(M_BH, SPIN, PROGRADE);
        R_DISK_INNER = kerr_isco_radius(M_BH, SPIN, PROGRADE);
    } else {
        R_HORIZON = 2.0 * M_BH;
        R_PHOTON_SPHERE = schwarzschild_photon_sphere_radius(M_BH);
        R_DISK_INNER = schwarzschild_isco_radius(M_BH);
    }
    R_CAPTURE = R_HORIZON * 1.15;
    R_ESCAPE = 40.0 * M_BH;

    R_DISK_OUTER *= M_BH;
    D_APPROACH_START *= M_BH;
    D_PEAK *= M_BH;
    D_DEPART_END *= M_BH;
}

double smoothstep01(double t) {
    t = std::clamp(t, 0.0, 1.0);
    return t * t * (3.0 - 2.0 * t);
}

// Shakura-Sunyaev-style thin-disk temperature: a power-law falling off
// with radius, tapered to exactly zero at the ISCO (the standard
// zero-torque inner boundary condition) so the disk fades to black at
// its own inner edge instead of an artificial hard cutoff.
double disk_temperature_profile(double r) {
    double taper = std::max(0.0, 1.0 - std::sqrt(R_DISK_INNER / r));
    return DISK_T0 * std::pow(r / R_DISK_INNER, -0.75) * std::pow(taper, 0.25);
}

// A first attempt at texture (sine-sum modulation of the smooth
// temperature profile, the same technique io_helpers.hpp uses for
// nebula/cloud texture) made the disk read as fragmented, disconnected
// patches rather than one coherent glowing structure -- backwards from
// professional references, which read as coherent BECAUSE they simplify
// (Double Negative's published account of rendering Gargantua for
// Interstellar describes deliberately dialing down the Doppler-driven
// brightness asymmetry a physically exact disk would have, for
// legibility). The real source of the reference's streaky, flowing
// look isn't texture on a smooth gradient at all -- it's structural:
// Eric Bruneton's open-source real-time black hole renderer
// (github.com/ebruneton/black_hole_shader, black_hole/model.glsl's
// DefaultDiscColor()) builds the disk from a SUM of many distinct
// "clumps" of material, each on its own slightly-elliptical, precessing
// orbit, each contributing a soft glow that fades away from its own
// centerline. Adapted here with the same structure: geometry stays this
// file's own (Schwarzschild/Kerr geodesics, redshift, blackbody color all
// untouched), only the disk's spatial density pattern is rebuilt on this
// reference's actual technique instead of a guess. 20 particles (up from
// an earlier 14) after a review complaint that the disk read too flat
// and uniform -- more overlapping streaks breaks that up.
struct DiscParticleParams { double u1, u2, phi0, dtheta_dphi; };

std::vector<DiscParticleParams> make_disc_particles() {
    // 20 was tried (up from 14) for more turbulence; combined with the
    // performance issue documented at AA_LEVEL's declaration, dialed
    // back to 16 as a middle ground -- most of the visual benefit of
    // more overlapping streaks, without scaling disc_density()'s
    // per-hit cost as far. A later attempt at 24 (this session) was
    // reverted -- reported worse, not better.
    constexpr std::size_t N = 16;
    std::vector<DiscParticleParams> particles(N);
    for (std::size_t i = 0; i < N; ++i) {
        double t = double(i) / double(N);
        double r_center = R_DISK_INNER + (R_DISK_OUTER - R_DISK_INNER) * (0.1 + 0.8 * t);
        double ecc = 0.10 + 0.20 * std::abs(std::sin(t * 11.0));
        double r_min = r_center * (1.0 - ecc);
        double r_max = r_center * (1.0 + ecc);
        particles[i] = DiscParticleParams{
            1.0 / r_max, 1.0 / r_min, t * 2.0 * std::numbers::pi * 5.0,
            1.5 + 3.0 * std::fmod(t * 7.0, 1.0)};
    }
    return particles;
}
// Filled by recompute_derived()'s caller in main(), once R_DISK_INNER/
// R_DISK_OUTER are final -- NOT a static initializer, since both of
// those depend on the (possibly CLI-overridden) mass/spin.
std::vector<DiscParticleParams> DISC_PARTICLES;

// Cheap value-noise substitute (sine-sum, matching io_helpers.hpp's own
// "no permutation tables" technique) for texture within one streak's
// falloff envelope -- not the source of the streak shape itself, which
// comes from the orbital geometry above. A 4-octave fBm variant (this
// session) was tried for a more organic/turbulent look and reverted --
// reported worse, not better (made the shadow-boundary spike render as
// blocky squares instead of a thin line, among other regressions).
double streak_noise(double x, double y) {
    double s = std::sin(x * 4.0 + y * 2.3) + 0.6 * std::sin(x * 7.3 - y * 5.1 + 1.3) +
               0.4 * std::sin(x * 13.0 + y * 9.0 - 0.7);
    return std::clamp(0.5 + 0.25 * s, 0.0, 1.0);
}

// Sum of the particles' contributions at disk-plane position (r, phi)
// and time p_t (p_t drives each particle's own Keplerian orbital
// advance). Mirrors DefaultDiscColor()'s density loop. The noise
// argument's radial scaling (7.0, up from an earlier 4.0) increases the
// effective spatial frequency of the streak texture with radius, so the
// outer disk reads as more scattered/grainy than the smoother inner
// region -- a review complaint asked for the halo to look "wider and
// more scattered the further out, like it's actually swirling."
double disc_density(double r, double phi, double p_t) {
    double density = 0.0;
    for (const auto& p : DISC_PARTICLES) {
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
        double noise = streak_noise(dx * (r / R_DISK_OUTER) * 7.0, dy * 4.0);
        density += falloff * noise;
    }
    // A small ambient floor: a literal zero here would commit a disk hit
    // to pure black between streaks, which reads as a hole punched in
    // the disk rather than sparse, wispy gas.
    return std::max(density, 0.05);
}

// Volumetric envelope: isothermal-hydrostatic-equilibrium vertical
// profile (density falls off as a Gaussian in height above the disk
// midplane, scale height growing linearly with r -- a constant-aspect-
// ratio "flared" disk) times a radial power law -- the standard
// simplified density model behind general-relativistic radiative-
// transfer renders (see e.g. the structured-torus GRRT literature this
// session researched directly, following the BHAC paper reference).
// Softly tapered to zero beyond ~1.5*R_DISK_OUTER (real disks have some
// genuine outer edge, wherever the accretion flow's specific angular
// momentum stops matching this simplified power law) -- an EARLIER,
// literally cutoff-free version of this function let the radial power
// law's slowly-decaying tail stay non-negligible all the way out to the
// camera's own orbit (measured: at r=45M with DISK_DENSITY_POWER=2.5,
// the tail was still ~0.0065, far above the >1e-6 "worth evaluating"
// gate below), so disc_density()'s expensive 20-particle texture sum was
// firing on ~74% of ALL geodesic steps for every ray, not just genuine
// disk crossings -- found via a standalone probe after a Kerr test
// render ran far slower than expected. The outer taper below drops that
// to ~7%, an order of magnitude, with no visible change to the disk
// itself (it only removes density values already far too faint to see).
double disc_density_envelope(double r, double height) {
    double scale_height = DISK_ASPECT_RATIO * r;
    double vertical = std::exp(-(height * height) / (2.0 * scale_height * scale_height));
    double radial = std::pow(r / R_DISK_INNER, -DISK_DENSITY_POWER);
    double inner_taper = std::clamp((r - R_HORIZON) / (R_DISK_INNER - R_HORIZON), 0.0, 1.0);
    double outer_taper = 1.0 - smoothstep01((r - R_DISK_OUTER) / (0.5 * R_DISK_OUTER));
    return vertical * radial * inner_taper * outer_taper;
}

// Soft HDR-style bloom: the Interstellar reference frame is diffuse and
// glowing throughout, not sharp-edged the way this renderer's raw,
// independent per-pixel colors are. A separable box blur (two passes,
// approximating a Gaussian) of the image's bright regions, added back
// additively, is the standard cheap post-process for this -- the same
// role a real bloom pass plays in the "Gargantua With HDR Bloom"
// ShaderToy reference found alongside the disc-density technique above.
// Operates on the unclamped HDR buffer (img may hold values well past
// 255 -- see trace_ray()'s own comment on why clamping earlier than
// this was a real bug, not a simplification): a pixel at 1500 and one
// at 300 both used to read as an identical clamped 255 by the time
// bloom saw them, so bloom could only blur an already-flat block. Here
// they stay genuinely different until the single tone-map pass in
// main(), after this function returns.
void apply_bloom(std::vector<double>& img, int w, int h) {
    // THRESHOLD raised (140->190) and STRENGTH lowered (0.6->0.35):
    // reported as reading "glossy"/plastic-smooth with no texture --
    // real cause, not a subjective nitpick: most of the disk body sits
    // above 140 HDR, so almost the whole disk (not just its brightest
    // peak) was getting box-blurred and added back onto itself, wiping
    // out disc_density()'s own streak texture under a soft, textureless
    // glow. Raising the threshold confines bloom to genuinely
    // saturated peaks (a real highlight halo), lowering strength keeps
    // it a subtle addition instead of the dominant contribution.
    constexpr double THRESHOLD = 190.0;
    constexpr int RADIUS = 6;
    constexpr double STRENGTH = 0.35;

    std::vector<double> bright(img.size(), 0.0);
    for (std::size_t i = 0; i < img.size() / 3; ++i) {
        double r = img[3 * i], g = img[3 * i + 1], b = img[3 * i + 2];
        double lum = 0.299 * r + 0.587 * g + 0.114 * b;
        double scale = std::max(0.0, lum - THRESHOLD) / (lum + 1e-6);
        bright[3 * i] = r * scale;
        bright[3 * i + 1] = g * scale;
        bright[3 * i + 2] = b * scale;
    }

    std::vector<double> tmp(bright.size());
    auto box_blur_pass = [&] {
        for (int y = 0; y < h; ++y) {
            for (int x = 0; x < w; ++x) {
                double sum[3] = {0.0, 0.0, 0.0};
                int count = 0;
                for (int dx = -RADIUS; dx <= RADIUS; ++dx) {
                    int xx = x + dx;
                    if (xx < 0 || xx >= w) continue;
                    std::size_t idx = 3 * (static_cast<std::size_t>(y) * w + xx);
                    sum[0] += bright[idx]; sum[1] += bright[idx + 1]; sum[2] += bright[idx + 2];
                    ++count;
                }
                std::size_t idx = 3 * (static_cast<std::size_t>(y) * w + x);
                tmp[idx] = sum[0] / count; tmp[idx + 1] = sum[1] / count; tmp[idx + 2] = sum[2] / count;
            }
        }
        for (int x = 0; x < w; ++x) {
            for (int y = 0; y < h; ++y) {
                double sum[3] = {0.0, 0.0, 0.0};
                int count = 0;
                for (int dy = -RADIUS; dy <= RADIUS; ++dy) {
                    int yy = y + dy;
                    if (yy < 0 || yy >= h) continue;
                    std::size_t idx = 3 * (static_cast<std::size_t>(yy) * w + x);
                    sum[0] += tmp[idx]; sum[1] += tmp[idx + 1]; sum[2] += tmp[idx + 2];
                    ++count;
                }
                std::size_t idx = 3 * (static_cast<std::size_t>(y) * w + x);
                bright[idx] = sum[0] / count; bright[idx + 1] = sum[1] / count; bright[idx + 2] = sum[2] / count;
            }
        }
    };
    box_blur_pass();
    box_blur_pass();

    for (std::size_t i = 0; i < img.size(); ++i) img[i] += STRENGTH * bright[i];
}

// Single tone-map + clamp, run once at the very end (after bloom), not
// per-ray -- see trace_ray()'s comment for why clamping earlier flattens
// distinct bright values together before bloom can tell them apart.
// Identity below KNEE (everything this session's disk/color tuning
// already targets stays untouched); above it, a smooth exponential
// rolls extreme HDR values off toward 255 instead of hard-clipping --
// C1-continuous at the knee (HEADROOM chosen so the rolloff's initial
// slope is exactly 1, matching identity's slope) so there's no visible
// seam where the curve switches on.
std::vector<std::uint8_t> tone_map_to_8bit(const std::vector<double>& hdr) {
    // KNEE fixes where the smooth rolloff starts (identity below it,
    // matching this session's whole disk/color calibration); HEADROOM
    // fixes where the curve asymptotes (always 255, so it never
    // over/undershoots the display range). SCALE is independent of both
    // -- it's how much INPUT range the curve takes to traverse that
    // headroom. A first version set SCALE=HEADROOM (forced for an exact
    // C1-continuous knee), which measured out to an effective ceiling
    // around c~400 -- but a real frame at close camera range measured
    // max=1235 with 36% of ALL channels above 300, so most of the
    // overexposed image was landing in that already-saturated tail,
    // indistinguishable from pure white. SCALE=400 spreads that same
    // 300-1235 range across roughly 213-251 in output instead, trading
    // a barely-visible derivative kink at the knee (imperceptible in
    // practice, the same tradeoff real filmic tone curves make) for
    // actually showing structure across the range that matters.
    constexpr double KNEE = 200.0;
    constexpr double HEADROOM = 255.0 - KNEE;
    constexpr double SCALE = 400.0;
    std::vector<std::uint8_t> out(hdr.size());
    for (std::size_t i = 0; i < hdr.size(); ++i) {
        double c = hdr[i];
        double mapped = c <= KNEE ? c : KNEE + HEADROOM * (1.0 - std::exp(-(c - KNEE) / SCALE));
        out[i] = static_cast<std::uint8_t>(std::clamp(mapped, 0.0, 255.0));
    }
    return out;
}

// ---------------------------------------------------------------------
// Oblate-spheroidal (Boyer-Lindquist <-> quasi-Cartesian) machinery --
// only needed by the Kerr path. r, theta, phi surfaces are mutually
// orthogonal in BOTH the ambient flat embedding below and under the
// actual curved metric (g_r_theta = g_r_phi = g_theta_phi = 0
// identically in BL form) -- a genuinely orthogonal coordinate system,
// not a coincidence, so normalizing these three ambient tangent vectors
// gives an exact orthonormal frame for decomposing a flat screen-space
// ray direction into (n_r, n_theta, n_phi) local components.
// ---------------------------------------------------------------------

struct SpheroidalTetrad { Vec<double, 3> e_r, e_theta, e_phi; };

SpheroidalTetrad spheroidal_tetrad(double r, double theta, double phi, double a) {
    double s = std::sin(theta), c = std::cos(theta);
    double cp = std::cos(phi), sp = std::sin(phi);
    double rho = std::sqrt(r * r + a * a);
    Vec<double, 3> e_r{(r / rho) * s * cp, (r / rho) * s * sp, c};
    Vec<double, 3> e_th{rho * c * cp, rho * c * sp, -r * s};
    Vec<double, 3> e_ph{-rho * s * sp, rho * s * cp, 0.0};
    return SpheroidalTetrad{e_r.normalized(), e_th.normalized(), e_ph.normalized()};
}

// Exact inverse of the oblate-spheroidal embedding X=rho*sin(theta)*
// cos(phi), Y=rho*sin(theta)*sin(phi), Z=r*cos(theta) (rho=sqrt(r^2+a^2))
// -- the defining relation between Boyer-Lindquist and quasi-Cartesian
// coordinates, not an approximation. Used only to place the CAMERA (a
// fixed 3-vector from the existing flat choreography) into (r,theta,phi);
// ray-direction decomposition uses the exact metric-orthonormal tetrad
// above, not this embedding.
struct SpheroidalPoint { double r, theta, phi; };

SpheroidalPoint to_spheroidal(const Vec<double, 3>& pos, double a) {
    double X = pos[0], Y = pos[1], Z = pos[2];
    double r2sum = X * X + Y * Y + Z * Z;
    double a2 = a * a;
    double r2 = 0.5 * (r2sum - a2) + std::sqrt(0.25 * (r2sum - a2) * (r2sum - a2) + a2 * Z * Z);
    double r = std::sqrt(std::max(r2, 1e-12));
    double theta = std::acos(std::clamp(Z / r, -1.0, 1.0));
    double phi = std::atan2(Y, X);
    return SpheroidalPoint{r, theta, phi};
}

// ---------------------------------------------------------------------
// Shared per-ray integration: volumetric emission-absorption along the
// geodesic, identical for both the Schwarzschild and Kerr paths. Only
// the initial/exit tetrad construction and the disk-height test differ
// between the two callers (ray_color_schwarzschild/ray_color_kerr) --
// injected here as small callables so the disk physics itself can't
// silently drift between the two paths. b_photon comes from the
// conserved Killing-vector quantities (killing_angular_momentum/
// killing_energy), which hold for ANY stationary axisymmetric metric --
// not a Schwarzschild-specific shortcut, so this same line already
// covers Kerr's g_t_phi cross term correctly.
// ---------------------------------------------------------------------

template<typename Metric, typename HeightFn, typename RedshiftFn, typename ExitDirFn>
Vec<double, 3> trace_ray(const Metric& metric, Vec<double, 8> state, double disk_time,
                          const Sky& sky, double D_cam, const Vec<double, 3>& dir_world,
                          HeightFn&& height_of, RedshiftFn&& redshift_of, ExitDirFn&& exit_dir_of) {
    double b_photon = killing_angular_momentum(metric, state) / killing_energy(metric, state);

    double r_min = state[1];
    Vec<double, 3> accumulated{0.0, 0.0, 0.0};
    double transmittance = 1.0;

    // No clamp here, deliberately -- this return value stays in
    // unclamped HDR range all the way through apply_bloom(); clamping
    // this early flattens every sufficiently bright pixel to the same
    // value BEFORE bloom sees it, which was the actual cause of a
    // solid-white blowout with no internal structure at close camera
    // range (found 2026-08-26). A single tone-map + clamp happens once,
    // in main(), after bloom.
    for (long step = 0; step < MAX_GEODESIC_STEPS; ++step) {
        double r = state[1];
        double theta = state[2];
        double phi = state[3];
        double dl = BASE_DLAMBDA * std::max(DLAMBDA_FLOOR_FRAC, r / D_cam);
        state = geodesic_step(metric, state, dl);
        double r_new = state[1];
        // Defensive: even with R_CAPTURE's margin, treat a non-finite
        // step (numerical blowup right at the horizon) as capture --
        // return whatever emission was already accumulated, add nothing
        // further.
        if (!std::isfinite(r_new)) return accumulated;
        if (r_new < r_min) r_min = r_new;
        double theta_new = state[2];
        double phi_new = state[3];

        // Volumetric emission-absorption integration, evaluated at the
        // step's midpoint: no separate "disk hit" test and no painted
        // glow. Material anywhere along the bent path -- including near
        // the photon sphere, where multi-orbit windings produce the
        // bright ring -- contributes exactly like the actual dust the
        // reference frame shows. dlambda (the affine-parameter step)
        // stands in for the path-length differential here -- a
        // practical, not rigorously calibrated, choice for a visually-
        // tuned integration.
        double r_mid = 0.5 * (r + r_new);
        double theta_mid = 0.5 * (theta + theta_new);
        double phi_mid = 0.5 * (phi + phi_new);
        double height = height_of(r_mid, theta_mid, phi_mid);
        double envelope = disc_density_envelope(r_mid, height);
        if (envelope > 1e-6) {
            double texture = disc_density(r_mid, phi_mid, disk_time);
            double density = envelope * texture;

            double z = redshift_of(r_mid, b_photon);
            double z_display = 1.0 + DOPPLER_DAMPING * (z - 1.0);
            double t_emit = disk_temperature_profile(r_mid);
            double t_obs = t_emit / std::max(z_display, 1e-3);
            double intensity = std::clamp(std::pow(t_obs / DISK_T_REF, 1.2), 0.0, 1.0);
            Vec<double, 3> emit_color = boost_saturation(
                Vec<double, 3>{blackbody_to_rgb255(t_obs) * intensity}, 1.15);

            double optical_depth = density * ABSORPTION * dl;
            double emit_amount = density * EMISSION_SCALE * dl;
            accumulated = Vec<double, 3>{accumulated + emit_color * (transmittance * emit_amount)};
            transmittance *= std::exp(-optical_depth);
        }

        if (r_new <= R_CAPTURE || transmittance < TRANSMITTANCE_CUTOFF) return accumulated;

        if (r_new >= R_ESCAPE && state[5] > 0.0) {
            Vec<double, 3> outdir = exit_dir_of(state);

            // Near the photon sphere, the true background is an
            // unresolvably fine, multiply-imaged Einstein-ring texture
            // (magnification formally diverges at b_crit) -- sampling
            // the discrete starfield there aliases into moire noise
            // under finite AA. Suppressed (not recolored -- no color is
            // added here, only the noisy sample's own contribution
            // reduced) in proportion to proximity to the critical curve.
            double crit_distance = r_min - R_PHOTON_SPHERE;
            double noise_suppress = std::exp(-(crit_distance * crit_distance) / (2.0 * 0.6 * 0.6));
            Vec<double, 3> sky_color{sample_sky_color(sky, outdir) * (1.0 - noise_suppress)};
            return Vec<double, 3>{accumulated + sky_color * transmittance};
        }
    }
    // Ran out of steps without resolving (deep near-critical grazer):
    // fold into captured (return accumulated with no sky term) rather
    // than faking an escape along the ORIGINAL camera-ray direction --
    // that fallback (this file's previous behavior, matching
    // blackhole_demo.cpp's own) produced a visible hard-edged, blocky
    // discontinuity exactly at the resolved/unresolved boundary
    // (reported as looking like a "censored" strip of squares in the
    // GPU port, where the same fallback existed). A ray still this
    // close to the photon sphere after max_steps is, for practical
    // purposes, indistinguishable from captured -- it's either about
    // to plunge or would need many more orbits than the step budget
    // allows either way; blending it into the shadow instead of
    // injecting an unrelated sky sample is the artifact-minimizing
    // choice, not just a fallback of convenience.
    return accumulated;
}

// Position only, parametrized by t in [0,1] -- factored out so the
// "weak focus" blend below can evaluate it at a frozen t as well as the
// current one.
Vec<double, 3> camera_position_for_t(double t) {
    double ease = smoothstep01(t);
    // Distance: approach to a close peak (the flyby), then pull back out
    // -- not a permanent close orbit. This is what produces the ~4/5-
    // frame-width dominant framing, only near the peak, the way an
    // actual flythrough would, rather than forcing that framing as the
    // whole shot's baseline.
    double D = t < 0.4 ? D_APPROACH_START + (D_PEAK - D_APPROACH_START) * smoothstep01(t / 0.4)
                        : D_PEAK + (D_DEPART_END - D_PEAK) * smoothstep01((t - 0.4) / 0.6);
    double elev = (ELEV_START_DEG + (ELEV_PEAK_DEG - ELEV_START_DEG) * ease) * std::numbers::pi / 180.0;
    double az = AZIMUTH_TOTAL_DEG * ease * std::numbers::pi / 180.0;
    return Vec<double, 3>{D * std::cos(elev) * std::cos(az), D * std::cos(elev) * std::sin(az),
                           D * std::sin(elev)};
}

// Direction of travel at t, via a centered finite difference of
// camera_position_for_t -- used only by the departure look-blend below.
Vec<double, 3> motion_direction_for_t(double t) {
    constexpr double eps = 1e-3;
    Vec<double, 3> p0 = camera_position_for_t(std::clamp(t - eps, 0.0, 1.0));
    Vec<double, 3> p1 = camera_position_for_t(std::clamp(t + eps, 0.0, 1.0));
    return Vec<double, 3>{(p1 - p0).normalized()};
}

// Where the departure look-blend (below) starts pulling the camera's
// gaze away from the hole and toward its own direction of travel.
double DEPARTURE_LOOK_START_T = 0.6;

Camera<double> camera_for_frame(int frame, int n_frames) {
    double t = n_frames > 1 ? double(frame) / double(n_frames - 1) : 0.0;
    Vec<double, 3> pos = camera_position_for_t(t);

    // Weak focus, blended continuously through the whole sequence (not a
    // hard switch from full tracking to none): mostly the direction to
    // the black hole as it was early on (frozen at FOCUS_FREEZE_T), with
    // only a gentle live correction. Because the camera doesn't actively
    // re-aim to keep the subject centered, its own arcing motion is what
    // sweeps the black hole across the frame -- small and off-center
    // early, dominant at closest approach, sliding out of frame as the
    // camera continues past. See file header's CAMERA CHOREOGRAPHY note.
    Vec<double, 3> frozen_pos = camera_position_for_t(FOCUS_FREEZE_T);
    Vec<double, 3> frozen_fwd{(-frozen_pos).normalized()};
    Vec<double, 3> live_fwd{(-pos).normalized()};
    Vec<double, 3> focus_fwd{
        Vec<double, 3>{live_fwd * FOCUS_BLEND + frozen_fwd * (1.0 - FOCUS_BLEND)}.normalized()};

    // Departure look-blend: past DEPARTURE_LOOK_START_T, the "weak
    // focus" above (still anchored toward the hole) is progressively
    // replaced by the camera's own direction of travel, reaching pure
    // motion-direction by t=1. Without this, D growing during departure
    // alone was never enough to actually lose the hole from frame --
    // the weak focus kept the gaze anchored toward it (and therefore
    // roughly centered) for the entire sequence, D or no D. A real
    // flyby stops looking back at what it passed and looks where it's
    // going; this is that, applied only late enough not to disturb the
    // approach/peak framing this session already tuned.
    double departure = smoothstep01((t - DEPARTURE_LOOK_START_T) / (1.0 - DEPARTURE_LOOK_START_T));
    Vec<double, 3> motion_fwd = motion_direction_for_t(t);
    Vec<double, 3> fwd{
        Vec<double, 3>{focus_fwd * (1.0 - departure) + motion_fwd * departure}.normalized()};
    Vec<double, 3> target{pos + fwd * pos.norm()};

    return Camera<double>{.position = pos, .target = target, .up = {0.0, 0.0, 1.0},
                           .fov_deg = FOV_DEG};
}

// "SPATIUM" watermark, spelled out as literal point stars on the
// celestial sphere rather than rendered text -- reuses Sky's own star
// primitive (add_star_to_sky mirrors make_starfield's own bucket-insert
// logic), so it costs nothing beyond a handful more entries in the
// existing star lookup. A hardcoded 5x7 dot-matrix table covers just the
// 7 distinct letters needed; there's no font dependency because there's
// no rasterization step -- each "on" cell IS one star.
struct Glyph5x7 { std::uint8_t rows[7]; };

const Glyph5x7* find_glyph(char c) {
    static const Glyph5x7 kS{{0b01111, 0b10000, 0b10000, 0b01110, 0b00001, 0b00001, 0b11110}};
    static const Glyph5x7 kP{{0b11110, 0b10001, 0b10001, 0b11110, 0b10000, 0b10000, 0b10000}};
    static const Glyph5x7 kA{{0b01110, 0b10001, 0b10001, 0b11111, 0b10001, 0b10001, 0b10001}};
    static const Glyph5x7 kT{{0b11111, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100}};
    static const Glyph5x7 kI{{0b11111, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100, 0b11111}};
    static const Glyph5x7 kU{{0b10001, 0b10001, 0b10001, 0b10001, 0b10001, 0b10001, 0b01110}};
    static const Glyph5x7 kM{{0b10001, 0b11011, 0b10101, 0b10001, 0b10001, 0b10001, 0b10001}};
    switch (c) {
        case 'S': return &kS;
        case 'P': return &kP;
        case 'A': return &kA;
        case 'T': return &kT;
        case 'I': return &kI;
        case 'U': return &kU;
        case 'M': return &kM;
        default: return nullptr;
    }
}

void add_star_to_sky(Sky& sky, const Vec<double, 3>& dir, double ang_size, double brightness) {
    sky.stars.push_back({dir, ang_size, brightness});
    int idx = static_cast<int>(sky.stars.size()) - 1;
    constexpr double pi = std::numbers::pi;
    double theta = std::acos(std::clamp(dir[1], -1.0, 1.0));
    double phi = std::atan2(dir[2], dir[0]);
    int tb = std::clamp(static_cast<int>(theta / pi * Sky::kThetaBuckets), 0, Sky::kThetaBuckets - 1);
    int pb = std::clamp(static_cast<int>((phi + pi) / (2.0 * pi) * Sky::kPhiBuckets), 0,
                         Sky::kPhiBuckets - 1);
    sky.star_buckets[static_cast<std::size_t>(tb) * Sky::kPhiBuckets + static_cast<std::size_t>(pb)]
        .push_back(idx);
}

// Places `text` as a block of point stars centered on `center_dir`, via
// small-angle tangent-plane projection (fine for the sub-radian block
// widths used here). `cell_ang` is the angular size of one font cell;
// `brightness` matches Star::brightness's own [0,1]-ish scale.
void add_text_stars(Sky& sky, const std::string& text, const Vec<double, 3>& center_dir,
                     double cell_ang, double brightness) {
    Vec<double, 3> c{center_dir.normalized()};
    Vec<double, 3> arbitrary =
        std::abs(c[1]) < 0.9 ? Vec<double, 3>{0.0, 1.0, 0.0} : Vec<double, 3>{1.0, 0.0, 0.0};
    // v must come from the UN-negated u_raw: negating u and then
    // deriving v = c.cross(u) from that already-flipped u flips v right
    // along with it, so the pair (u,v) -> (-u,-v) together, which is a
    // 180-degree ROTATION (determinant +1) -- confirmed by trying exactly
    // that first: the word moved to a new spot but still read as the same
    // "MUITA92". Only ONE axis may flip for a genuine mirror.
    //
    // The mirror itself is real, not a guess: this tangent basis spells
    // the word correctly for a viewer OUTSIDE the celestial sphere
    // looking in at its surface, but the camera sits at the sphere's
    // center looking OUT -- the same flip as text on the inside of a
    // dome. Confirmed by the original reported artifact: without this
    // fix "SPATIUM" rendered as reversed-and-mirrored "MUITA92" (M/U/I/T/A
    // are left-right symmetric so looked unchanged; P and S are not, and
    // their mirror images read as 9 and 2).
    Vec<double, 3> u_raw{arbitrary.cross(c).normalized()};
    Vec<double, 3> v{c.cross(u_raw).normalized()};
    Vec<double, 3> u{-u_raw};

    double total_w = static_cast<double>(text.size()) * 6.0 - 1.0;  // 5 cols + 1 gap/letter
    double x0 = -total_w * 0.5;
    for (std::size_t li = 0; li < text.size(); ++li) {
        const Glyph5x7* g = find_glyph(text[li]);
        if (!g) continue;
        for (int row = 0; row < 7; ++row) {
            for (int col = 0; col < 5; ++col) {
                if (!((g->rows[row] >> (4 - col)) & 1)) continue;
                double dx = (x0 + static_cast<double>(li) * 6.0 + col) * cell_ang;
                double dy = (3.0 - row) * cell_ang;
                Vec<double, 3> raw{c + u * dx + v * dy};
                add_star_to_sky(sky, raw.normalized(), cell_ang * 0.55, brightness);
            }
        }
    }
}

void print_usage() {
    std::print(
        "Usage: blackhole_gr_demo [options]\n"
        "  Rigorous-GR ray-tracer (Schwarzschild or Kerr): full 4-coordinate\n"
        "  geodesic integration (Dual<T>-exact Christoffel symbols), real\n"
        "  tilted/lensed volumetric accretion disk, camera choreography.\n"
        "  All render constants below are overridable so iteration doesn't\n"
        "  require recompiling.\n"
        "\n"
        "  Sequence:\n"
        "    --frames N            total sequence length (default {})\n"
        "    --start-frame N       render one frame at this index (spot-check\n"
        "                          without rendering everything before it)\n"
        "    --force               overwrite existing output files\n"
        "    --width N / --height N   output resolution (default {}x{})\n"
        "    --aa N                supersample grid, N*N rays/pixel (default {})\n"
        "\n"
        "  Black hole:\n"
        "    --mass M              mass in geometric units (default {})\n"
        "    --spin A              Kerr spin parameter, 0 <= |A| < mass\n"
        "                          (default 0 = Schwarzschild fast path)\n"
        "    --retrograde          disk rotates opposite the spin (default prograde)\n"
        "\n"
        "  Camera:\n"
        "    --fov DEG             field of view (default {})\n"
        "    --d-start / --d-peak / --d-end M   flyby distances, in units of\n"
        "                          mass (defaults {}/{}/{})\n"
        "    --elev-start / --elev-peak DEG     elevation above the disk plane\n"
        "                          (defaults {}/{})\n"
        "    --azimuth DEG         total azimuthal sweep (default {})\n"
        "    --focus-blend F       0=fully frozen look direction, 1=fully live\n"
        "                          tracking (default {})\n"
        "    --focus-freeze-t T    sequence position the frozen direction is\n"
        "                          taken from, in [0,1] (default {})\n"
        "\n"
        "  Disk:\n"
        "    --disk-outer M        outer disk radius, in units of mass (default {})\n"
        "    --disk-aspect R       scale-height / radius ratio (default {})\n"
        "    --disk-density-power P   radial falloff exponent (default {})\n"
        "    --disk-t0 K / --disk-t-ref K   blackbody temperature scale/reference\n"
        "                          (defaults {}/{})\n"
        "    --disk-time-per-frame T   rotation speed (default {})\n"
        "    --doppler-damping F   0=no asymmetry, 1=full physical (default {})\n"
        "    --absorption A        optical thickness per unit density (default {})\n"
        "    --emission E          emission strength per unit density (default {})\n"
        "\n"
        "  --help                  show this message\n"
        "  Output:                 blackhole_gr_frames/frame_%04d.png\n",
        N_FRAMES, W, H, AA_LEVEL, M_BH, FOV_DEG, D_APPROACH_START, D_PEAK, D_DEPART_END,
        ELEV_START_DEG, ELEV_PEAK_DEG, AZIMUTH_TOTAL_DEG, FOCUS_BLEND, FOCUS_FREEZE_T, R_DISK_OUTER,
        DISK_ASPECT_RATIO, DISK_DENSITY_POWER, DISK_T0, DISK_T_REF, DISK_TIME_PER_FRAME,
        DOPPLER_DAMPING, ABSORPTION, EMISSION_SCALE);
}

} // namespace

int main(int argc, char** argv) {
    bool force = false;
    int start_frame = -1;
    for (int i = 1; i < argc; ++i) {
        std::string_view a = argv[i];
        auto next_double = [&] { return std::atof(argv[++i]); };
        auto next_int = [&] { return std::atoi(argv[++i]); };
        if (a == "--help" || a == "-h") { print_usage(); return 0; }
        if (a == "--force") { force = true; continue; }
        if (a == "--retrograde") { PROGRADE = false; continue; }
        if (a == "--frames" && i + 1 < argc) { N_FRAMES = next_int(); continue; }
        if (a == "--start-frame" && i + 1 < argc) { start_frame = next_int(); continue; }
        if (a == "--width" && i + 1 < argc) { W = next_int(); continue; }
        if (a == "--height" && i + 1 < argc) { H = next_int(); continue; }
        if (a == "--aa" && i + 1 < argc) { AA_LEVEL = next_int(); continue; }
        if (a == "--mass" && i + 1 < argc) { M_BH = next_double(); continue; }
        if (a == "--spin" && i + 1 < argc) { SPIN = next_double(); continue; }
        if (a == "--fov" && i + 1 < argc) { FOV_DEG = next_double(); continue; }
        if (a == "--d-start" && i + 1 < argc) { D_APPROACH_START = next_double(); continue; }
        if (a == "--d-peak" && i + 1 < argc) { D_PEAK = next_double(); continue; }
        if (a == "--d-end" && i + 1 < argc) { D_DEPART_END = next_double(); continue; }
        if (a == "--elev-start" && i + 1 < argc) { ELEV_START_DEG = next_double(); continue; }
        if (a == "--elev-peak" && i + 1 < argc) { ELEV_PEAK_DEG = next_double(); continue; }
        if (a == "--azimuth" && i + 1 < argc) { AZIMUTH_TOTAL_DEG = next_double(); continue; }
        if (a == "--focus-blend" && i + 1 < argc) { FOCUS_BLEND = next_double(); continue; }
        if (a == "--focus-freeze-t" && i + 1 < argc) { FOCUS_FREEZE_T = next_double(); continue; }
        if (a == "--disk-outer" && i + 1 < argc) { R_DISK_OUTER = next_double(); continue; }
        if (a == "--disk-aspect" && i + 1 < argc) { DISK_ASPECT_RATIO = next_double(); continue; }
        if (a == "--disk-density-power" && i + 1 < argc) { DISK_DENSITY_POWER = next_double(); continue; }
        if (a == "--disk-t0" && i + 1 < argc) { DISK_T0 = next_double(); continue; }
        if (a == "--disk-t-ref" && i + 1 < argc) { DISK_T_REF = next_double(); continue; }
        if (a == "--disk-time-per-frame" && i + 1 < argc) { DISK_TIME_PER_FRAME = next_double(); continue; }
        if (a == "--doppler-damping" && i + 1 < argc) { DOPPLER_DAMPING = next_double(); continue; }
        if (a == "--absorption" && i + 1 < argc) { ABSORPTION = next_double(); continue; }
        if (a == "--emission" && i + 1 < argc) { EMISSION_SCALE = next_double(); continue; }
        std::print(stderr, "unknown option: {}\n", a);
        return 1;
    }
    recompute_derived();
    DISC_PARTICLES = make_disc_particles();

    int render_first = start_frame >= 0 ? start_frame : 0;
    int render_count = start_frame >= 0 ? 1 : N_FRAMES;

    SchwarzschildMetric<double> schw_metric{M_BH};
    KerrMetric<double> kerr_metric{M_BH, SPIN};
    Sky sky = make_starfield(10000, /*seed=*/7, /*tint=*/{8.0, 6.0, 18.0});

    // SPATIUM watermark: a plainly-legible copy off to one side (reads
    // normally in the wide establishing/departure shots), plus a second
    // copy centered on the exact antipodal sky direction behind the hole
    // AT THE PEAK FRAME. Since the camera always looks toward the hole,
    // that direction is where a b=0 ray would land, so stars placed
    // there straddle the photon sphere's capture cone: the ones nearer
    // center fall in (swallowed, texturing the shadow edge), the rest
    // get bent around it into arcs at closest approach -- real lensing,
    // not a compositing trick, exactly where the user asked for proof
    // the render isn't "merging with the library" (i.e. isn't secretly
    // a flat overlay).
    {
        // Legible copy: offset from frame 0's actual forward direction
        // (not a guessed world vector) so it lands in-frame during the
        // wide establishing shot, beside the hole's small far-away
        // silhouette rather than on top of it.
        Camera<double> cam0 = camera_for_frame(0, N_FRAMES);
        Vec<double, 3> fwd0{Vec<double, 3>{cam0.target - cam0.position}.normalized()};
        Vec<double, 3> arbitrary0 =
            std::abs(fwd0[1]) < 0.9 ? Vec<double, 3>{0.0, 1.0, 0.0} : Vec<double, 3>{1.0, 0.0, 0.0};
        Vec<double, 3> u0{arbitrary0.cross(fwd0).normalized()};
        Vec<double, 3> v0{fwd0.cross(u0).normalized()};
        Vec<double, 3> side_dir{Vec<double, 3>{fwd0 + u0 * 0.22 + v0 * 0.05}.normalized()};
        add_text_stars(sky, "SPATIUM", side_dir, 0.0008, 1.5);
    }
    {
        // Hero copy: centered on the exact antipodal sky direction
        // behind the hole AT THE PEAK FRAME -- see the comment above
        // add_text_stars' call site in the header block for why.
        Vec<double, 3> peak_pos = camera_position_for_t(0.4);
        Vec<double, 3> dir_behind_hole{(-peak_pos).normalized()};
        add_text_stars(sky, "SPATIUM", dir_behind_hole, 0.003, 2.5);
    }
    std::error_code ec;
    std::filesystem::create_directories("blackhole_gr_frames", ec);

    auto t0 = std::chrono::steady_clock::now();

    for (int frame = render_first; frame < render_first + render_count; ++frame) {
        Camera<double> cam = camera_for_frame(frame, N_FRAMES);
        const auto basis = make_camera_basis(cam);
        const double D_cam = cam.position.norm();
        // disc_density()'s p_t drives each particle's own Keplerian
        // orbital advance -- DISK_TIME_PER_FRAME is chosen so the
        // fastest (innermost) particles complete roughly half an orbit
        // over the full sequence: dphi_dt = u_avg*sqrt(0.5*u_avg) ~ 0.037
        // rad per unit p_t for a particle near R_DISK_INNER, so
        // pi/0.037 ~ 85 units of total p_t over N_FRAMES frames.
        const double disk_time = frame * DISK_TIME_PER_FRAME;

        // HDR buffer -- unclamped, may hold values well past 255. See
        // trace_ray()'s and apply_bloom()'s comments: tone_map_to_8bit()
        // below is the only place this gets clamped, once, after bloom.
        std::vector<double> img(3 * static_cast<std::size_t>(W) * H, 0.0);

        if (SPIN == 0.0) {
            // --- Schwarzschild fast path -------------------------------
            const Vec<double, 3> r_hat = cam.position / D_cam;
            const Vec<double, 3> n_disk{0.0, 0.0, 1.0};
            const double f_cam = 1.0 - 2.0 * M_BH / D_cam;

            parallel_for_rows(H, [&](int y) {
                for (int x = 0; x < W; ++x) {
                    double* px = &img[3 * (static_cast<std::size_t>(y) * W + x)];

                    auto ray_color = [&](double sx, double sy) -> Vec<double, 3> {
                        Vec<double, 3> dir_world = camera_ray_dir(basis, sx, sy);

                        // Exact tetrad-corrected initial state: a camera
                        // at coordinate radius D is a local static
                        // observer with e_r_hat=sqrt(f)*d/dr -- so a
                        // photon whose LOCAL direction has radial/
                        // tangential components (n_r,n_t) (the flat dot-
                        // product decomposition against the radial
                        // direction, exact by construction) has 4-
                        // velocity u^t=1/sqrt(f), u^r=n_r*sqrt(f),
                        // u^phi=n_t/D.
                        double dir_dot_r = dir_world.dot(r_hat);
                        Vec<double, 3> t_raw = dir_world - r_hat * dir_dot_r;
                        double t_norm = t_raw.norm();
                        if (t_norm < 1e-9) {
                            if (dir_dot_r < 0.0) return Vec<double, 3>{0.0, 0.0, 0.0};
                            return sample_sky_color(sky, dir_world);
                        }
                        Vec<double, 3> t_hat = t_raw / t_norm;

                        double ut = 1.0 / std::sqrt(f_cam);
                        double ur = dir_dot_r * std::sqrt(f_cam);
                        double uphi = t_norm / D_cam;

                        Vec<double, 3> n_ray = r_hat.cross(t_hat);
                        Vec<double, 3> line_dir = n_ray.cross(n_disk);
                        double line_norm = line_dir.norm();
                        bool disk_plane_ok = line_norm > 1e-9;
                        double phi_disk = 0.0;
                        if (disk_plane_ok) {
                            line_dir = line_dir / line_norm;
                            phi_disk = std::atan2(line_dir.dot(t_hat), line_dir.dot(r_hat));
                        }

                        Vec<double, 8> state{0.0, D_cam, std::numbers::pi / 2.0, 0.0, ut, ur, 0.0, uphi};

                        auto height_of = [&](double r_mid, double /*theta_mid*/, double phi_mid) {
                            return disk_plane_ok ? r_mid * std::sin(phi_mid - phi_disk) * line_norm : 0.0;
                        };
                        auto redshift_of = [&](double r_mid, double b) {
                            return disk_redshift_factor(M_BH, r_mid, b);
                        };
                        auto exit_dir_of = [&](const Vec<double, 8>& st) {
                            // Exact exit direction from the photon's
                            // actual outgoing 4-velocity, not the
                            // rotated radial vector (only exact as
                            // r->infinity).
                            double r_new = st[1], phi_new = st[3];
                            double f_esc = 1.0 - 2.0 * M_BH / r_new;
                            double nr = st[5] / std::sqrt(f_esc);
                            double nphi = r_new * st[7];
                            double mag = std::sqrt(nr * nr + nphi * nphi);
                            nr /= mag; nphi /= mag;
                            double cph = std::cos(phi_new), sph = std::sin(phi_new);
                            Vec<double, 3> r_hat_exit{r_hat * cph + t_hat * sph};
                            Vec<double, 3> t_hat_exit{t_hat * cph - r_hat * sph};
                            return Vec<double, 3>{r_hat_exit * nr + t_hat_exit * nphi};
                        };

                        return trace_ray(schw_metric, state, disk_time, sky, D_cam, dir_world,
                                          height_of, redshift_of, exit_dir_of);
                    };

                    Vec<double, 3> c = supersample_pixel_hdr(x, y, W, H, basis.tan_half,
                                                              static_cast<double>(W) / H, ray_color,
                                                              AA_LEVEL);
                    px[0] = c[0]; px[1] = c[1]; px[2] = c[2];
                }
            });
        } else {
            // --- Kerr path: genuine 3D (t,r,theta,phi) integration -----
            // See file header's KERR PATH note for the oblate-spheroidal
            // setup and the static-observer tetrad approximation at the
            // (always-far) camera.
            SpheroidalPoint cam_sp = to_spheroidal(cam.position, SPIN);
            SpheroidalTetrad tet_cam = spheroidal_tetrad(cam_sp.r, cam_sp.theta, cam_sp.phi, SPIN);
            Matrix<double, 4, 4> g_cam =
                kerr_metric(Vec<double, 4>{0.0, cam_sp.r, cam_sp.theta, cam_sp.phi});
            double ut_cam = 1.0 / std::sqrt(-g_cam(0, 0));
            double sqrt_grr_cam = std::sqrt(g_cam(1, 1));
            double sqrt_gthth_cam = std::sqrt(g_cam(2, 2));
            double sqrt_gphph_cam = std::sqrt(g_cam(3, 3));

            parallel_for_rows(H, [&](int y) {
                for (int x = 0; x < W; ++x) {
                    double* px = &img[3 * (static_cast<std::size_t>(y) * W + x)];

                    auto ray_color = [&](double sx, double sy) -> Vec<double, 3> {
                        Vec<double, 3> dir_world = camera_ray_dir(basis, sx, sy);

                        double n_r = dir_world.dot(tet_cam.e_r);
                        double n_theta = dir_world.dot(tet_cam.e_theta);
                        double n_phi = dir_world.dot(tet_cam.e_phi);

                        // up=(0,0,1) is always the spin axis, and
                        // e_phi's z-component is identically zero by
                        // construction -- so up.e_phi==0 always,
                        // meaning there's always exactly one screen
                        // column where n_phi (and hence uphi) is
                        // exactly 0 for every ray in it: a genuine
                        // degenerate-orbit separatrix (not a bug) that
                        // renders badly under finite step budgets/
                        // precision (found this session, reported as a
                        // "censored"-looking blocky column). Nudged on
                        // n_phi, not uphi itself: n_phi is a unit-
                        // vector dot product (always O(1)), while
                        // uphi's own scale shrinks with camera distance
                        // -- a fixed epsilon directly on uphi clamped
                        // far more than just the degenerate column at
                        // wide shots (a real regression, caught before
                        // being committed to the long render).
                        constexpr double kDegeneracyEpsilon = 1e-4;
                        if (std::abs(n_phi) < kDegeneracyEpsilon) n_phi = kDegeneracyEpsilon;

                        double ur = n_r / sqrt_grr_cam;
                        double utheta = n_theta / sqrt_gthth_cam;
                        double uphi = n_phi / sqrt_gphph_cam;

                        Vec<double, 8> state{0.0, cam_sp.r, cam_sp.theta, cam_sp.phi,
                                              ut_cam, ur, utheta, uphi};

                        // The disk sits in the spin-aligned equatorial
                        // plane (theta=pi/2) by construction -- height
                        // above it is just r*cos(theta), directly from
                        // the state, no per-ray plane trick needed since
                        // theta is now a real dynamical variable.
                        auto height_of = [](double r_mid, double theta_mid, double /*phi_mid*/) {
                            return r_mid * std::cos(theta_mid);
                        };
                        auto redshift_of = [&](double r_mid, double b) {
                            return kerr_disk_redshift_factor(M_BH, SPIN, r_mid, b, PROGRADE);
                        };
                        auto exit_dir_of = [&](const Vec<double, 8>& st) {
                            double r_e = st[1], theta_e = st[2], phi_e = st[3];
                            Matrix<double, 4, 4> g_e =
                                kerr_metric(Vec<double, 4>{0.0, r_e, theta_e, phi_e});
                            double nr_e = st[5] * std::sqrt(g_e(1, 1));
                            double nth_e = st[6] * std::sqrt(g_e(2, 2));
                            double nph_e = st[7] * std::sqrt(g_e(3, 3));
                            double mag = std::sqrt(nr_e * nr_e + nth_e * nth_e + nph_e * nph_e);
                            nr_e /= mag; nth_e /= mag; nph_e /= mag;
                            SpheroidalTetrad tet_e = spheroidal_tetrad(r_e, theta_e, phi_e, SPIN);
                            return Vec<double, 3>{tet_e.e_r * nr_e + tet_e.e_theta * nth_e +
                                                   tet_e.e_phi * nph_e};
                        };

                        return trace_ray(kerr_metric, state, disk_time, sky, D_cam, dir_world,
                                          height_of, redshift_of, exit_dir_of);
                    };

                    Vec<double, 3> c = supersample_pixel_hdr(x, y, W, H, basis.tan_half,
                                                              static_cast<double>(W) / H, ray_color,
                                                              AA_LEVEL);
                    px[0] = c[0]; px[1] = c[1]; px[2] = c[2];
                }
            });
        }

        apply_bloom(img, W, H);
        std::vector<std::uint8_t> img8 = tone_map_to_8bit(img);

        char path[64];
        std::snprintf(path, sizeof(path), "blackhole_gr_frames/frame_%04d.png", frame);
        if (spatium::examples::confirm_overwrite(path, force)) write_png_rgb(path, W, H, img8);

        std::print("\r  frame {}/{}", frame + 1, N_FRAMES);
        std::fflush(stdout);
    }

    auto t1 = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::print("\rblackhole_gr_demo: {} frames at {}x{}, {:.0f} ms ({:.0f} ms/frame)\n", render_count, W,
               H, ms, ms / render_count);
    std::print("Assemble with:\n"
               "  ffmpeg -framerate 30 -i blackhole_gr_frames/frame_%04d.png "
               "-pix_fmt yuv420p blackhole_gr.mp4\n");
    return 0;
}
