#pragma once
// The compute shader that renders a SpacetimeScene: one light ray per
// pixel, traced backwards from an observer through the scene's metric
// (physics/relativity/metric_glsl.hpp), collecting the disk's emission on
// the way and the sky where it escapes.
//
// What each part is:
//
//   Camera     an observer at rest in the scene's coordinates, with a
//              tetrad orthonormal in the metric there (built on the host,
//              `observer_tetrad`); a pixel's ray is k = -e0 + n, null by
//              construction, integrated towards the past.
//   Step       the affine step scales with the distance to the nearest
//              hole, from a floor near the horizon to a cap far out.
//   Capture    closer to a hole than its horizon radius is black.
//   Disk       emission and absorption in a thick disk about the centre of
//              mass, temperature ~ r^(-3/4), matter on Keplerian circles
//              about the total mass. The shift g = (p . u_obs)/(p . u_emit)
//              brightens it by g^4 and moves its temperature by g -- the
//              Doppler and gravitational shift together, from the metric.
//              Colour from render/spectral.hpp's blackbody_to_rgb255 as a
//              table, so device and host share one curve.
//   Sky        stars as point sources: a jittered star per cell of a cube
//              map, brightness by a power law, colour by temperature, drawn
//              as a Gaussian a pixel wide so a star is a point rather than
//              a disc; a galactic band with dust lanes from value noise.
//
// Not part of any module: a header of strings, like render/gpu_trace_glsl.hpp.

#include <spatium/io/field_glsl.hpp>
#include <spatium/physics/relativity/metric_glsl.hpp>
#include <spatium/physics/relativity/spacetime_scene.hpp>
#include <spatium/render/spectral.hpp>
#include <cmath>
#include <format>
#include <string>
#include <vector>

namespace spatium::render {

// Push constants, 128 bytes: the most every device guarantees.
struct BlackholePush {
    float e0[4], e1[4], e2[4], e3[4];   // the observer's tetrad: time, right, up, forward
    float cam[4];                        // the observer's event
    float view[4];                       // width, height, tan(fov/2), exposure
    float disk[4];                       // inner, outer, Doppler damping, temperature scale
    std::uint32_t flags[4];              // bits (see kFlag*), max steps, sample index, sky seed
};
static_assert(sizeof(BlackholePush) == 128);
// flags[0]: the target stores B,G,R; the disk is on; accumulate samples
// (flags[2] is then which sample this is, 0 starting afresh).
inline constexpr std::uint32_t kFlagBgr = 1u, kFlagDisk = 2u, kFlagAccumulate = 4u;

// 256 entries of blackbody_to_rgb255 over 1000..40000 K, log-spaced, in 0..1.
inline std::vector<float> blackbody_table() {
    std::vector<float> t(4 * 256);
    for (int i = 0; i < 256; ++i) {
        const double k = 1000.0 * std::pow(40.0, i / 255.0);
        const auto c = blackbody_to_rgb255(k);
        for (int j = 0; j < 3; ++j) t[4 * i + j] = static_cast<float>(c[static_cast<std::size_t>(j)] / 255.0);
        t[4 * i + 3] = 1.0f;
    }
    return t;
}

// The observer at rest at `x`, looking at the origin with `up` roughly up:
// e0 along d/dt normalised in the metric, then forward, up and right by
// Gram-Schmidt in the metric, each orthogonal to e0 and to the ones before.
template<typename Metric>
std::array<Vec<double, 4>, 4> observer_tetrad(const Metric& metric, const Vec<double, 4>& x,
                                              const Vec<double, 3>& up_hint) {
    const auto g = metric(x);
    auto dot = [&](const Vec<double, 4>& a, const Vec<double, 4>& b) {
        double s = 0;
        for (std::size_t i = 0; i < 4; ++i)
            for (std::size_t j = 0; j < 4; ++j) s += g(i, j) * a[i] * b[j];
        return s;
    };
    Vec<double, 4> e0{1.0, 0.0, 0.0, 0.0};
    e0 = Vec<double, 4>{e0 * (1.0 / std::sqrt(-dot(e0, e0)))};
    auto orth = [&](Vec<double, 4> v, std::initializer_list<const Vec<double, 4>*> against, bool timelike_first) {
        for (auto* a : against) {
            const double s = dot(*a, *a);
            v = Vec<double, 4>{v - *a * (dot(v, *a) / s)};
        }
        (void)timelike_first;
        return Vec<double, 4>{v * (1.0 / std::sqrt(dot(v, v)))};
    };
    const Vec<double, 3> f3{Vec<double, 3>{-x[1], -x[2], -x[3]} * (1.0 / std::sqrt(x[1] * x[1] + x[2] * x[2] + x[3] * x[3]))};
    const Vec<double, 4> fwd = orth({0.0, f3[0], f3[1], f3[2]}, {&e0}, false);
    const Vec<double, 4> up = orth({0.0, up_hint[0], up_hint[1], up_hint[2]}, {&e0, &fwd}, false);
    // right = fwd x up in the coordinates, then made orthonormal too.
    const Vec<double, 3> r3{f3[1] * up_hint[2] - f3[2] * up_hint[1], f3[2] * up_hint[0] - f3[0] * up_hint[2],
                            f3[0] * up_hint[1] - f3[1] * up_hint[0]};
    const Vec<double, 4> right = orth({0.0, r3[0], r3[1], r3[2]}, {&e0, &fwd, &up}, false);
    return {e0, right, up, fwd};
}

// What every shader of a scene shares: the hole paths, the metric and its
// geodesic step, the horizon test, integer hashes. `form` is the Kerr-
// Schild form traced in; rays and dust both use the outgoing one, so their
// coordinates agree.
inline Result<std::string> scene_common_glsl(const physics::relativity::SpacetimeScene<double>& scene) {
    // Rays are traced backwards, towards the past horizon: the outgoing
    // form is the one regular there (see SpacetimeScene::metric).
    auto metric = scene.metric(physics::relativity::SpacetimeScene<double>::Form::outgoing);
    if (!metric) return std::unexpected(metric.error());

    io::build::GlslModule mod;
    mod.code += io::build::glsl_prelude(2, 3);
    const auto& holes = scene.holes();
    std::string hole_data;
    for (std::size_t i = 0; i < holes.size(); ++i) {
        for (const char* c : {"x", "y", "z"}) {
            const auto& f = std::string(c) == "x" ? holes[i].x : std::string(c) == "y" ? holes[i].y : holes[i].z;
            auto r = io::build::emit_scalar(mod, f, std::format("hole{}_{}", i, c));
            if (!r) return std::unexpected(r.error());
        }
        const double m = holes[i].mass, a = holes[i].spin;
        // The horizon is r = r+ in the Boyer-Lindquist radius, which in
        // Kerr-Schild coordinates is not a sphere: on the equator it sits
        // at |x| = sqrt(r+^2 + a^2). So the test takes r from |x| and z,
        // as the metric does, rather than the coordinate distance.
        hole_data += std::format(
            "    {{ FieldIn fi; fi.u = t; vec3 q = p - vec3(hole{0}_x(fi), hole{0}_y(fi), hole{0}_z(fi));"
            " float w = dot(q, q) - {2}; float r = sqrt(0.5 * w + sqrt(0.25 * w * w + {2} * q.z * q.z));"
            " float d = r / {1}; if (d < best) best = d; }}\n",
            i, io::build::detail::glsl_float(m + std::sqrt(std::max(0.0, m * m - a * a))),
            io::build::detail::glsl_float(a * a));
    }
    const double M = scene.total_mass();
    std::string src = mod.code;
    src += physics::relativity::emit_metric_glsl(*metric, "metric");
    src += physics::relativity::geodesic_glsl("metric");
    src += std::format(R"GLSL(
const float M_TOTAL = {0};

// Distance to the nearest hole in units of its horizon radius: < 1 is inside.
float horizon_distance(vec3 p, float t) {{
    float best = 1e30;
{1}    return best;
}}

uint hash_u(uvec3 v) {{
    v = v * 1664525u + 1013904223u;
    v.x += v.y * v.z; v.y += v.z * v.x; v.z += v.x * v.y;
    v ^= v >> 16u;
    v.x += v.y * v.z; v.y += v.z * v.x; v.z += v.x * v.y;
    return v.x ^ v.y ^ v.z;
}}
float hash_f(uvec3 v) {{ return float(hash_u(v) & 0xffffffu) / 16777216.0; }}
)GLSL", io::build::detail::glsl_float(M), hole_data);
    return src;
}

// The whole shader for `scene`. Bindings: 0 pixels (uint, RGBA8), 1 the
// blackbody table, 2 and 3 the field tables the hole paths may read, 4 the
// unused, 5 the accumulation, 6 and 7 the matter's tree and instances.
inline Result<std::string> blackhole_shader(const physics::relativity::SpacetimeScene<double>& scene) {
    auto common = scene_common_glsl(scene);
    if (!common) return std::unexpected(common.error());
    const auto& dust = scene.dust();

    std::string src = "#version 450\nlayout(local_size_x = 8, local_size_y = 8) in;\n";
    src += R"GLSL(
layout(std430, binding = 0) writeonly buffer Image { uint pixels[]; };
layout(std430, binding = 1) readonly buffer Blackbody { vec4 bb[256]; };
layout(std430, binding = 5) buffer Accum { vec4 accum[]; };   // rgb sum, samples
struct LNode    { vec3 lo; uint left; vec3 hi; uint right; };
struct Instance { vec4 r0, r1, r2, scale_quadric, color_rough, emissive_opacity; };
layout(std430, binding = 6) readonly buffer Nodes { LNode nodes[]; };   // the matter's tree
layout(std430, binding = 7) readonly buffer Insts { Instance insts[]; };
layout(push_constant) uniform Push {
    vec4 e0, e1, e2, e3, cam, view, disk; uvec4 flags;
} pc;
)GLSL";
    src += *common;
    src += std::format(R"GLSL(
const bool MATTER_ON = {0};
const uint SKY_SEED = {2}u;
const float MATTER_GAIN = {1};
)GLSL", dust.count > 0 ? "true" : "false", io::build::detail::glsl_float(120000.0 / std::max<double>(1.0, dust.count)), scene.sky().seed);
    src += R"GLSL(
// ── Noise and stars, from integer hashes only ────────────────────
float value_noise(vec3 p) {
    vec3 i = floor(p), f = fract(p);
    vec3 w = f * f * (3.0 - 2.0 * f);
    float n = 0.0;
    for (int dz = 0; dz < 2; ++dz) for (int dy = 0; dy < 2; ++dy) for (int dx = 0; dx < 2; ++dx) {
        vec3 o = vec3(dx, dy, dz);
        float h = hash_f(uvec3(ivec3(i + o) + 65536));
        vec3 k = mix(1.0 - w, w, o);
        n += h * k.x * k.y * k.z;
    }
    return n;
}
float fbm(vec3 p) {
    float s = 0.0, a = 0.5;
    for (int o = 0; o < 5; ++o) { s += a * value_noise(p); p *= 2.07; a *= 0.5; }
    return s;
}
vec3 bb_color(float kelvin) {
    float u = clamp(log(kelvin / 1000.0) / log(40.0), 0.0, 1.0) * 255.0;
    int i = int(u);
    return mix(bb[i].rgb, bb[min(i + 1, 255)].rgb, u - float(i));
}

// A cube-map cell per direction; each cell may hold one star, jittered
// inside it. A star is drawn as a Gaussian of `sigma` radians, so a star
// is a point a pixel wide, not a disc.
vec3 stars(vec3 d, float sigma, uint seed) {
    vec3 a = abs(d);
    int face; vec2 uv;
    if (a.x >= a.y && a.x >= a.z) { face = d.x > 0.0 ? 0 : 1; uv = d.yz / a.x; }
    else if (a.y >= a.z)          { face = d.y > 0.0 ? 2 : 3; uv = d.xz / a.y; }
    else                          { face = d.z > 0.0 ? 4 : 5; uv = d.xy / a.z; }
    const float N = 360.0;
    vec2 c = (uv * 0.5 + 0.5) * N;
    vec3 sum = vec3(0.0);
    for (int dy = -1; dy <= 1; ++dy) for (int dx = -1; dx <= 1; ++dx) {
        ivec2 cell = ivec2(floor(c)) + ivec2(dx, dy);
        uvec3 key = uvec3(uint(cell.x + 4096), uint(cell.y + 4096), uint(face) + seed * 8u);
        if (hash_f(key) > 0.02) continue;                         // about 15 000 stars in all
        vec2 j = vec2(hash_f(key + 11u), hash_f(key + 23u));
        vec2 suv = ((vec2(cell) + j) / N) * 2.0 - 1.0;
        vec3 sd;
        if (face == 0) sd = vec3(1.0, suv); else if (face == 1) sd = vec3(-1.0, suv);
        else if (face == 2) sd = vec3(suv.x, 1.0, suv.y); else if (face == 3) sd = vec3(suv.x, -1.0, suv.y);
        else if (face == 4) sd = vec3(suv, 1.0); else sd = vec3(suv, -1.0);
        sd = normalize(sd);
        // The chord, not acos of the dot: in float, acos near 1 steps in
        // units of about a pixel's angle, which drew every star as rings.
        float ang = length(sd - d);
        // brightness: power law, most stars faint
        float b = pow(hash_f(key + 37u), 12.0) * 40.0 + 0.15;
        float kelvin = mix(3000.0, 12000.0, pow(hash_f(key + 53u), 1.5));
        sum += bb_color(kelvin) * b * exp(-0.5 * ang * ang / (sigma * sigma));
    }
    return sum;
}

vec3 sky(vec3 d, float sigma, uint seed) {
    // The galactic plane, tilted against the scene's.
    vec3 n = normalize(vec3(0.3, 0.85, 0.42));
    float b = dot(d, n);
    float band = exp(-b * b / 0.018);
    float body = fbm(d * 3.0 + 7.0);
    float lanes = smoothstep(0.35, 0.65, fbm(d * 9.0 + 3.0));
    vec3 glow = mix(vec3(0.55, 0.45, 0.35), vec3(0.35, 0.42, 0.6), fbm(d * 1.5)) * band * body * (1.0 - 0.8 * lanes * band);
    vec3 s = stars(d, sigma, seed) * (1.0 - 0.6 * lanes * band);
    return 0.18 * glow + s;
}

// ── The disk ─────────────────────────────────────────────────────
// Emission and absorption at event (t, p) for a photon with covector
// p_mu = g_mu_nu k^nu; returns (emission rgb, absorption).
vec4 disk_at(vec4 x, vec4 k, float g[10]) {
    vec3 p = x.yzw;
    float R = length(p.xy);
    float inner = pc.disk.x, outer = pc.disk.y;
    if (R < inner * 0.8 || R > outer * 1.3) return vec4(0.0);
    float h = 0.018 * R;
    float vert = exp(-0.5 * p.z * p.z / (h * h));
    if (vert < 1e-3) return vec4(0.0);
    float edge = smoothstep(inner * 0.8, inner, R) * (1.0 - smoothstep(outer, outer * 1.3, R));
    // Keplerian matter about the total mass; the pattern turns with it.
    float omega = sqrt(M_TOTAL / (R * R * R));
    float phi = atan(p.y, p.x) - omega * x.x;
    float swirl = fbm(vec3(cos(phi) * 3.0, sin(phi) * 3.0, log(R) * 6.0));
    // Streaks: the noise sharpened, so the disk has lanes rather than haze.
    float lanes = smoothstep(0.35, 0.75, swirl);
    float rho = vert * edge * (0.15 + 1.6 * lanes);
    // Emitter 4-velocity, and the shift to the observer.
    vec4 v = vec4(1.0, -omega * p.y, omega * p.x, 0.0);
    mat4 gm;
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) gm[i][j] = g[metric_entry(i, j)];
    float vv = dot(v, gm * v);
    if (vv >= 0.0) return vec4(0.0);          // no timelike circular orbit here
    vec4 u = v / sqrt(-vv);
    vec4 pl = gm * k;                           // lowered
    float shift = dot(pl, pc.e0) / dot(pl, u);  // nu_obs / nu_emit
    // Damped towards 1 before it reaches the picture, as Double Negative
    // did for Gargantua: the exact asymmetry, tens to one, reads as a
    // one-sided disk. pc.disk.z is how much of it is kept, 1 the physics.
    shift = 1.0 + pc.disk.z * (shift - 1.0);
    // T ~ r^(-3/4), so what a blackbody puts out, ~ T^4, falls as r^(-3):
    // the inner edge outshines the rim by orders of magnitude.
    float T = pc.disk.w * 7000.0 * pow(inner / R, 0.75) * shift;
    float I = rho * pow(inner / R, 3.0) * pow(shift, 4.0) * 3.0;
    // Absorption kept low: the thin ring at the shadow's edge is the disk
    // seen after one, two, three windings about the photon orbit, and a
    // ray that has lost its transmittance by then never shows it.
    return vec4(bb_color(T) * I, rho * 0.8);
}

// ── The matter, drawn through its tree ───────────────────────────
// The light the particles near the chord a -> b give off: each a Gaussian
// blob of radius sigma about its centre, its temperature from where it is
// (T ~ r^(-3/4)) and its shift from its own 4-velocity, so the Doppler is
// the particle's, not a formula's. The chord is one RK4 step, short near
// the matter, so a straight segment stands for the curved path there.
bool seg_box(vec3 a, vec3 inv, float len, vec3 lo, vec3 hi) {
    vec3 t0 = (lo - a) * inv, t1 = (hi - a) * inv;
    vec3 tn = min(t0, t1), tf = max(t0, t1);
    float n = max(max(tn.x, tn.y), max(tn.z, 0.0));
    float f = min(min(tf.x, tf.y), min(tf.z, len));
    return n <= f;
}
vec3 matter_along(vec3 a, vec3 b, vec4 k, float g[10], uint node_count) {
    vec3 light = vec3(0.0);
    if (!MATTER_ON || node_count == 0u) return light;
    vec3 d = b - a;
    float len = length(d);
    if (len <= 0.0) return light;
    vec3 dir = d / len;
    vec3 inv = 1.0 / (dir + vec3(1e-12));
    if (!seg_box(a, inv, len, nodes[0].lo, nodes[0].hi)) return light;
    mat4 gm;
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) gm[i][j] = g[metric_entry(i, j)];
    vec4 pl = gm * k;
    float inner = pc.disk.x;
    uint stack[48];
    int sp = 0;
    stack[sp++] = 0u;
    while (sp > 0) {
        LNode n = nodes[stack[--sp]];
        if ((n.left & 0x80000000u) != 0u) {
            Instance q = insts[n.left & 0x7fffffffu];
            vec3 c = vec3(q.r0.w, q.r1.w, q.r2.w);
            float sigma = q.scale_quadric.x;
            float t = clamp(dot(c - a, dir), 0.0, len);
            vec3 off = a + dir * t - c;
            float w = exp(-0.5 * dot(off, off) / (sigma * sigma));
            if (w < 1e-3) continue;
            float R = max(length(c.xy), inner * 0.5);
            float shift = dot(pl, pc.e0) / dot(pl, q.color_rough);
            shift = 1.0 + pc.disk.z * (shift - 1.0);
            float T = pc.disk.w * 7000.0 * pow(inner / R, 0.75) * shift;
            float I = pow(inner / R, 3.0) * pow(shift, 4.0);
            light += bb_color(T) * I * w * min(len, 2.5 * sigma) / sigma * MATTER_GAIN;
            continue;
        }
        LNode l = nodes[n.left], r = nodes[n.right];
        if (sp < 46 && seg_box(a, inv, len, l.lo, l.hi)) stack[sp++] = n.left;
        if (sp < 46 && seg_box(a, inv, len, r.lo, r.hi)) stack[sp++] = n.right;
    }
    return light;
}

void main() {
    uvec2 id = gl_GlobalInvocationID.xy;
    uint W = uint(pc.view.x), H = uint(pc.view.y);
    if (id.x >= W || id.y >= H) return;
    float aspect = pc.view.x / pc.view.y;
    // Accumulating, each sample lands somewhere else in the pixel.
    vec2 jitter = vec2(0.5);
    if ((pc.flags.x & 4u) != 0u && pc.flags.z > 0u)
        jitter = vec2(hash_f(uvec3(id, pc.flags.z)), hash_f(uvec3(id, pc.flags.z + 7919u)));
    vec2 ndc = (vec2(id) + jitter) / pc.view.xy * 2.0 - 1.0;
    ndc.y = -ndc.y;
    vec3 n = normalize(vec3(ndc.x * pc.view.z * aspect, ndc.y * pc.view.z, 1.0));
    vec4 x = pc.cam;
    vec4 k = -pc.e0 + n.x * pc.e1 + n.y * pc.e2 + n.z * pc.e3;
    float sigma = 0.7 * 2.0 * pc.view.z / pc.view.y;   // a pixel, in radians

    vec3 light = vec3(0.0);
    float trans = 1.0;
    uint max_steps = pc.flags.y;
    float escape = 3.0 * length(pc.cam.yzw);
    bool captured = false, escaped = false;
    for (uint s = 0u; s < max_steps && trans > 0.001; ++s) {
        float hd = horizon_distance(x.yzw, x.x);
        if (hd < 1.02) { captured = true; break; }
        if (length(x.yzw) > escape) { escaped = true; break; }
        // The step grows with the distance to the nearest horizon, in its
        // radii: fine near a hole, where the path bends, coarse far out.
        float dl = clamp(0.05 * hd, 0.01, 1.5);
        // Finer near the photon orbit, where a ray winds and the rings
        // of higher order are exponentially thin.
        if (hd < 2.5) dl *= 0.35;
        float ds = dl * length(k.yzw);
        vec4 x0 = x, k0 = k;
        geodesic_step(x, k, dl);
        if ((pc.flags.x & 2u) != 0u) {
            float g[10]; vec4 dg[10];
            metric(x0, g, dg);
            if (MATTER_ON) {
                // The disk is its particles: one substance, not a field
                // with dust laid over it.
                light += trans * matter_along(x0.yzw, x.yzw, k0, g, pc.flags.w);
            } else {
                vec4 e = disk_at(x0, k0, g);
                light += trans * e.rgb * ds * 0.6;
                trans *= exp(-e.a * ds);
            }
        }
    }
    // A ray that neither escaped nor fell in within the step budget is
    // circling close to the photon orbit; it is drawn as captured.
    if (escaped) light += trans * sky(normalize(k.yzw), sigma, SKY_SEED);

    // Accumulate in HDR, before the curve, so the average is of light.
    uint pix = id.y * W + id.x;
    if ((pc.flags.x & 4u) != 0u) {
        vec4 sum = pc.flags.z == 0u ? vec4(0.0) : accum[pix];
        sum += vec4(light, 1.0);
        accum[pix] = sum;
        light = sum.rgb / sum.a;
    }
    // A curve on luminance that keeps the colour: the filmic curve on each
    // channel pressed every bright pixel to the same white, and the
    // brighter, bluer side of the disk with it.
    vec3 c = light * pc.view.w;
    float L = dot(c, vec3(0.2126, 0.7152, 0.0722));
    float Lm = L * (1.0 + L / 64.0) / (1.0 + L);
    c *= Lm / max(L, 1e-6);
    c = clamp(c, 0.0, 1.0);
    c = pow(c, vec3(1.0 / 2.2));
    uvec3 b = uvec3(c * 255.0 + 0.5);
    if ((pc.flags.x & 1u) != 0u) b = b.bgr;
    // Alpha marks a ray that ended at a horizon (0) from one that did not
    // (255), so a check can measure the shadow rather than guess it from
    // colour; a window's copy ignores alpha.
    uint alpha = escaped ? 255u : 0u;
    pixels[pix] = b.r | (b.g << 8) | (b.b << 16) | (alpha << 24);
}
)GLSL";
    return src;
}

// ── Dust: test particles on the scene's geodesics ────────────────
//
// Each particle is (x, u), a timelike geodesic of the same metric the rays
// are traced through -- in the same outgoing form, so a particle and a ray
// agree on where it is. That form is singular on the future horizon,
// which a falling particle approaches; a particle within 1.5 horizon radii
// of a hole is born again at the dust's outer edge on a circular orbit, and
// the ones lost that way are behind the shadow anyway. No gas: the dust
// feels gravity and nothing else.

// Advances every particle until its coordinate time reaches `target`, in
// at most 32 RK4 steps. Binding 0 the particles (x, u per particle), 2 and
// 3 the field tables. Push: target, count, frame, outer radius.
inline Result<std::string> dust_move_shader(const physics::relativity::SpacetimeScene<double>& scene) {
    auto common = scene_common_glsl(scene);
    if (!common) return std::unexpected(common.error());
    std::string src = "#version 450\nlayout(local_size_x = 64) in;\n";
    src += R"GLSL(
layout(std430, binding = 0) buffer Particles { vec4 s[]; };
layout(push_constant) uniform Push { float target; uint count; uint frame; float outer; } pc;
)GLSL";
    src += *common;
    src += R"GLSL(
// A circular orbit about the total mass at radius R, angle phi, height z:
// the Keplerian direction, normalised in the metric.
void born(uint i, out vec4 x, out vec4 u) {
    float a = hash_f(uvec3(i, pc.frame, 7u)), b = hash_f(uvec3(i, pc.frame, 13u)), c = hash_f(uvec3(i, pc.frame, 29u));
    float R = pc.outer * sqrt(mix(0.5, 1.0, a));
    float phi = 6.2831853 * b;
    float z = (c - 0.5) * 0.04 * R;
    x = vec4(pc.target, R * cos(phi), R * sin(phi), z);
    float omega = sqrt(M_TOTAL / (R * R * R));
    vec4 v = vec4(1.0, -omega * x.z, omega * x.y, 0.0);
    float g[10]; vec4 dg[10];
    metric(x, g, dg);
    mat4 gm;
    for (int p = 0; p < 4; ++p) for (int q = 0; q < 4; ++q) gm[p][q] = g[metric_entry(p, q)];
    u = v / sqrt(max(1e-6, -dot(v, gm * v)));
}

void main() {
    uint i = gl_GlobalInvocationID.x;
    if (i >= pc.count) return;
    vec4 x = s[2u * i], u = s[2u * i + 1u];
    for (int k = 0; k < 32 && x.x < pc.target; ++k) {
        float hd = horizon_distance(x.yzw, x.x);
        if (hd < 1.5 || length(x.yzw) > 3.0 * pc.outer) { born(i, x, u); break; }
        // The proper-time step, capped so the coordinate time lands on target.
        float dl = min(clamp(0.05 * hd, 0.01, 1.5) * 4.0, (pc.target - x.x) / max(u.x, 1e-3));
        geodesic_step(x, u, dl);
    }
    s[2u * i] = x;
    s[2u * i + 1u] = u;
}
)GLSL";
    return src;
}

// Each particle as an instance of the unit sphere scaled to sigma, for
// viewer/gpu_lbvh.hpp's tree: translation its position, color_rough its
// 4-velocity (what the shader reads the Doppler from). Binding 0 the
// particles, 1 the instances. Push: count, sigma.
inline const char* kMatterPackGlsl = R"GLSL(
#version 450
layout(local_size_x = 64) in;
struct Instance { vec4 r0, r1, r2, scale_quadric, color_rough, emissive_opacity; };
layout(std430, binding = 0) readonly buffer Particles { vec4 s[]; };
layout(std430, binding = 1) writeonly buffer Insts { Instance insts[]; };
layout(push_constant) uniform Push { uint count; float sigma; uint pad0; uint pad1; } pc;
void main() {
    uint i = gl_GlobalInvocationID.x;
    if (i >= pc.count) return;
    vec4 x = s[2u * i], u = s[2u * i + 1u];
    Instance q;
    q.r0 = vec4(1.0, 0.0, 0.0, x.y);
    q.r1 = vec4(0.0, 1.0, 0.0, x.z);
    q.r2 = vec4(0.0, 0.0, 1.0, x.w);
    q.scale_quadric = vec4(pc.sigma, uintBitsToFloat(0u), 0.0, 0.0);
    q.color_rough = u;
    q.emissive_opacity = vec4(x.x, 0.0, 0.0, 1.0);
    insts[i] = q;
}
)GLSL";

// The particles' first state: clumps on circular orbits about the total
// mass, each already sheared as Keplerian rotation would have sheared it
// over its age -- a particle at radius R has turned by Omega(R) * age from
// where the clump began -- so the disk starts as the streams differential
// rotation makes, and the geodesics carry on from there. Normalised in the
// metric; no gas and no magnetic field, only the shear.
inline std::vector<float> dust_initial(const physics::relativity::SpacetimeScene<double>& scene) {
    const auto& d = scene.dust();
    const auto metric = *scene.metric(physics::relativity::SpacetimeScene<double>::Form::outgoing);
    const double M = scene.total_mass(), outer = d.outer * M, inner = d.inner * M;
    std::vector<float> s(8 * static_cast<std::size_t>(d.count));
    std::uint64_t state = d.seed * 0x9E3779B97F4A7C15ull + 1;
    auto uni = [&state] {
        state = state * 6364136223846793005ull + 1442695040888963407ull;
        return static_cast<double>(state >> 11) / static_cast<double>(1ull << 53);
    };
    auto gauss = [&] {
        const double a = std::max(1e-12, uni()), b = uni();
        return std::sqrt(-2.0 * std::log(a)) * std::cos(2 * std::numbers::pi * b);
    };
    const std::uint32_t clumps = std::max(1u, d.count / 400u);
    for (std::uint32_t i = 0; i < d.count; ++i) {
        // Which clump: its centre radius (area-weighted), phase and age.
        std::uint64_t cs = (i % clumps) * 0x9E3779B97F4A7C15ull + d.seed;
        auto cuni = [&cs] {
            cs = cs * 6364136223846793005ull + 1442695040888963407ull;
            return static_cast<double>(cs >> 11) / static_cast<double>(1ull << 53);
        };
        const double Rc = std::sqrt(inner * inner + (outer * outer - inner * inner) * cuni());
        const double phic = 2 * std::numbers::pi * cuni();
        const double age = 300.0 * M * cuni();
        const double width = 0.06 * Rc;
        const double R = std::max(inner * 0.9, Rc + width * gauss());
        const double omega = std::sqrt(M / (R * R * R)), omegac = std::sqrt(M / (Rc * Rc * Rc));
        const double phi = phic + (omega - omegac) * age + 0.02 * gauss();
        const double z = 0.01 * R * gauss();
        const Vec<double, 4> x{0.0, R * std::cos(phi), R * std::sin(phi), z};
        const Vec<double, 4> v{1.0, -omega * x[2], omega * x[1], 0.0};
        const auto g = metric(x);
        double vv = 0;
        for (std::size_t a = 0; a < 4; ++a)
            for (std::size_t b2 = 0; b2 < 4; ++b2) vv += g(a, b2) * v[a] * v[b2];
        const double n = 1.0 / std::sqrt(std::max(1e-12, -vv));
        for (int c = 0; c < 4; ++c) {
            s[8 * i + c] = static_cast<float>(x[static_cast<std::size_t>(c)]);
            s[8 * i + 4 + c] = static_cast<float>(v[static_cast<std::size_t>(c)] * n);
        }
    }
    return s;
}

} // namespace spatium::render
