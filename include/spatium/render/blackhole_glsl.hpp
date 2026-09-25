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
    float disk[4];                       // inner, outer, aspect h/r, temperature scale
    std::uint32_t flags[4];              // bgr, max steps, disk on, sky seed
};
static_assert(sizeof(BlackholePush) == 128);

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

// The whole shader for `scene`. Bindings: 0 pixels (uint, RGBA8), 1 the
// blackbody table, 2 and 3 the field tables the hole paths may read.
inline Result<std::string> blackhole_shader(const physics::relativity::SpacetimeScene<double>& scene) {
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

    std::string src = "#version 450\nlayout(local_size_x = 8, local_size_y = 8) in;\n";
    src += R"GLSL(
layout(std430, binding = 0) writeonly buffer Image { uint pixels[]; };
layout(std430, binding = 1) readonly buffer Blackbody { vec4 bb[256]; };
layout(push_constant) uniform Push {
    vec4 e0, e1, e2, e3, cam, view, disk; uvec4 flags;
} pc;
)GLSL";
    src += mod.code;
    src += physics::relativity::emit_metric_glsl(*metric, "metric");
    src += physics::relativity::geodesic_glsl("metric");
    src += std::format(R"GLSL(
const float M_TOTAL = {0};

// Distance to the nearest hole in units of its horizon radius: < 1 is inside.
float horizon_distance(vec3 p, float t) {{
    float best = 1e30;
{1}    return best;
}}
)GLSL", io::build::detail::glsl_float(M), hole_data);
    src += R"GLSL(
// ── Noise and stars, from integer hashes only ────────────────────
uint hash_u(uvec3 v) {
    v = v * 1664525u + 1013904223u;
    v.x += v.y * v.z; v.y += v.z * v.x; v.z += v.x * v.y;
    v ^= v >> 16u;
    v.x += v.y * v.z; v.y += v.z * v.x; v.z += v.x * v.y;
    return v.x ^ v.y ^ v.z;
}
float hash_f(uvec3 v) { return float(hash_u(v) & 0xffffffu) / 16777216.0; }
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
    float h = pc.disk.z * R;
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
    // T ~ r^(-3/4), so what a blackbody puts out, ~ T^4, falls as r^(-3):
    // the inner edge outshines the rim by orders of magnitude.
    float T = pc.disk.w * 7000.0 * pow(inner / R, 0.75) * shift;
    float I = rho * pow(inner / R, 3.0) * pow(shift, 4.0) * 3.0;
    return vec4(bb_color(T) * I, rho * 3.0);
}

void main() {
    uvec2 id = gl_GlobalInvocationID.xy;
    uint W = uint(pc.view.x), H = uint(pc.view.y);
    if (id.x >= W || id.y >= H) return;
    float aspect = pc.view.x / pc.view.y;
    vec2 ndc = (vec2(id) + 0.5) / pc.view.xy * 2.0 - 1.0;
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
    for (uint s = 0u; s < max_steps && trans > 0.01; ++s) {
        float hd = horizon_distance(x.yzw, x.x);
        if (hd < 1.02) { captured = true; break; }
        if (length(x.yzw) > escape) { escaped = true; break; }
        // The step grows with the distance to the nearest horizon, in its
        // radii: fine near a hole, where the path bends, coarse far out.
        float dl = clamp(0.05 * hd, 0.01, 1.5);
        if (pc.flags.z != 0u) {
            float g[10]; vec4 dg[10];
            metric(x, g, dg);
            vec4 e = disk_at(x, k, g);
            float ds = dl * length(k.yzw);
            light += trans * e.rgb * ds * 0.6;
            trans *= exp(-e.a * ds);
        }
        geodesic_step(x, k, dl);
    }
    // A ray that neither escaped nor fell in within the step budget is
    // circling close to the photon orbit; it is drawn as captured.
    if (escaped) light += trans * sky(normalize(k.yzw), sigma, pc.flags.w);

    // ACES-style filmic curve, then sRGB.
    vec3 c = light * pc.view.w;
    c = clamp((c * (2.51 * c + 0.03)) / (c * (2.43 * c + 0.59) + 0.14), 0.0, 1.0);
    c = pow(c, vec3(1.0 / 2.2));
    uvec3 b = uvec3(c * 255.0 + 0.5);
    if (pc.flags.x != 0u) b = b.bgr;
    // Alpha marks a ray that ended at a horizon (0) from one that did not
    // (255), so a check can measure the shadow rather than guess it from
    // colour; a window's copy ignores alpha.
    uint alpha = escaped ? 255u : 0u;
    pixels[id.y * W + id.x] = b.r | (b.g << 8) | (b.b << 16) | (alpha << 24);
}
)GLSL";
    return src;
}

} // namespace spatium::render
