#pragma once
// Particles smaller than a pixel, drawn by projection rather than traced.
//
// Measured before this was written: through two million sub-pixel specks
// a primary ray opened 239 interior nodes of the instance tree and tested
// 0.6 instances, and 85% of rays hit nothing -- the cost of a ray through a
// haze is the number of clusters it crosses, not the number of specks it
// hits, and no amount of tuning the tree changes that order. Projecting a
// speck costs one invocation per speck, whatever the rays do.
//
// So the renderer chooses per instance per frame, by projected size: at
// or above a threshold in pixels a speck stays in the instance tree and is
// traced with everything else; below it, it is left out of the tree
// (`kLbvhCullGlsl`) and splatted here. What a splatted speck gives up is
// what a sub-pixel speck could never show anyway: shadows, reflections, its
// outline.
//
// Three passes after the trace:
//
//   kSplat*Glsl     each small speck into the 2x2 pixels around its centre,
//                   bilinearly -- a speck moving across a pixel boundary
//                   fades across it instead of jumping -- hidden where the
//                   traced depth is nearer. Coverage (projected area times
//                   opacity) and coverage-weighted colour accumulate in
//                   16.16 fixed point, because integer atomics are the ones
//                   every device has.
//   kCompositeGlsl  coverage summed as optical depth, transmittance
//                   exp(-sum): traced colour times transmittance, plus the
//                   specks' mean colour times its complement. Order-free,
//                   which is what makes the accumulation valid at all.
//
// Not part of any module: a header of strings, like `gpu_trace_glsl.hpp`.

namespace spatium::render::gpu {

// The size test both the tree's keys and the splat pass apply, so the two
// agree on which specks are whose. A speck behind the camera is never
// splatted -- the tree decides whether anything of it is visible.
inline constexpr const char* kSplatClassifyGlsl = R"GLSL(
// cam.xyz: camera position; cam.w: pixels per unit of tan at distance 1,
// i.e. (H / 2) / tan(fov / 2). fwd.xyz: forward; fwd.w: the threshold in
// pixels, zero to splat nothing.
bool splat_small(vec3 lo, vec3 hi, vec4 cam, vec4 fwd, out float radius_px) {
    vec3 c = 0.5 * (lo + hi);
    float z = dot(c - cam.xyz, fwd.xyz);
    radius_px = 0.0;
    if (fwd.w <= 0.0 || z <= 0.0) return false;
    radius_px = 0.5 * length(hi - lo) / z * cam.w;
    return radius_px < fwd.w;
}
)GLSL";

// kSplatHeadGlsl + kSplatClassifyGlsl + kSplatBodyGlsl is the kernel.
inline constexpr const char* kSplatHeadGlsl = R"GLSL(
#version 450
layout(local_size_x = 256) in;
struct Instance { vec4 r0, r1, r2, scale_quadric, color_rough, emissive_opacity; };
layout(std430, binding = 0) readonly buffer Insts { Instance insts[]; };
layout(std430, binding = 1) readonly buffer Boxes { vec4 boxes[]; };   // lo (w = live), hi
layout(std430, binding = 2) readonly buffer Depth { float depth[]; };
layout(std430, binding = 3) buffer Accum { uint accum[]; };            // coverage, r, g, b; 16.16

layout(push_constant) uniform Push {
    vec4 cam;       // position; w = (H/2)/tan(fov/2)
    vec4 fwd;       // forward; w = threshold in pixels
    vec4 right;     // w = tan(fov/2)
    vec4 up;        // w = aspect
    vec4 key;       // w = fill strength
    vec4 fill;
    uvec4 size;     // W, H, count, 0
} pc;
)GLSL";

inline constexpr const char* kSplatBodyGlsl = R"GLSL(
void add(uint pix, float a, vec3 c) {
    if (a <= 0.0) return;
    atomicAdd(accum[4u * pix + 0u], uint(a * 65536.0));
    atomicAdd(accum[4u * pix + 1u], uint(a * c.r * 65536.0));
    atomicAdd(accum[4u * pix + 2u], uint(a * c.g * 65536.0));
    atomicAdd(accum[4u * pix + 3u], uint(a * c.b * 65536.0));
}

void main() {
    uint i = gl_GlobalInvocationID.x;
    if (i >= pc.size.z) return;
    vec4 lo = boxes[2u * i], hi = boxes[2u * i + 1u];
    if (lo.w == 0.0) return;
    float r_px;
    if (!splat_small(lo.xyz, hi.xyz, pc.cam, pc.fwd, r_px)) return;

    vec3 c = 0.5 * (lo.xyz + hi.xyz);
    vec3 v = c - pc.cam.xyz;
    float z = dot(v, pc.fwd.xyz);
    float tan_half = pc.right.w, aspect = pc.up.w;
    float W = float(pc.size.x), H = float(pc.size.y);
    // The inverse of camera_pixel_dir: a direction's right and up
    // components over its forward one are the pixel's nx and ny.
    float nx = dot(v, pc.right.xyz) / z / (aspect * tan_half);
    float ny = dot(v, pc.up.xyz) / z / tan_half;
    float x = (nx + 1.0) * 0.5 * W - 0.5;
    float y = (1.0 - ny) * 0.5 * H - 0.5;
    if (x < -1.0 || y < -1.0 || x > W || y > H) return;
    float dist = length(v);

    // shade() in gpu_trace_glsl.hpp, with the flake's own axis as its
    // normal: the sheet a flake is lies across its local z.
    Instance g = insts[i];
    vec3 n = normalize(vec3(g.r0.z, g.r1.z, g.r2.z));
    vec3 dir = v / dist;
    if (dot(n, dir) > 0.0) n = -n;
    float diff = max(0.0, dot(n, pc.key.xyz));
    diff += pc.key.w * max(0.0, dot(n, pc.fill.xyz));
    diff = min(diff, 1.0);
    vec3 col = g.color_rough.xyz * (0.22 + 0.78 * diff) + g.emissive_opacity.xyz;

    // Coverage: the speck's disc over one pixel, times its opacity.
    float a = min(1.0, 3.14159265 * r_px * r_px) * g.emissive_opacity.w;

    int x0 = int(floor(x)), y0 = int(floor(y));
    float fx = x - float(x0), fy = y - float(y0);
    for (int k = 0; k < 4; ++k) {
        int px = x0 + (k & 1), py = y0 + (k >> 1);
        if (px < 0 || py < 0 || px >= int(pc.size.x) || py >= int(pc.size.y)) continue;
        uint pix = uint(py) * pc.size.x + uint(px);
        if (dist >= depth[pix]) continue;   // behind what the trace found
        float w = ((k & 1) != 0 ? fx : 1.0 - fx) * ((k >> 1) != 0 ? fy : 1.0 - fy);
        add(pix, a * w, col);
    }
}
)GLSL";

inline constexpr const char* kCompositeGlsl = R"GLSL(
#version 450
layout(local_size_x = 8, local_size_y = 8) in;
layout(std430, binding = 0) buffer Image { uint pixels[]; };
layout(std430, binding = 1) readonly buffer Accum { uint accum[]; };
layout(push_constant) uniform Push { uvec4 size; } pc;   // W, H, bgra, 0

void main() {
    uvec2 p = gl_GlobalInvocationID.xy;
    if (p.x >= pc.size.x || p.y >= pc.size.y) return;
    uint i = p.y * pc.size.x + p.x;
    float A = float(accum[4u * i]) / 65536.0;
    if (A <= 0.0) return;
    vec3 mean = vec3(float(accum[4u * i + 1u]), float(accum[4u * i + 2u]), float(accum[4u * i + 3u]))
              / float(accum[4u * i]);
    uint px = pixels[i];
    vec3 traced = vec3(float(px & 255u), float((px >> 8) & 255u), float((px >> 16) & 255u)) / 255.0;
    if (pc.size.z != 0u) mean = mean.bgr;   // the image stores B,G,R
    float T = exp(-A);
    vec3 c = clamp(traced * T + mean * (1.0 - T), 0.0, 1.0);
    uvec3 b = uvec3(c * 255.0 + 0.5);
    pixels[i] = b.r | (b.g << 8) | (b.b << 16) | (255u << 24);
}
)GLSL";

}  // namespace spatium::render::gpu
