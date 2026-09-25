#pragma once
// The linear BVH's build, as compute shaders -- `render/lbvh.hpp`'s
// `build_lbvh()` step for step, less the sort, which stays on the host.
//
// Four kernels, each one of build_lbvh()'s phases:
//
//   kBoxesGlsl   an instance's world box, and the bounds of the centres
//   kKeysGlsl    the Morton code of each centre, with the index beside it --
//                after render/gpu_splat_glsl.hpp's kSplatClassifyGlsl, since
//                an instance drawn by projection is left out
//   (host)       the sort, over 8 bytes a key read from shared memory
//   kTreeGlsl    a leaf and an internal node per invocation, Karras 2012
//   kBoundsGlsl  boxes climbing from the leaves on per-node counters
//
// The sort is the one phase that is not a map over independent items -- a
// parallel radix sort on a device is its own subsystem -- and on memory
// shared with the host it is also the cheapest phase to leave where it is:
// the host reads 8 bytes per instance instead of the 96 of a slot.
//
// Not part of any module: a header of strings, like `gpu_trace_glsl.hpp`.

namespace spatium::render::gpu {

inline constexpr const char* kLbvhCommonGlsl = R"GLSL(
#version 450
layout(local_size_x = 256) in;
struct Quadric  { mat4 q; vec4 lo; vec4 hi; };
struct Instance { vec4 r0, r1, r2, scale_quadric, color_rough, emissive_opacity; };
struct LNode    { vec3 lo; uint left; vec3 hi; uint right; };
const uint LEAF = 0x80000000u;
const uint NONE = 0xffffffffu;

// A float as a uint whose unsigned order is the float's order, so an
// atomic min or max on the uint is one on the float.
uint ordered(float f) {
    uint u = floatBitsToUint(f);
    return (u & 0x80000000u) != 0u ? ~u : (u | 0x80000000u);
}
float unordered(uint u) {
    return uintBitsToFloat((u & 0x80000000u) != 0u ? (u & 0x7fffffffu) : ~u);
}
)GLSL";

// detail::instance_box and the centre bounds of build_lbvh().
inline constexpr const char* kBoxesGlsl = R"GLSL(
layout(std430, binding = 0) readonly buffer Insts { Instance insts[]; };
layout(std430, binding = 1) readonly buffer Quads { Quadric quads[]; };
layout(std430, binding = 2) writeonly buffer Boxes { vec4 boxes[]; };   // lo (w = live), hi
layout(std430, binding = 3) coherent buffer Bounds { uint bounds[]; };  // min xyz, max xyz, ordered
layout(push_constant) uniform Push { uvec4 info; } pc;                   // count

void main() {
    uint i = gl_GlobalInvocationID.x;
    if (i >= pc.info.x) return;
    Instance g = insts[i];
    float sc = g.scale_quadric.x;
    if (sc == 0.0) {
        boxes[2u * i] = vec4(0.0, 0.0, 0.0, 0.0);
        boxes[2u * i + 1u] = vec4(0.0);
        return;
    }
    Quadric q = quads[floatBitsToUint(g.scale_quadric.y)];
    float s = abs(sc);
    vec3 c = 0.5 * (q.lo.xyz + q.hi.xyz) * sc;
    vec3 hl = 0.5 * (q.hi.xyz - q.lo.xyz) * s;
    vec4 rows[3] = vec4[3](g.r0, g.r1, g.r2);
    vec3 lo, hi;
    for (int r = 0; r < 3; ++r) {
        float wc = rows[r].w + dot(rows[r].xyz, c);
        float wh = dot(abs(rows[r].xyz), hl);
        float pad = wh * 1e-6 + abs(wc) * 1e-7;
        lo[r] = wc - wh - pad;
        hi[r] = wc + wh + pad;
    }
    boxes[2u * i] = vec4(lo, 1.0);
    boxes[2u * i + 1u] = vec4(hi, 0.0);
    vec3 ctr = 0.5 * (lo + hi);
    for (int k = 0; k < 3; ++k) {
        atomicMin(bounds[k], ordered(ctr[k]));
        atomicMax(bounds[3 + k], ordered(ctr[k]));
    }
}
)GLSL";

// build_lbvh()'s keys: the 30-bit Morton code of the centre above the
// instance index -- as a uvec2 (index, code), which is the little-endian
// layout of the host's `(code << 32) | index`. Dead and splatted instances
// get no key at all; the live ones are packed at the front.
inline constexpr const char* kKeysGlsl = R"GLSL(
layout(std430, binding = 0) readonly buffer Boxes { vec4 boxes[]; };
layout(std430, binding = 1) coherent buffer Bounds { uint bounds[]; };   // [6]: live key count
layout(std430, binding = 2) writeonly buffer Keys { uvec2 keys[]; };
// cam and fwd are render/gpu_splat_glsl.hpp's: an instance splat_small()
// accepts is drawn by projection and left out of the tree. A zero
// threshold (fwd.w) keeps every instance.
layout(push_constant) uniform Push { uvec4 info; vec4 cam; vec4 fwd; } pc;

uint spread10(uint v) {
    v &= 0x3ffu;
    v = (v | (v << 16)) & 0x030000FFu;
    v = (v | (v << 8)) & 0x0300F00Fu;
    v = (v | (v << 4)) & 0x030C30C3u;
    v = (v | (v << 2)) & 0x09249249u;
    return v;
}

void main() {
    uint i = gl_GlobalInvocationID.x;
    if (i >= pc.info.x) return;
    vec4 lo = boxes[2u * i], hi = boxes[2u * i + 1u];
    // Only live keys are written, packed at the front by a counter: the
    // host sorts what the tree holds, not every instance. Their order here
    // is whatever the atomics give -- the sort is what orders them.
    float r_px;
    if (lo.w == 0.0 || splat_small(lo.xyz, hi.xyz, pc.cam, pc.fwd, r_px)) return;
    uint code = 0u;
    for (int k = 0; k < 3; ++k) {
        float a = unordered(bounds[k]), z = unordered(bounds[3 + k]);
        float ext = z > a ? z - a : 1.0;
        float c = 0.5 * (lo[k] + hi[k]);
        float u = clamp((c - a) / ext, 0.0, 1.0);
        code |= spread10(uint(u * 1023.0)) << (2 - k);
    }
    keys[atomicAdd(bounds[6], 1u)] = uvec2(i, code);
}
)GLSL";

// build_lbvh()'s leaves and internal nodes. `info.x` is the count of live
// keys after the sort, n; nodes [0, n-1) are internal and [n-1, 2n-1) the
// leaves, as on the host.
inline constexpr const char* kTreeGlsl = R"GLSL(
layout(std430, binding = 0) readonly buffer Keys { uvec2 keys[]; };
layout(std430, binding = 1) readonly buffer Boxes { vec4 boxes[]; };
layout(std430, binding = 2) writeonly buffer Nodes { LNode nodes[]; };
layout(std430, binding = 3) writeonly buffer Parent { uint parent[]; };
layout(push_constant) uniform Push { uvec4 info; } pc;

// Leading zeros of the 64-bit xor of two keys; -1 past either end.
int delta(int i, int j) {
    int n = int(pc.info.x);
    if (j < 0 || j >= n) return -1;
    uvec2 a = keys[i], b = keys[j];
    uint xh = a.y ^ b.y, xl = a.x ^ b.x;
    if (xh != 0u) return 31 - findMSB(xh);
    return 32 + 31 - findMSB(xl);
}

void main() {
    int i = int(gl_GlobalInvocationID.x);
    int n = int(pc.info.x);
    if (i >= n) return;

    uint inst = keys[i].x;
    LNode leaf;
    leaf.lo = boxes[2u * inst].xyz;
    leaf.hi = boxes[2u * inst + 1u].xyz;
    leaf.left = inst | LEAF;
    leaf.right = 0u;
    nodes[n - 1 + i] = leaf;
    if (i >= n - 1) return;

    int d = delta(i, i + 1) - delta(i, i - 1) > 0 ? 1 : -1;
    int dmin = delta(i, i - d);
    int lmax = 2;
    while (delta(i, i + lmax * d) > dmin) lmax *= 2;
    int l = 0;
    for (int t = lmax / 2; t >= 1; t /= 2)
        if (delta(i, i + (l + t) * d) > dmin) l += t;
    int j = i + l * d;
    int dnode = delta(i, j);
    int s = 0;
    for (int div = 2;; div *= 2) {
        int t = (l + div - 1) / div;
        if (delta(i, i + (s + t) * d) > dnode) s += t;
        if (t <= 1) break;
    }
    int g = i + s * d + min(d, 0);
    uint left = min(i, j) == g ? uint(n - 1 + g) : uint(g);
    uint right = max(i, j) == g + 1 ? uint(n - 1 + g + 1) : uint(g + 1);
    LNode node;
    node.lo = vec3(0.0); node.hi = vec3(0.0);
    node.left = left;
    node.right = right;
    nodes[i] = node;
    parent[left] = uint(i);
    parent[right] = uint(i);
}
)GLSL";

// build_lbvh()'s bounds: each leaf climbs; at each parent the first child
// to arrive stops and the second, whose sibling's box is by then written,
// takes the union and climbs on. The barrier before each counter makes
// the box written below it visible to whichever invocation arrives second.
inline constexpr const char* kBoundsGlsl = R"GLSL(
layout(std430, binding = 0) coherent buffer Nodes { LNode nodes[]; };
layout(std430, binding = 1) readonly buffer Parent { uint parent[]; };
layout(std430, binding = 2) coherent buffer Arrived { uint arrived[]; };
layout(push_constant) uniform Push { uvec4 info; } pc;

void main() {
    uint i = gl_GlobalInvocationID.x;
    uint n = pc.info.x;
    if (i >= n) return;
    uint at = parent[n - 1u + i];
    memoryBarrierBuffer();
    while (at != NONE) {
        if (atomicAdd(arrived[at], 1u) == 0u) return;
        memoryBarrierBuffer();
        LNode a = nodes[nodes[at].left];
        LNode b = nodes[nodes[at].right];
        nodes[at].lo = min(a.lo, b.lo);
        nodes[at].hi = max(a.hi, b.hi);
        memoryBarrierBuffer();
        at = parent[at];
    }
}
)GLSL";

}  // namespace spatium::render::gpu
