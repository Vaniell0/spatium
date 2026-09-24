#pragma once
// The compute shader for a packed scene, as a string compiled at run time.
//
// It is `render/gpu_scene.hpp`'s `trace()` and `shade()` transcribed, and
// each function here keeps the name of the one it transcribes so the two
// can be read side by side. A change to either is a change to both; the
// comparison in `examples/donut_live.cpp` is what notices when they drift.
//
// A string in a header rather than a .comp file next to the binary: a
// file found through a relative path is how `donut_demo` came to divide by
// zero when started from the wrong directory.

namespace spatium::render::gpu {

inline constexpr const char* kTraceGlsl = R"GLSL(
#version 450
layout(local_size_x = 8, local_size_y = 8) in;

struct Node     { vec3 lo; uint first; vec3 hi; uint count; };
struct Triangle { vec4 v0, v1, v2, n0, n1, n2, color_rough, emissive; };
struct Quadric  { mat4 q; vec4 lo; vec4 hi; };
struct Instance { vec4 r0, r1, r2, scale_quadric, color_rough, emissive_opacity; };

layout(std430, binding = 0) readonly buffer TriNodes  { Node tri_nodes[]; };
layout(std430, binding = 1) readonly buffer Tris      { Triangle tris[]; };
layout(std430, binding = 2) readonly buffer InstNodes { Node inst_nodes[]; };
layout(std430, binding = 3) readonly buffer Quads     { Quadric quads[]; };
layout(std430, binding = 4) readonly buffer Insts     { Instance insts[]; };
layout(std430, binding = 5) writeonly buffer Image    { uint pixels[]; };
layout(std430, binding = 6) writeonly buffer HitIds   { uint hit_ids[]; };

layout(push_constant) uniform Push {
    vec4 cam_pos;      // .w = tan(fov/2)
    vec4 fwd;          // .w = aspect
    vec4 right;
    vec4 up;
    vec4 key;          // .w = fill strength
    vec4 fill;
    vec4 background;
    uvec4 size;        // width, height, tri node count, inst node count
} pc;

const float EPS = 1.1920929e-7 * 128.0;   // epsilon<float>()
const float INF = 1.0 / 0.0;
const float FLT_MAX = 3.4028235e38;

struct Hit { uint kind; uint index; float t; float u; float v; vec3 normal; };

bool slab(vec3 o, vec3 d, Node n, float t_max, out float t_enter) {
    float lo = 0.0, hi = FLT_MAX;
    for (int i = 0; i < 3; ++i) {
        if (abs(d[i]) < EPS) {
            if (o[i] < n.lo[i] || o[i] > n.hi[i]) { t_enter = 0.0; return false; }
        } else {
            float inv = 1.0 / d[i];
            float t1 = (n.lo[i] - o[i]) * inv, t2 = (n.hi[i] - o[i]) * inv;
            if (t1 > t2) { float s = t1; t1 = t2; t2 = s; }
            lo = max(lo, t1);
            hi = min(hi, t2);
            if (lo > hi) { t_enter = 0.0; return false; }
        }
    }
    t_enter = lo;
    return lo <= t_max;
}

bool hit_triangle(vec3 o, vec3 d, Triangle tri, inout Hit h) {
    vec3 a = tri.v0.xyz, e1 = tri.v1.xyz - a, e2 = tri.v2.xyz - a;
    vec3 p = cross(d, e2);
    float det = dot(e1, p);
    if (abs(det) < EPS) return false;
    float f = 1.0 / det;
    vec3 s = o - a;
    float u = f * dot(s, p);
    if (u < 0.0 || u > 1.0) return false;
    vec3 q = cross(s, e1);
    float v = f * dot(d, q);
    if (v < 0.0 || u + v > 1.0) return false;
    float t = f * dot(e2, q);
    if (t < 0.0 || !(t < h.t)) return false;
    h.t = t; h.u = u; h.v = v;
    vec3 n = cross(e1, e2);
    h.normal = n * (1.0 / sqrt(dot(n, n)));
    return true;
}

bool hit_instance(vec3 ro, vec3 rd, Instance g, Quadric qd, inout Hit h) {
    float s = g.scale_quadric.x;
    if (s == 0.0) return false;
    vec3 r0 = g.r0.xyz, r1 = g.r1.xyz, r2 = g.r2.xyz;
    vec3 tr = vec3(g.r0.w, g.r1.w, g.r2.w);
    // The rows of R, taken as columns, are the columns of R-transpose.
    mat3 Rt = mat3(r0, r1, r2);
    vec3 o = Rt * ((ro - tr) * (1.0 / s));
    vec3 d = Rt * (rd * (1.0 / s));
    vec3 centre = (qd.lo.xyz + qd.hi.xyz) * 0.5;
    float t_shift = dot(centre - o, d) / dot(d, d);
    o = o + d * t_shift;

    mat4 Q = qd.q;                            // column-major, Q[c][r]
    vec3 Qd = mat3(Q) * d;
    vec4 oh = vec4(o, 1.0);
    vec4 Qo = Q * oh;
    float a = dot(d, Qd);
    float b = 2.0 * dot(d, Qo.xyz);
    float c = dot(oh, Qo);

    float roots[2];
    int n_roots = 0;
    if (a == 0.0) {
        if (b != 0.0) { roots[0] = -c / b; n_roots = 1; }
    } else {
        float disc = b * b - 4.0 * a * c;
        if (!(disc >= 0.0)) return false;
        float sq = sqrt(disc);
        roots[0] = (-b + sq) / (2.0 * a);
        roots[1] = (-b - sq) / (2.0 * a);
        n_roots = 2;
    }
    bool found = false;
    for (int k = 0; k < n_roots; ++k) {
        float t = roots[k] + t_shift;
        if (!(t >= 0.0) || !(t < h.t)) continue;
        vec3 p = o + d * roots[k];
        bool inside = true;
        for (int i = 0; i < 3 && inside; ++i) {
            float e = qd.hi[i] - qd.lo[i];
            float tol = EPS * max(abs(e), 1.0) * 8.0;
            inside = !(p[i] < qd.lo[i] - tol) && !(p[i] > qd.hi[i] + tol);
        }
        if (!inside) continue;
        vec3 grad = (Q * vec4(p, 1.0)).xyz;
        float len = sqrt(dot(grad, grad));
        vec3 nl = len > EPS ? grad * (1.0 / len) : vec3(0.0);
        h.t = t; h.u = 0.0; h.v = 0.0;
        h.normal = vec3(dot(r0, nl), dot(r1, nl), dot(r2, nl));
        found = true;
    }
    return found;
}

// `walk` twice rather than once over a template: GLSL has no templates,
// and a flag switching the leaf test would be a branch in the inner loop.
void walk_triangles(vec3 o, vec3 d, inout Hit h) {
    if (pc.size.z == 0u) return;
    uint stack[64];
    float enter[64];
    int sp = 0;
    stack[sp++] = 0u; enter[0] = 0.0;
    while (sp > 0) {
        --sp;
        uint self = stack[sp];
        // Tested when it was pushed; only a hit found since can rule it out.
        if (enter[sp] > h.t) continue;
        Node n = tri_nodes[self];
        if (self == 0u) { float te; if (!slab(o, d, n, h.t, te)) continue; }
        if (n.count > 0u) {
            for (uint i = 0u; i < n.count; ++i)
                if (hit_triangle(o, d, tris[n.first + i], h)) { h.kind = 1u; h.index = n.first + i; }
            continue;
        }
        uint left = self + 1u, right = n.first;
        float tl, trr;
        bool l_ok = slab(o, d, tri_nodes[left], h.t, tl);
        bool r_ok = slab(o, d, tri_nodes[right], h.t, trr);
        if (l_ok && r_ok) {
            if (tl > trr) { enter[sp] = tl; stack[sp++] = left; enter[sp] = trr; stack[sp++] = right; }
            else          { enter[sp] = trr; stack[sp++] = right; enter[sp] = tl; stack[sp++] = left; }
        } else if (l_ok) { enter[sp] = tl; stack[sp++] = left; }
        else if (r_ok)   { enter[sp] = trr; stack[sp++] = right; }
    }
}

void walk_instances(vec3 o, vec3 d, inout Hit h) {
    if (pc.size.w == 0u) return;
    uint stack[64];
    float enter[64];
    int sp = 0;
    stack[sp++] = 0u; enter[0] = 0.0;
    while (sp > 0) {
        --sp;
        uint self = stack[sp];
        // Tested when it was pushed; only a hit found since can rule it out.
        if (enter[sp] > h.t) continue;
        Node n = inst_nodes[self];
        if (self == 0u) { float te; if (!slab(o, d, n, h.t, te)) continue; }
        if (n.count > 0u) {
            for (uint i = 0u; i < n.count; ++i) {
                Instance g = insts[n.first + i];
                uint qi = floatBitsToUint(g.scale_quadric.y);
                if (hit_instance(o, d, g, quads[qi], h)) { h.kind = 2u; h.index = n.first + i; }
            }
            continue;
        }
        uint left = self + 1u, right = n.first;
        float tl, trr;
        bool l_ok = slab(o, d, inst_nodes[left], h.t, tl);
        bool r_ok = slab(o, d, inst_nodes[right], h.t, trr);
        if (l_ok && r_ok) {
            if (tl > trr) { enter[sp] = tl; stack[sp++] = left; enter[sp] = trr; stack[sp++] = right; }
            else          { enter[sp] = trr; stack[sp++] = right; enter[sp] = tl; stack[sp++] = left; }
        } else if (l_ok) { enter[sp] = tl; stack[sp++] = left; }
        else if (r_ok)   { enter[sp] = trr; stack[sp++] = right; }
    }
}

vec3 shade(vec3 d, Hit h) {
    if (h.kind == 0u) return pc.background.xyz;
    vec3 n, base, emit;
    if (h.kind == 1u) {
        Triangle t = tris[h.index];
        float w0 = 1.0 - h.u - h.v;
        n = t.n0.xyz * w0 + t.n1.xyz * h.u + t.n2.xyz * h.v;
        n = n * (1.0 / sqrt(dot(n, n)));
        base = t.color_rough.xyz;
        emit = t.emissive.xyz;
    } else {
        Instance g = insts[h.index];
        n = h.normal;
        base = g.color_rough.xyz;
        emit = g.emissive_opacity.xyz;
    }
    if (dot(n, d) > 0.0) n = -n;
    float diff = max(0.0, dot(n, pc.key.xyz));
    diff += pc.key.w * max(0.0, dot(n, pc.fill.xyz));
    diff = min(diff, 1.0);
    return base * (0.22 + 0.78 * diff) + emit;
}

void main() {
    uvec2 px = gl_GlobalInvocationID.xy;
    uint W = pc.size.x, H = pc.size.y;
    if (px.x >= W || px.y >= H) return;
    // camera_pixel_dir, at the pixel centre.
    float nx = (2.0 * (float(px.x) + 0.5) / float(W) - 1.0) * pc.fwd.w * pc.cam_pos.w;
    float ny = (1.0 - 2.0 * (float(px.y) + 0.5) / float(H)) * pc.cam_pos.w;
    vec3 d = pc.fwd.xyz + pc.right.xyz * nx + pc.up.xyz * ny;
    d = d * (1.0 / sqrt(dot(d, d)));
    vec3 o = pc.cam_pos.xyz;

    Hit h;
    h.kind = 0u; h.index = 0u; h.t = INF; h.u = 0.0; h.v = 0.0; h.normal = vec3(0.0);
    walk_triangles(o, d, h);
    walk_instances(o, d, h);

    vec3 c = clamp(shade(d, h), 0.0, 1.0);
    uvec3 b = uvec3(c * 255.0 + 0.5);
    uint i = px.y * W + px.x;
    // background.w says the target image stores B,G,R,A; writing in its
    // order here keeps the copy into a swapchain a plain copy.
    if (pc.background.w > 0.5) b = b.bgr;
    pixels[i] = b.r | (b.g << 8) | (b.b << 16) | (255u << 24);
    hit_ids[i] = (h.kind << 30) | h.index;
}
)GLSL";

}  // namespace spatium::render::gpu
