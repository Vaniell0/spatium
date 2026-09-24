#pragma once
// A laid-out scene as the arrays a compute shader reads, and the same
// traversal written once more on the host in the shader's own precision.
//
// Two decisions carry the rest.
//
// **The trees are converted, not rebuilt.** `BVH`'s node array is already
// what a device wants -- flat, index-linked, left child next, walked with a
// stack -- so a device tree is that array in fp32 with the primitives put in
// leaf order. A second build on the device side could split differently,
// and then a disagreement between the two renderers would be a question
// about two trees rather than about arithmetic.
//
// **`trace()` below is the shader's specification, not a convenience.** It
// is written the way GLSL has to be -- fp32 throughout, no recursion, no
// library templates -- so it can be compared against the fp64 path on the
// host, where a disagreement can be inspected, and the shader can then be
// compared against *it*. Without it the only check on the device would be
// looking at the picture.
//
// Everything is fp32 because the device is: Iris Xe has no `shaderFloat64`.
// That is not a loss for the black hole, whose precision gate passed in
// fp32, but it is not automatically free here either -- see the quadric.
//
// Not part of the `spatium.render` module, for the reason
// `render/cooked_scene.hpp` is not.
#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/render/cooked_scene.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <algorithm>
#  include <array>
#  include <cmath>
#  include <cstdint>
#  include <cstring>
#  include <limits>
#  include <optional>
#  include <vector>
#endif

namespace spatium::render::gpu {

// std430 layouts, and each one is written so that the GLSL struct of the
// same name has the same offsets with no explicit padding: a `vec3`
// followed by a 4-byte scalar packs into one 16-byte slot, and a `vec4`
// is used everywhere else. Sizes are asserted rather than trusted.

// GLSL: struct Node { vec3 lo; uint first; vec3 hi; uint count; };
// Same meaning as `BVH::Node`: a leaf has count > 0 and `first` indexes the
// primitive array; an internal node has count == 0, its left child is the
// next node and `first` is the right child.
struct Node {
    float lo[3];
    std::uint32_t first;
    float hi[3];
    std::uint32_t count;
};
static_assert(sizeof(Node) == 32);

// GLSL: struct Triangle { vec4 v0, v1, v2, n0, n1, n2, color_rough, emissive; };
// Vertices and normals in world space; `.w` unused except where named.
struct Triangle {
    float v0[4], v1[4], v2[4];
    float n0[4], n1[4], n2[4];
    float color_rough[4];   // rgb, roughness
    float emissive[4];      // rgb, unused
};
static_assert(sizeof(Triangle) == 128);

// GLSL: struct Quadric { mat4 q; vec4 lo; vec4 hi; };  (q column-major)
// One per shape, shared by every instance of it -- the whole point of an
// instance.
struct Quadric {
    float q[16];      // column-major, so a GLSL mat4 reads it unchanged
    float lo[4];      // clip box
    float hi[4];
};
static_assert(sizeof(Quadric) == 96);

// GLSL: struct Instance { vec4 r0, r1, r2; vec4 scale_quadric; vec4 color_rough; vec4 emissive_opacity; };
// Rotation rows with the translation in `.w`, so a row is one dot product
// away from a world coordinate. The quadric index travels as the bits of a
// float (`floatBitsToUint` on the device) to keep the struct one layout.
struct Instance {
    float r0[4], r1[4], r2[4];   // rows of R; .w = translation
    float scale_quadric[4];      // scale, quadric index (as bits), unused, unused
    float color_rough[4];
    float emissive_opacity[4];
};
static_assert(sizeof(Instance) == 96);

struct Scene {
    std::vector<Node> tri_nodes;
    std::vector<Triangle> triangles;          // in leaf order
    std::vector<std::uint32_t> triangle_source;  // leaf order -> index into CookedScene::prims
    std::vector<Node> inst_nodes;
    std::vector<Quadric> quadrics;
    std::vector<Instance> instances;          // in leaf order
    std::vector<std::uint32_t> instance_source;  // leaf order -> index into CookedScene::looks

    std::size_t bytes() const {
        return tri_nodes.size() * sizeof(Node) + triangles.size() * sizeof(Triangle) +
               inst_nodes.size() * sizeof(Node) + quadrics.size() * sizeof(Quadric) +
               instances.size() * sizeof(Instance);
    }
};

namespace detail {

inline std::uint32_t bits_of(float f) {
    std::uint32_t u;
    std::memcpy(&u, &f, sizeof u);
    return u;
}
inline float float_of(std::uint32_t u) {
    float f;
    std::memcpy(&f, &u, sizeof f);
    return f;
}

template<typename BvhNode>
std::vector<Node> convert_nodes(const std::vector<BvhNode>& in) {
    std::vector<Node> out(in.size());
    for (std::size_t i = 0; i < in.size(); ++i) {
        for (std::size_t k = 0; k < 3; ++k) {
            // Rounded outward, so a box never shrinks past a primitive
            // it holds when it drops to fp32 -- a box that did would cull
            // a hit the fp64 tree keeps, and that disagreement would be
            // the tree's rather than the primitive's.
            out[i].lo[k] = std::nextafter(static_cast<float>(in[i].bounds.min_corner[k]),
                                          -std::numeric_limits<float>::infinity());
            out[i].hi[k] = std::nextafter(static_cast<float>(in[i].bounds.max_corner[k]),
                                          std::numeric_limits<float>::infinity());
        }
        out[i].first = in[i].first;
        out[i].count = in[i].count;
    }
    return out;
}

inline void put3(float* dst, const Vec<double, 3>& v, float w = 0.0f) {
    dst[0] = static_cast<float>(v[0]);
    dst[1] = static_cast<float>(v[1]);
    dst[2] = static_cast<float>(v[2]);
    dst[3] = w;
}

}  // namespace detail

inline Scene pack(const CookedScene<double>& laid) {
    Scene out;

    out.tri_nodes = detail::convert_nodes(laid.triangles.nodes());
    const auto& tri_order = laid.triangles.prim_indices();
    out.triangles.resize(tri_order.size());
    out.triangle_source.resize(tri_order.size());
    for (std::size_t k = 0; k < tri_order.size(); ++k) {
        const auto& p = laid.prims[tri_order[k]];
        auto& t = out.triangles[k];
        detail::put3(t.v0, Vec<double, 3>{p.tri[0]});
        detail::put3(t.v1, Vec<double, 3>{p.tri[1]});
        detail::put3(t.v2, Vec<double, 3>{p.tri[2]});
        detail::put3(t.n0, p.vertex_normals[0]);
        detail::put3(t.n1, p.vertex_normals[1]);
        detail::put3(t.n2, p.vertex_normals[2]);
        detail::put3(t.color_rough, p.color, static_cast<float>(p.roughness));
        detail::put3(t.emissive, p.emissive);
        out.triangle_source[k] = static_cast<std::uint32_t>(tri_order[k]);
    }

    out.quadrics.resize(laid.quadrics.size());
    for (std::size_t i = 0; i < laid.quadrics.size(); ++i) {
        const auto& bq = laid.quadrics[i];
        for (std::size_t c = 0; c < 4; ++c)
            for (std::size_t r = 0; r < 4; ++r)
                out.quadrics[i].q[c * 4 + r] = static_cast<float>(bq.surface.Q(r, c));
        detail::put3(out.quadrics[i].lo, Vec<double, 3>{bq.clip.min_corner});
        detail::put3(out.quadrics[i].hi, Vec<double, 3>{bq.clip.max_corner});
    }

    out.inst_nodes = detail::convert_nodes(laid.instances.nodes());
    const auto& inst_order = laid.instances.prim_indices();
    const auto& shapes = laid.instances.shapes();
    out.instances.resize(inst_order.size());
    out.instance_source.resize(inst_order.size());
    for (std::size_t k = 0; k < inst_order.size(); ++k) {
        const auto& in = shapes[inst_order[k]];
        const auto& look = laid.looks[inst_order[k]];
        auto& g = out.instances[k];
        float* rows[3] = {g.r0, g.r1, g.r2};
        for (std::size_t r = 0; r < 3; ++r) {
            for (std::size_t c = 0; c < 3; ++c) rows[r][c] = static_cast<float>(in.rotation(r, c));
            rows[r][3] = static_cast<float>(in.translation[r]);
        }
        const auto qi = static_cast<std::uint32_t>(in.shape - laid.quadrics.data());
        g.scale_quadric[0] = static_cast<float>(in.scale);
        g.scale_quadric[1] = detail::float_of(qi);
        g.scale_quadric[2] = g.scale_quadric[3] = 0.0f;
        detail::put3(g.color_rough, look.color, static_cast<float>(look.roughness));
        detail::put3(g.emissive_opacity, look.emissive, static_cast<float>(look.opacity));
        out.instance_source[k] = static_cast<std::uint32_t>(inst_order[k]);
    }
    return out;
}

// ── The shader's specification ───────────────────────────────────────
//
// Plain fp32 and small structs, so that each function below has a line
// for line GLSL counterpart. Where the fp64 library makes a numerical
// choice -- the parallel test's epsilon, the clip tolerance, the textbook
// quadratic formula -- the same choice is made here in fp32, so that a
// disagreement measures precision and nothing else.

struct V3 { float x, y, z; };
inline V3 operator+(V3 a, V3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline V3 operator-(V3 a, V3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline V3 operator*(V3 a, float s) { return {a.x * s, a.y * s, a.z * s}; }
inline float dot(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline V3 cross(V3 a, V3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline V3 v3(const float* p) { return {p[0], p[1], p[2]}; }

struct Ray32 { V3 o, d; };

enum class HitKind : std::uint8_t { None, Triangle, Instance };

struct Hit {
    HitKind kind = HitKind::None;
    std::uint32_t index = 0;   // leaf-order index into the packed array
    float t = std::numeric_limits<float>::infinity();
    float u = 0, v = 0;        // barycentrics, triangles only
    V3 normal{};               // geometric for a triangle, exact for a quadric
};

// The slab test, as `intersect_parameters` does it: a component below
// epsilon is treated as parallel rather than divided by.
inline bool slab(const Ray32& r, const Node& n, float t_max, float& t_enter) {
    const float eps = epsilon<float>();
    float lo = 0.0f, hi = std::numeric_limits<float>::max();
    const float o[3] = {r.o.x, r.o.y, r.o.z}, d[3] = {r.d.x, r.d.y, r.d.z};
    for (int i = 0; i < 3; ++i) {
        if (std::abs(d[i]) < eps) {
            if (o[i] < n.lo[i] || o[i] > n.hi[i]) return false;
        } else {
            const float inv = 1.0f / d[i];
            float t1 = (n.lo[i] - o[i]) * inv, t2 = (n.hi[i] - o[i]) * inv;
            if (t1 > t2) std::swap(t1, t2);
            lo = std::max(lo, t1);
            hi = std::min(hi, t2);
            if (lo > hi) return false;
        }
    }
    t_enter = lo;
    return lo <= t_max;
}

// Möller-Trumbore, as `ray_triangle` does it.
inline bool hit_triangle(const Ray32& r, const Triangle& tri, Hit& h) {
    const V3 a = v3(tri.v0), e1 = v3(tri.v1) - a, e2 = v3(tri.v2) - a;
    const V3 p = cross(r.d, e2);
    const float det = dot(e1, p);
    if (std::abs(det) < epsilon<float>()) return false;
    const float f = 1.0f / det;
    const V3 s = r.o - a;
    const float u = f * dot(s, p);
    if (u < 0.0f || u > 1.0f) return false;
    const V3 q = cross(s, e1);
    const float v = f * dot(r.d, q);
    if (v < 0.0f || u + v > 1.0f) return false;
    const float t = f * dot(e2, q);
    if (t < 0.0f || !(t < h.t)) return false;
    h.t = t; h.u = u; h.v = v;
    const V3 n = cross(e1, e2);
    h.normal = n * (1.0f / std::sqrt(dot(n, n)));
    return true;
}

// Into the instance's local space and a quadric there, as `ray_hit` on an
// `Instanced<BoundedQuadric>` does it. The direction is divided by the
// scale and not renormalised, which is what keeps `t` a world distance.
//
// One departure from the fp64 path, and it is measured rather than
// assumed. The roots come from the textbook formula, as in
// `solve_quadratic`, and for a small quadric far from the ray origin b^2
// and 4ac are large and nearly equal, so their difference is where fp32's
// precision goes: on the donut at t=1.5, a sprinkle's distance came back
// with a p99 relative error of 9.0e-5. So the origin is first moved along
// the ray to the point nearest the clip box's centre, the quadratic is
// solved there, and the shift is added back. Same formula, operands of
// the size of the object instead of the size of the scene: p99 2.3e-7,
// and half the instance pixels that disagreed with fp64 stop disagreeing.
// fp64 does not need it and is left alone, because it is the reference.
inline bool hit_instance(const Ray32& r, const Instance& g, const Quadric& qd, Hit& h) {
    const float s = g.scale_quadric[0];
    if (s == 0.0f) return false;
    const V3 r0 = v3(g.r0), r1 = v3(g.r1), r2 = v3(g.r2);
    const V3 tr{g.r0[3], g.r1[3], g.r2[3]};
    // Rᵀ applied to a vector is the columns of R dotted with it, and the
    // columns of R are the rows' components.
    auto apply_rt = [&](V3 w) {
        return V3{r0.x * w.x + r1.x * w.y + r2.x * w.z,
                  r0.y * w.x + r1.y * w.y + r2.y * w.z,
                  r0.z * w.x + r1.z * w.y + r2.z * w.z};
    };
    V3 o = apply_rt((r.o - tr) * (1.0f / s));
    const V3 d = apply_rt(r.d * (1.0f / s));
    const V3 centre{(qd.lo[0] + qd.hi[0]) * 0.5f, (qd.lo[1] + qd.hi[1]) * 0.5f,
                    (qd.lo[2] + qd.hi[2]) * 0.5f};
    const float t_shift = dot(centre - o, d) / dot(d, d);
    o = o + d * t_shift;

    const float* Q = qd.q;   // column-major: Q(r, c) = Q[c * 4 + r]
    auto Qm = [Q](int row, int col) { return Q[col * 4 + row]; };
    const V3 Qd{Qm(0, 0) * d.x + Qm(0, 1) * d.y + Qm(0, 2) * d.z,
                Qm(1, 0) * d.x + Qm(1, 1) * d.y + Qm(1, 2) * d.z,
                Qm(2, 0) * d.x + Qm(2, 1) * d.y + Qm(2, 2) * d.z};
    const float oh[4] = {o.x, o.y, o.z, 1.0f};
    float Qo[4];
    for (int row = 0; row < 4; ++row)
        Qo[row] = Qm(row, 0) * oh[0] + Qm(row, 1) * oh[1] + Qm(row, 2) * oh[2] + Qm(row, 3) * oh[3];
    const float a = dot(d, Qd);
    const float b = 2.0f * dot(d, V3{Qo[0], Qo[1], Qo[2]});
    const float c = oh[0] * Qo[0] + oh[1] * Qo[1] + oh[2] * Qo[2] + oh[3] * Qo[3];

    float roots[2];
    int n_roots = 0;
    if (a == 0.0f) {
        if (b != 0.0f) roots[n_roots++] = -c / b;
    } else {
        const float disc = b * b - 4.0f * a * c;
        if (!(disc >= 0.0f)) return false;
        const float sq = std::sqrt(disc);
        roots[n_roots++] = (-b + sq) / (2.0f * a);
        roots[n_roots++] = (-b - sq) / (2.0f * a);
    }

    bool found = false;
    for (int k = 0; k < n_roots; ++k) {
        const float t = roots[k] + t_shift;
        if (!(t >= 0.0f) || !(t < h.t)) continue;
        const V3 p = o + d * roots[k];
        bool inside = true;
        const float pc[3] = {p.x, p.y, p.z};
        for (int i = 0; i < 3 && inside; ++i) {
            const float e = qd.hi[i] - qd.lo[i];
            const float tol = epsilon<float>() * std::max(std::abs(e), 1.0f) * 8.0f;
            inside = !(pc[i] < qd.lo[i] - tol) && !(pc[i] > qd.hi[i] + tol);
        }
        if (!inside) continue;
        // Gradient of the quadric at p, turned back into world space.
        const float ph[4] = {p.x, p.y, p.z, 1.0f};
        V3 grad{};
        float* gp[3] = {&grad.x, &grad.y, &grad.z};
        for (int row = 0; row < 3; ++row)
            *gp[row] = Qm(row, 0) * ph[0] + Qm(row, 1) * ph[1] + Qm(row, 2) * ph[2] + Qm(row, 3) * ph[3];
        const float len = std::sqrt(dot(grad, grad));
        const V3 nl = len > epsilon<float>() ? grad * (1.0f / len) : V3{};
        h.t = t; h.u = h.v = 0.0f;
        h.normal = V3{dot(r0, nl), dot(r1, nl), dot(r2, nl)};
        found = true;
    }
    return found;
}

// Stack traversal, near child first, as `BVH::ray_cast` walks it.
template<typename Leaf>
inline bool walk(const std::vector<Node>& nodes, const Ray32& r, Hit& h, Leaf&& leaf) {
    if (nodes.empty()) return false;
    std::uint32_t stack[64];
    int sp = 0;
    stack[sp++] = 0;
    bool any = false;
    while (sp > 0) {
        const Node& n = nodes[stack[--sp]];
        float t_enter;
        if (!slab(r, n, h.t, t_enter)) continue;
        if (n.count > 0) {
            for (std::uint32_t i = 0; i < n.count; ++i)
                if (leaf(n.first + i)) any = true;
            continue;
        }
        const std::uint32_t self = static_cast<std::uint32_t>(&n - nodes.data());
        const std::uint32_t left = self + 1, right = n.first;
        float tl = 0, tr = 0;
        const bool l_ok = slab(r, nodes[left], h.t, tl);
        const bool r_ok = slab(r, nodes[right], h.t, tr);
        if (l_ok && r_ok) {
            // Push the far child first so the near one is popped next.
            if (tl > tr) { stack[sp++] = left; stack[sp++] = right; }
            else         { stack[sp++] = right; stack[sp++] = left; }
        } else if (l_ok) {
            stack[sp++] = left;
        } else if (r_ok) {
            stack[sp++] = right;
        }
    }
    return any;
}

// The nearer of the two trees' hits, as the renderer takes it.
inline Hit trace(const Scene& sc, const Ray32& r) {
    Hit h;
    walk(sc.tri_nodes, r, h, [&](std::uint32_t i) {
        if (!hit_triangle(r, sc.triangles[i], h)) return false;
        h.kind = HitKind::Triangle;
        h.index = i;
        return true;
    });
    walk(sc.inst_nodes, r, h, [&](std::uint32_t i) {
        const auto& g = sc.instances[i];
        const auto qi = detail::bits_of(g.scale_quadric[1]);
        if (!hit_instance(r, g, sc.quadrics[qi], h)) return false;
        h.kind = HitKind::Instance;
        h.index = i;
        return true;
    });
    return h;
}

}  // namespace spatium::render::gpu
