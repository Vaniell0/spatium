#pragma once
// A Scatter's instances moved on the device: one kernel invocation per
// site, writing the `gpu::Instance` slot `cook()` + `lay_out()` + `pack()`
// would have written on the host.
//
// This is what the cook() cost was blocking. Re-cooking two million
// particles is seconds a frame on the host; a particle's motion, once it
// is an expression, is a few hundred arithmetic ops, and the device has
// the lanes to run two million of them. Everything that does not change
// with time -- the sites, their frames, the shape, the material's
// constants -- is uploaded once, and time arrives as one push constant.
//
// The kernel is generated per node from the node's own fields
// (`io/field_glsl.hpp`), so nothing in it is specific to the dust: any
// Scatter whose motion is a structural placement and whose colours are
// structural can be moved this way.
//
// Not part of the `spatium.render` module, for the reason
// `render/cooked_scene.hpp` is not.
#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/io/build.hpp>
#  include <spatium/io/field_glsl.hpp>
#  include <spatium/render/gpu_scene.hpp>
#  include <cstdint>
#  include <format>
#  include <string>
#  include <vector>
#endif

namespace spatium::render::gpu {

// A site as the kernel reads it: the seat, and the site frame's columns.
// GLSL: struct Site { vec4 pos; vec4 c0, c1, c2; };
struct Site {
    float pos[4];
    float c0[4], c1[4], c2[4];
};
static_assert(sizeof(Site) == 64);

// The push constants. GLSL: see `InstancePush` in the generated source.
struct InstancePush {
    float rest_centroid_t[4];     // the shape's rest centroid; .w = time
    float color_rough[4];         // material constants when no colour field
    float emissive_opacity[4];
    std::uint32_t info[4];        // site count, quadric index, has colour field, has glow field
    std::uint32_t base[4];        // where in the instance array this node's slots start
};
static_assert(sizeof(InstancePush) == 80);

struct InstanceKernel {
    std::string source;           // a complete compute shader
    io::build::GlslModule module; // its noise and gather tables, bindings 0 and 1
    std::vector<Site> sites;      // binding 2; the output is binding 3
    InstancePush push{};          // everything but the time
    std::size_t node = 0;
};

// The kernel for Scatter node `idx`, whose instances share `shape`'s exact
// form (`quadric` is its index in the packed quadric array). Fails, naming
// the reason, when the node cannot be moved this way: a motion that
// deforms, or any field with a closure in it.
inline Result<InstanceKernel> make_instance_kernel(const io::build::Trace<double>& trace,
                                                   std::size_t idx,
                                                   const Vec<double, 3>& rest_centroid,
                                                   std::uint32_t quadric) {
    namespace bd = io::build;
    const auto& n = trace.node(idx);
    if (n.kind != bd::Kind::Scatter)
        return std::unexpected(Error(ErrorCode::InvalidArgument,
                                     std::format("node {} is not a Scatter", idx)));

    InstanceKernel k;
    k.node = idx;
    k.module.code = bd::glsl_prelude(0, 1);
    if (auto r = bd::emit_placement(k.module, n.transform, "motion"); !r)
        return std::unexpected(Error(r.error().code, std::format("node {}: {}", idx, r.error().message)));
    const bool has_colour = static_cast<bool>(n.color_fn);
    const bool has_glow = static_cast<bool>(n.emissive_fn);
    if (has_colour)
        if (auto r = bd::emit_vector(k.module, n.color_fn, "colour"); !r)
            return std::unexpected(Error(r.error().code, std::format("node {} colour: {}", idx, r.error().message)));
    if (has_glow)
        if (auto r = bd::emit_vector(k.module, n.emissive_fn, "glow"); !r)
            return std::unexpected(Error(r.error().code, std::format("node {} glow: {}", idx, r.error().message)));

    for (const auto& s : bd::scatter_spots(trace, idx, 0.0)) {
        Site site{};
        detail::put3(site.pos, s.position);
        detail::put3(site.c0, Vec<double, 3>{s.frame(0, 0), s.frame(1, 0), s.frame(2, 0)});
        detail::put3(site.c1, Vec<double, 3>{s.frame(0, 1), s.frame(1, 1), s.frame(2, 1)});
        detail::put3(site.c2, Vec<double, 3>{s.frame(0, 2), s.frame(1, 2), s.frame(2, 2)});
        k.sites.push_back(site);
    }

    const auto& mat = n.material;
    detail::put3(k.push.rest_centroid_t, rest_centroid);
    detail::put3(k.push.color_rough, mat.base_color, static_cast<float>(mat.roughness));
    detail::put3(k.push.emissive_opacity, mat.emissive, static_cast<float>(mat.opacity));
    k.push.info[0] = static_cast<std::uint32_t>(k.sites.size());
    k.push.info[1] = quadric;
    k.push.info[2] = has_colour ? 1u : 0u;
    k.push.info[3] = has_glow ? 1u : 0u;

    // cook()'s composition, written out: the placement `p -> R s p + tr`
    // applied to a site seated at `pos` and turned by `F` gives the
    // instance translation `tr + R (s pos)` and rotation `R F`, and the
    // material is resolved where the shape's centroid lands.
    k.source = std::string("#version 450\nlayout(local_size_x = 64) in;\n") + k.module.code +
               (has_colour ? "" : "vec3 colour(FieldIn in_) { return vec3(0.0); }\n") +
               (has_glow ? "" : "vec3 glow(FieldIn in_) { return vec3(0.0); }\n") + R"GLSL(
struct Site { vec4 pos; vec4 c0; vec4 c1; vec4 c2; };
struct Instance { vec4 r0, r1, r2, scale_quadric, color_rough, emissive_opacity; };
layout(std430, binding = 2) readonly buffer Sites { Site sites[]; };
layout(std430, binding = 3) writeonly buffer Out { Instance insts[]; };
layout(push_constant) uniform InstancePush {
    vec4 rest_centroid_t;
    vec4 color_rough;
    vec4 emissive_opacity;
    uvec4 info;
    uvec4 base;
} pc;

void main() {
    uint i = gl_GlobalInvocationID.x;
    if (i >= pc.info.x) return;
    Site st = sites[i];
    FieldIn fin;
    fin.u = pc.rest_centroid_t.w;
    fin.v = 0.0;
    fin.id = i;
    fin.origin = st.pos.xyz;
    fin.p = vec3(0.0);

    vec3 tr; mat3 R; float s;
    motion_place(fin, tr, R, s);
    vec3 T = tr + R * (st.pos.xyz * s);
    mat3 W = R * mat3(st.c0.xyz, st.c1.xyz, st.c2.xyz);

    FieldIn cin = fin;
    cin.p = W * (pc.rest_centroid_t.xyz * s) + T;
    vec3 col = pc.info.z != 0u ? colour(cin) : pc.color_rough.xyz;
    vec3 emit = pc.info.w != 0u ? glow(cin) : pc.emissive_opacity.xyz;

    Instance o;
    // Rows of W, translation in .w -- gpu::Instance's layout. GLSL
    // indexes a matrix by column first.
    o.r0 = vec4(W[0][0], W[1][0], W[2][0], T.x);
    o.r1 = vec4(W[0][1], W[1][1], W[2][1], T.y);
    o.r2 = vec4(W[0][2], W[1][2], W[2][2], T.z);
    o.scale_quadric = vec4(s, uintBitsToFloat(pc.info.y), 0.0, 0.0);
    o.color_rough = vec4(col, pc.color_rough.w);
    o.emissive_opacity = vec4(emit, pc.emissive_opacity.w);
    insts[pc.base.x + i] = o;
}
)GLSL";
    return k;
}

}  // namespace spatium::render::gpu
