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
    const bool has_colour = static_cast<bool>(n.color_fn);
    const bool has_glow = static_cast<bool>(n.emissive_fn);

    // The whole of main() in one value-numbered scope: the placement, the
    // composition with the site, then colour and glow at the point the
    // shape's centroid lands on. A hash, a gathered target or a pull the
    // motion already computed is not computed again for the colour.
    bd::GlslScope sc;
    const bd::GlslInputs at_rest{"fin", "fin.p"};
    auto pl = bd::emit_placement_into(sc, k.module, n.transform, at_rest);
    if (!pl)
        return std::unexpected(Error(pl.error().code, std::format("node {}: {}", idx, pl.error().message)));
    sc.line(std::format("vec3 T = {} + {} * (st.pos.xyz * {});", pl->translation, pl->rotation, pl->scale));
    sc.line(std::format("mat3 W = {} * mat3(st.c0.xyz, st.c1.xyz, st.c2.xyz);", pl->rotation));
    sc.line(std::format("float S = {};", pl->scale));
    sc.line("vec3 at = W * (pc.rest_centroid_t.xyz * S) + T;");
    const bd::GlslInputs at_point{"fin", "at"};
    std::string colour = "pc.color_rough.xyz", glow = "pc.emissive_opacity.xyz";
    auto vector_at = [&](const bd::VecField<double>& f, const char* what) -> Result<std::string> {
        auto pod = bd::lower(f);
        if (!pod)
            return std::unexpected(Error(pod.error().code, std::format("node {} {}: {}", idx, what, pod.error().message)));
        const auto b = bd::detail::append_data(k.module, pod->scalars.tables, pod->scalars.points);
        return bd::emit_vector_into(sc, *pod, b, at_point, "at");
    };
    if (has_colour) {
        auto r = vector_at(n.color_fn, "colour");
        if (!r) return std::unexpected(std::move(r.error()));
        colour = *r;
    }
    if (has_glow) {
        auto r = vector_at(n.emissive_fn, "glow");
        if (!r) return std::unexpected(std::move(r.error()));
        glow = *r;
    }
    sc.line(std::format("vec3 col = {};", colour));
    sc.line(std::format("vec3 emit = {};", glow));

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
    // The noise tables into workgroup memory when they are small enough to
    // sit there comfortably -- see field_perm_at in the prelude.
    const std::size_t perm = k.module.perm.size();
    const std::string shared_perm =
        perm > 0 && perm <= 4096 ? std::format("#define FIELD_PERM_SHARED {}\n", perm) : "";
    k.source = std::string("#version 450\n") + shared_perm + "layout(local_size_x = 64) in;\n" +
               k.module.code + R"GLSL(
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
    // Before the early return: the copy ends in a barrier every invocation
    // of the workgroup must reach.
    field_load_perm();
    uint i = gl_GlobalInvocationID.x;
    if (i >= pc.info.x) return;
    Site st = sites[i];
    FieldIn fin;
    fin.u = pc.rest_centroid_t.w;
    fin.v = 0.0;
    fin.id = i;
    fin.origin = st.pos.xyz;
    fin.p = vec3(0.0);
)GLSL" + sc.body + R"GLSL(
    Instance o;
    // Rows of W, translation in .w -- gpu::Instance's layout. GLSL
    // indexes a matrix by column first.
    o.r0 = vec4(W[0][0], W[1][0], W[2][0], T.x);
    o.r1 = vec4(W[0][1], W[1][1], W[2][1], T.y);
    o.r2 = vec4(W[0][2], W[1][2], W[2][2], T.z);
    o.scale_quadric = vec4(S, uintBitsToFloat(pc.info.y), 0.0, 0.0);
    o.color_rough = vec4(col, pc.color_rough.w);
    o.emissive_opacity = vec4(emit, pc.emissive_opacity.w);
    insts[pc.base.x + i] = o;
}
)GLSL";
    return k;
}

}  // namespace spatium::render::gpu
