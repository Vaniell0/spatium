#pragma once
// A lowered field written out as GLSL: straight-line code, one statement
// per op, for a driver's compiler to allocate registers for and fold.
//
// Emitted rather than interpreted. An interpreter on a GPU keeps a register
// array the size of the whole program per lane -- far past what fits in
// registers for the dust's hundred-odd ops -- and branches on every op. The
// same program as code has neither cost, and one more thing comes free:
// the three components of a `Make` repeat the hash, the target and the
// pull, and a compiler merges repeated subexpressions where an interpreter
// could only evaluate them three times.
//
// **What it must agree with.** The POD interpreter, run in fp32 on the
// host, is the reference; the device is compared against it within a
// tolerance rather than to the bit, because a GPU compiler may reassociate
// and contract arithmetic that the host keeps in order. The definitions
// that carry real content -- Perlin noise, the instance hash, the SO(3)
// exponential, the gather's clamp -- are transcribed from the C++ ones
// line for line below, and each names what it transcribes.
//
// Not part of any module: it sits beside `field_pod.hpp`, which sits on
// `io/build.hpp`.
#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/io/field_pod.hpp>
#  include <cmath>
#  include <cstdint>
#  include <format>
#  include <string>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::io::build {

// Everything a set of emitted fields reads, shared by all of them: the
// noise tables (one uint per permutation byte, so a lookup is one load)
// and the gather tables, x y z per point. The kernel binds these two
// buffers; each emitted program knows its own offsets into them.
struct GlslModule {
    std::string code;
    std::vector<std::uint32_t> perm;
    std::vector<float> points;
};

// The helpers every emitted program calls, and the input record. Binding
// numbers are the caller's: `perm_binding` and `points_binding` name the
// two storage buffers holding `GlslModule::perm` and `::points`.
inline std::string glsl_prelude(int perm_binding, int points_binding) {
    return std::format(R"GLSL(
layout(std430, binding = {0}) readonly buffer FieldPerm   {{ uint field_perm[]; }};
layout(std430, binding = {1}) readonly buffer FieldPoints {{ float field_points[]; }};

// FieldInputs: time as the first parameter, then the instance.
struct FieldIn {{ float u; float v; uint id; vec3 origin; vec3 p; }};

// instance_hash / instance_unit in io/field.hpp.
uint field_instance_hash(uint id, uint salt) {{
    uint h = 2166136261u ^ salt;
    h ^= id;
    h *= 16777619u;
    h ^= h >> 13;
    h *= 0x85ebca6bu;
    h ^= h >> 16;
    return h;
}}
float field_instance_unit(uint id, uint salt) {{
    return float(field_instance_hash(id, salt) % 1000000u) / 1000000.0;
}}

// PerlinNoise::sample in algebra/noise.hpp, over a table at `base`.
int field_perm_at(uint base, int i) {{ return int(field_perm[base + uint(i & 511)]); }}
float field_fade(float t) {{ return t * t * t * (t * (t * 6.0 - 15.0) + 10.0); }}
float field_lerp(float t, float a, float b) {{ return a + t * (b - a); }}
float field_grad(int hash, float x, float y, float z) {{
    int h = hash & 15;
    float u = h < 8 ? x : y;
    float v = h < 4 ? y : ((h == 12 || h == 14) ? x : z);
    return ((h & 1) == 0 ? u : -u) + ((h & 2) == 0 ? v : -v);
}}
float field_noise(uint base, float x, float y, float z) {{
    int X = int(floor(x)) & 255;
    int Y = int(floor(y)) & 255;
    int Z = int(floor(z)) & 255;
    x -= floor(x);
    y -= floor(y);
    z -= floor(z);
    float u = field_fade(x), v = field_fade(y), w = field_fade(z);
    int A = field_perm_at(base, X) + Y, AA = field_perm_at(base, A) + Z, AB = field_perm_at(base, A + 1) + Z;
    int B = field_perm_at(base, X + 1) + Y, BA = field_perm_at(base, B) + Z, BB = field_perm_at(base, B + 1) + Z;
    return field_lerp(w,
        field_lerp(v, field_lerp(u, field_grad(field_perm_at(base, AA), x, y, z),
                                    field_grad(field_perm_at(base, BA), x - 1.0, y, z)),
                      field_lerp(u, field_grad(field_perm_at(base, AB), x, y - 1.0, z),
                                    field_grad(field_perm_at(base, BB), x - 1.0, y - 1.0, z))),
        field_lerp(v, field_lerp(u, field_grad(field_perm_at(base, AA + 1), x, y, z - 1.0),
                                    field_grad(field_perm_at(base, BA + 1), x - 1.0, y, z - 1.0)),
                      field_lerp(u, field_grad(field_perm_at(base, AB + 1), x, y - 1.0, z - 1.0),
                                    field_grad(field_perm_at(base, BB + 1), x - 1.0, y - 1.0, z - 1.0))));
}}

// Field::gather_at in io/field.hpp: truncate a non-negative index, clamp.
float field_gather(uint start, uint extent, float index, uint k) {{
    if (extent == 0u) return 0.0;
    float lo = index < 0.0 ? 0.0 : index;
    uint at = uint(lo);
    if (at >= extent) at = extent - 1u;
    return field_points[(start + at) * 3u + k];
}}

// SO3::exp in algebra/groups/so3.hpp: Rodrigues from theta^2, with the
// Taylor branch below epsilon<float>().
mat3 field_so3_exp(vec3 w) {{
    float theta2 = dot(w, w);
    float angle = sqrt(theta2);
    // Columns of skew(w): K = [[0,-z,y],[z,0,-x],[-y,x,0]].
    mat3 K = mat3(vec3(0.0, w.z, -w.y), vec3(-w.z, 0.0, w.x), vec3(w.y, -w.x, 0.0));
    float a, b;
    if (angle < 1.1920929e-7 * 128.0) {{
        a = 1.0 - theta2 / 6.0;
        b = 0.5 - theta2 / 24.0;
    }} else {{
        a = sin(angle) / angle;
        b = (1.0 - cos(angle)) / theta2;
    }}
    return mat3(1.0) + K * a + (K * K) * b;
}}
)GLSL",
                       perm_binding, points_binding);
}

namespace detail {

// A float literal GLSL reads back as the same float: nine significant
// digits round-trip any fp32, and a bare integer needs its point.
inline std::string glsl_float(double v) {
    const float f = static_cast<float>(v);
    if (std::isinf(f)) return f > 0 ? "(1.0/0.0)" : "(-1.0/0.0)";
    if (std::isnan(f)) return "(0.0/0.0)";
    std::string s = std::format("{:.9g}", f);
    if (s.find_first_of(".eE") == std::string::npos) s += ".0";
    return s;
}

// The offsets a program's own table and point indices shift by once its
// data is appended to the module.
struct GlslBases {
    std::uint32_t perm = 0;
    std::uint32_t points = 0;   // in points, not floats
};

inline GlslBases append_data(GlslModule& m, const std::vector<std::uint8_t>& tables,
                             const std::vector<double>& points) {
    GlslBases b{static_cast<std::uint32_t>(m.perm.size()),
                static_cast<std::uint32_t>(m.points.size() / 3)};
    for (auto byte : tables) m.perm.push_back(byte);
    for (auto x : points) m.points.push_back(static_cast<float>(x));
    return b;
}

// One scalar program as a function body: `float r<i> = ...;` per op.
inline std::string scalar_body(const PodOp<double>* ops, std::uint32_t count, GlslBases b) {
    std::string out;
    auto r = [](std::uint32_t i) { return std::format("r{}", i); };
    for (std::uint32_t i = 0; i < count; ++i) {
        const auto& n = ops[i];
        std::string e;
        switch (static_cast<Op>(n.code)) {
            case Op::Const:  e = glsl_float(n.value); break;
            case Op::U:      e = "in_.u"; break;
            case Op::V:      e = "in_.v"; break;
            case Op::Add:    e = std::format("{} + {}", r(n.a), r(n.b)); break;
            case Op::Sub:    e = std::format("{} - {}", r(n.a), r(n.b)); break;
            case Op::Mul:    e = std::format("{} * {}", r(n.a), r(n.b)); break;
            case Op::Div:    e = std::format("{} / {}", r(n.a), r(n.b)); break;
            case Op::Min:    e = std::format("min({}, {})", r(n.a), r(n.b)); break;
            case Op::Max:    e = std::format("max({}, {})", r(n.a), r(n.b)); break;
            case Op::Sin:    e = std::format("sin({})", r(n.a)); break;
            case Op::Cos:    e = std::format("cos({})", r(n.a)); break;
            case Op::Less:   e = std::format("({} < {} ? 1.0 : 0.0)", r(n.a), r(n.b)); break;
            case Op::Noise:
                e = std::format("field_noise({}u, {}, {}, {})",
                                b.perm + n.table * static_cast<std::uint32_t>(kNoiseTableBytes),
                                r(n.a), r(n.b), r(n.c));
                break;
            case Op::Id:     e = "float(in_.id)"; break;
            case Op::Origin: e = std::format("in_.origin[{}]", n.k); break;
            case Op::Point:  e = std::format("in_.p[{}]", n.k); break;
            case Op::Hash:   e = std::format("field_instance_unit(in_.id, {}u)", n.k); break;
            case Op::Sqrt:   e = std::format("sqrt({})", r(n.a)); break;
            case Op::Gather:
                e = std::format("field_gather({}u, {}u, {}, {}u)", b.points + n.table, n.extent,
                                r(n.a), n.k);
                break;
            case Op::Opaque: e = "0.0"; break;   // unreachable: lower() refuses it
        }
        out += std::format("    float {} = {};\n", r(i), e);
    }
    out += std::format("    return {};\n", r(count - 1));
    return out;
}

}  // namespace detail

// `float <name>(FieldIn in_)` for a scalar field.
inline Result<void> emit_scalar(GlslModule& m, const Field<double>& f, const std::string& name) {
    auto pod = lower(f);
    if (!pod) return std::unexpected(std::move(pod.error()));
    const auto b = detail::append_data(m, pod->tables, pod->points);
    m.code += std::format("float {}(FieldIn in_) {{\n{}}}\n", name,
                          detail::scalar_body(pod->ops.data(),
                                              static_cast<std::uint32_t>(pod->ops.size()), b));
    return {};
}

// `vec3 <name>(FieldIn in_)` for a point field, with each factor and each
// Make component emitted as its own function `<name>_f<op>_<k>` ahead of
// it. `in_.u` is time and `in_.p` the point, as `VecField::inputs_of` has
// them.
inline Result<void> emit_vector(GlslModule& m, const VecField<double>& f, const std::string& name) {
    auto pod = lower(f);
    if (!pod) return std::unexpected(std::move(pod.error()));
    const auto b = detail::append_data(m, pod->scalars.tables, pod->scalars.points);
    auto fname = [&](std::size_t i, int k) { return std::format("{}_f{}_{}", name, i, k); };

    for (std::size_t i = 0; i < pod->ops.size(); ++i) {
        const auto& n = pod->ops[i];
        const auto op = static_cast<VecOp>(n.code);
        const int count = op == VecOp::Scale ? 1 : (op == VecOp::Rotate || op == VecOp::Make) ? 3 : 0;
        for (int k = 0; k < count; ++k)
            m.code += std::format(
                "float {}(FieldIn in_) {{\n{}}}\n", fname(i, k),
                detail::scalar_body(pod->scalars.ops.data() + n.factor[k].begin, n.factor[k].count, b));
    }

    std::string body;
    auto v = [](std::uint32_t i) { return std::format("v{}", i); };
    for (std::size_t i = 0; i < pod->ops.size(); ++i) {
        const auto& n = pod->ops[i];
        std::string e;
        switch (static_cast<VecOp>(n.code)) {
            case VecOp::Point:  e = "in_.p"; break;
            case VecOp::Const:
                e = std::format("vec3({}, {}, {})", detail::glsl_float(n.value[0]),
                                detail::glsl_float(n.value[1]), detail::glsl_float(n.value[2]));
                break;
            case VecOp::Add:    e = std::format("{} + {}", v(n.a), v(n.b)); break;
            case VecOp::Sub:    e = std::format("{} - {}", v(n.a), v(n.b)); break;
            case VecOp::Scale:  e = std::format("{} * {}(in_)", v(n.a), fname(i, 0)); break;
            case VecOp::Rotate:
                e = std::format("field_so3_exp(vec3({}(in_), {}(in_), {}(in_))) * {}", fname(i, 0),
                                fname(i, 1), fname(i, 2), v(n.a));
                break;
            case VecOp::Make:
                e = std::format("vec3({}(in_), {}(in_), {}(in_))", fname(i, 0), fname(i, 1),
                                fname(i, 2));
                break;
            case VecOp::Opaque: e = "vec3(0.0)"; break;   // unreachable: lower() refuses it
        }
        body += std::format("    vec3 {} = {};\n", v(static_cast<std::uint32_t>(i)), e);
    }
    body += std::format("    return {};\n", v(static_cast<std::uint32_t>(pod->ops.size() - 1)));
    m.code += std::format("vec3 {}(FieldIn in_) {{\n{}}}\n", name, body);
    return {};
}

// `void <name>_place(FieldIn in_, out vec3 tr, out mat3 R, out float s)`:
// `VecField::placement_at`, for a motion that is a placement. Translation
// is the motion at the point's origin; the linear part is the Scale and
// Rotate factors along the one path that reaches the point, innermost
// first, as `linear_on_path` composes them. Emits the motion itself as
// `<name>` too.
inline Result<void> emit_placement(GlslModule& m, const VecField<double>& f,
                                   const std::string& name) {
    if (!f.is_placement())
        return std::unexpected(Error(ErrorCode::InvalidArgument,
            "emit_placement: the motion is a deformation, and a deformation has no placement"));
    if (auto r = emit_vector(m, f, name); !r) return r;

    // The path from the root to the Point, then walked back out.
    std::vector<std::string> steps;
    auto walk = [&](auto&& self, std::size_t i) -> void {
        const auto& n = f.op(i);
        switch (n.op) {
            case VecOp::Point: return;
            case VecOp::Scale:
                self(self, n.a);
                steps.push_back(std::format("    s *= {}_f{}_0(in_);\n", name, i));
                return;
            case VecOp::Rotate:
                self(self, n.a);
                steps.push_back(std::format(
                    "    R = field_so3_exp(vec3({0}_f{1}_0(in_), {0}_f{1}_1(in_), {0}_f{1}_2(in_))) * R;\n",
                    name, i));
                return;
            case VecOp::Add:
            case VecOp::Sub:
                if (f.reads_point(n.a)) { self(self, n.a); return; }
                if (f.reads_point(n.b)) {
                    self(self, n.b);
                    if (n.op == VecOp::Sub) steps.push_back("    s = -s;\n");
                    return;
                }
                steps.push_back("    s = 0.0;\n");   // the point does not appear
                return;
            default:
                steps.push_back("    s = 0.0;\n");
                return;
        }
    };
    walk(walk, f.size() - 1);

    std::string body = std::format(
        "    FieldIn at0 = in_;\n    at0.p = vec3(0.0);\n    tr = {}(at0);\n"
        "    R = mat3(1.0);\n    s = 1.0;\n", name);
    for (const auto& st : steps) body += st;
    m.code += std::format("void {}_place(FieldIn in_, out vec3 tr, out mat3 R, out float s) {{\n{}}}\n",
                          name, body);
    return {};
}

}  // namespace spatium::io::build
