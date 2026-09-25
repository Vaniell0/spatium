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
#  include <unordered_map>
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
struct FieldIn {{ float u; float v; uint id; vec3 origin; vec3 p; vec4 x; }};

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
//
// A kernel that defines FIELD_PERM_SHARED to the table's length reads it
// from workgroup memory instead, after field_load_perm() has copied it
// there: sixteen dependent random loads per noise sample are what a
// particle's motion spends most of its time on, and the table is 2 KB.
#ifdef FIELD_PERM_SHARED
shared uint field_perm_s[FIELD_PERM_SHARED];
int field_perm_at(uint base, int i) {{ return int(field_perm_s[base + uint(i & 511)]); }}
void field_load_perm() {{
    for (uint k = gl_LocalInvocationIndex; k < uint(FIELD_PERM_SHARED); k += gl_WorkGroupSize.x * gl_WorkGroupSize.y)
        field_perm_s[k] = field_perm[k];
    barrier();
}}
#else
int field_perm_at(uint base, int i) {{ return int(field_perm[base + uint(i & 511)]); }}
void field_load_perm() {{}}
#endif
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
    // A block already in the module is reused whole, so two programs that
    // gather from one table name the same offsets -- and a gather in a
    // motion and the same gather in its colour become one value in a scope.
    auto find_block = [](const auto& hay, const auto& needle, std::size_t stride) -> std::ptrdiff_t {
        if (needle.empty()) return 0;
        for (std::size_t at = 0; at + needle.size() <= hay.size(); at += stride) {
            bool same = true;
            for (std::size_t k = 0; k < needle.size() && same; ++k)
                same = hay[at + k] == needle[k];
            if (same) return static_cast<std::ptrdiff_t>(at);
        }
        return -1;
    };
    std::vector<std::uint32_t> perm(tables.begin(), tables.end());
    std::vector<float> pts;
    for (auto x : points) pts.push_back(static_cast<float>(x));
    GlslBases b;
    auto pa = find_block(m.perm, perm, kNoiseTableBytes);
    if (pa < 0) { pa = static_cast<std::ptrdiff_t>(m.perm.size()); m.perm.insert(m.perm.end(), perm.begin(), perm.end()); }
    auto qa = find_block(m.points, pts, 3);
    if (qa < 0) { qa = static_cast<std::ptrdiff_t>(m.points.size()); m.points.insert(m.points.end(), pts.begin(), pts.end()); }
    b.perm = static_cast<std::uint32_t>(pa);
    b.points = static_cast<std::uint32_t>(qa / 3);
    return b;
}

}  // namespace detail

// A run of straight-line code in which an expression is written once and
// named, however many times it is asked for: value numbering, keyed by the
// expression's text with its operands already replaced by the names of
// their values. Two equal subtrees therefore produce equal text and one
// variable -- across the three components of a Make, across a motion and
// its colour, across everything emitted into one scope. That is what the
// flat IR's cheap structural identity was for, and it is also where the
// dust's cost was: each component computed the hash, the target and the
// pull again.
struct GlslScope {
    std::string body;
    std::unordered_map<std::string, std::string> named;
    int next = 0;

    // The name of `expr`'s value, emitting it on first use.
    std::string value(const std::string& type, const std::string& expr) {
        const std::string key = type + '|' + expr;
        if (auto it = named.find(key); it != named.end()) return it->second;
        std::string name = std::format("c{}", next++);
        body += std::format("    {} {} = {};\n", type, name, expr);
        named.emplace(key, name);
        return name;
    }
    void line(const std::string& text) { body += "    " + text + "\n"; }
};

// Where an emitted program reads its inputs from: `rec` for time, instance
// and origin, and `point` for the point, which is the one input a colour
// and a motion evaluated for the same instance do not share.
struct GlslInputs {
    std::string rec = "in_";
    std::string point = "in_.p";
};

// One scalar program into `sc`; returns the name of its value.
inline std::string emit_scalar_into(GlslScope& sc, const PodOp<double>* ops, std::uint32_t count,
                                    detail::GlslBases b, const GlslInputs& in) {
    std::vector<std::string> v(count);
    for (std::uint32_t i = 0; i < count; ++i) {
        const auto& n = ops[i];
        const auto& A = n.a < i ? v[n.a] : v[0];
        const auto& B = n.b < i ? v[n.b] : v[0];
        const auto& C = n.c < i ? v[n.c] : v[0];
        std::string e;
        switch (static_cast<Op>(n.code)) {
            // Literals are not named: they are free to repeat, and inline
            // they let the compiler fold.
            case Op::Const:  v[i] = detail::glsl_float(n.value); continue;
            case Op::U:      e = in.rec + ".u"; break;
            case Op::V:      e = in.rec + ".v"; break;
            case Op::Add:    e = std::format("{} + {}", A, B); break;
            case Op::Sub:    e = std::format("{} - {}", A, B); break;
            case Op::Mul:    e = std::format("{} * {}", A, B); break;
            case Op::Div:    e = std::format("{} / {}", A, B); break;
            case Op::Min:    e = std::format("min({}, {})", A, B); break;
            case Op::Max:    e = std::format("max({}, {})", A, B); break;
            case Op::Sin:    e = std::format("sin({})", A); break;
            case Op::Cos:    e = std::format("cos({})", A); break;
            case Op::Less:   e = std::format("({} < {} ? 1.0 : 0.0)", A, B); break;
            case Op::Noise:
                e = std::format("field_noise({}u, {}, {}, {})",
                                b.perm + n.table * static_cast<std::uint32_t>(kNoiseTableBytes), A, B, C);
                break;
            case Op::Id:     e = "float(" + in.rec + ".id)"; break;
            case Op::Origin: e = std::format("{}.origin[{}]", in.rec, n.k); break;
            case Op::Point:  e = std::format("{}[{}]", in.point, n.k); break;
            case Op::Hash:   e = std::format("field_instance_unit({}.id, {}u)", in.rec, n.k); break;
            case Op::Sqrt:   e = std::format("sqrt({})", A); break;
            case Op::Gather:
                e = std::format("field_gather({}u, {}u, {}, {}u)", b.points + n.table, n.extent, A, n.k);
                break;
            case Op::Opaque: v[i] = "0.0"; continue;   // unreachable: lower() refuses it
            case Op::Coord:  e = std::format("{}.x[{}]", in.rec, n.k); break;
        }
        v[i] = sc.value("float", e);
    }
    return v[count - 1];
}

// One lowered point field into `sc`; returns the name of its value.
// `point_expr` is what the Point op reads -- the input point, or `vec3(0.0)`
// for a placement's translation.
inline std::string emit_vector_into(GlslScope& sc, const PodVecField<double>& pod, detail::GlslBases b,
                                    const GlslInputs& in, const std::string& point_expr) {
    auto factor = [&](const PodRange& r) {
        return emit_scalar_into(sc, pod.scalars.ops.data() + r.begin, r.count, b, in);
    };
    std::vector<std::string> v(pod.ops.size());
    for (std::size_t i = 0; i < pod.ops.size(); ++i) {
        const auto& n = pod.ops[i];
        std::string e;
        switch (static_cast<VecOp>(n.code)) {
            case VecOp::Point:  v[i] = point_expr; continue;
            case VecOp::Const:
                v[i] = std::format("vec3({}, {}, {})", detail::glsl_float(n.value[0]),
                                   detail::glsl_float(n.value[1]), detail::glsl_float(n.value[2]));
                continue;
            case VecOp::Add:    e = std::format("{} + {}", v[n.a], v[n.b]); break;
            case VecOp::Sub:    e = std::format("{} - {}", v[n.a], v[n.b]); break;
            case VecOp::Scale:  e = std::format("{} * {}", v[n.a], factor(n.factor[0])); break;
            case VecOp::Rotate: {
                const auto x = factor(n.factor[0]), y = factor(n.factor[1]), z = factor(n.factor[2]);
                const auto rot = sc.value("mat3", std::format("field_so3_exp(vec3({}, {}, {}))", x, y, z));
                e = std::format("{} * {}", rot, v[n.a]);
                break;
            }
            case VecOp::Make: {
                const auto x = factor(n.factor[0]), y = factor(n.factor[1]), z = factor(n.factor[2]);
                e = std::format("vec3({}, {}, {})", x, y, z);
                break;
            }
            case VecOp::Opaque: v[i] = "vec3(0.0)"; continue;   // unreachable
        }
        v[i] = sc.value("vec3", e);
    }
    return v[pod.ops.size() - 1];
}

// The names a placement leaves in its scope: `p -> R s p + tr`.
struct GlslPlacement {
    std::string translation, rotation, scale;
};

// `VecField::placement_at`, into `sc`. Translation is the motion at the
// point's origin; the linear part is the Scale and Rotate factors along
// the one path that reaches the point, innermost first, as
// `linear_on_path` composes them. The factors are the same values the
// translation's evaluation already named, so they cost nothing twice.
inline Result<GlslPlacement> emit_placement_into(GlslScope& sc, GlslModule& m,
                                                 const VecField<double>& f, const GlslInputs& in) {
    if (!f.is_placement())
        return std::unexpected(Error(ErrorCode::InvalidArgument,
            "emit_placement: the motion is a deformation, and a deformation has no placement"));
    auto pod = lower(f);
    if (!pod) return std::unexpected(std::move(pod.error()));
    const auto b = detail::append_data(m, pod->scalars.tables, pod->scalars.points);

    GlslPlacement out;
    out.translation = emit_vector_into(sc, *pod, b, in, "vec3(0.0)");
    auto factor = [&](const PodRange& r) {
        return emit_scalar_into(sc, pod->scalars.ops.data() + r.begin, r.count, b, in);
    };

    std::string R = "mat3(1.0)", s = "1.0";
    auto walk = [&](auto&& self, std::size_t i) -> void {
        const auto& n = f.op(i);
        const auto& p = pod->ops[i];
        switch (n.op) {
            case VecOp::Point: return;
            case VecOp::Scale:
                self(self, n.a);
                s = sc.value("float", std::format("{} * {}", s, factor(p.factor[0])));
                return;
            case VecOp::Rotate: {
                self(self, n.a);
                const auto x = factor(p.factor[0]), y = factor(p.factor[1]), z = factor(p.factor[2]);
                const auto rot = sc.value("mat3", std::format("field_so3_exp(vec3({}, {}, {}))", x, y, z));
                R = sc.value("mat3", std::format("{} * {}", rot, R));
                return;
            }
            case VecOp::Add:
            case VecOp::Sub:
                if (f.reads_point(n.a)) { self(self, n.a); return; }
                if (f.reads_point(n.b)) {
                    self(self, n.b);
                    if (n.op == VecOp::Sub) s = sc.value("float", "-" + s);
                    return;
                }
                s = "0.0";   // the point does not appear
                return;
            default:
                s = "0.0";
                return;
        }
    };
    walk(walk, f.size() - 1);
    out.rotation = R;
    out.scale = s;
    return out;
}

// `float <name>(FieldIn in_)` for a scalar field.
inline Result<void> emit_scalar(GlslModule& m, const Field<double>& f, const std::string& name) {
    auto pod = lower(f);
    if (!pod) return std::unexpected(std::move(pod.error()));
    const auto b = detail::append_data(m, pod->tables, pod->points);
    GlslScope sc;
    const auto r = emit_scalar_into(sc, pod->ops.data(), static_cast<std::uint32_t>(pod->ops.size()), b, {});
    m.code += std::format("float {}(FieldIn in_) {{\n{}    return {};\n}}\n", name, sc.body, r);
    return {};
}

// `vec3 <name>(FieldIn in_)` for a point field: `in_.u` is time and `in_.p`
// the point, as `VecField::inputs_of` has them.
inline Result<void> emit_vector(GlslModule& m, const VecField<double>& f, const std::string& name) {
    auto pod = lower(f);
    if (!pod) return std::unexpected(std::move(pod.error()));
    const auto b = detail::append_data(m, pod->scalars.tables, pod->scalars.points);
    GlslScope sc;
    const auto r = emit_vector_into(sc, *pod, b, {}, "in_.p");
    m.code += std::format("vec3 {}(FieldIn in_) {{\n{}    return {};\n}}\n", name, sc.body, r);
    return {};
}

// `void <name>_place(FieldIn in_, out vec3 tr, out mat3 R, out float s)`.
inline Result<void> emit_placement(GlslModule& m, const VecField<double>& f, const std::string& name) {
    GlslScope sc;
    auto pl = emit_placement_into(sc, m, f, {});
    if (!pl) return std::unexpected(std::move(pl.error()));
    m.code += std::format(
        "void {}_place(FieldIn in_, out vec3 tr, out mat3 R, out float s) {{\n{}"
        "    tr = {};\n    R = {};\n    s = {};\n}}\n",
        name, sc.body, pl->translation, pl->rotation, pl->scale);
    return {};
}

}  // namespace spatium::io::build
