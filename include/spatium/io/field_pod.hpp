#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/error.hpp>
#  include <spatium/io/field.hpp>
#  include <spatium/algebra/noise.hpp>
#  include <spatium/algebra/groups/so3.hpp>
#  include <cmath>
#  include <cstddef>
#  include <cstdint>
#  include <cstring>
#  include <format>
#  include <string>
#  include <type_traits>
#  include <vector>
#endif

// A field as plain data, and the interpreter that runs it.
//
// This is the step `docs/gpu-abi-design.md` names as the first one before
// any backend: lower a `Field` or `VecField` to arrays of op codes, child
// indices and numbers -- nothing erased, nothing templated over a
// callable, nothing pointing back into the field it came from -- and
// evaluate that, alongside `eval_into`, until the two agree **bit for
// bit** over a real scene. Device code cannot run a `std::function`, so
// what cannot be written this way cannot be evaluated on a GPU at all;
// what can be written this way and agrees exactly has shown it crosses.
//
// **Bit for bit, never within a tolerance.** A tolerance passes under a
// different order of operations, and order of operations is precisely
// what an export can get wrong without anything else noticing -- so the
// interpreter below is written switch-case for switch-case against
// `Field::eval_into` and `VecField::eval_into`, down to evaluating a
// factor as its own field at (t, 0) exactly where a Scale or a Rotate
// reads it. It is a second implementation on purpose: sharing the loop
// with the first would make agreement a tautology. What it does share are
// definitions rather than ways of walking a tree -- `std::sin`, SO3's
// `exp`, and `PerlinNoise::sample`, which is the noise itself run over a
// table that is now plain data.
//
// **Lowering fails, it does not approximate.** An opaque leaf is a C++
// closure and has no plain-data form, so `lower()` returns an error
// naming the op and where it sits rather than a program that silently
// omits it. The count of fields that lower is then a fact about the
// scene, and it must equal what `field_report()` calls structural -- if
// the two ever differ, one of them is describing something it is not.

SPATIUM_EXPORT namespace spatium::io::build {

// One scalar op. Every member is a number; `code` is an `Op` by value.
//
// Child indices are relative to the start of the program the op belongs
// to, not to any larger buffer, so a lowered field can be placed anywhere
// in a pool -- which is how a VecField's factors share one -- without
// rewriting a single index.
template<Scalar T = double>
struct PodOp {
    std::uint8_t  code  = 0;
    std::uint32_t a = 0, b = 0, c = 0;
    std::uint32_t table = 0;   // Noise: which 512-byte table in `tables`
    T             value{};     // Const
};

// A lowered scalar field: its ops in the original topological order, and
// the noise tables they name, 512 bytes each, back to back. Two noise ops
// built from the same seed share one table.
template<Scalar T = double>
struct PodField {
    std::vector<PodOp<T>>     ops;
    std::vector<std::uint8_t> tables;
};

// Where one factor's program sits inside a VecField's shared scalar pool.
struct PodRange {
    std::uint32_t begin = 0, count = 0;
};

// One vector op. A Scale reads `factor[0]`; a Rotate reads all three, as
// the components of the axis-angle vector in order -- the same layout
// `VecFieldOp::factor_fields` has.
template<Scalar T = double>
struct PodVecOp {
    std::uint8_t  code = 0;
    std::uint32_t a = 0, b = 0;
    T             value[3]{};
    PodRange      factor[3]{};
};

// A lowered point field. Every factor's scalar program lives in
// `scalars.ops`, each at its own `PodRange`, over one table buffer.
template<Scalar T = double>
struct PodVecField {
    std::vector<PodVecOp<T>> ops;
    PodField<T>              scalars;
};

// Checked rather than hoped: these are what crosses a boundary as bytes,
// and a member that stopped being trivially copyable -- a string, an
// erased callable -- would stop that without failing to compile anywhere
// else.
static_assert(std::is_trivially_copyable_v<PodOp<double>> &&
              std::is_standard_layout_v<PodOp<double>>);
static_assert(std::is_trivially_copyable_v<PodVecOp<double>> &&
              std::is_standard_layout_v<PodVecOp<double>>);

inline constexpr std::size_t kNoiseTableBytes = 512;

namespace detail {

// A table already in the buffer is reused, compared by content: the
// shared_ptr a Field holds is an ownership detail, and two fields built
// from equal seeds should not carry two copies just because they were
// built apart.
inline std::uint32_t intern_table(std::vector<std::uint8_t>& tables,
                                  const algebra::PerlinNoise& n) {
    const auto& perm = n.permutation();
    const std::size_t count = tables.size() / kNoiseTableBytes;
    for (std::size_t k = 0; k < count; ++k)
        if (std::memcmp(tables.data() + k * kNoiseTableBytes, perm.data(), kNoiseTableBytes) == 0)
            return static_cast<std::uint32_t>(k);
    tables.insert(tables.end(), perm.begin(), perm.end());
    return static_cast<std::uint32_t>(count);
}

// Append `f`'s ops to `out`, returning where they landed. `where` names
// the field in the error, so a failure deep inside a VecField's factor
// still says which one.
template<Scalar T>
Result<PodRange> lower_into(PodField<T>& out, const Field<T>& f, const std::string& where) {
    PodRange r{static_cast<std::uint32_t>(out.ops.size()),
               static_cast<std::uint32_t>(f.size())};
    for (std::size_t i = 0; i < f.size(); ++i) {
        const auto& n = f.op(i);
        if (n.op == Op::Opaque)
            return std::unexpected(Error(ErrorCode::InvalidArgument,
                std::format("lower: {}op {} is {} -- a closure has no plain-data form",
                            where, i, op_name(n.op))));
        PodOp<T> p{};
        p.code  = static_cast<std::uint8_t>(n.op);
        p.a     = n.a;
        p.b     = n.b;
        p.c     = n.c;
        p.value = n.value;
        if (n.op == Op::Noise) p.table = intern_table(out.tables, *n.noise);
        out.ops.push_back(p);
    }
    return r;
}

} // namespace detail

template<Scalar T>
Result<PodField<T>> lower(const Field<T>& f) {
    PodField<T> out;
    auto r = detail::lower_into(out, f, "");
    if (!r) return std::unexpected(std::move(r.error()));
    return out;
}

template<Scalar T>
Result<PodVecField<T>> lower(const VecField<T>& f) {
    PodVecField<T> out;
    for (std::size_t i = 0; i < f.size(); ++i) {
        const auto& n = f.op(i);
        PodVecOp<T> p{};
        p.code = static_cast<std::uint8_t>(n.op);
        p.a    = n.a;
        p.b    = n.b;
        for (std::size_t k = 0; k < 3; ++k) p.value[k] = n.value[k];

        if (n.op == VecOp::Opaque)
            return std::unexpected(Error(ErrorCode::InvalidArgument,
                std::format("lower: op {} is {} -- a closure has no plain-data form",
                            i, vec_op_name(n.op))));

        // A factor given as a callable is a closure exactly as an Opaque
        // leaf is, and the report counts it as one; lowering has to refuse
        // it for the same reason, or the two counts would part company.
        if (n.op == VecOp::Scale || n.op == VecOp::Rotate) {
            const std::size_t want = n.op == VecOp::Scale ? 1 : 3;
            if (n.factor_fields.size() != want)
                return std::unexpected(Error(ErrorCode::InvalidArgument,
                    std::format("lower: op {} is {} with a closure factor -- a closure has "
                                "no plain-data form", i, vec_op_name(n.op))));
            for (std::size_t k = 0; k < want; ++k) {
                auto r = detail::lower_into(out.scalars, n.factor_fields[k],
                                            std::format("op {} ({}) factor {}: ",
                                                        i, vec_op_name(n.op), k));
                if (!r) return std::unexpected(std::move(r.error()));
                p.factor[k] = *r;
            }
        }
        out.ops.push_back(p);
    }
    return out;
}

// ── The interpreter ──────────────────────────────────────────────
//
// Reads nothing but the arrays. `scratch` must hold `count` values; it is
// the caller's so that a consumer evaluating many points can keep one,
// the way a GPU lane keeps its registers.

template<Scalar T>
T interpret(const PodOp<T>* ops, std::uint32_t count, const std::uint8_t* tables,
            T* s, T u, T v) {
    for (std::uint32_t i = 0; i < count; ++i) {
        const auto& n = ops[i];
        switch (static_cast<Op>(n.code)) {
            case Op::Const: s[i] = n.value; break;
            case Op::U:     s[i] = u; break;
            case Op::V:     s[i] = v; break;
            case Op::Add:   s[i] = s[n.a] + s[n.b]; break;
            case Op::Sub:   s[i] = s[n.a] - s[n.b]; break;
            case Op::Mul:   s[i] = s[n.a] * s[n.b]; break;
            case Op::Div:   s[i] = s[n.a] / s[n.b]; break;
            case Op::Min:   { using std::min; s[i] = min(s[n.a], s[n.b]); break; }
            case Op::Max:   { using std::max; s[i] = max(s[n.a], s[n.b]); break; }
            case Op::Sin:   { using std::sin; s[i] = sin(s[n.a]); break; }
            case Op::Cos:   { using std::cos; s[i] = cos(s[n.a]); break; }
            case Op::Less:  s[i] = s[n.a] < s[n.b] ? T{1} : T{0}; break;
            case Op::Noise:
                s[i] = algebra::PerlinNoise::sample(tables + std::size_t{n.table} * kNoiseTableBytes,
                                                    s[n.a], s[n.b], s[n.c]);
                break;
            // Unreachable from `lower()`, which refuses it. A zero rather
            // than an assert, so a hand-built program with one in it gives
            // a wrong number a bit-exact comparison will catch, not a
            // crash in a lane that has no way to report one.
            case Op::Opaque: s[i] = T{0}; break;
        }
    }
    return s[count - 1];
}

template<Scalar T>
T interpret(const PodField<T>& f, T u, T v) {
    std::vector<T> scratch(f.ops.size());
    return interpret(f.ops.data(), static_cast<std::uint32_t>(f.ops.size()),
                     f.tables.data(), scratch.data(), u, v);
}

template<Scalar T>
Vec<T, 3> interpret(const PodVecField<T>& f, const MotionEnv<T>& env) {
    std::vector<Vec<T, 3>> s(f.ops.size());
    std::vector<T> scratch(f.scalars.ops.size());

    // A factor is its own field evaluated at (t, 0) -- the convention
    // `scaled` and `rotated` fix -- run from its own start in the pool.
    auto factor = [&](const PodRange& r) {
        return interpret(f.scalars.ops.data() + r.begin, r.count, f.scalars.tables.data(),
                         scratch.data(), env.t, T{0});
    };

    for (std::size_t i = 0; i < f.ops.size(); ++i) {
        const auto& n = f.ops[i];
        switch (static_cast<VecOp>(n.code)) {
            case VecOp::Point:  s[i] = env.p; break;
            case VecOp::Const:  s[i] = Vec<T, 3>{n.value[0], n.value[1], n.value[2]}; break;
            case VecOp::Add:    s[i] = Vec<T, 3>{s[n.a] + s[n.b]}; break;
            case VecOp::Sub:    s[i] = Vec<T, 3>{s[n.a] - s[n.b]}; break;
            case VecOp::Scale:  s[i] = Vec<T, 3>{s[n.a] * factor(n.factor[0])}; break;
            case VecOp::Rotate: {
                // Three evaluations in component order, then one exp, as
                // the expression form of `rotated` does it.
                const T ax = factor(n.factor[0]);
                const T ay = factor(n.factor[1]);
                const T az = factor(n.factor[2]);
                s[i] = Vec<T, 3>{algebra::SO3<T>{}.exp(Vec<T, 3>{ax, ay, az}) * s[n.a]};
                break;
            }
            case VecOp::Opaque: s[i] = Vec<T, 3>{}; break;   // see the scalar case
        }
    }
    return s[f.ops.size() - 1];
}

} // namespace spatium::io::build
