#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/error.hpp>
#  include <spatium/io/field.hpp>
#  include <spatium/algebra/noise.hpp>
#  include <spatium/algebra/groups/so3.hpp>
#  include <bit>
#  include <charconv>
#  include <concepts>
#  include <optional>
#  include <cmath>
#  include <cstddef>
#  include <cstdint>
#  include <cstring>
#  include <algorithm>
#  include <format>
#  include <limits>
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
// factor as its own field -- time first, then the instance, its origin and
// the point -- exactly where a Scale, a Rotate or a Make reads it. It is a second implementation on purpose: sharing the loop
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
//
// **A symbol is a slot.** A named symbol lowers to a read of the slot its
// name is given in `symbols`, and its value is passed in when the program
// runs: the program is the field's structure, the slot its number. So a
// parameter moved by a slider, or time, changes data and not the program
// -- on a device, a uniform rather than a new shader. Run with no value for
// a slot, a symbol reads NaN, as the unbound field it came from does.

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
    std::uint32_t table = 0;   // Noise: which 512-byte table in `tables`;
                               // Gather: the first point of its table in `points`;
                               // Net: its first weight in `weights`
    std::uint32_t extent = 0;  // Gather: how many points that table has; Net: its inputs, F
    std::uint32_t k = 0;       // Origin, Point, Gather: component; Hash: salt; Net: its units, H
    T             value{};     // Const
};

// A lowered scalar field: its ops in the original topological order, and
// the noise tables they name, 512 bytes each, back to back. Two noise ops
// built from the same seed share one table.
template<Scalar T = double>
struct PodField {
    std::vector<PodOp<T>>     ops;
    std::vector<std::uint8_t> tables;
    // Every Gather's points, x y z per point, tables back to back and
    // shared by content the way noise tables are.
    std::vector<T>            points;
    // The name of each slot a Symbol op reads (its `k`), in the order the
    // ops first read them; one slot per name, however many ops read it.
    std::vector<std::string>  symbols;
    // Every Net's weights, in `TanhNet`'s order, nets back to back and
    // shared by content. Numbers the program reads, so a net trained
    // further changes these and not the ops.
    std::vector<T>            weights;
};

// Where one factor's program sits inside a VecField's shared scalar pool.
struct PodRange {
    std::uint32_t begin = 0, count = 0;
};

// One vector op. A Scale reads `factor[0]`; a Rotate reads all three, as
// the components of the axis-angle vector in order, and a Make as the
// components of the vector it makes -- the layout `VecFieldOp::factor_fields`
// has.
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

// Equal as numbers would join 0.0 with -0.0, and a table that differs from
// another only in the sign of a zero is a different table (1/x tells them
// apart), so floating-point tables are compared by bit pattern.
template<Scalar T>
bool same_number(const T& x, const T& y) {
    if constexpr (std::is_same_v<T, float>)
        return std::bit_cast<std::uint32_t>(x) == std::bit_cast<std::uint32_t>(y);
    else if constexpr (std::is_same_v<T, double>)
        return std::bit_cast<std::uint64_t>(x) == std::bit_cast<std::uint64_t>(y);
    else
        return x == y;
}

template<Scalar T>
std::uint32_t intern_points(std::vector<T>& points, const std::vector<Vec<T, 3>>& table,
                            std::uint32_t& start) {
    // Search for an equal table already present, at any point boundary --
    // two gathers from one table are the common case, and only that case
    // is worth the search.
    const std::size_t n = table.size(), have = points.size() / 3;
    for (std::size_t s0 = 0; s0 + n <= have; ++s0) {
        bool same = true;
        for (std::size_t i = 0; i < n && same; ++i)
            for (std::size_t k = 0; k < 3 && same; ++k)
                same = same_number(points[(s0 + i) * 3 + k], table[i][k]);
        if (same) { start = static_cast<std::uint32_t>(s0); return static_cast<std::uint32_t>(n); }
    }
    start = static_cast<std::uint32_t>(have);
    for (const auto& p : table)
        for (std::size_t k = 0; k < 3; ++k) points.push_back(p[k]);
    return static_cast<std::uint32_t>(n);
}

// Append `f`'s ops to `out`, returning where they landed. `where` names
// the field in the error, so a failure deep inside a VecField's factor
// still says which one.
template<Scalar T>
Result<PodRange> lower_into(PodField<T>& out, const Field<T>& whole, const std::string& where) {
    // Only what the root reads: an op nothing reads would be emitted, sent
    // to a device and evaluated for nothing. (The builders never leave
    // one; a pool assembled elsewhere can.) Op numbers in an error below
    // are positions in this culled field.
    const Field<T> f = whole.reachable_from_root();
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
        p.k     = n.k;
        if (n.op == Op::Noise) p.table = intern_table(out.tables, *n.noise);
        if (n.op == Op::Gather) p.extent = intern_points(out.points, *n.table, p.table);
        if (n.op == Op::Net) {
            const auto& w = *n.weights;
            std::size_t at = 0;
            while (at + w.size() <= out.weights.size() && !std::equal(w.begin(), w.end(), out.weights.begin() + at)) ++at;
            if (at + w.size() > out.weights.size()) at = out.weights.size(), out.weights.insert(out.weights.end(), w.begin(), w.end());
            p.table = static_cast<std::uint32_t>(at);
            p.extent = static_cast<std::uint32_t>(Field<T>::net_inputs(w.size(), n.k));
        }
        if (n.op == Op::Symbol) {
            const auto at = std::ranges::find(out.symbols, n.name);
            p.k = static_cast<std::uint32_t>(at - out.symbols.begin());
            if (at == out.symbols.end()) out.symbols.push_back(n.name);
        }
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
        if (n.op == VecOp::Scale || n.op == VecOp::Rotate || n.op == VecOp::Make) {
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

// ── Other arithmetics ────────────────────────────────────────────
//
// The interpreter below is one program for every arithmetic it is run in: a float, as a device
// runs it; a double; Real50; a fraction (core/rational.hpp) or a decimal (core/decimal.hpp), in
// which 0.1 is exact. An arithmetic with no square root has no sine either: such an arithmetic
// computes only the ops `exact_op` names, and a program reading any other is refused by that op's
// name (`beyond`) before it runs.

namespace detail {
using std::sqrt, std::sin, std::cos, std::log, std::exp;
template<class T>
concept transcendental = requires(T a) { sqrt(a); sin(a); cos(a); log(a); exp(a); };
} // namespace detail

// The ops an arithmetic with nothing but + - * / and order computes.
inline bool exact_op(Op o) {
    switch (o) {
        case Op::Const: case Op::U: case Op::V: case Op::Add: case Op::Sub: case Op::Mul: case Op::Div:
        case Op::Min: case Op::Max: case Op::Less: case Op::Abs: case Op::Symbol: case Op::Coord:
        case Op::Origin: case Op::Point: case Op::Opaque:
            return true;
        default:
            return false;
    }
}

// The first op of `f` that arithmetic T cannot compute, if any.
template<Scalar T>
std::optional<Op> beyond(const PodField<double>& f) {
    if constexpr (!detail::transcendental<T>)
        for (const auto& n : f.ops)
            if (!exact_op(static_cast<Op>(n.code))) return static_cast<Op>(n.code);
    return std::nullopt;
}

// A double as arithmetic T: rounded for a float, and read through its shortest decimal form for
// an arithmetic that holds decimals (Real50 from its text, a fraction or a decimal by their own
// constructor), so a literal 0.1 is the 0.1 that was typed.
template<Scalar T>
T from_double(double x) {
    if constexpr (std::is_floating_point_v<T>) {
        return static_cast<T>(x);
    } else if constexpr (std::constructible_from<T, const char*>) {
        char buf[32];
        const auto [end, ec] = std::to_chars(buf, buf + sizeof buf - 1, x);
        *(ec == std::errc{} ? end : buf) = '\0';
        return T(buf);
    } else {
        return T(x);
    }
}

// The same program in arithmetic T: each number as T holds it, the ops and tables as they are.
template<Scalar T>
PodField<T> converted(const PodField<double>& f) {
    PodField<T> out;
    out.ops.reserve(f.ops.size());
    for (const auto& n : f.ops)
        out.ops.push_back(PodOp<T>{n.code, n.a, n.b, n.c, n.table, n.extent, n.k, from_double<T>(n.value)});
    out.tables = f.tables;
    for (double p : f.points) out.points.push_back(from_double<T>(p));
    out.symbols = f.symbols;
    for (double w : f.weights) out.weights.push_back(from_double<T>(w));
    return out;
}

// The same program at fp32, as a device runs it: each number rounded to
// the nearest float, the ops and tables as they are. What a shader is
// compared against on the host.
inline PodField<float> narrowed(const PodField<double>& f) { return converted<float>(f); }

namespace detail {

// The ops past + - * / and order, for an arithmetic that has them.
template<Scalar T>
T inexact_step(const PodOp<T>& n, const T* s, const std::uint8_t* tables, const T* points, const T* weights,
               const FieldInputs<T>& in) {
    switch (static_cast<Op>(n.code)) {
        case Op::Sin:  { using std::sin; return sin(s[n.a]); }
        case Op::Cos:  { using std::cos; return cos(s[n.a]); }
        case Op::Sqrt: { using std::sqrt; return sqrt(s[n.a]); }
        case Op::Log:  { using std::log; return log(s[n.a]); }
        case Op::Exp:  { using std::exp; return exp(s[n.a]); }
        case Op::Noise:
            return algebra::PerlinNoise::sample(tables + std::size_t{n.table} * kNoiseTableBytes, s[n.a], s[n.b], s[n.c]);
        case Op::Id:   return static_cast<T>(in.id);
        case Op::Hash: return instance_unit<T>(in.id, n.k);
        case Op::Gather: {
            // Field::gather_at's rounding, over the flat buffer.
            if (n.extent == 0) return T{0};
            const T lo = s[n.a] < T{0} ? T{0} : s[n.a];
            auto at = static_cast<std::size_t>(lo);
            if (at >= n.extent) at = n.extent - 1;
            return points[(std::size_t{n.table} + at) * 3 + n.k];
        }
        case Op::Net: {
            const T x[3] = {s[n.a], s[n.b], s[n.c]};
            return Field<T>::net_at(weights + n.table, n.extent, n.k, x);
        }
        default: return T{0};
    }
}

} // namespace detail

// ── The interpreter ──────────────────────────────────────────────
//
// Reads nothing but the arrays. `scratch` must hold `count` values; it is
// the caller's so that a consumer evaluating many points can keep one,
// the way a GPU lane keeps its registers. `slots` holds `slot_count`
// values, a symbol's by its slot; one past them reads NaN.

template<Scalar T>
T interpret(const PodOp<T>* ops, std::uint32_t count, const std::uint8_t* tables,
            const T* points, T* s, const FieldInputs<T>& in,
            const T* slots = nullptr, std::size_t slot_count = 0, const T* weights = nullptr) {
    for (std::uint32_t i = 0; i < count; ++i) {
        const auto& n = ops[i];
        switch (static_cast<Op>(n.code)) {
            case Op::Const: s[i] = n.value; break;
            case Op::U:     s[i] = in.u; break;
            case Op::V:     s[i] = in.v; break;
            case Op::Add:   s[i] = s[n.a] + s[n.b]; break;
            case Op::Sub:   s[i] = s[n.a] - s[n.b]; break;
            case Op::Mul:   s[i] = s[n.a] * s[n.b]; break;
            case Op::Div:   s[i] = s[n.a] / s[n.b]; break;
            case Op::Min:   { using std::min; s[i] = min(s[n.a], s[n.b]); break; }
            case Op::Max:   { using std::max; s[i] = max(s[n.a], s[n.b]); break; }
            case Op::Less:  s[i] = s[n.a] < s[n.b] ? T{1} : T{0}; break;
            case Op::Origin: s[i] = in.origin[n.k]; break;
            case Op::Point:  s[i] = in.p[n.k]; break;
            case Op::Abs:    { using std::abs; s[i] = abs(s[n.a]); break; }
            default:
                // the rest are computed in an arithmetic that has them; `beyond` refuses them elsewhere
                if constexpr (detail::transcendental<T>) s[i] = detail::inexact_step(n, s, tables, points, weights, in);
                else s[i] = T(std::numeric_limits<double>::quiet_NaN());
                break;
            // Unreachable from `lower()`, which refuses it. A zero rather
            // than an assert, so a hand-built program with one in it gives
            // a wrong number a bit-exact comparison will catch, not a
            // crash in a lane that has no way to report one.
            case Op::Opaque: s[i] = T{0}; break;
            case Op::Coord:  s[i] = in.x[n.k]; break;
            case Op::Symbol:
                s[i] = n.k < slot_count ? slots[n.k] : T(std::numeric_limits<double>::quiet_NaN());
                break;
        }
    }
    return s[count - 1];
}

// `slots` are the values of `f.symbols`, in that order.
template<Scalar T>
T interpret(const PodField<T>& f, const FieldInputs<T>& in, const std::vector<T>& slots = {}) {
    std::vector<T> scratch(f.ops.size());
    return interpret(f.ops.data(), static_cast<std::uint32_t>(f.ops.size()),
                     f.tables.data(), f.points.data(), scratch.data(), in, slots.data(), slots.size(), f.weights.data());
}

template<Scalar T>
T interpret(const PodField<T>& f, T u, T v, const std::vector<T>& slots = {}) {
    return interpret(f, FieldInputs<T>{u, v}, slots);
}

// `slots` are the values of `f.scalars.symbols`, in that order.
template<Scalar T>
Vec<T, 3> interpret(const PodVecField<T>& f, const MotionEnv<T>& env, const std::vector<T>& slots = {}) {
    std::vector<Vec<T, 3>> s(f.ops.size());
    std::vector<T> scratch(f.scalars.ops.size());

    // A factor is its own field evaluated at (t, 0) -- the convention
    // `scaled` and `rotated` fix -- run from its own start in the pool.
    const FieldInputs<T> in = VecField<T>::inputs_of(env);
    auto factor = [&](const PodRange& r) {
        return interpret(f.scalars.ops.data() + r.begin, r.count, f.scalars.tables.data(),
                         f.scalars.points.data(), scratch.data(), in, slots.data(), slots.size(), f.scalars.weights.data());
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
            case VecOp::Make: {
                const T x = factor(n.factor[0]);
                const T y = factor(n.factor[1]);
                const T z = factor(n.factor[2]);
                s[i] = Vec<T, 3>{x, y, z};
                break;
            }
            case VecOp::Opaque: s[i] = Vec<T, 3>{}; break;   // see the scalar case
        }
    }
    return s[f.ops.size() - 1];
}

} // namespace spatium::io::build
