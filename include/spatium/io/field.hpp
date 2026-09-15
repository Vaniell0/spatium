#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <cassert>
#  include <algorithm>
#  include <cstddef>
#  include <cstdint>
#  include <functional>
#  include <string>
#  include <type_traits>
#  include <typeindex>
#  include <utility>
#  include <vector>
#endif

// A Field is an expression over (u, v) in which an opaque callable is a
// *leaf*, not an alternative to the expression.
//
// The distinction is the whole design. Under "a field is either a tree or
// a callable", a user picks a representation where they declare the field
// and silently forfeits caching, lowering and hoisting for the entire
// field. Under "a tree whose leaves may be opaque", composing a
// structural field with an opaque one keeps the outer structure and loses
// exactly one leaf. "Is this lowerable" is then *derived* by looking for
// opaque leaves, never declared and never kept in sync by hand.
//
// ── Representation: a flat pool in topological order ─────────────
//
// A Field is logically a tree. It is represented as a vector of tagged
// ops with children stored before their parents, the root last.
//
// The cheap structural hash and the cheap deduplication are consequences,
// not the reason. The reason is that a topologically ordered array
// evaluates in ONE LINEAR PASS: a loop over indices where each op reads
// children already computed. No recursion, no explicit stack, no
// branching on depth. That is the shape a GPU can run, and the shape a
// linked tree cannot provide without first being flattened -- so building
// the linked form would mean writing the flattening anyway.
//
// **APPEND-ONLY. Inserting into the middle is forbidden.** Topological
// order is an obligation this type owes its own evaluator, not a property
// it gets for free. It holds today for a structural reason: a node can
// only name indices that already exist, and larger indices do not exist
// yet. An operation that inserted in the middle would renumber nothing
// and break it *silently* -- the first evaluation would read an
// uncomputed child and return whatever was in the scratch slot. Every
// builder below therefore appends, and the invariant is asserted per
// appended node in Debug rather than trusted -- O(1) where it can break,
// not an O(n) re-derivation afterwards that it still holds.
//
// **A Field owns its pool. This is a value type.** Copying copies the
// ops; `+=` appends into this field's own vector and cannot be observed
// by any other field, the same way `std::string::operator+=` cannot.
// `c = a + b` leaves `a` and `b` untouched, which is why it must copy and
// is therefore O(n) -- and why `+=` exists as the linear way to build a
// field term by term.
//
// The price of that, stated so nobody has to discover it: **copying a
// Field copies its pool.** `f1 = f2` is O(ops), like `std::string`
// without COW, not a pointer bump. Fields are small (single digits of
// ops in every use so far) so this is cheap, but it is linear rather than
// free, and passing one to something that stores a
// `std::function<T(T,T)>` -- `offset_surface`, for instance -- makes a
// copy there too. The node keeps the Field, so structure stays visible to
// the report and to lowering; it is the evaluation path that wraps.
//
// Saying that precisely matters, because the obvious phrase for what the
// builders do -- "append, never mutate" -- is borrowed from a *shared*
// pool design (ROADMAP's, where a Field is a range into a pool other
// fields also index) and would be describing a guarantee this code does
// not need and does not give. Here append-only protects exactly one
// thing: **topological order**, because a node can only name indices that
// already exist. It protects nothing about aliasing, since there is no
// aliasing to protect. If Field ever becomes (shared_ptr<Pool>, range) --
// which is what substitution f(g) across many fields would want -- then
// the aliasing guarantee has to be added *and* proved, not inherited from
// this sentence.
//
// **Known boundary: this leaks under interactive editing.** Append-never-
// mutate means superseded ops stay. For a trace built once and rendered
// many times that is free; for a live-reloaded scene it grows without
// bound, and this pool is the first thing that would have to change.
// Named here so it is not rediscovered as a bug.
//
// **RTTI is required**, and not newly so: `Trace` already stores the
// exact analytic form in a `std::any`, which cannot exist without it.
// `typeid(F)` on an opaque leaf adds no requirement that the DSL did not
// already have. Spatium does not build with `-fno-rtti` and does not
// pretend to.

SPATIUM_EXPORT namespace spatium::io::build {

enum class Op : std::uint8_t {
    Const,   // a literal; `value`
    U,       // the first parameter
    V,       // the second parameter
    Add, Sub, Mul, Div,
    Opaque,  // a callable leaf -- the escape hatch, kept first-class
};

inline const char* op_name(Op o) {
    switch (o) {
        case Op::Const:  return "Const";
        case Op::U:      return "U";
        case Op::V:      return "V";
        case Op::Add:    return "Add";
        case Op::Sub:    return "Sub";
        case Op::Mul:    return "Mul";
        case Op::Div:    return "Div";
        case Op::Opaque: return "Opaque";
    }
    return "?";
}

inline bool is_binary(Op o) {
    return o == Op::Add || o == Op::Sub || o == Op::Mul || o == Op::Div;
}

template<Scalar T = double>
struct FieldOp {
    Op op = Op::Const;

    T value{};                          // Const
    std::uint32_t a = 0, b = 0;         // child indices, both < this index

    // Opaque leaf. The callable itself, plus what can be known about it
    // without asking the user for a name.
    std::function<T(T, T)> fn;

    // `typeid` of the closure type. Compiler-generated, so it costs the
    // user nothing and opens no door to registering function *names*:
    // identical for every instance of one lambda (including 19 800
    // instances of a lambda written once inside a loop), different for
    // two lambdas written separately. std::move_only_function has no
    // target_type(), which is why identity is captured here at
    // construction rather than asked of the erased callable later.
    std::type_index type = std::type_index(typeid(void));

    // sizeof the closure, so a heavy capture is visible rather than
    // mysterious. A lambda capturing a 512-byte PerlinNoise *by value*
    // costs 9.67 MB across 19 800 fields; one capturing it by reference
    // costs 8 bytes. The report prints the total so the difference shows
    // up as a number instead of as memory pressure nobody can source.
    std::size_t payload = 0;

    // Stateless closures are interchangeable between instances, which is
    // what makes deduplication by type alone sound for them and unsound
    // for capturing ones. Decided here, not asked of the user.
    bool stateless = false;
};

template<Scalar T = double>
class Field {
public:
    Field() { push_const(T{0}); }
    Field(T constant) { push_const(constant); }   // NOLINT: implicit on purpose

    static Field u() { Field f; f.ops_.clear(); f.push(Op::U); return f; }
    static Field v() { Field f; f.ops_.clear(); f.push(Op::V); return f; }

    // The escape hatch, and deliberately as ordinary to write as the
    // structural builders. An IR whose hatch is second-class becomes a
    // language you have to learn.
    template<typename F>
        requires std::is_invocable_r_v<T, const F&, T, T>
    static Field opaque(F fn) {
        Field f;
        f.ops_.clear();
        FieldOp<T> n{};
        n.op        = Op::Opaque;
        n.type      = std::type_index(typeid(F));
        n.payload   = sizeof(F);
        n.stateless = std::is_empty_v<F>;
        n.fn        = std::move(fn);
        f.ops_.push_back(std::move(n));
        f.structural_ = false;
        return f;
    }

    // Any callable becomes an opaque leaf. Implicit on purpose: every
    // place that used to take a std::function still reads the same, and
    // what makes the opacity visible is the report, not a spelling the
    // user has to remember. `Field::opaque` stays for when saying it
    // explicitly is clearer.
    template<typename F>
        requires std::is_invocable_r_v<T, const F&, T, T> &&
                 (!std::is_same_v<std::remove_cvref_t<F>, Field>) &&
                 (!std::is_convertible_v<F, T>)
    Field(F fn) : Field(opaque(std::move(fn))) {}   // NOLINT: implicit on purpose

    // ── Evaluation: one pass, children already computed ──────────

    T operator()(T u, T v) const {
        constexpr std::size_t kInline = 32;
        if (ops_.size() <= kInline) {
            T scratch[kInline];
            return eval_into(scratch, u, v);
        }
        std::vector<T> scratch(ops_.size());
        return eval_into(scratch.data(), u, v);
    }

    // Derived by construction, never declared: a field is structural iff
    // no op in it is opaque. Cached rather than walked, because this is a
    // method someone will call inside a loop over vertices, where an O(ops)
    // walk per vertex is a different cost from an O(fields) report per
    // frame.
    bool is_structural() const { return structural_; }

    std::size_t size() const { return ops_.size(); }
    const FieldOp<T>& op(std::size_t i) const { return ops_[i]; }
    const std::vector<FieldOp<T>>& ops() const { return ops_; }

    // ── Composition: append, never mutate ────────────────────────

    // Taken by value: an rvalue operand is moved in and appended to, an
    // lvalue is copied because the original has to survive the
    // expression. That copy is not a defect to optimise away -- `c = a + b`
    // leaving `a` unchanged is what the operator means -- but it does make
    // `acc = acc + term` quadratic in a loop, which is why the
    // accumulating form below exists and is the one to reach for.
    friend Field operator+(Field x, const Field& y) { x.append(Op::Add, y); return x; }
    friend Field operator-(Field x, const Field& y) { x.append(Op::Sub, y); return x; }
    friend Field operator*(Field x, const Field& y) { x.append(Op::Mul, y); return x; }
    friend Field operator/(Field x, const Field& y) { x.append(Op::Div, y); return x; }

    // The accumulating form. Appends into this field's own pool, so
    // building a field term by term is linear rather than quadratic:
    // nothing else can observe this pool mid-expression, so appending to
    // it is not mutation of anything anyone holds.
    Field& operator+=(const Field& y) { append(Op::Add, y); return *this; }
    Field& operator-=(const Field& y) { append(Op::Sub, y); return *this; }
    Field& operator*=(const Field& y) { append(Op::Mul, y); return *this; }
    Field& operator/=(const Field& y) { append(Op::Div, y); return *this; }

    // Kept public so a test can assert the invariant independently of the
    // per-append assert that maintains it -- a check nothing executes
    // reads exactly like a check that passes.
    bool topologically_ordered() const {
        for (std::size_t i = 0; i < ops_.size(); ++i) {
            const auto& n = ops_[i];
            if (!is_binary(n.op)) continue;
            if (n.a >= i || n.b >= i) return false;
        }
        return true;
    }
private:
    std::vector<FieldOp<T>> ops_;
    bool structural_ = true;

    void push_const(T value) {
        FieldOp<T> n{};
        n.op = Op::Const;
        n.value = value;
        ops_.push_back(std::move(n));
    }

    void push(Op o) {
        FieldOp<T> n{};
        n.op = o;
        ops_.push_back(std::move(n));
    }

    T eval_into(T* s, T u, T v) const {
        for (std::size_t i = 0; i < ops_.size(); ++i) {
            const auto& n = ops_[i];
            switch (n.op) {
                case Op::Const:  s[i] = n.value; break;
                case Op::U:      s[i] = u; break;
                case Op::V:      s[i] = v; break;
                case Op::Add:    s[i] = s[n.a] + s[n.b]; break;
                case Op::Sub:    s[i] = s[n.a] - s[n.b]; break;
                case Op::Mul:    s[i] = s[n.a] * s[n.b]; break;
                case Op::Div:    s[i] = s[n.a] / s[n.b]; break;
                case Op::Opaque: s[i] = n.fn(u, v); break;
            }
        }
        return s[ops_.size() - 1];
    }

    // Append y's ops, with their child indices shifted, then the parent
    // naming both roots. This is the whole of composition: the parent can
    // only name indices that already exist, so topological order holds by
    // construction rather than being repaired afterwards.
    //
    // The order check is per node and O(1), not a walk of the array. That
    // is both cheaper and a stronger statement: it verifies that the
    // operation which is supposed to preserve the invariant actually does,
    // at the moment it could break it, rather than re-deriving afterwards
    // that the invariant still holds. The O(ops) version it replaced cost
    // 21% of the build time of a 4 000-term field -- measured, not guessed,
    // and the other 79% was the copying that `operator+=` now avoids.
    void append(Op o, const Field& y) {
        const auto shift  = static_cast<std::uint32_t>(ops_.size());
        const auto x_root = shift - 1;

        // Reserve only past the current capacity, and then geometrically.
        // A bare `reserve(size() + k)` on every append reallocates every
        // time -- it asks for exactly what is needed, which defeats the
        // vector's own doubling and turns a linear accumulation into a
        // quadratic one. Measured: the accumulating form stayed quadratic
        // (12.5 s for 19 800 terms) until this line was written this way.
        const std::size_t need = ops_.size() + y.ops_.size() + 1;
        if (need > ops_.capacity())
            ops_.reserve(std::max(need, ops_.capacity() * 2));

        for (const auto& n : y.ops_) {
            FieldOp<T> c = n;
            if (is_binary(c.op)) { c.a += shift; c.b += shift; }
            assert((!is_binary(c.op) ||
                    (c.a < ops_.size() && c.b < ops_.size())) &&
                   "Field: a child must precede its parent");
            ops_.push_back(std::move(c));
        }

        FieldOp<T> parent{};
        parent.op = o;
        parent.a  = x_root;
        parent.b  = static_cast<std::uint32_t>(ops_.size() - 1);
        assert(parent.a < ops_.size() && parent.b < ops_.size() &&
               "Field: a child must precede its parent");
        ops_.push_back(std::move(parent));

        structural_ = structural_ && y.structural_;
    }

};

// ── What a field report has to say ───────────────────────────────
//
// Three rules, each of which exists because of a specific way this went
// wrong or could go wrong.
//
// 1. **Separate the units.** An earlier version of this struct counted
//    fields and leaves in the same tally, so `fields=5, recognized=4,
//    structural=1` read as three alternatives that happened to add up.
//    They are not alternatives and they do not add up in general: a
//    field can be non-structural and still hold several opaque leaves,
//    and then nothing sums. Fields and leaves are now separate groups,
//    each with an identity that must hold, and `check()` asserts both.
//
// 2. **Count `unknown` separately and print it even when zero.** A tally
//    reading "39 600 fields, 1 distinct type" cannot be told apart from
//    "all 39 600 fell into one unknown bucket" -- the same shape of
//    defect as a NaN root reporting itself as real, which
//    `conventions.md` names. Zero unknowns is information; a missing
//    line is not.
//
// 3. **Give `payload_bytes` its comparison.** A number with no threshold
//    is a number with no answer: the reader cannot tell whether 1 544 B
//    is fine and 9.67 MB is not. The threshold is not arbitrary and not
//    a size -- it is the cache. `benchmarks/bench_trace.cpp:84` measured
//    19 800 closures each capturing a 512-byte PerlinNoise by value as
//    "~9.7 MB of identical tables against a 12 MB L3", costing ~3% of
//    frame time. So the comparison to make is against L3, and
//    `payload_verdict()` below states it rather than leaving the reader
//    to know it.

struct FieldStats {
    // Fields examined. structural_fields + opaque_fields == fields.
    std::size_t fields = 0;
    std::size_t structural_fields = 0;  // no opaque leaf anywhere in it
    std::size_t opaque_fields = 0;      // at least one opaque leaf

    // Leaves found. recognized_leaves + unknown_leaves == opaque_leaves.
    std::size_t opaque_leaves = 0;
    std::size_t recognized_leaves = 0;  // closure type known
    std::size_t unknown_leaves = 0;     // could not be classified

    std::size_t distinct_types = 0;
    std::size_t payload_bytes = 0;      // captured state in opaque leaves

    // The two identities, checkable. They are the reason the groups are
    // split: if either fails, the report is describing something other
    // than what it says.
    bool consistent() const {
        return structural_fields + opaque_fields == fields &&
               recognized_leaves + unknown_leaves == opaque_leaves;
    }
};

// Typical L3 on the machine bench_trace.cpp was measured on. Not a
// portable fact, which is why it is named rather than hidden in a
// comparison -- a reader on different hardware should substitute theirs.
inline constexpr std::size_t kReferenceL3Bytes = 12u * 1024 * 1024;

// The comparison a bare byte count cannot make for itself.
inline const char* payload_verdict(std::size_t payload_bytes) {
    if (payload_bytes * 4 >= kReferenceL3Bytes)
        return "comparable to L3 -- captures are evicting the scene; capture by reference or share";
    if (payload_bytes * 32 >= kReferenceL3Bytes)
        return "a visible fraction of L3 -- worth looking at what is captured by value";
    return "small against L3";
}

template<Scalar T>
void accumulate(FieldStats& stats, const Field<T>& f,
                std::vector<std::type_index>& seen) {
    ++stats.fields;
    if (f.is_structural()) ++stats.structural_fields;
    else                   ++stats.opaque_fields;

    for (std::size_t i = 0; i < f.size(); ++i) {
        const auto& n = f.op(i);
        if (n.op != Op::Opaque) continue;
        ++stats.opaque_leaves;
        if (n.type == std::type_index(typeid(void))) {
            ++stats.unknown_leaves;
            continue;
        }
        ++stats.recognized_leaves;
        stats.payload_bytes += n.payload;
        if (std::find(seen.begin(), seen.end(), n.type) == seen.end())
            seen.push_back(n.type);
    }
    stats.distinct_types = seen.size();
}

} // namespace spatium::io::build
