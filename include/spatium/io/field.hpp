#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <spatium/algebra/matrix.hpp>
#  include <spatium/algebra/groups/so3.hpp>
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

// ── A field over points, and the environment it reads ────────────
//
// The same design as Field one dimension up: a flat pool, children before
// parents, an opaque callable as a leaf rather than as an alternative to
// the expression.
//
// What differs is the input, and that difference is the point. A motion
// today is a function of a point and a time, and writing `(p, t)` into
// every signature would mean that giving an object its own clock, or
// letting a solver write its pose, changes every signature and every call
// site in every scene anyone has written. So the input is a **named
// environment**. Today it carries exactly what it carried before and
// behaves identically; tomorrow a field is added to this struct and
// nothing else moves.
//
// That is the cheapest possible hedge against questions that are open on
// purpose -- where mutable state lives, whether a timer is a function or
// a state, whether physics bodies live in the trace or beside it. None of
// them has a consumer yet, so none is being answered here. This costs
// nothing and stops the answers from being expensive.
template<Scalar T = double>
struct MotionEnv {
    Vec<T, 3> p{};   // the point being moved
    T t{};           // scene time

    // Where this *instance* sits at rest -- one value per object, not one
    // per vertex, and that difference is the whole reason it can be here
    // at all.
    //
    // A motion that tells instances apart has to read something that
    // differs between them. Before this the only such thing was `p`, and
    // reading `p` makes a motion a deformation: the affine test refuses
    // it, correctly, because a term that varies per vertex cannot be a
    // shared transform. So a `Scatter` had one motion for all its
    // instances and they all flew the same way.
    //
    // An origin varies per object and is constant across an object's
    // vertices, so a field reading it stays affine in the point and stays
    // instanceable. The particle's trajectory comes from where it
    // started; its shape does not move.
    //
    // This is also what a plain radial burst needs, which is less obvious.
    // A placement is `b + R*s*p` applied to every vertex, and a scattered
    // vertex is `seat + F*local` -- so one scalar `s` scales the distance
    // flown *and* the speck, and expanding the cloud inflates every fleck
    // with it. `seat` and `local` arrive already summed into `p`;
    // separating them is exactly this field.
    Vec<T, 3> origin{};
};

enum class VecOp : std::uint8_t {
    Point,     // the environment's p
    Const,     // a literal vector
    Add, Sub,
    Scale,     // child a, times a scalar read from the environment
    Rotate,    // child a, turned by a rotation read from the environment
    Opaque,    // a callable leaf; `reads_point` says whether it uses p
};

inline const char* vec_op_name(VecOp o) {
    switch (o) {
        case VecOp::Point:  return "Point";
        case VecOp::Const:  return "Const";
        case VecOp::Add:    return "Add";
        case VecOp::Sub:    return "Sub";
        case VecOp::Scale:  return "Scale";
        case VecOp::Rotate: return "Rotate";
        case VecOp::Opaque: return "Opaque";
    }
    return "?";
}

inline bool is_binary(VecOp o) { return o == VecOp::Add || o == VecOp::Sub; }

// One child, named once rather than spelled out at each of the four
// places that shift indices, test affinity or walk a path. Adding a
// second unary op is then one entry here instead of four edits that must
// agree.
inline bool is_unary(VecOp o) { return o == VecOp::Scale || o == VecOp::Rotate; }

template<Scalar T = double>
struct VecFieldOp {
    VecOp op = VecOp::Point;
    Vec<T, 3> value{};
    std::uint32_t a = 0, b = 0;

    // move_only_function, not std::function, and deliberately not the
    // same choice the scalar field makes. A motion hook is allowed to own
    // move-only state -- tests/test_build_dsl.cpp has one owning a
    // unique_ptr -- and that capability was added on purpose in PR #26.
    // A thickness cannot make the same choice, because it flows into
    // ParametricSurface's ParamFn, which is a std::function.
    //
    // This costs copyability: a VecField is move-only, so a TraceNode and
    // a Trace are too, so cook() consumes a trace rather than copying it.
    // That is the status quo and it is fine. What it does NOT cost is the
    // report: identity is captured as typeid(F) at construction, below, so
    // a motion is classifiable whether or not the erasure that holds it
    // can be asked. The report's blind spot was never move_only_function's
    // missing target_type() -- it was that the type did not belong to us.
    std::move_only_function<Vec<T, 3>(const MotionEnv<T>&) const> fn;

    // The scalar a Scale multiplies by. Read from the environment, never
    // from the point -- a Scale whose factor depended on the point would
    // be a deformation again.
    std::move_only_function<T(const MotionEnv<T>&) const> scale_fn;

    // The rotation a Rotate turns by, under exactly the same rule: read
    // from the environment, never from the point. A rotation that varied
    // per vertex would bend the object rather than orient it.
    //
    // Stored as the matrix rather than as an axis-angle vector, because
    // this is the form both consumers want -- evaluation multiplies a
    // point by it, and a placement hands it to a ray test that needs the
    // transpose. Callers who think in axis-angle get the `rotated()`
    // overload that exponentiates through SO3, so nothing is lost at the
    // authoring end and no exp() runs per vertex.
    std::move_only_function<Matrix<T, 3, 3>(const MotionEnv<T>&) const> rot_fn;

    // **The distinction that decides whether a motion can be instanced**,
    // and it is not "structural versus opaque". A leaf may be as opaque as
    // it likes -- the donut's particle position is Perlin noise and always
    // will be -- as long as it does not read the *point*.
    //
    // Cost is per vertex. A term that ignores p is evaluated once per
    // object per frame; a term that reads p is evaluated once per vertex.
    // So a motion of the form A(t) + p * s(t), with A and s opaque but
    // p-blind, still places a shared geometry with a transform: A and s
    // are 19 800 evaluations, not 19 800 x 8. Requiring A and s to be
    // *structural* would have meant making PerlinNoise an expression,
    // which is both enormous and unnecessary.
    //
    // Defaults to true, because assuming a leaf ignores p when it does not
    // would silently collapse a deformation onto one shape. Saying it
    // ignores p is a promise the author makes, via `opaque_of_time`.
    bool reads_point = true;
    std::type_index type = std::type_index(typeid(void));
    std::size_t payload = 0;
    bool stateless = false;
};

template<Scalar T = double>
class VecField {
public:
    // Default is the identity motion, structurally: not an opaque lambda
    // that happens to return its argument, but a Point node a reader can
    // see is the identity.
    VecField() { push(VecOp::Point); }

    static VecField point() { return VecField{}; }

    static VecField constant(const Vec<T, 3>& v) {
        VecField f;
        f.ops_.clear();
        VecFieldOp<T> n{};
        n.op = VecOp::Const;
        n.value = v;
        f.ops_.push_back(std::move(n));
        return f;
    }

    // A pure translation, structurally -- the case that matters most,
    // because a translation is an isometry and so a node carrying one
    // could in principle keep its exact analytic form instead of dropping
    // it. Recognising that is what an opaque callable can never support
    // and this representation can; the check itself is the next step, not
    // this one.
    static VecField translation(const Vec<T, 3>& by) {
        return point() + constant(by);
    }

    // An opaque leaf reading where this instance started, and the time.
    //
    // Same promise as `opaque_of_time` and kept the same way: the callable
    // is handed an origin and a `t` and has no route to the point, so the
    // *signature* enforces what a comment could only ask for. That is what
    // lets `reads_point` stay false and the motion stay a placement while
    // every instance still flies its own way.
    template<typename F>
        requires std::is_invocable_r_v<Vec<T, 3>, const F&, const Vec<T, 3>&, T>
    static VecField opaque_per_instance(F fn) {
        VecField f;
        f.ops_.clear();
        VecFieldOp<T> n{};
        n.op          = VecOp::Opaque;
        n.type        = std::type_index(typeid(F));
        n.payload     = sizeof(F);
        n.stateless   = std::is_empty_v<F>;
        n.reads_point = false;
        n.fn = [fn = std::move(fn)](const MotionEnv<T>& e) { return fn(e.origin, e.t); };
        f.ops_.push_back(std::move(n));
        f.structural_ = false;
        return f;
    }

    // An opaque leaf that promises not to read the point. The promise is
    // the whole content of this factory -- the callable takes an
    // environment and cannot reach a point through it, so the signature
    // enforces what the comment on `reads_point` describes.
    template<typename F>
        requires std::is_invocable_r_v<Vec<T, 3>, const F&, T>
    static VecField opaque_of_time(F fn) {
        VecField f;
        f.ops_.clear();
        VecFieldOp<T> n{};
        n.op          = VecOp::Opaque;
        n.type        = std::type_index(typeid(F));
        n.payload     = sizeof(F);
        n.stateless   = std::is_empty_v<F>;
        n.reads_point = false;
        n.fn = [fn = std::move(fn)](const MotionEnv<T>& e) { return fn(e.t); };
        f.ops_.push_back(std::move(n));
        f.structural_ = false;
        return f;
    }

    // `field * s(t)`. The factor reads time, never the point: a factor
    // that varied per vertex would make this a deformation again, which
    // is the thing the whole distinction exists to keep out.
    //
    // Consumes the field rather than copying it, because a pool of
    // move-only leaves cannot be copied -- the same reason `operator+`
    // takes its left operand by value.
    // A factor that varies per instance: same rule as the motion itself,
    // reading the origin rather than the point, so a scale can differ
    // between instances without becoming a deformation.
    template<typename S>
        requires (!std::is_invocable_r_v<T, const S&, T>) &&
                 std::is_invocable_r_v<T, const S&, const Vec<T, 3>&, T>
    friend VecField scaled(VecField x, S s) {
        VecFieldOp<T> n{};
        n.op       = VecOp::Scale;
        n.a        = static_cast<std::uint32_t>(x.ops_.size() - 1);
        n.scale_fn = [s = std::move(s)](const MotionEnv<T>& e) { return s(e.origin, e.t); };
        assert(n.a < x.ops_.size() && "VecField: a child must precede its parent");
        x.ops_.push_back(std::move(n));
        return x;
    }

    template<typename S>
        requires std::is_invocable_r_v<T, const S&, T>
    friend VecField scaled(VecField x, S s) {
        VecFieldOp<T> n{};
        n.op       = VecOp::Scale;
        n.a        = static_cast<std::uint32_t>(x.ops_.size() - 1);
        n.scale_fn = [s = std::move(s)](const MotionEnv<T>& e) { return s(e.t); };
        assert(n.a < x.ops_.size() && "VecField: a child must precede its parent");
        x.ops_.push_back(std::move(n));
        return x;
    }

    // `R(t) * field`. The companion to `scaled`, and the op that turns a
    // placement from "somewhere else" into "somewhere else, facing
    // somewhere else" -- without it every instance of a shared geometry
    // is not merely in the same pose but in the *same orientation*, which
    // is invisible on a sphere and impossible to miss on anything with a
    // side.
    //
    // Two overloads, because the two ways of saying a rotation are both
    // natural and neither should have to be converted by hand. This one
    // takes the matrix; the next takes an axis-angle vector and
    // exponentiates it through SO3, which is what a constant orientation
    // (`[o](T) { return o; }`) and a spin (`[w](T t) { return w * t; }`)
    // both want to write.
    // Per-instance turns, for the same reason `scaled` has them: a
    // scatter has one motion field, so anything that must differ between
    // its instances has to read the origin. Without this every flake in a
    // cloud tumbles in lockstep.
    template<typename R>
        requires (!std::is_invocable_r_v<Matrix<T, 3, 3>, const R&, T>) &&
                 (!std::is_invocable_r_v<Vec<T, 3>, const R&, T>) &&
                 std::is_invocable_r_v<Vec<T, 3>, const R&, const Vec<T, 3>&, T>
    friend VecField rotated(VecField x, R r) {
        return rotated(std::move(x), [r = std::move(r)](const Vec<T, 3>& o, T t) {
            return algebra::SO3<T>{}.exp(r(o, t));
        });
    }

    template<typename R>
        requires (!std::is_invocable_r_v<Matrix<T, 3, 3>, const R&, T>) &&
                 std::is_invocable_r_v<Matrix<T, 3, 3>, const R&, const Vec<T, 3>&, T>
    friend VecField rotated(VecField x, R r) {
        VecFieldOp<T> n{};
        n.op     = VecOp::Rotate;
        n.a      = static_cast<std::uint32_t>(x.ops_.size() - 1);
        n.rot_fn = [r = std::move(r)](const MotionEnv<T>& e) { return r(e.origin, e.t); };
        assert(n.a < x.ops_.size() && "VecField: a child must precede its parent");
        x.ops_.push_back(std::move(n));
        return x;
    }

    template<typename R>
        requires std::is_invocable_r_v<Matrix<T, 3, 3>, const R&, T>
    friend VecField rotated(VecField x, R r) {
        VecFieldOp<T> n{};
        n.op     = VecOp::Rotate;
        n.a      = static_cast<std::uint32_t>(x.ops_.size() - 1);
        n.rot_fn = [r = std::move(r)](const MotionEnv<T>& e) { return r(e.t); };
        assert(n.a < x.ops_.size() && "VecField: a child must precede its parent");
        x.ops_.push_back(std::move(n));
        return x;
    }

    template<typename R>
        requires (!std::is_invocable_r_v<Matrix<T, 3, 3>, const R&, T>) &&
                 std::is_invocable_r_v<Vec<T, 3>, const R&, T>
    friend VecField rotated(VecField x, R r) {
        return rotated(std::move(x), [r = std::move(r)](T t) {
            return algebra::SO3<T>{}.exp(r(t));
        });
    }

    // A callable leaf. Accepts either shape: a function of the whole
    // environment, or the historical (p, t) -- the latter wrapped, so
    // every scene written against the old signature still compiles while
    // new code can read whatever the environment grows.
    template<typename F>
        requires std::is_invocable_r_v<Vec<T, 3>, const F&, const MotionEnv<T>&>
    static VecField opaque(F fn) {
        return make_opaque<F>(std::move(fn));
    }

    template<typename F>
        requires (!std::is_invocable_r_v<Vec<T, 3>, const F&, const MotionEnv<T>&>) &&
                 std::is_invocable_r_v<Vec<T, 3>, const F&, const Vec<T, 3>&, T>
    static VecField opaque(F fn) {
        return make_opaque<F>(
            [fn = std::move(fn)](const MotionEnv<T>& e) { return fn(e.p, e.t); },
            sizeof(F), std::type_index(typeid(F)), std::is_empty_v<F>);
    }

    template<typename F>
        requires (!std::is_same_v<std::remove_cvref_t<F>, VecField>) &&
                 (std::is_invocable_r_v<Vec<T, 3>, const F&, const MotionEnv<T>&> ||
                  std::is_invocable_r_v<Vec<T, 3>, const F&, const Vec<T, 3>&, T>)
    VecField(F fn) : VecField(opaque(std::move(fn))) {}   // NOLINT: implicit on purpose

    Vec<T, 3> operator()(const MotionEnv<T>& env) const {
        constexpr std::size_t kInline = 32;
        if (ops_.size() <= kInline) {
            Vec<T, 3> scratch[kInline];
            return eval_into(scratch, env);
        }
        std::vector<Vec<T, 3>> scratch(ops_.size());
        return eval_into(scratch.data(), env);
    }

    // The historical call shape, kept so call sites read unchanged.
    Vec<T, 3> operator()(const Vec<T, 3>& p, T t) const {
        return (*this)(MotionEnv<T>{p, t});
    }

    bool is_structural() const { return structural_; }
    bool is_identity() const {
        return ops_.size() == 1 && ops_[0].op == VecOp::Point;
    }

    // ── Placement: what makes a motion instanceable ──────────────
    //
    // Does this op's value depend on the point being moved? Derived by
    // walking, not declared -- the same call as `is_structural()`, and for
    // the same reason: a flag kept in sync by hand goes stale.
    bool reads_point(std::size_t i) const {
        const auto& n = ops_[i];
        switch (n.op) {
            case VecOp::Point:  return true;
            case VecOp::Const:  return false;
            case VecOp::Opaque: return n.reads_point;
            case VecOp::Scale:
            case VecOp::Rotate: return reads_point(n.a);
            case VecOp::Add:
            case VecOp::Sub:    return reads_point(n.a) || reads_point(n.b);
        }
        return true;
    }

    // A *placement* is a motion that is affine in the point: it moves the
    // whole object without changing its shape, so N instances can share
    // one geometry and differ only by a transform. Anything else is a
    // deformation, and instances that deform differently have nothing to
    // share.
    //
    // The test is structural and cheap: exactly one path through the
    // expression reaches the Point, and along it the Point is either bare
    // or under Scale factors -- never under an Opaque that reads it, and
    // never on both sides of one operator. Everything else in the tree may
    // be as opaque as it likes, because it does not touch the point and so
    // is evaluated once per object rather than once per vertex.
    bool is_placement() const { return affine_in_point(ops_.size() - 1); }

    // The transform this placement is, at one moment: `p -> R*s*p + b`.
    //
    // Translation is the whole motion evaluated with the point at the
    // origin, and that definition needed no change when rotation arrived,
    // which is worth saying because the obvious worry is the opposite. A
    // translation introduced *below* a rotation on the path is itself
    // rotated, and evaluating at p = 0 already reports where the origin
    // ends up -- through every rotation on the way out. Only the linear
    // part had to learn anything new.
    //
    // Rotation and scale are the product of the Rotate and Scale factors
    // along the path to the Point. Order matters for the matrices and does
    // not for the scalars, and the two separate cleanly for exactly that
    // reason: R1*s1*R2*s2 == (R1*R2)*(s1*s2), because a scalar commutes
    // with everything. So the pair stays a pair rather than collapsing
    // into one general 3x3 -- and that is not a saving, it is the
    // contract: an orthonormal R and a uniform s are what let a ray test
    // transform into local space without rescaling `t` or fixing up a
    // normal. A general linear map would lose both.
    struct Placement {
        Vec<T, 3> translation{};
        Matrix<T, 3, 3> rotation = Matrix<T, 3, 3>::identity();
        T scale = T{1};
    };

    Placement placement_at(const MotionEnv<T>& env) const {
        assert(is_placement() && "placement_at: this motion is a deformation");
        MotionEnv<T> at_origin = env;
        at_origin.p = Vec<T, 3>{};
        auto [rotation, scale] = linear_on_path(ops_.size() - 1, env);
        return Placement{(*this)(at_origin), rotation, scale};
    }

    // "Is there a motion here" -- what `if (node.transform)` meant when
    // the slot was an erased function that could be empty. A VecField is
    // never empty: the default *is* the identity, and saying so
    // structurally is the improvement. So the question becomes whether it
    // is anything other than the identity, and it is answered by looking
    // at one op rather than by a null check.
    explicit operator bool() const { return !is_identity(); }

    std::size_t size() const { return ops_.size(); }
    const VecFieldOp<T>& op(std::size_t i) const { return ops_[i]; }

    // Both operands are consumed. A pool of move-only leaves cannot be
    // copied, so `c = a + b` leaving `a` and `b` usable is not available
    // here -- unlike the scalar Field, whose leaves are copyable. This is
    // the price of letting a motion hook own move-only state, and it is
    // stated rather than discovered: composition takes rvalues.
    friend VecField operator+(VecField x, VecField y) { x.append(VecOp::Add, y); return x; }
    friend VecField operator-(VecField x, VecField y) { x.append(VecOp::Sub, y); return x; }
    VecField& operator+=(VecField y) { append(VecOp::Add, y); return *this; }
    VecField& operator-=(VecField y) { append(VecOp::Sub, y); return *this; }

    bool topologically_ordered() const {
        for (std::size_t i = 0; i < ops_.size(); ++i) {
            if (!is_binary(ops_[i].op)) continue;
            if (ops_[i].a >= i || ops_[i].b >= i) return false;
        }
        return true;
    }

private:
    std::vector<VecFieldOp<T>> ops_;
    bool structural_ = true;

    // Affine in the point: the Point is reached along exactly one branch,
    // through Add/Sub/Scale only. An Opaque that reads the point fails
    // here even though it might in truth be affine -- we cannot ask a
    // closure what it does, which is the whole reason the structural form
    // exists.
    bool affine_in_point(std::size_t i) const {
        const auto& n = ops_[i];
        switch (n.op) {
            case VecOp::Point:  return true;
            case VecOp::Const:  return true;
            case VecOp::Opaque: return !n.reads_point;
            case VecOp::Scale:
            case VecOp::Rotate: return affine_in_point(n.a);
            case VecOp::Add:
            case VecOp::Sub:
                // Both sides affine, and at most one of them touching the
                // point: p + p would be affine arithmetically but is not a
                // rigid placement of one geometry.
                if (!affine_in_point(n.a) || !affine_in_point(n.b)) return false;
                return !(reads_point(n.a) && reads_point(n.b));
        }
        return false;
    }

    // The linear part of the placement, walked along the one path that
    // reaches the Point. Returns the pair rather than a single matrix so
    // that the orthonormal factor and the uniform factor stay
    // distinguishable at the far end -- see `Placement` for why that
    // distinction is the contract and not an optimization.
    //
    // A `Sub` whose point-bearing side is the right one contributes -1 to
    // the *scale*, not to the rotation: a negated uniform scale is a
    // reflection and has no place in SO(3), while as a scale it is the
    // behaviour that was already there and already correct -- a ray test
    // dividing origin and direction by the same negative number still
    // describes the same points.
    struct Linear {
        Matrix<T, 3, 3> rotation = Matrix<T, 3, 3>::identity();
        T scale = T{1};
    };

    Linear linear_on_path(std::size_t i, const MotionEnv<T>& env) const {
        const auto& n = ops_[i];
        switch (n.op) {
            case VecOp::Point:  return Linear{};
            case VecOp::Scale: {
                auto below = linear_on_path(n.a, env);
                below.scale *= n.scale_fn(env);
                return below;
            }
            case VecOp::Rotate: {
                auto below = linear_on_path(n.a, env);
                below.rotation = n.rot_fn(env) * below.rotation;
                return below;
            }
            case VecOp::Add:
            case VecOp::Sub: {
                if (reads_point(n.a)) return linear_on_path(n.a, env);
                if (reads_point(n.b)) {
                    auto below = linear_on_path(n.b, env);
                    if (n.op == VecOp::Sub) below.scale = -below.scale;
                    return below;
                }
                // The point does not appear: the shape collapses.
                return Linear{Matrix<T, 3, 3>::identity(), T{0}};
            }
            default: return Linear{Matrix<T, 3, 3>::identity(), T{0}};
        }
    }

    void push(VecOp o) {
        VecFieldOp<T> n{};
        n.op = o;
        ops_.push_back(std::move(n));
    }

    template<typename F, typename Fn>
    static VecField make_opaque(Fn fn, std::size_t payload, std::type_index type, bool stateless) {
        VecField f;
        f.ops_.clear();
        VecFieldOp<T> n{};
        n.op = VecOp::Opaque;
        n.type = type;
        n.payload = payload;
        n.stateless = stateless;
        n.fn = std::move(fn);
        f.ops_.push_back(std::move(n));
        f.structural_ = false;
        return f;
    }

    template<typename F>
    static VecField make_opaque(F fn) {
        return make_opaque<F>(std::move(fn), sizeof(F),
                              std::type_index(typeid(F)), std::is_empty_v<F>);
    }

    Vec<T, 3> eval_into(Vec<T, 3>* s, const MotionEnv<T>& env) const {
        for (std::size_t i = 0; i < ops_.size(); ++i) {
            const auto& n = ops_[i];
            switch (n.op) {
                case VecOp::Point:  s[i] = env.p; break;
                case VecOp::Const:  s[i] = n.value; break;
                case VecOp::Add:    s[i] = Vec<T, 3>{s[n.a] + s[n.b]}; break;
                case VecOp::Sub:    s[i] = Vec<T, 3>{s[n.a] - s[n.b]}; break;
                case VecOp::Scale:  s[i] = Vec<T, 3>{s[n.a] * n.scale_fn(env)}; break;
                case VecOp::Rotate: s[i] = Vec<T, 3>{n.rot_fn(env) * s[n.a]}; break;
                case VecOp::Opaque: s[i] = n.fn(env); break;
            }
        }
        return s[ops_.size() - 1];
    }

    // Same rule and the same reason as Field::append -- see there for why
    // reserving exactly what is needed on every append is quadratic.
    void append(VecOp o, VecField& y) {
        const auto shift  = static_cast<std::uint32_t>(ops_.size());
        const auto x_root = shift - 1;
        const std::size_t need = ops_.size() + y.ops_.size() + 1;
        if (need > ops_.capacity())
            ops_.reserve(std::max(need, ops_.capacity() * 2));

        for (auto& n : y.ops_) {
            VecFieldOp<T> c = std::move(n);
            if (is_binary(c.op) || is_unary(c.op)) {
                c.a += shift;
                if (is_binary(c.op)) c.b += shift;
            }
            assert((!(is_binary(c.op) || is_unary(c.op)) ||
                    (c.a < ops_.size() && (!is_binary(c.op) || c.b < ops_.size()))) &&
                   "VecField: a child must precede its parent");
            ops_.push_back(std::move(c));
        }

        VecFieldOp<T> parent{};
        parent.op = o;
        parent.a  = x_root;
        parent.b  = static_cast<std::uint32_t>(ops_.size() - 1);
        assert(parent.a < ops_.size() && parent.b < ops_.size() &&
               "VecField: a child must precede its parent");
        ops_.push_back(std::move(parent));

        structural_ = structural_ && y.structural_;
    }
};

// A point field counts into the same report. Before this existed, a
// motion slot could not be classified at all -- std::move_only_function
// has no target_type() -- so a report covering only thicknesses would say
// "unknown = 0" while every motion in the scene sat in a slot it could
// not see. An honest zero and an invisible population are exactly the two
// states a single number cannot distinguish.
template<Scalar T>
void accumulate(FieldStats& stats, const VecField<T>& f,
                std::vector<std::type_index>& seen) {
    ++stats.fields;
    if (f.is_structural()) ++stats.structural_fields;
    else                   ++stats.opaque_fields;

    for (std::size_t i = 0; i < f.size(); ++i) {
        const auto& n = f.op(i);
        if (n.op != VecOp::Opaque) continue;
        ++stats.opaque_leaves;
        if (n.type == std::type_index(typeid(void))) { ++stats.unknown_leaves; continue; }
        ++stats.recognized_leaves;
        stats.payload_bytes += n.payload;
        if (std::find(seen.begin(), seen.end(), n.type) == seen.end())
            seen.push_back(n.type);
    }
    stats.distinct_types = seen.size();
}

} // namespace spatium::io::build
