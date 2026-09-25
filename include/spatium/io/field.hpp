#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <spatium/algebra/matrix.hpp>
#  include <spatium/algebra/groups/so3.hpp>
#  include <spatium/algebra/noise.hpp>
#  include <cassert>
#  include <algorithm>
#  include <cmath>
#  include <cstddef>
#  include <cstdint>
#  include <functional>
#  include <memory>
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

    // Min and Max, and they earn their place rather than rounding out a
    // set. Without them the vocabulary cannot clamp, and without a clamp
    // it cannot write a smoothstep -- which is what twenty of the donut's
    // twenty-four motion factors are. The arithmetic half of a smoothstep
    // (`e * e * (3 - 2 * e)`) was always expressible; `std::clamp` in
    // front of it was the entire reason those factors had to be closures.
    //
    // They are also the right kind of addition: total on every input,
    // needing no branch a consumer has to model, and a single native
    // instruction on every target an export would aim at.
    Min, Max,

    // The rest were each found by reading a leaf that could not be written
    // without them, not by filling out a list. See `docs/gpu-abi-design.md`
    // for the vocabulary they were checked against.
    //
    // Sin and Cos, because the icing's drip is noise sampled around the
    // ring at (cos u, sin u) -- a circle, so that the field closes on
    // itself where u wraps. One child each.
    Sin, Cos,

    // 1 where the first child is less than the second, 0 otherwise.
    //
    // This is the decision about booleans the cube's hard step was waiting
    // on, and the decision is to not have them: a comparison yields a
    // scalar that is exactly 0 or exactly 1, the way GLSL's `step` does,
    // so the IR keeps one value type and a select is a multiply. The
    // cube's `time < 0.12 ? 1.0 : 0.0` is then `less(t, 0.12)`, the same
    // two values bit for bit.
    Less,

    // Perlin noise of three children, from the table carried in `noise`.
    // Every field on the donut's surfaces is built on it, and a noise
    // table is data -- 512 bytes fixed by a seed -- so it does not need to
    // be a closure to be carried along.
    Noise,

    // Inputs that are not a surface parameter, for fields evaluated per
    // *instance*: which instance (`Id`), where it started (`Origin`, one
    // component, `k`), and the point being coloured (`Point`, `k`). A
    // field evaluated over (u, v) reads them as zero.
    //
    // They exist because the donut's dust is a million particles sharing
    // one field, and every one of its four remaining closures began by
    // telling the particles apart -- which, without these, only a closure
    // could do.
    Id, Origin, Point,

    // A number in [0, 1) that is fixed per instance and per `salt` (in
    // `k`): the instance index through FNV and an avalanche, reduced mod
    // 10^6. Integer arithmetic only, so a host in fp64 and a device in fp32
    // hash the same integers to the same bits before either converts.
    Hash,

    // One child. For the distance a colour fades over.
    Sqrt,

    // Component `k` of the point at `floor(child)` in `table`, the index
    // clamped into range. The dust's targets are points on the BOOM
    // letterforms, a fixed table like a noise permutation is; a gather is
    // how a field reads one without a closure holding the table.
    Gather,

    Opaque,  // a callable leaf -- the escape hatch, kept first-class

    // Coordinate `k` of a point of spacetime, 0..3: the input a metric
    // reads, so a metric is a field like any other -- lowered, compiled,
    // and differentiated through the same pool. Appended after Opaque
    // rather than beside the other inputs, because an op's number is its
    // code in lowered data and inserting would renumber what is saved.
    Coord,
};

// The FNV-and-avalanche hash `Op::Hash` is defined by, spelled once so the
// field, a POD interpreter and any shader written against it cannot drift.
inline std::uint32_t instance_hash(std::uint32_t id, std::uint32_t salt) {
    std::uint32_t h = 2166136261u ^ salt;
    h ^= id;
    h *= 16777619u;
    h ^= h >> 13;
    h *= 0x85ebca6bu;
    h ^= h >> 16;
    return h;
}

template<Scalar T>
inline T instance_unit(std::uint32_t id, std::uint32_t salt) {
    return static_cast<T>(instance_hash(id, salt) % 1000000u) / static_cast<T>(1000000);
}

inline const char* op_name(Op o) {
    switch (o) {
        case Op::Const:  return "Const";
        case Op::U:      return "U";
        case Op::V:      return "V";
        case Op::Add:    return "Add";
        case Op::Sub:    return "Sub";
        case Op::Mul:    return "Mul";
        case Op::Div:    return "Div";
        case Op::Min:    return "Min";
        case Op::Max:    return "Max";
        case Op::Sin:    return "Sin";
        case Op::Cos:    return "Cos";
        case Op::Less:   return "Less";
        case Op::Noise:  return "Noise";
        case Op::Id:     return "Id";
        case Op::Origin: return "Origin";
        case Op::Point:  return "Point";
        case Op::Hash:   return "Hash";
        case Op::Sqrt:   return "Sqrt";
        case Op::Gather: return "Gather";
        case Op::Opaque: return "Opaque";
        case Op::Coord:  return "Coord";
    }
    return "?";
}

inline bool is_binary(Op o) {
    return o == Op::Add || o == Op::Sub || o == Op::Mul || o == Op::Div ||
           o == Op::Min || o == Op::Max || o == Op::Less;
}

// How many children an op reads, named once for the three places that
// shift indices, check order and lower -- each of which would otherwise
// have to learn separately that Sin has one child and Noise three.
inline std::uint32_t arity(Op o) {
    if (o == Op::Sin || o == Op::Cos || o == Op::Sqrt || o == Op::Gather) return 1;
    if (o == Op::Noise)               return 3;
    return is_binary(o) ? 2 : 0;
}

template<Scalar T = double>
struct FieldOp {
    Op op = Op::Const;

    T value{};                          // Const
    std::uint32_t a = 0, b = 0, c = 0;  // child indices, all < this index

    // Noise: the table. Shared rather than held by value because a Field
    // copies its ops whenever it is composed, and 512 bytes per copy is
    // what the capture accounting below exists to catch. Never mutated
    // after construction, so sharing it is not aliasing anything.
    std::shared_ptr<const algebra::PerlinNoise> noise;

    // An integer immediate: the component for Origin, Point and Gather, the
    // salt for Hash. Separate from the child indices so that an op with no
    // children never has a number in a slot that means "child".
    std::uint32_t k = 0;

    // Gather: the table. Shared for the reason `noise` is.
    std::shared_ptr<const std::vector<Vec<T, 3>>> table;

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

// What a field can read. A surface field reads (u, v) and nothing else;
// a field evaluated per instance -- a motion factor, a colour -- also has
// an instance, the point it started from and the point being coloured.
// One record rather than a growing argument list, for the reason
// `MotionEnv` is one.
template<Scalar T = double>
struct FieldInputs {
    T u{}, v{};
    std::uint32_t id = 0;
    Vec<T, 3> origin{};
    Vec<T, 3> p{};
    Vec<T, 4> x{};   // a point of spacetime, for a metric; see Op::Coord
};

template<Scalar T = double>
class Field {
public:
    Field() { push_const(T{0}); }
    Field(T constant) { push_const(constant); }   // NOLINT: implicit on purpose

    static Field u() { Field f; f.ops_.clear(); f.push(Op::U); return f; }
    static Field v() { Field f; f.ops_.clear(); f.push(Op::V); return f; }

    // The same slot as `u()`, named for what a motion factor binds it to.
    //
    // `scaled` and `rotated` read the first parameter as time, and that is
    // a convention rather than something a reader can derive -- so a call
    // site that says `t()` states it where the reader meets it, instead of
    // relying on them having read the overload. `docs/gpu-abi-design.md`
    // wrote this spelling down before the overloads existed, which is a
    // decent sign it is the one people reach for.
    static Field t() { return u(); }

    // Per-instance inputs; see Op::Id, Op::Origin, Op::Point.
    static Field id() { Field f; f.ops_.clear(); f.push(Op::Id); return f; }
    static Field origin(std::uint32_t k) { return input(Op::Origin, k); }
    static Field point(std::uint32_t k) {
        Field f = input(Op::Point, k);
        f.reads_point_ = true;
        return f;
    }

    // A fixed number in [0, 1) per instance; `salt` picks an independent
    // one. See Op::Hash.
    static Field hash(std::uint32_t salt) { return input(Op::Hash, salt); }

    // Coordinate k of the spacetime point a metric is evaluated at.
    static Field coord(std::uint32_t k) { return input(Op::Coord, k); }

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

    T operator()(T u, T v) const { return (*this)(FieldInputs<T>{u, v}); }

    T operator()(const FieldInputs<T>& in) const {
        constexpr std::size_t kInline = 32;
        if (ops_.size() <= kInline) {
            T scratch[kInline];
            return eval_into(scratch, in);
        }
        std::vector<T> scratch(ops_.size());
        return eval_into(scratch.data(), in);
    }

    // Whether any op reads the point being coloured or moved. Derived like
    // `is_structural()`: a factor that reads it varies per vertex, and a
    // motion scaled by such a factor is a deformation however it is spelled.
    bool reads_point() const { return reads_point_; }

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

    // Hidden friends, like the operators above, so they are found by ADL
    // on a Field and do not sit in the enclosing namespace waiting to be
    // confused with `std::min`.
    friend Field min(Field x, const Field& y) { x.append(Op::Min, y); return x; }
    friend Field max(Field x, const Field& y) { x.append(Op::Max, y); return x; }

    // Written once here rather than at each call site, because `min(max(x,
    // lo), hi)` is the kind of thing that is spelled backwards eventually.
    // Not a new op: it lowers to the two above, so a consumer that knows
    // Min and Max needs to learn nothing to read a clamp.
    friend Field clamp(Field x, const Field& lo, const Field& hi) {
        return min(max(std::move(x), lo), hi);
    }

    // The smoothstep the scene's motion factors are made of, as one
    // expression: `e = clamp(x, 0, 1)`, then `e * e * (3 - 2 * e)`.
    //
    // It lives here rather than in a demo because it is the shape a growth
    // curve takes every time, and because a factor written by hand each
    // time is a factor that becomes a lambda the moment it gets long.
    friend Field smoothstep(Field x) {
        Field e = clamp(std::move(x), Field{T{0}}, Field{T{1}});
        return e * e * (Field{T{3}} - Field{T{2}} * e);
    }

    friend Field sin(Field x) { x.push_unary(Op::Sin); return x; }
    friend Field cos(Field x) { x.push_unary(Op::Cos); return x; }

    // `less(a, b)` is 1 where a < b and 0 elsewhere -- a scalar, not a
    // boolean, so the result composes with everything else here by
    // multiplication. See Op::Less for why there is no boolean type.
    friend Field less(Field x, const Field& y) { x.append(Op::Less, y); return x; }

    friend Field sqrt(Field x) { x.push_unary(Op::Sqrt); return x; }

    // Component `k` of `table[floor(index)]`. The table is shared, not
    // copied into every field that gathers from it.
    friend Field gather(std::shared_ptr<const std::vector<Vec<T, 3>>> table, Field index,
                        std::uint32_t k) {
        index.push_unary(Op::Gather);
        index.ops_.back().table = std::move(table);
        index.ops_.back().k = k;
        return index;
    }

    // Noise over three fields. The table is copied once into shared
    // storage here, so the caller's PerlinNoise can go out of scope and
    // every later copy of this field costs a reference count, not 512
    // bytes.
    friend Field noise(const algebra::PerlinNoise& n, Field x, const Field& y, const Field& z) {
        const auto xr = static_cast<std::uint32_t>(x.ops_.size() - 1);
        const auto yr = x.append_operand(y);
        const auto zr = x.append_operand(z);
        FieldOp<T> parent{};
        parent.op    = Op::Noise;
        parent.a     = xr;
        parent.b     = yr;
        parent.c     = zr;
        parent.noise = std::make_shared<const algebra::PerlinNoise>(n);
        x.ops_.push_back(std::move(parent));
        return x;
    }

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
            const auto k = arity(n.op);
            if ((k > 0 && n.a >= i) || (k > 1 && n.b >= i) || (k > 2 && n.c >= i))
                return false;
        }
        return true;
    }
private:
    std::vector<FieldOp<T>> ops_;
    bool structural_ = true;
    bool reads_point_ = false;

    static Field input(Op o, std::uint32_t k) {
        Field f;
        f.ops_.clear();
        f.push(o);
        f.ops_.back().k = k;
        return f;
    }

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

    T eval_into(T* s, const FieldInputs<T>& in) const {
        for (std::size_t i = 0; i < ops_.size(); ++i) {
            const auto& n = ops_[i];
            switch (n.op) {
                case Op::Const:  s[i] = n.value; break;
                case Op::U:      s[i] = in.u; break;
                case Op::V:      s[i] = in.v; break;
                case Op::Add:    s[i] = s[n.a] + s[n.b]; break;
                case Op::Sub:    s[i] = s[n.a] - s[n.b]; break;
                case Op::Mul:    s[i] = s[n.a] * s[n.b]; break;
                case Op::Div:    s[i] = s[n.a] / s[n.b]; break;
                // ADL-friendly like the rest of the library's math, so a
                // Field over a user scalar picks up that type's own min.
                case Op::Min:    { using std::min; s[i] = min(s[n.a], s[n.b]); break; }
                case Op::Max:    { using std::max; s[i] = max(s[n.a], s[n.b]); break; }
                case Op::Sin:    { using std::sin; s[i] = sin(s[n.a]); break; }
                case Op::Cos:    { using std::cos; s[i] = cos(s[n.a]); break; }
                case Op::Less:   s[i] = s[n.a] < s[n.b] ? T{1} : T{0}; break;
                case Op::Noise:  s[i] = (*n.noise)(s[n.a], s[n.b], s[n.c]); break;
                case Op::Id:     s[i] = static_cast<T>(in.id); break;
                case Op::Origin: s[i] = in.origin[n.k]; break;
                case Op::Point:  s[i] = in.p[n.k]; break;
                case Op::Hash:   s[i] = instance_unit<T>(in.id, n.k); break;
                case Op::Sqrt:   { using std::sqrt; s[i] = sqrt(s[n.a]); break; }
                case Op::Gather: s[i] = gather_at(*n.table, s[n.a], n.k); break;
                case Op::Opaque: s[i] = n.fn(in.u, in.v); break;
                case Op::Coord:  s[i] = in.x[n.k]; break;
            }
        }
        return s[ops_.size() - 1];
    }

public:
    // The gather's arithmetic, public so an interpreter of the lowered form
    // shares the rounding rather than repeating it: truncation toward zero
    // is `static_cast` to an unsigned index, which for a non-negative index
    // is a floor, then clamped into the table.
    static T gather_at(const std::vector<Vec<T, 3>>& table, T index, std::uint32_t k) {
        if (table.empty()) return T{0};
        const T lo = index < T{0} ? T{0} : index;
        auto i = static_cast<std::size_t>(lo);
        if (i >= table.size()) i = table.size() - 1;
        return table[i][k];
    }

private:

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
        const auto x_root = static_cast<std::uint32_t>(ops_.size() - 1);
        const auto y_root = append_operand(y);

        FieldOp<T> parent{};
        parent.op = o;
        parent.a  = x_root;
        parent.b  = y_root;
        assert(parent.a < ops_.size() && parent.b < ops_.size() &&
               "Field: a child must precede its parent");
        ops_.push_back(std::move(parent));
    }

    // One child, which is always the current root -- so there is nothing
    // to shift and nothing that could land out of order.
    void push_unary(Op o) {
        FieldOp<T> n{};
        n.op = o;
        n.a  = static_cast<std::uint32_t>(ops_.size() - 1);
        ops_.push_back(std::move(n));
    }

    // y's ops, shifted, and the index its root landed at. The half of
    // `append` that does not know how many operands the parent will have,
    // split out so a three-child op appends two of them the same way.
    std::uint32_t append_operand(const Field& y) {
        const auto shift = static_cast<std::uint32_t>(ops_.size());

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
            const auto k = arity(c.op);
            if (k > 0) c.a += shift;
            if (k > 1) c.b += shift;
            if (k > 2) c.c += shift;
            assert((k < 1 || c.a < ops_.size()) && (k < 2 || c.b < ops_.size()) &&
                   (k < 3 || c.c < ops_.size()) &&
                   "Field: a child must precede its parent");
            ops_.push_back(std::move(c));
        }

        structural_ = structural_ && y.structural_;
        reads_point_ = reads_point_ || y.reads_point_;
        return static_cast<std::uint32_t>(ops_.size() - 1);
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

// The leaf half of a field's accounting, on its own. Split out because a
// factor written as an expression has leaves to report -- it may itself
// hold an opaque one -- but is not a field in its own right: it is part of
// the motion it scales, and counting it separately would make `fields`
// grow when a scene is rewritten more structurally, which is backwards.
template<Scalar T>
void accumulate_leaves(FieldStats& stats, const Field<T>& f,
                       std::vector<std::type_index>& seen) {
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
}

template<Scalar T>
void accumulate(FieldStats& stats, const Field<T>& f,
                std::vector<std::type_index>& seen) {
    ++stats.fields;
    if (f.is_structural()) ++stats.structural_fields;
    else                   ++stats.opaque_fields;

    accumulate_leaves(stats, f, seen);
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

    // Which instance this is, counting a Scatter's sites in order: the
    // integer a per-instance quantity is hashed from.
    //
    // The origin used to be that input, hashed through the bit patterns of
    // its coordinates, and it cannot be on a GPU that has no fp64: the bits
    // being hashed do not exist there, so the device would give every
    // particle different parameters from the host. An index is the same
    // integer on both. A node that is its own object is instance 0.
    std::uint32_t id = 0;
};

enum class VecOp : std::uint8_t {
    Point,     // the environment's p
    Const,     // a literal vector
    Add, Sub,
    Scale,     // child a, times a scalar read from the environment
    Rotate,    // child a, turned by a rotation read from the environment
    // A vector made of three scalar fields, one per component, held in
    // `factor_fields`. The way a scalar expression -- a hash, a gathered
    // target, a noise -- becomes a position or a colour without a closure.
    Make,
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
        case VecOp::Make:   return "Make";
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

    // The same factor, written as an expression instead of as a closure --
    // empty when the caller passed a callable, which is what the two
    // erasures above have always held.
    //
    // One entry for a Scale (the scalar), three for a Rotate (the
    // axis-angle components, in order). The closures above are still the
    // evaluation path even when this is populated: they are then thin
    // wrappers around these, so nothing in eval, in `is_placement()` or in
    // `placement_at()` has to learn a second way to read a factor. What
    // this buys is not speed but *readability by the library* -- a factor
    // in here can be inspected, exported and rewritten, and a factor in a
    // closure cannot be any of the three.
    //
    // Twenty of the donut's twenty-four factors are one stateless
    // `smoothstep` of time. They were opaque for want of somewhere to
    // write them down, not for want of structure.
    std::vector<Field<T>> factor_fields;

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

    // What a factor or a component field reads, from a motion's
    // environment: time as the first parameter under the convention the
    // factors have always used, then the instance, its origin and the point.
    static FieldInputs<T> inputs_of(const MotionEnv<T>& e) {
        return FieldInputs<T>{e.t, T{0}, e.id, e.origin, e.p};
    }

    // A vector from three scalar fields. Structural exactly when all three
    // are, and reading the point exactly when one of them does -- which is
    // what decides whether a motion built from it can still be a placement.
    static VecField make(Field<T> x, Field<T> y, Field<T> z) {
        VecField f;
        f.ops_.clear();
        VecFieldOp<T> n{};
        n.op = VecOp::Make;
        const bool structural = x.is_structural() && y.is_structural() && z.is_structural();
        n.factor_fields.push_back(std::move(x));
        n.factor_fields.push_back(std::move(y));
        n.factor_fields.push_back(std::move(z));
        f.ops_.push_back(std::move(n));
        f.structural_ = structural;
        return f;
    }

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

    // A Scale or a Rotate carries its factor as a closure, and that
    // closure is every bit as opaque as a leaf built by `make_opaque`.
    // It has to be recorded as such, in the same three fields, or the
    // report cannot see it.
    //
    // It could not, and said so in our favour. `accumulate` only looks at
    // nodes whose op is `Opaque`, and these are `Scale` and `Rotate`, so
    // twenty-four closures in the donut scene were invisible to both the
    // report and `is_structural()` -- which kept answering "structural"
    // for a field with a lambda inside it. The printed figure was 66 of
    // 73 fields structural, reading as "nearly all of it"; the truth was
    // under half. A measurement that is wrong towards "everything is
    // fine" is the one that never gets questioned, and this one is the
    // instrument the export question is measured with.
    template<typename F>
    static void note_factor(VecFieldOp<T>& n) {
        n.type      = std::type_index(typeid(F));
        n.payload   = sizeof(F);
        n.stateless = std::is_empty_v<F>;
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
        note_factor<S>(n);
        assert(n.a < x.ops_.size() && "VecField: a child must precede its parent");
        x.ops_.push_back(std::move(n));
        x.structural_ = false;
        return x;
    }

    template<typename S>
        requires std::is_invocable_r_v<T, const S&, T>
    friend VecField scaled(VecField x, S s) {
        VecFieldOp<T> n{};
        n.op       = VecOp::Scale;
        n.a        = static_cast<std::uint32_t>(x.ops_.size() - 1);
        n.scale_fn = [s = std::move(s)](const MotionEnv<T>& e) { return s(e.t); };
        note_factor<S>(n);
        assert(n.a < x.ops_.size() && "VecField: a child must precede its parent");
        x.ops_.push_back(std::move(n));
        x.structural_ = false;
        return x;
    }

    // The same factor, written as an expression rather than as a callable,
    // and the overload that makes a growing object structural.
    //
    // It exists because the two above cannot be fixed by writing a scene
    // more carefully. With only callable factors, any scene that grows or
    // spins holds a closure *by construction*, so "the scene exports" is
    // not a claim that is currently false -- it is one that cannot become
    // true. That is the difference between a gap and a locked door.
    //
    // **The convention, and it is a choice rather than a consequence: the
    // factor reads time as the field's first parameter, and its second is
    // zero.** `Field` names its inputs "first" and "second" instead of u
    // and v exactly so they can be bound to something that is not a
    // surface domain, and this is that something. The second is zero and
    // not also t, so that nothing reads as meaningful by accident.
    //
    // A factor that must differ between instances stays a callable, and
    // that is not a temporary shortfall: the only such factor in the donut
    // scene picks its axis through an integer hash and a modulo, which
    // Add/Sub/Mul/Div/Const cannot express and which would therefore stay
    // an opaque leaf under any widening of this signature.
    friend VecField scaled(VecField x, Field<T> s) {
        VecFieldOp<T> n{};
        n.op = VecOp::Scale;
        n.a  = static_cast<std::uint32_t>(x.ops_.size() - 1);

        const bool structural = s.is_structural();
        n.factor_fields.push_back(s);
        n.scale_fn = [s = std::move(s)](const MotionEnv<T>& e) { return s(inputs_of(e)); };

        // No note_factor here, deliberately. The identity of a factor that
        // has an expression behind it *is* the expression, and the report
        // reads it out of `factor_fields` rather than out of a typeid.
        assert(n.a < x.ops_.size() && "VecField: a child must precede its parent");
        x.ops_.push_back(std::move(n));

        // Structural unless the expression itself holds an opaque leaf --
        // which it can, since any (u, v) callable converts to a Field. So
        // this overload does not promise structurality, it *propagates* it.
        if (!structural) x.structural_ = false;
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
        // Recorded after the fact so the report names the type the caller
        // actually wrote, not the axis-angle wrapper put around it here.
        auto f = rotated(std::move(x), [r = std::move(r)](const Vec<T, 3>& o, T t) {
            return algebra::SO3<T>{}.exp(r(o, t));
        });
        note_factor<R>(f.ops_.back());
        return f;
    }

    template<typename R>
        requires (!std::is_invocable_r_v<Matrix<T, 3, 3>, const R&, T>) &&
                 std::is_invocable_r_v<Matrix<T, 3, 3>, const R&, const Vec<T, 3>&, T>
    friend VecField rotated(VecField x, R r) {
        VecFieldOp<T> n{};
        n.op     = VecOp::Rotate;
        n.a      = static_cast<std::uint32_t>(x.ops_.size() - 1);
        n.rot_fn = [r = std::move(r)](const MotionEnv<T>& e) { return r(e.origin, e.t); };
        note_factor<R>(n);
        assert(n.a < x.ops_.size() && "VecField: a child must precede its parent");
        x.ops_.push_back(std::move(n));
        x.structural_ = false;
        return x;
    }

    template<typename R>
        requires std::is_invocable_r_v<Matrix<T, 3, 3>, const R&, T>
    friend VecField rotated(VecField x, R r) {
        VecFieldOp<T> n{};
        n.op     = VecOp::Rotate;
        n.a      = static_cast<std::uint32_t>(x.ops_.size() - 1);
        n.rot_fn = [r = std::move(r)](const MotionEnv<T>& e) { return r(e.t); };
        note_factor<R>(n);
        assert(n.a < x.ops_.size() && "VecField: a child must precede its parent");
        x.ops_.push_back(std::move(n));
        x.structural_ = false;
        return x;
    }

    template<typename R>
        requires (!std::is_invocable_r_v<Matrix<T, 3, 3>, const R&, T>) &&
                 std::is_invocable_r_v<Vec<T, 3>, const R&, T>
    friend VecField rotated(VecField x, R r) {
        auto f = rotated(std::move(x), [r = std::move(r)](T t) {
            return algebra::SO3<T>{}.exp(r(t));
        });
        note_factor<R>(f.ops_.back());
        return f;
    }

    // The turn as an expression: the axis-angle vector, one field per
    // component, under the same convention as `scaled` -- time is the
    // first parameter and the second is zero.
    //
    // Three fields rather than one vector-valued field because `Field` is
    // scalar-valued and `VecField` is the wrong shape here: a VecField
    // reads a point, and a factor that read the point would be a
    // deformation, which is the one thing this whole distinction exists to
    // keep out. Three scalars cannot express that mistake.
    //
    // Exponentiated through SO3 exactly as the callable axis-angle
    // overload above does, so a constant orientation is
    // `rotated(f, 0, 0, yaw)` and a spin is `rotated(f, 0, 0, w *
    // Field<T>::u())` -- both structural, neither needing a lambda.
    friend VecField rotated(VecField x, Field<T> ax, Field<T> ay, Field<T> az) {
        VecFieldOp<T> n{};
        n.op = VecOp::Rotate;
        n.a  = static_cast<std::uint32_t>(x.ops_.size() - 1);

        const bool structural =
            ax.is_structural() && ay.is_structural() && az.is_structural();
        n.factor_fields.push_back(ax);
        n.factor_fields.push_back(ay);
        n.factor_fields.push_back(az);
        n.rot_fn = [ax = std::move(ax), ay = std::move(ay), az = std::move(az)](
                       const MotionEnv<T>& e) {
            const auto in = inputs_of(e);
            return algebra::SO3<T>{}.exp(Vec<T, 3>{ax(in), ay(in), az(in)});
        };

        assert(n.a < x.ops_.size() && "VecField: a child must precede its parent");
        x.ops_.push_back(std::move(n));
        if (!structural) x.structural_ = false;
        return x;
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
            case VecOp::Make:   return factor_reads_point(n);
            case VecOp::Scale:
            case VecOp::Rotate: return reads_point(n.a) || factor_reads_point(n);
            case VecOp::Add:
            case VecOp::Sub:    return reads_point(n.a) || reads_point(n.b);
        }
        return true;
    }

    // A factor written as an expression can now read the point, and one
    // that does varies per vertex -- so a Scale by it is a deformation.
    // A closure factor cannot reach the point at all: its signatures take
    // an origin and a time.
    static bool factor_reads_point(const VecFieldOp<T>& n) {
        for (const auto& g : n.factor_fields)
            if (g.reads_point()) return true;
        return false;
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
            case VecOp::Make:   return !factor_reads_point(n);
            case VecOp::Scale:
            case VecOp::Rotate: return !factor_reads_point(n) && affine_in_point(n.a);
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
                case VecOp::Make: {
                    const auto in = inputs_of(env);
                    s[i] = Vec<T, 3>{n.factor_fields[0](in), n.factor_fields[1](in),
                                     n.factor_fields[2](in)};
                    break;
                }
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
        // A Scale's or a Rotate's factor is a closure too. Counting only
        // `Opaque` was the blind spot that let a field with a lambda in
        // it report as structural.
        const bool factor = (n.op == VecOp::Scale  && n.scale_fn) ||
                            (n.op == VecOp::Rotate && n.rot_fn) ||
                            n.op == VecOp::Make;

        // A factor with an expression behind it is not a leaf, and the
        // closure holding it is a wrapper rather than a subject. Report
        // what the expression itself contains -- usually nothing, which is
        // the entire point of writing it as an expression, but not always:
        // any (u, v) callable converts to a Field, so an opaque leaf can
        // still arrive this way and must still be counted.
        if (factor && !n.factor_fields.empty()) {
            for (const auto& g : n.factor_fields) accumulate_leaves(stats, g, seen);
            continue;
        }

        if (n.op != VecOp::Opaque && !factor) continue;
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
