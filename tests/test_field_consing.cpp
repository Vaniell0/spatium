#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_message.hpp>
#include <spatium/algebra/noise.hpp>
#include <spatium/io/build.hpp>
#include <spatium/io/field_glsl.hpp>
#include <spatium/io/field_pod.hpp>

#include "donut_scene.hpp"

#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <random>
#include <string>
#include <tuple>
#include <vector>

using namespace spatium;
namespace bd = spatium::io::build;

// A pool with an op nothing reads, assembled by hand: the public builders
// cannot make one, and culling has to be shown on something that has one.
namespace spatium::io::build {
template<>
struct FieldPoolAccess<double> {
    static Field<double> make(std::vector<FieldOp<double>> ops) {
        return Field<double>::adopt(std::move(ops));
    }
};
} // namespace spatium::io::build

namespace {

using F = bd::Field<double>;

std::uint64_t bits(double x) { return std::bit_cast<std::uint64_t>(x); }

double nan_with(std::uint64_t payload) {
    return std::bit_cast<double>(0x7ff8000000000000ull | payload);
}

bd::FieldOp<double> op_of(bd::Op o, std::uint32_t a = 0, std::uint32_t b = 0, std::uint32_t c = 0) {
    bd::FieldOp<double> n{};
    n.op = o;
    n.a = a;
    n.b = b;
    n.c = c;
    return n;
}

bd::FieldOp<double> const_op(double v) {
    auto n = op_of(bd::Op::Const);
    n.value = v;
    return n;
}

std::size_t count_op(const F& f, bd::Op o) {
    std::size_t n = 0;
    for (const auto& x : f.ops()) n += x.op == o;
    return n;
}

// ── A random expression, built twice ─────────────────────────────

enum class K {
    Const, U, V, Id, Origin, Point, Hash, Coord, OpaqueA, OpaqueB, OpaqueStateless,
    Sin, Cos, Sqrt, Smooth,
    Add, Sub, Mul, Div, Min, Max, Less,
    Noise, Gather,
};

struct Expr {
    K kind{};
    double num = 0;
    std::uint32_t k = 0;
    int pick = 0;
    std::vector<std::shared_ptr<Expr>> kids;
};
using ExprP = std::shared_ptr<Expr>;

const std::vector<double>& constants() {
    static const std::vector<double> v{
        0.0, -0.0, 1.0, 2.0, 0.5, -3.0, 3.0, 1e-300,
        std::numeric_limits<double>::infinity(), nan_with(1), nan_with(2), nan_with(1)};
    return v;
}

// Two tables with the same content held apart, and a third that differs
// from them only in the sign of one zero.
struct Shared {
    algebra::PerlinNoise n1a{1}, n1b{1}, n2{2};
    std::shared_ptr<const std::vector<Vec<double, 3>>> ta, tb, tc;
    Shared() {
        std::vector<Vec<double, 3>> a{{1.0, 2.0, 3.0}, {0.0, 5.5, 0.25}, {7.0, -8.0, 9.5}};
        ta = std::make_shared<const std::vector<Vec<double, 3>>>(a);
        tb = std::make_shared<const std::vector<Vec<double, 3>>>(a);
        a[1][0] = -0.0;
        tc = std::make_shared<const std::vector<Vec<double, 3>>>(a);
    }
};

ExprP random_expr(std::mt19937_64& rng, int depth, std::vector<ExprP>& seen) {
    auto pick = [&](int n) { return static_cast<int>(rng() % static_cast<std::uint64_t>(n)); };
    // Repeated subexpressions are the point: reuse one built earlier.
    if (!seen.empty() && pick(4) == 0) return seen[static_cast<std::size_t>(pick(static_cast<int>(seen.size())))];

    auto e = std::make_shared<Expr>();
    const bool leaf = depth <= 0 || pick(6) == 0;
    if (leaf) {
        switch (pick(11)) {
            case 0: case 1: case 2:
                e->kind = K::Const;
                e->num = constants()[static_cast<std::size_t>(pick(static_cast<int>(constants().size())))];
                break;
            case 3: e->kind = K::U; break;
            case 4: e->kind = K::V; break;
            case 5: e->kind = K::Id; break;
            case 6: e->kind = K::Origin; e->k = static_cast<std::uint32_t>(pick(3)); break;
            case 7: e->kind = K::Point;  e->k = static_cast<std::uint32_t>(pick(3)); break;
            case 8: e->kind = K::Hash;   e->k = static_cast<std::uint32_t>(pick(3)); break;
            case 9: e->kind = K::Coord;  e->k = static_cast<std::uint32_t>(pick(4)); break;
            default:
                e->kind = pick(3) == 0 ? K::OpaqueStateless : (pick(2) ? K::OpaqueA : K::OpaqueB);
                e->pick = pick(3);
                break;
        }
    } else {
        static const K ops[] = {K::Sin, K::Cos, K::Sqrt, K::Smooth, K::Add, K::Add, K::Sub, K::Mul,
                                K::Mul, K::Div, K::Min, K::Max, K::Less, K::Noise, K::Gather};
        e->kind = ops[pick(static_cast<int>(sizeof ops / sizeof ops[0]))];
        e->pick = pick(3);
        e->k = static_cast<std::uint32_t>(pick(3));
        std::size_t n = 2;
        if (e->kind == K::Sin || e->kind == K::Cos || e->kind == K::Sqrt || e->kind == K::Smooth ||
            e->kind == K::Gather) n = 1;
        if (e->kind == K::Noise) n = 3;
        for (std::size_t i = 0; i < n; ++i) e->kids.push_back(random_expr(rng, depth - 1, seen));
    }
    seen.push_back(e);
    return e;
}

// The same expression with merging on or off. Every leaf of the unmerged
// build is made `without_consing()`, and a composite takes its left
// operand's mode, so the whole pool is appended without a lookup.
F build(const Expr& e, bool consed, const Shared& sh) {
    auto leaf = [&](F f) { return consed ? f : f.without_consing(); };
    auto kid = [&](std::size_t i) { return build(*e.kids[i], consed, sh); };
    const double cs[3] = {0.5, -2.0, 1.25};
    switch (e.kind) {
        case K::Const: return leaf(F{e.num});
        case K::U:     return leaf(F::u());
        case K::V:     return leaf(F::v());
        case K::Id:    return leaf(F::id());
        case K::Origin: return leaf(F::origin(e.k));
        case K::Point: return leaf(F::point(e.k));
        case K::Hash:  return leaf(F::hash(e.k));
        case K::Coord: return leaf(F::coord(e.k));
        case K::OpaqueA: {
            const double c = cs[e.pick];
            return leaf(F::opaque([c](double u, double v) { return std::sin(c * u) + v; }));
        }
        case K::OpaqueB: {
            const double c = cs[e.pick];
            return leaf(F::opaque([c](double u, double v) { return c * u - v; }));
        }
        case K::OpaqueStateless:
            return leaf(F::opaque([](double u, double) { return u * u; }));
        case K::Sin:    return sin(kid(0));
        case K::Cos:    return cos(kid(0));
        case K::Sqrt:   return sqrt(kid(0));
        case K::Smooth: return smoothstep(kid(0));
        case K::Add:    return kid(0) + kid(1);
        case K::Sub:    return kid(0) - kid(1);
        case K::Mul:    return kid(0) * kid(1);
        case K::Div:    return kid(0) / kid(1);
        case K::Min:    return min(kid(0), kid(1));
        case K::Max:    return max(kid(0), kid(1));
        case K::Less:   return less(kid(0), kid(1));
        case K::Noise: {
            const algebra::PerlinNoise* ns[3] = {&sh.n1a, &sh.n1b, &sh.n2};
            return noise(*ns[e.pick], kid(0), kid(1), kid(2));
        }
        case K::Gather: {
            const auto* ts = e.pick == 0 ? &sh.ta : e.pick == 1 ? &sh.tb : &sh.tc;
            return gather(*ts, kid(0), e.k);
        }
    }
    return F{};
}

std::vector<bd::FieldInputs<double>> sample_inputs(std::mt19937_64& rng) {
    const double xs[] = {-1.0, 0.0, -0.0, 0.37, 2.5, 100.0, -7.25};
    auto x = [&] { return xs[rng() % 7]; };
    std::vector<bd::FieldInputs<double>> v;
    for (int i = 0; i < 8; ++i)
        v.push_back({x(), x(), static_cast<std::uint32_t>(rng() % 5),
                     {x(), x(), x()}, {x(), x(), x()}, {x(), x(), x(), x()}});
    return v;
}

// Whether any two ops that could have been merged were not: the pool a
// merging build leaves must have none.
bool has_mergeable_pair(const F& f) {
    using Key = std::tuple<int, std::uint32_t, std::uint32_t, std::uint32_t, std::uint32_t,
                           std::uint64_t, std::string, std::vector<std::uint64_t>>;
    std::map<Key, int> seen;
    for (const auto& n : f.ops()) {
        if (n.op == bd::Op::Opaque) continue;
        const auto k = bd::arity(n.op);
        std::string perm;
        std::vector<std::uint64_t> pts;
        if (n.op == bd::Op::Noise) {
            const auto& p = n.noise->permutation();
            perm.assign(p.begin(), p.end());
        }
        if (n.op == bd::Op::Gather)
            for (const auto& p : *n.table)
                for (int j = 0; j < 3; ++j) pts.push_back(bits(p[static_cast<std::size_t>(j)]));
        Key key{static_cast<int>(n.op), k > 0 ? n.a : 0u, k > 1 ? n.b : 0u, k > 2 ? n.c : 0u, n.k,
                n.op == bd::Op::Const ? bits(n.value) : 0u, perm, pts};
        if (++seen[key] > 1) return true;
    }
    return false;
}

} // namespace

// The oracle. A merged field must be the same function as the same
// expression built with no merging, on every bit, for all of: repeated
// subexpressions, constants that differ only in the sign of a zero or the
// payload of a NaN, noise tables and gather tables equal by content but
// held apart, and opaque leaves of one closure type with different state.
TEST_CASE("A merged field evaluates bit for bit like the unmerged build",
          "[field][consing]") {
    const Shared sh;
    std::mt19937_64 rng(20261022);
    std::size_t cases = 0, shrank = 0, ops_merged = 0, ops_unmerged = 0, evaluations = 0;

    for (int t = 0; t < 600; ++t) {
        std::vector<ExprP> seen;
        const auto e = random_expr(rng, 2 + t % 5, seen);
        const F merged = build(*e, true, sh);
        const F plain = build(*e, false, sh);
        ++cases;

        CAPTURE(t, merged.size(), plain.size());
        REQUIRE(merged.topologically_ordered());
        REQUIRE(plain.topologically_ordered());
        CHECK(merged.size() <= plain.size());
        CHECK_FALSE(has_mergeable_pair(merged));
        CHECK(merged.is_structural() == plain.is_structural());
        CHECK(merged.reads_point() == plain.reads_point());
        // Opaque leaves are copied, never joined.
        CHECK(count_op(merged, bd::Op::Opaque) == count_op(plain, bd::Op::Opaque));
        shrank += merged.size() < plain.size();
        ops_merged += merged.size();
        ops_unmerged += plain.size();

        // Nothing in a builder's pool is dead, so culling changes nothing.
        const F culled = merged.reachable_from_root();
        REQUIRE(culled.size() == merged.size());
        for (std::size_t i = 0; i < merged.size(); ++i) {
            CHECK(culled.op(i).op == merged.op(i).op);
            CHECK(culled.op(i).a == merged.op(i).a);
            CHECK(culled.op(i).b == merged.op(i).b);
            CHECK(culled.op(i).c == merged.op(i).c);
        }

        std::optional<bd::PodField<double>> pod;
        if (merged.is_structural()) {
            auto r = bd::lower(merged);
            REQUIRE(r.has_value());
            pod = std::move(*r);
        } else {
            CHECK_FALSE(bd::lower(merged).has_value());
        }

        for (const auto& in : sample_inputs(rng)) {
            ++evaluations;
            const auto want = bits(plain(in));
            CHECK(bits(merged(in)) == want);
            CHECK(bits(culled(in)) == want);
            if (pod) CHECK(bits(bd::interpret(*pod, in)) == want);
        }
    }
    INFO(cases << " expressions, " << shrank << " shrank, ops " << ops_unmerged << " -> "
               << ops_merged << ", " << evaluations << " evaluations");
    CHECK(shrank > cases / 4);
    CHECK(ops_merged < ops_unmerged);
}

TEST_CASE("Equal ops are merged, and the root stays last", "[field][consing]") {
    // (u + v) * (u + v): U, V, Add, Mul -- not seven.
    const F s = F::u() + F::v();
    const F sq = s * s;
    CHECK(sq.size() == 4);
    CHECK(sq.op(sq.size() - 1).op == bd::Op::Mul);
    CHECK(sq.topologically_ordered());
    CHECK(s.without_consing().size() == 3);
    CHECK((s.without_consing() * s).size() == 7);

    // x += x used to append to the pool it was reading.
    F x = F::u() + F{2.0};
    const double before = x(3.0, 0.0);
    x += x;
    CHECK(x(3.0, 0.0) == before + before);
    CHECK(x.size() == 4);

    // The accumulating form agrees with the copying form.
    F acc{0.0}, by_copy{0.0};
    for (int i = 0; i < 50; ++i) {
        const F term = (F::u() * F{static_cast<double>(i)} + F::hash(static_cast<std::uint32_t>(i % 4))) *
                       sin(F::v());
        acc += term;
        by_copy = by_copy + term;
    }
    CHECK(acc.size() == by_copy.size());
    CHECK(bits(acc(0.3, 0.7)) == bits(by_copy(0.3, 0.7)));
}

TEST_CASE("Constants merge by bit pattern, not by value", "[field][consing]") {
    // 0.0 and -0.0 compare equal and are different numbers (1/x).
    const F z = F{0.0} + F{-0.0};
    CHECK(count_op(z, bd::Op::Const) == 2);
    CHECK(std::signbit((F{1.0} / F{-0.0})(0, 0)));
    CHECK_FALSE(std::signbit((F{1.0} / F{0.0})(0, 0)));
    // Two NaNs with different payloads stay two; the same NaN twice is one.
    CHECK(count_op(F{nan_with(1)} + F{nan_with(2)}, bd::Op::Const) == 2);
    CHECK(count_op(F{nan_with(1)} + F{nan_with(1)}, bd::Op::Const) == 1);
    CHECK(count_op(F{0.5} + F{0.5}, bd::Op::Const) == 1);
    // And the payload survives evaluation on the merged side.
    CHECK(bits((F{nan_with(2)} + F{nan_with(1)})(0, 0)) == bits(nan_with(2) + nan_with(1)));
}

TEST_CASE("An opaque leaf is never merged", "[field][consing]") {
    auto fn = [](double u, double v) { return u - 2.0 * v; };
    const F leaf = F::opaque(fn);
    const F twice = leaf + leaf;
    CHECK(count_op(twice, bd::Op::Opaque) == 2);
    CHECK(twice.size() == 3);

    // One closure type, two states: each keeps its own value.
    auto make = [](double c) { return F::opaque([c](double u, double) { return u * c; }); };
    const F both = make(2.0) + make(5.0) + make(2.0);
    CHECK(count_op(both, bd::Op::Opaque) == 3);
    CHECK(both(3.0, 0.0) == 3.0 * 2.0 + 3.0 * 5.0 + 3.0 * 2.0);

    // A parent of two different leaves is two different ops.
    const F p = (make(2.0) + F::u()) * (make(5.0) + F::u());
    CHECK(count_op(p, bd::Op::Add) == 2);
    CHECK_FALSE(bd::lower(p).has_value());

    // Squaring a field that holds a leaf copies the leaf.
    const F sq = (make(2.0) + F::u());
    CHECK(count_op(sq * sq, bd::Op::Opaque) == 2);
}

TEST_CASE("Noise and gather merge only when their tables are equal in full",
          "[field][consing]") {
    const Shared sh;
    auto nz = [&](const algebra::PerlinNoise& n) { return noise(n, F::u(), F::v(), F{0.25}); };
    // Two PerlinNoise objects from one seed: equal by content, merged.
    CHECK(count_op(nz(sh.n1a) + nz(sh.n1b), bd::Op::Noise) == 1);
    CHECK(count_op(nz(sh.n1a) + nz(sh.n2), bd::Op::Noise) == 2);

    auto ga = [&](const auto& t) { return gather(t, F::u(), 0); };
    CHECK(count_op(ga(sh.ta) + ga(sh.ta), bd::Op::Gather) == 1);   // same pointer
    CHECK(count_op(ga(sh.ta) + ga(sh.tb), bd::Op::Gather) == 1);   // equal content
    CHECK(count_op(ga(sh.ta) + ga(sh.tc), bd::Op::Gather) == 2);   // 0.0 against -0.0
    CHECK(count_op(gather(sh.ta, F::u(), 0) + gather(sh.ta, F::u(), 1), bd::Op::Gather) == 2);

    // The one place the two tables differ is visible in the value.
    const F g = F{1.0} / (gather(sh.ta, F{1.0}, 0) * F{0.0}) + F{1.0} / (gather(sh.tc, F{1.0}, 0) * F{0.0});
    const F h = F{1.0} / (gather(sh.ta, F{1.0}, 0) * F{0.0}) + F{1.0} / (gather(sh.ta, F{1.0}, 0) * F{0.0});
    CHECK(bits(g(0, 0)) != bits(h(0, 0)));

    // Lowering keeps the two tables apart as well: a lowered table is looked
    // for by bit pattern, so 0.0 does not stand in for -0.0.
    auto pod = bd::lower(g);
    REQUIRE(pod.has_value());
    CHECK(pod->points.size() == 18);
    CHECK(bits(bd::interpret(*pod, 0.0, 0.0)) == bits(g(0, 0)));
    bd::GlslModule m;
    REQUIRE(bd::emit_scalar(m, g, "g").has_value());
    CHECK(m.points.size() == 18);
}

// ── Culling ──────────────────────────────────────────────────────

TEST_CASE("Culling keeps what the root reads, in order, and lowering emits nothing else",
          "[field][consing][pod]") {
    using bd::Op;
    // 0 Const 1 (dead)  1 U  2 Sin(1) (dead)  3 Const 2.5  4 V  5 Add(3,4)  6 Cos(dead)
    // is not allowed (the root must be last), so the root is 7.
    std::vector<bd::FieldOp<double>> ops;
    ops.push_back(const_op(1.0));
    ops.push_back(op_of(Op::U));
    ops.push_back(op_of(Op::Sin, 1));
    ops.push_back(const_op(2.5));
    ops.push_back(op_of(Op::V));
    ops.push_back(op_of(Op::Add, 3, 4));
    ops.push_back(op_of(Op::Cos, 2));
    ops.push_back(op_of(Op::Mul, 5, 3));
    const F f = bd::FieldPoolAccess<double>::make(ops);
    REQUIRE(f.topologically_ordered());
    REQUIRE(f.size() == 8);

    const F c = f.reachable_from_root();
    REQUIRE(c.size() == 4);
    CHECK(c.topologically_ordered());
    CHECK(c.op(0).op == Op::Const);
    CHECK(c.op(1).op == Op::V);
    CHECK(c.op(2).op == Op::Add);
    CHECK(c.op(3).op == Op::Mul);
    CHECK(c.op(2).a == 0);
    CHECK(c.op(2).b == 1);
    CHECK(c.op(3).a == 2);
    CHECK(c.op(3).b == 0);
    for (double u : {0.0, -0.0, 1.5}) for (double v : {-2.0, 0.25})
        CHECK(bits(c(u, v)) == bits(f(u, v)));
    CHECK(c.reachable_from_root().size() == c.size());

    auto pod = bd::lower(f);
    REQUIRE(pod.has_value());
    CHECK(pod->ops.size() == 4);
    for (const auto& p : pod->ops) {
        CHECK(p.code != static_cast<std::uint8_t>(Op::Sin));
        CHECK(p.code != static_cast<std::uint8_t>(Op::Cos));
        CHECK(p.code != static_cast<std::uint8_t>(Op::U));
    }
    CHECK(bits(bd::interpret(*pod, 0.7, 1.9)) == bits(f(0.7, 1.9)));

    bd::GlslModule m;
    REQUIRE(bd::emit_scalar(m, f, "g").has_value());
    CHECK(m.code.find("sin(") == std::string::npos);
    CHECK(m.code.find("cos(") == std::string::npos);
    CHECK(m.code.find("in_.u") == std::string::npos);
}

TEST_CASE("A dead opaque leaf does not stop a field from lowering", "[field][consing][pod]") {
    using bd::Op;
    std::vector<bd::FieldOp<double>> ops;
    auto dead = op_of(Op::Opaque);
    dead.fn = [](double u, double) { return u; };
    ops.push_back(dead);
    ops.push_back(op_of(Op::U));
    ops.push_back(op_of(Op::Sin, 1));
    const F f = bd::FieldPoolAccess<double>::make(ops);
    CHECK_FALSE(f.is_structural());
    const F c = f.reachable_from_root();
    CHECK(c.is_structural());
    CHECK(c.size() == 2);
    auto pod = bd::lower(f);
    REQUIRE(pod.has_value());
    CHECK(bits(bd::interpret(*pod, 0.4, 0.0)) == bits(f(0.4, 0.0)));
}

TEST_CASE("Culling a random pool with dead ops interleaved restores the original",
          "[field][consing]") {
    const Shared sh;
    std::mt19937_64 rng(77);
    for (int t = 0; t < 200; ++t) {
        std::vector<ExprP> seen;
        const auto e = random_expr(rng, 2 + t % 4, seen);
        const F f = build(*e, true, sh);

        // Before each op, maybe a dead op reading earlier live ones.
        std::vector<bd::FieldOp<double>> ops;
        std::vector<std::uint32_t> at(f.size());
        for (std::size_t i = 0; i < f.size(); ++i) {
            if (i > 0 && rng() % 3 == 0) {
                const auto a = at[rng() % i], b = at[rng() % i];
                ops.push_back(op_of(rng() % 2 ? bd::Op::Add : bd::Op::Mul, a, b));
            }
            auto n = f.op(i);
            const auto k = bd::arity(n.op);
            if (k > 0) n.a = at[n.a];
            if (k > 1) n.b = at[n.b];
            if (k > 2) n.c = at[n.c];
            at[i] = static_cast<std::uint32_t>(ops.size());
            ops.push_back(std::move(n));
        }
        const F padded = bd::FieldPoolAccess<double>::make(ops);
        REQUIRE(padded.topologically_ordered());
        const F back = padded.reachable_from_root();
        REQUIRE(back.size() == f.size());
        for (std::size_t i = 0; i < f.size(); ++i) {
            CHECK(back.op(i).op == f.op(i).op);
            CHECK(back.op(i).a == f.op(i).a);
            CHECK(back.op(i).b == f.op(i).b);
            CHECK(back.op(i).c == f.op(i).c);
        }
        for (const auto& in : sample_inputs(rng))
            CHECK(bits(padded(in)) == bits(f(in)));
        if (f.is_structural()) {
            auto pod = bd::lower(padded);
            REQUIRE(pod.has_value());
            CHECK(pod->ops.size() == f.size());
        }
    }
}

// ── The real scene ───────────────────────────────────────────────

namespace {

struct OpTally {
    std::size_t fields = 0, scalar_ops = 0, tree_ops = 0, vec_ops = 0;
};

// What an append-only pool with no merging holds for `f`: the size of the
// tree it denotes, every shared subexpression counted at each use.
std::uint64_t tree_size(const F& f) {
    std::vector<std::uint64_t> t(f.size());
    for (std::size_t i = 0; i < f.size(); ++i) {
        const auto& n = f.op(i);
        const auto k = bd::arity(n.op);
        t[i] = 1 + (k > 0 ? t[n.a] : 0) + (k > 1 ? t[n.b] : 0) + (k > 2 ? t[n.c] : 0);
    }
    return t.back();
}

void tally(OpTally& o, const F& f) {
    ++o.fields;
    o.scalar_ops += f.size();
    o.tree_ops += tree_size(f);
}

void tally(OpTally& o, const bd::VecField<double>& v) {
    o.vec_ops += v.size();
    for (std::size_t i = 0; i < v.size(); ++i)
        for (const auto& f : v.op(i).factor_fields) tally(o, f);
}

} // namespace

TEST_CASE("The donut scene's fields get smaller, and still agree with the interpreter",
          "[field][consing][donut]") {
    auto boom = donut::load_boom_points(SPATIUM_SOURCE_DIR "/examples/data/boom_points.txt");
    REQUIRE_FALSE(boom.empty());
    bd::Trace<double> scene;
    donut::build_scene(scene, 11000, 35200, boom);

    OpTally o;
    std::size_t culled_smaller = 0, structural = 0, mismatches = 0, evaluations = 0;
    for (std::size_t i = 0; i < scene.size(); ++i) {
        const auto& n = scene.node(i);
        if (n.kind == bd::Kind::Offset) tally(o, n.thickness);
        tally(o, n.transform);
        tally(o, n.color_fn);
        tally(o, n.emissive_fn);

        // Every factor of every point field: culled equals as built, and
        // lowering gives the same bits.
        for (const auto* v : {&n.transform, &n.color_fn, &n.emissive_fn})
            for (std::size_t k = 0; k < v->size(); ++k)
                for (const auto& f : v->op(k).factor_fields) {
                    const F c = f.reachable_from_root();
                    culled_smaller += c.size() < f.size();
                    if (!f.is_structural()) continue;
                    ++structural;
                    auto pod = bd::lower(f);
                    REQUIRE(pod.has_value());
                    CHECK(pod->ops.size() == c.size());
                    for (double t : {-0.5, 0.0, 0.12, 0.7, 2.2, 5.9}) {
                        for (std::uint32_t id : {0u, 17u, 35199u}) {
                            bd::FieldInputs<double> in{t, 0.0, id, {0.1, -0.2, 0.3}, {1.0, 2.0, -0.5}};
                            ++evaluations;
                            if (bits(bd::interpret(*pod, in)) != bits(f(in))) ++mismatches;
                            if (bits(c(in)) != bits(f(in))) ++mismatches;
                        }
                    }
                }
    }
    WARN("donut scene scalar ops: " << o.tree_ops << " unmerged -> " << o.scalar_ops
         << " merged over " << o.fields << " fields (point-field ops " << o.vec_ops << ")");
    CHECK(culled_smaller == 0);   // the builders leave nothing to cull
    CHECK(structural > 0);
    CHECK(evaluations > 0);
    CHECK(mismatches == 0);
    CHECK(o.scalar_ops < o.tree_ops);
}
