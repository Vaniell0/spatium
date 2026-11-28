// Laws of the library's spaces, found by enumeration with no concept told.
//
// Terms are built over the operations of core/access.hpp and a little arithmetic, typed by sort:
//
//   R a number, P a point, V@x a tangent at the point term x
//   log(x, y): V@x   exp(x, v@x): P   g_x(u@x, w@x): R   d(x, y): R
//   2v, v/2, -v, u + w (one base)     a + b, a * b, the constants 0, 1/2, 2
//
// and by a number's degree in length (d and log 1, g the sum of its tangents', a constant 0):
// numbers of different degree are never added or compared, so "d <= 4", a fact about where a
// sampler puts its points, cannot be written.
//
// Every term is evaluated on every space of the connectivity matrix (tests/connectivity/probes.hpp)
// on the same tuples (p, q, r, u, w): drawn by the space's sampler, then moved by a breaker so that
// what holds only in general position or at one scale does not survive -- each subset of q, r, u, w
// drawn toward p, coincident and aligned inputs, and tangents made long enough to pass any cut
// locus here. Built bottom up by size, a term equal to an older one wherever both are defined is
// merged into it, and the merge is a law of every space; of the terms kept, a pair equal or
// ordered on some spaces is a law of those spaces. A law that is an instance of a smaller one, or
// follows from smaller ones by equality, congruence or a constant factor, is not reported.
//
// The enumeration and the merge, the laws and their pruning are rsc/core/law_engine.hpp; this file is what
// spatium gives it: the sorts and operations above, the spaces of the matrix, the breaker, the terms of each size.
//
// What it prints is JSON, which rsc_laws.py reads: the spaces, the laws with the spaces each holds
// on and is defined on, and the operations a space's types do not have.
//
//   rsc_laws [max size] [derivations.jsonl]   the second writes each law the closure derived, with its chain
#include "probes.hpp"

#include "../core/law_engine.hpp"

#include <algorithm>
#include <bit>
#include <cstdint>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <tuple>
#include <variant>
#include <vector>

namespace laws {

using namespace connectivity;

enum class Op { VarP, VarV, Log, Exp, G, Dist, Scale, Neg, VAdd, RAdd, RMul, RConst };
enum class Sort { R, P, V };

bool space_ops(const std::string& t);

struct Sig {
    using Op = laws::Op;
    using Sort = laws::Sort;
    static bool is_var(Op o) { return o == Op::VarP || o == Op::VarV; }
    static bool is_const(Op o) { return o == Op::RConst; }
    static bool is_mul(Op o) { return o == Op::RMul; }
    static bool is_add(Op o) { return o == Op::RAdd; }
    static bool ordered(Sort s) { return s == Sort::R; }
    static bool interesting(const rsc::laws::Term<Sig>& t) { return space_ops(t.text); }
};
using Term = rsc::laws::Term<Sig>;
using Rows = rsc::laws::Rows;
using World = rsc::laws::World<Sig>;

constexpr int K_WIDE = 192;    // tuples (p, q, r, u, w) per space at the sampler's scale
constexpr int GROUP = 12;      // then the breaker's, a group of each kind:
constexpr int SHRINKS = 15;    //   each nonempty subset of {q, r, u, w} drawn toward p, to a fiftieth
constexpr int DEGENERATE = 5;  //   r = q, q = p, u = log(p, q), w = log(p, r), r = exp(p, u)
constexpr int ALIGNED = 22;    //   u or w along log(p, q) or log(p, r), at -2, -1, -1/2, 1/2, 1, 2 (less two above)
constexpr int LONG = 3;        //   u, w, or the way to q made 8 long in the space's own metric
constexpr int K = K_WIDE + GROUP * (SHRINKS + DEGENERATE + ALIGNED + LONG);
constexpr double TOL = 1e-6;

template<std::size_t N> void flat(const Vec<double, N>& v, std::vector<double>& out) {
    for (std::size_t i = 0; i < N; ++i) out.push_back(v[i]);
}
template<std::size_t R_, std::size_t C_> void flat(const Matrix<double, R_, C_>& m, std::vector<double>& out) {
    for (std::size_t i = 0; i < R_; ++i)
        for (std::size_t j = 0; j < C_; ++j) out.push_back(m(i, j));
}
template<class A, class B> void flat(const spatium::product_detail::Pair<A, B>& p, std::vector<double>& out) {
    flat(p.first, out);
    flat(p.second, out);
}

bool rows_eq(const Rows& x, const Rows& y) { return rsc::laws::rows_eq(x, y, TOL); }

template<class Spc>
struct Impl final : rsc::laws::Space<Sig> {
    using S = decltype(Spc::template make<double>());
    using PT = P<S>;
    using VT = V<S>;
    struct PW { PT x; };
    struct VW { VT x; };
    using Val = std::variant<std::monostate, PW, VW, double>;
    S s = Spc::template make<double>();
    std::vector<std::vector<Val>> vals, staged_vals;
    std::vector<std::array<PT, 3>> pts;
    std::vector<std::array<VT, 2>> tans;

    explicit Impl(std::string n) {
        name = std::move(n);
        Rand r;
        std::vector<std::decay_t<decltype(Spc::points(r)[0])>> pd;
        while (pd.size() < 3 * K) for (const auto& x : Spc::points(r)) pd.push_back(x);
        for (int t = 0; t < K; ++t) {
            std::array<PT, 3> P3{PT(pd[3 * t]), PT(pd[3 * t + 1]), PT(pd[3 * t + 2])};
            std::array<VT, 2> V2{VT(Spc::tangent(r, pd[3 * t])), VT(Spc::tangent(r, pd[3 * t]))};
            if (t >= K_WIDE) brake(P3, V2, (t - K_WIDE) / GROUP);
            pts.push_back(P3);
            tans.push_back(V2);
        }
    }

    // The breaker: group g of the tuples past the sampler's.
    void brake(std::array<PT, 3>& P3, std::array<VT, 2>& V2, int g) const {
        auto toward_p = [&](int m, double f) {
            try { P3[m] = PT(spaces::exp_map(s, P3[0], VT(spaces::log_map(s, P3[0], P3[m])), f)); } catch (...) {}
        };
        auto lengthen = [&](const VT& v) {
            if constexpr (requires { spaces::metric(s, P3[0], v, v); }) {
                const double n = std::sqrt(double(spaces::metric(s, P3[0], v, v)));
                return n > 0 ? VT(v * (8.0 / n)) : v;
            } else return v;
        };
        try {
            if (g < SHRINKS) {
                const int subset = g + 1;   // bits: q, r, u, w
                if (subset & 1) toward_p(1, 0.02);
                if (subset & 2) toward_p(2, 0.02);
                if (subset & 4) V2[0] = VT(V2[0] * 0.02);
                if (subset & 8) V2[1] = VT(V2[1] * 0.02);
            } else if (g < SHRINKS + DEGENERATE) {
                switch (g - SHRINKS) {
                case 0: P3[2] = P3[1]; break;
                case 1: P3[1] = P3[0]; break;
                case 2: V2[0] = VT(spaces::log_map(s, P3[0], P3[1])); break;
                case 3: V2[1] = VT(spaces::log_map(s, P3[0], P3[2])); break;
                default: P3[2] = PT(spaces::exp_map(s, P3[0], V2[0], 1.0)); break;
                }
            } else if (g < SHRINKS + DEGENERATE + ALIGNED) {
                static const double alphas[] = {-2, -1, -0.5, 0.5, 1, 2};
                std::vector<std::tuple<int, int, double>> kinds;   // (tangent, point, alpha)
                for (int tv = 0; tv < 2; ++tv)
                    for (int pt = 1; pt < 3; ++pt)
                        for (double al : alphas)
                            if (!(al == 1 && tv + 1 == pt)) kinds.emplace_back(tv, pt, al);
                const auto [tv, pt, al] = kinds[g - SHRINKS - DEGENERATE];
                V2[tv] = VT(VT(spaces::log_map(s, P3[0], P3[pt])) * al);
            } else {
                switch (g - SHRINKS - DEGENERATE - ALIGNED) {
                case 0: V2[0] = lengthen(V2[0]); break;
                case 1: V2[1] = lengthen(V2[1]); break;
                default: P3[1] = PT(spaces::exp_map(s, P3[0], lengthen(VT(spaces::log_map(s, P3[0], P3[1]))), 1.0)); break;
                }
            }
        } catch (...) {}
    }

    const Val& arg(int a, int x) const { return vals[a][x]; }
    const PT& p(int a, int x) const { return std::get<PW>(vals[a][x]).x; }
    const VT& v(int a, int x) const { return std::get<VW>(vals[a][x]).x; }
    double rv(int a, int x) const { return std::get<double>(vals[a][x]); }

    void evaluate(const Term& t, std::vector<Val>& out, Rows& rows) const {
        out.assign(K, Val{});
        rows.at.assign(K, {});
        rows.count = 0;
        auto undef = [&](int a, int x) { return a >= 0 && arg(a, x).index() == 0; };
        for (int x = 0; x < K; ++x) {
            if (undef(t.a, x) || undef(t.b, x) || undef(t.c, x)) continue;
            try {
                switch (t.op) {
                case Op::VarP: out[x] = PW{pts[x][int(t.k)]}; break;
                case Op::VarV: out[x] = VW{tans[x][int(t.k)]}; break;
                case Op::Log: out[x] = VW{VT(spaces::log_map(s, p(t.a, x), p(t.b, x)))}; break;
                case Op::Exp: out[x] = PW{PT(spaces::exp_map(s, p(t.a, x), v(t.b, x), 1.0))}; break;
                case Op::G:
                    if constexpr (requires { spaces::metric(s, p(t.a, x), v(t.b, x), v(t.c, x)); })
                        out[x] = double(spaces::metric(s, p(t.a, x), v(t.b, x), v(t.c, x)));
                    break;
                case Op::Dist:
                    if constexpr (requires { spaces::distance(s, p(t.a, x), p(t.b, x)); })
                        out[x] = double(spaces::distance(s, p(t.a, x), p(t.b, x)));
                    break;
                case Op::Scale: out[x] = VW{VT(v(t.a, x) * t.k)}; break;
                // -v as v * -1: a matrix tangent has no unary minus, which missing() reports
                case Op::Neg: out[x] = VW{VT(v(t.a, x) * -1.0)}; break;
                case Op::VAdd: out[x] = VW{VT(v(t.a, x) + v(t.b, x))}; break;
                case Op::RAdd: out[x] = rv(t.a, x) + rv(t.b, x); break;
                case Op::RMul: out[x] = rv(t.a, x) * rv(t.b, x); break;
                case Op::RConst: out[x] = t.k; break;
                }
            } catch (...) { out[x] = Val{}; }
            if (out[x].index() == 0) continue;
            std::vector<double> c;
            std::visit([&](const auto& z) {
                using Z = std::decay_t<decltype(z)>;
                if constexpr (std::is_same_v<Z, double>) c.push_back(z);
                else if constexpr (!std::is_same_v<Z, std::monostate>) flat(z.x, c);
            }, out[x]);
            if (!std::all_of(c.begin(), c.end(), [](double y) { return std::isfinite(y); })) { out[x] = Val{}; continue; }
            rows.at[x] = std::move(c);
            ++rows.count;
        }
    }

    void stage(const std::vector<Term>& batch) override {
        staged_vals.assign(batch.size(), {});
        staged.assign(batch.size(), {});
        for (std::size_t i = 0; i < batch.size(); ++i) evaluate(batch[i], staged_vals[i], staged[i]);
    }
    void accept(std::size_t i) override {
        vals.push_back(std::move(staged_vals[i]));
        pool.push_back(std::move(staged[i]));
    }
    void clear_stage() override { staged_vals.clear(); staged.clear(); }

    std::vector<std::string> missing() const override {
        std::vector<std::string> m;
        if constexpr (!requires(const VT& a) { -a; }) m.push_back("unary minus on a tangent");
        if constexpr (!requires(const VT& a, const VT& b) { a - b; }) m.push_back("difference of two tangents");
        if constexpr (!requires(const VT& a, const VT& b) { a + b; }) m.push_back("sum of two tangents");
        if constexpr (!requires(const VT& a) { a * 2.0; }) m.push_back("a tangent times a number");
        if constexpr (!requires(const VT& a) { 2.0 * a; }) m.push_back("a number times a tangent");
        return m;
    }
};

std::string number(double k) { return k == 0.5 ? "1/2" : k == 2 ? "2" : k == -1 ? "-1" : "0"; }

// Prefix form: a variable or a constant itself, else (op args...).
std::string sexp_of(const World& w, const Term& t) {
    auto arg = [&](int i) { return " " + w.pool[i].sexp; };
    switch (t.op) {
    case Op::VarP: case Op::VarV: return t.text;
    case Op::RConst: return number(t.k);
    case Op::Log: return "(log" + arg(t.a) + arg(t.b) + ")";
    case Op::Exp: return "(exp" + arg(t.a) + arg(t.b) + ")";
    case Op::Dist: return "(d" + arg(t.a) + arg(t.b) + ")";
    case Op::G: return "(g" + arg(t.a) + arg(t.b) + arg(t.c) + ")";
    case Op::Scale: return "(scale " + number(t.k) + arg(t.a) + ")";
    case Op::Neg: return "(neg" + arg(t.a) + ")";
    case Op::VAdd: return "(vadd" + arg(t.a) + arg(t.b) + ")";
    case Op::RAdd: return "(add" + arg(t.a) + arg(t.b) + ")";
    case Op::RMul: return "(mul" + arg(t.a) + arg(t.b) + ")";
    }
    return "?";
}

// A term whose degree checks; nullopt where a sum or a comparison would mix degrees.
std::optional<Term> typed(const World& w, Term t) {
    auto D = [&](int i) { return w.pool[i].deg; };
    switch (t.op) {
    case Op::VarV: case Op::Log: case Op::Dist: t.deg = 1; break;
    case Op::Scale: case Op::Neg: t.deg = D(t.a); break;
    case Op::VAdd: if (D(t.a) != D(t.b)) return std::nullopt; t.deg = D(t.a); break;
    case Op::G: t.deg = D(t.b) + D(t.c); break;
    case Op::RAdd: if (D(t.a) != D(t.b)) return std::nullopt; t.deg = D(t.a); break;
    case Op::RMul: t.deg = D(t.a) + D(t.b); break;
    default: t.deg = 0;
    }
    t.sexp = sexp_of(w, t);
    return t;
}

// The terms of one size, from the pool's smaller ones.
std::vector<Term> level(const World& w, int size) {
    std::vector<Term> out;
    const int n = int(w.pool.size());
    auto of = [&](Sort s, int sz) {
        std::vector<int> r;
        for (int i = 0; i < n; ++i)
            if (w.pool[i].sort == s && w.pool[i].size == sz) r.push_back(i);
        return r;
    };
    auto txt = [&](int i) { return w.pool[i].text; };
    auto push = [&](Term t) { if (auto u = typed(w, std::move(t))) out.push_back(std::move(*u)); };
    for (int a : of(Sort::V, size - 1)) {
        const auto& A = w.pool[a];
        if (A.op == Op::Scale || A.op == Op::Neg) continue;   // one layer of scaling
        push({Op::Scale, a, -1, -1, 2.0, Sort::V, A.base, size, "2" + txt(a)});
        push({Op::Scale, a, -1, -1, 0.5, Sort::V, A.base, size, txt(a) + "/2"});
        push({Op::Neg, a, -1, -1, 0, Sort::V, A.base, size, "-" + txt(a)});
    }
    for (int s1 = 1; s1 <= size - 2; ++s1) {
        const int s2 = size - 1 - s1;
        for (int a : of(Sort::P, s1))
            for (int b : of(Sort::P, s2)) {
                push({Op::Dist, a, b, -1, 0, Sort::R, -1, size, "d(" + txt(a) + "," + txt(b) + ")"});
                if (a != b) push({Op::Log, a, b, -1, 0, Sort::V, a, size, "log(" + txt(a) + "," + txt(b) + ")"});
            }
        for (int a : of(Sort::P, s1))
            for (int v : of(Sort::V, s2))
                if (w.pool[v].base == a) push({Op::Exp, a, v, -1, 0, Sort::P, -1, size, "exp(" + txt(a) + "," + txt(v) + ")"});
        for (int a : of(Sort::V, s1))
            for (int b : of(Sort::V, s2))
                if (a < b && w.pool[a].base == w.pool[b].base)
                    push({Op::VAdd, a, b, -1, 0, Sort::V, w.pool[a].base, size, "(" + txt(a) + "+" + txt(b) + ")"});
        for (int a : of(Sort::R, s1))
            for (int b : of(Sort::R, s2))
                if (a <= b) {
                    const bool ca = w.pool[a].op == Op::RConst, cb = w.pool[b].op == Op::RConst;
                    const bool zero = (ca && w.pool[a].k == 0) || (cb && w.pool[b].k == 0);
                    if (!ca && !cb) push({Op::RAdd, a, b, -1, 0, Sort::R, -1, size, "(" + txt(a) + "+" + txt(b) + ")"});
                    if (((ca != cb) || a == b) && !zero) push({Op::RMul, a, b, -1, 0, Sort::R, -1, size, txt(a) + "*" + txt(b)});
                }
    }
    for (int s1 = 1; s1 <= size - 3; ++s1)
        for (int s2 = 1; s1 + s2 <= size - 2; ++s2) {
            const int s3 = size - 1 - s1 - s2;
            for (int x : of(Sort::P, s1))
                for (int u : of(Sort::V, s2))
                    for (int v : of(Sort::V, s3))
                        if (u <= v && w.pool[u].base == x && w.pool[v].base == x)
                            push({Op::G, x, u, v, 0, Sort::R, -1, size, "g_" + txt(x) + "(" + txt(u) + "," + txt(v) + ")"});
        }
    return out;
}

bool space_ops(const std::string& t) {
    return t.find("log") != std::string::npos || t.find("exp") != std::string::npos ||
           t.find("d(") != std::string::npos || t.find("g_") != std::string::npos;
}

template<class... Spcs> void add_spaces(World& w, const std::vector<std::string>& names) {
    std::size_t i = 0;
    (w.spaces.push_back(std::make_unique<Impl<Spcs>>(names[i++])), ...);
}

void all_spaces(World& w) {
    add_spaces<connectivity::E1, connectivity::E3, connectivity::E4, connectivity::S1, connectivity::S2,
               connectivity::S3, connectivity::S7, connectivity::H1, connectivity::H2, connectivity::H4,
               connectivity::S2xE1, connectivity::SPDLogE, connectivity::SPDAff, connectivity::TorusChart,
               connectivity::TorusTyped, connectivity::SphereLevelSet, connectivity::SphereLevelSetTyped,
               connectivity::S2Derived, connectivity::H2Derived, connectivity::CylinderADL, connectivity::S2xH2,
               connectivity::H2xE3, connectivity::S2xS2, connectivity::SPDAffxS2, connectivity::CylxE3>(
        w, {"E1", "E3", "E4", "S1", "S2", "S3", "S7", "H1", "H2", "H4", "S2xE1", "SPD2-logE", "SPD2-affine",
            "Torus-chart", "Torus-typed", "Sphere-implicit", "Sphere-implicit-typed", "S2-derived", "H2-derived",
            "Cylinder-ADL", "S2xH2", "H2xE3", "S2xS2", "SPD2-affinexS2", "CylinderxE3"});
    w.key_spaces = {1, 4, 8, 12, 20};   // E3, S2, H2, SPD2-affine, S2xH2: closed forms
}

// Size 1: the variables and constants. p is term 0, where the tangents live.
void atoms(World& w, bool merge) {
    std::vector<Term> first;
    const char* pn[] = {"p", "q", "r"};
    for (int i = 0; i < 3; ++i) first.push_back({Op::VarP, -1, -1, -1, double(i), Sort::P, -1, 1, pn[i]});
    first.push_back({Op::RConst, -1, -1, -1, 0.0, Sort::R, -1, 1, "0"});
    first.push_back({Op::RConst, -1, -1, -1, 0.5, Sort::R, -1, 1, "1/2"});
    first.push_back({Op::RConst, -1, -1, -1, 2.0, Sort::R, -1, 1, "2"});
    for (auto& t : first) t.sexp = t.op == Op::RConst ? number(t.k) : t.text;
    w.add(first, merge);
    w.add({{Op::VarV, -1, -1, -1, 0.0, Sort::V, 0, 1, "u", 1, "u"}, {Op::VarV, -1, -1, -1, 1.0, Sort::V, 0, 1, "w", 1, "w"}}, merge);
}

int run(int max_size, const char* derivations = nullptr) {
    World w;
    w.tol = TOL;
    all_spaces(w);
    atoms(w, true);
    for (int size = 2; size <= max_size; ++size) {
        w.add(level(w, size));
        std::fprintf(stderr, "size %d: %zu terms kept, %zu merged\n", size, w.pool.size(), w.ghosts.size());
    }
    return rsc::laws::report(w, K, max_size, derivations);
}

// Laws read back, one a line on stdin -- an id, the left side, the relation, the right side, in
// prefix form, tab-separated -- and evaluated on the named spaces alone: for each, an id and on
// which of those spaces it holds and is defined. What `discover` found once is checked again column
// by column, so a change pays for the spaces it touches, not for the whole enumeration.
int check(const std::vector<std::string>& only) {
    World w;
    w.tol = TOL;
    all_spaces(w);
    std::erase_if(w.spaces, [&](const auto& s) { return std::find(only.begin(), only.end(), s->name) == only.end(); });
    if (w.spaces.size() != only.size()) {
        std::fprintf(stderr, "rsc_laws check: a space named is not one of the matrix's\n");
        return 2;
    }
    w.key_spaces.clear();
    atoms(w, false);
    std::map<std::string, int> known;
    for (int i = 0; i < int(w.pool.size()); ++i) known[w.pool[i].sexp] = i;

    // (op args...) read into the pool, each subterm once
    std::function<int(const std::string&, std::size_t&)> parse = [&](const std::string& src, std::size_t& at) -> int {
        while (at < src.size() && src[at] == ' ') ++at;
        std::size_t start = at;
        if (src[at] != '(') {
            while (at < src.size() && src[at] != ' ' && src[at] != ')') ++at;
            const auto it = known.find(src.substr(start, at - start));
            if (it == known.end()) throw std::runtime_error("unknown atom " + src.substr(start, at - start));
            return it->second;
        }
        ++at;
        std::size_t e = at;
        while (src[e] != ' ') ++e;
        const std::string op = src.substr(at, e - at);
        at = e;
        double k = 0;
        if (op == "scale") {
            while (src[at] == ' ') ++at;
            std::size_t f = at;
            while (src[f] != ' ') ++f;
            const std::string num = src.substr(at, f - at);
            k = num == "1/2" ? 0.5 : std::stod(num);
            at = f;
        }
        std::vector<int> args;
        for (;;) {
            while (src[at] == ' ') ++at;
            if (src[at] == ')') { ++at; break; }
            args.push_back(parse(src, at));
        }
        const std::string form = src.substr(start, at - start);
        if (auto it = known.find(form); it != known.end()) return it->second;
        Term t{};
        t.k = k;
        t.text = form;
        t.sexp = form;
        auto A = [&](int i) { return i < int(args.size()) ? args[i] : -1; };
        t.a = A(0); t.b = A(1); t.c = A(2);
        if (op == "log") { t.op = Op::Log; t.sort = Sort::V; t.base = t.a; }
        else if (op == "exp") { t.op = Op::Exp; t.sort = Sort::P; }
        else if (op == "d") { t.op = Op::Dist; t.sort = Sort::R; }
        else if (op == "g") { t.op = Op::G; t.sort = Sort::R; }
        else if (op == "scale") { t.op = Op::Scale; t.sort = Sort::V; t.base = w.pool[t.a].base; }
        else if (op == "neg") { t.op = Op::Neg; t.sort = Sort::V; t.base = w.pool[t.a].base; }
        else if (op == "vadd") { t.op = Op::VAdd; t.sort = Sort::V; t.base = w.pool[t.a].base; }
        else if (op == "add") { t.op = Op::RAdd; t.sort = Sort::R; }
        else if (op == "mul") { t.op = Op::RMul; t.sort = Sort::R; }
        else throw std::runtime_error("unknown operation " + op);
        w.add({t}, false);
        known[form] = int(w.pool.size()) - 1;
        return int(w.pool.size()) - 1;
    };

    return rsc::laws::check(w, [&](const std::string& side) {
        std::size_t at = 0;
        return parse(side, at);
    });
}

}  // namespace laws

// rsc_laws [max size]            discover: enumerate, print the laws as JSON
// rsc_laws check space...        check: laws on stdin, evaluated on those spaces alone
// rsc_laws missing               the arithmetic each space's tangent lacks, a line each, no enumeration
int main(int argc, char** argv) {
    if (argc > 1 && std::string(argv[1]) == "check") return laws::check(std::vector<std::string>(argv + 2, argv + argc));
    if (argc > 1 && std::string(argv[1]) == "missing") {
        laws::World w;
        laws::all_spaces(w);
        for (auto& s : w.spaces)
            for (auto& m : s->missing()) std::printf("%s\t%s\n", s->name.c_str(), m.c_str());
        return 0;
    }
    return laws::run(argc > 1 ? std::atoi(argv[1]) : 6, argc > 2 ? argv[2] : nullptr);
}
