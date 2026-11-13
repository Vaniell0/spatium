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
// What it prints is JSON, which rsc_laws.py reads: the spaces, the laws with the spaces each holds
// on and is defined on, and the operations a space's types do not have.
//
//   rsc_laws [max size]
#include "probes.hpp"

#include <algorithm>
#include <atomic>
#include <bit>
#include <cstdint>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <map>
#include <memory>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <variant>
#include <vector>

namespace laws {

using namespace connectivity;

enum class Op { VarP, VarV, Log, Exp, G, Dist, Scale, Neg, VAdd, RAdd, RMul, RConst };
enum class Sort { R, P, V };

struct Term {
    Op op;
    int a = -1, b = -1, c = -1;
    double k = 0;
    Sort sort;
    int base = -1;   // for V: the P term it is tangent at
    int size = 1;
    std::string text;
    int deg = 0;
    std::string sexp;   // the same term in prefix form, which `check` reads back
};

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

// A term's values on one space: its coordinates on each tuple, empty where it is undefined there.
struct Rows {
    std::vector<std::vector<double>> at;
    int count = 0;
    // Defined on a space: on at least half of its tuples.
    bool defined() const { return count * 2 >= K; }
};

// Two terms compared on the tuples where both are defined, a quarter of them at least.
template<class Cmp> bool all_common(const Rows& x, const Rows& y, Cmp cmp) {
    if (!x.defined() || !y.defined()) return false;
    int common = 0;
    for (int t = 0; t < K; ++t) {
        const auto &a = x.at[t], &b = y.at[t];
        if (a.empty() || b.empty()) continue;
        if (a.size() != b.size() || !cmp(a, b)) return false;
        ++common;
    }
    return common * 4 >= K;
}
bool rows_eq(const Rows& x, const Rows& y) {
    return all_common(x, y, [](const auto& a, const auto& b) {
        for (std::size_t m = 0; m < a.size(); ++m)
            if (std::abs(a[m] - b[m]) > TOL * (1 + std::abs(a[m]))) return false;
        return true;
    });
}
bool rows_le(const Rows& x, const Rows& y) {
    return all_common(x, y, [](const auto& a, const auto& b) { return a[0] <= b[0] + TOL * (1 + std::abs(b[0])); });
}
// What merge candidates share, read on the first two tuples.
double rows_key(const Rows& x) {
    if (!x.defined() || x.at[0].empty() || x.at[1].empty()) return std::numeric_limits<double>::quiet_NaN();
    double h = 0;
    for (int t = 0; t < 2; ++t)
        for (double y : x.at[t]) h = h * 1.618 + y;
    return h;
}

struct Space {
    std::string name;
    std::vector<Rows> pool, staged;
    virtual ~Space() = default;
    // Evaluates a batch of terms whose arguments are in the pool, into `staged`.
    virtual void stage(const std::vector<Term>& batch) = 0;
    // Moves one staged term to the end of the pool.
    virtual void accept(std::size_t i) = 0;
    virtual void clear_stage() = 0;
    // The arithmetic a tangent of this space does not have.
    virtual std::vector<std::string> missing() const = 0;
};

template<class Spc>
struct Impl final : Space {
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

struct World {
    std::vector<std::unique_ptr<Space>> spaces;
    std::vector<Term> pool;
    std::vector<std::pair<Term, int>> ghosts;   // each merged term and the kept term it equals
    std::unordered_map<std::string, std::vector<int>> buckets;

    // Merge candidates share a key read on spaces whose operations are closed forms, where a
    // term is defined on every tuple or on none.
    std::vector<int> key_spaces;
    std::string key(const std::vector<const Rows*>& rows, const Term& t) const {
        std::string k = std::to_string(int(t.sort)) + ":" + std::to_string(t.base);
        for (int i : key_spaces) {
            const double h = rows_key(*rows[i]);
            char buf[32];
            if (std::isnan(h)) std::snprintf(buf, sizeof buf, "|u");
            else std::snprintf(buf, sizeof buf, "|%.3g", h);
            k += buf;
        }
        return k;
    }

    // Adds a batch of terms built from the pool: evaluated on every space at once, then each kept
    // unless it equals, wherever both are defined, a term already kept.
    void add(std::vector<Term> batch, bool merge = true) {
        std::vector<std::thread> workers;
        std::atomic<std::size_t> next{0};
        const unsigned n = std::max(1u, std::min<unsigned>(std::thread::hardware_concurrency(), unsigned(spaces.size())));
        for (unsigned w = 0; w < n; ++w)
            workers.emplace_back([&] {
                for (std::size_t i; (i = next++) < spaces.size();) spaces[i]->stage(batch);
            });
        for (auto& w : workers) w.join();
        for (std::size_t i = 0; i < batch.size(); ++i) {
            std::vector<const Rows*> rows;
            bool any = false;
            for (auto& s : spaces) { rows.push_back(&s->staged[i]); any = any || s->staged[i].defined(); }
            if (!any && merge) continue;
            auto& bucket = buckets[key(rows, batch[i])];
            int same = -1;
            for (int o : bucket) {
                bool met = false, differ = false;
                for (std::size_t k = 0; k < spaces.size() && !differ; ++k) {
                    const Rows &a = spaces[k]->pool[o], &b = *rows[k];
                    if (!(a.defined() && b.defined())) continue;
                    if (!rows_eq(a, b)) differ = true;
                    met = true;
                }
                if (met && !differ) { same = o; break; }
            }
            if (same >= 0 && merge) { ghosts.emplace_back(batch[i], same); continue; }
            for (auto& s : spaces) s->accept(i);
            bucket.push_back(int(pool.size()));
            pool.push_back(std::move(batch[i]));
        }
        for (auto& s : spaces) s->clear_stage();
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

std::string json_string(const std::string& s) {
    std::string out = "\"";
    for (char ch : s) {
        if (ch == '"' || ch == '\\') out += '\\';
        out += ch;
    }
    return out + "\"";
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

int run(int max_size) {
    World w;
    all_spaces(w);
    const int ns = int(w.spaces.size());
    atoms(w, true);
    for (int size = 2; size <= max_size; ++size) {
        w.add(level(w, size));
        std::fprintf(stderr, "size %d: %zu terms kept, %zu merged\n", size, w.pool.size(), w.ghosts.size());
    }
    const int n = int(w.pool.size());

    struct Law { int a, b; bool le; std::uint64_t hold, def; int size; };
    using Mask = std::uint64_t;
    std::vector<Law> found;
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) {
            const auto &A = w.pool[i], &B = w.pool[j];
            if (A.sort != B.sort || A.base != B.base) continue;
            if (!space_ops(A.text) && !space_ops(B.text)) continue;
            const bool zero = (A.op == Op::RConst && A.k == 0) || (B.op == Op::RConst && B.k == 0);
            if (A.sort == Sort::R && A.deg != B.deg && !zero) continue;
            Mask def = 0, eq = 0, l = 0, g = 0;
            for (int k = 0; k < ns; ++k) {
                const Rows &x = w.spaces[k]->pool[i], &y = w.spaces[k]->pool[j];
                if (!(x.defined() && y.defined())) continue;
                def |= Mask{1} << k;
                if (rows_eq(x, y)) eq |= Mask{1} << k;
                if (A.sort == Sort::R) {
                    if (rows_le(x, y)) l |= Mask{1} << k;
                    if (rows_le(y, x)) g |= Mask{1} << k;
                }
            }
            const int sz = A.size + B.size;
            if (eq) found.push_back({i, j, false, eq, def, sz});
            if (l & ~eq) found.push_back({i, j, true, l, def, sz});
            if (g & ~eq) found.push_back({j, i, true, g, def, sz});
        }
    std::stable_sort(found.begin(), found.end(), [](const Law& x, const Law& y) { return x.size < y.size; });

    // A kept law matched onto a candidate: its variables bound to terms, a tangent's point with it.
    std::function<bool(int, int, std::map<int, int>&)> match = [&](int pat, int tgt, std::map<int, int>& sub) -> bool {
        const auto &Pt = w.pool[pat], &Tt = w.pool[tgt];
        if (Pt.sort != Tt.sort) return false;
        if (Pt.op == Op::VarP || Pt.op == Op::VarV) {
            if (Pt.op == Op::VarV) {
                auto it = sub.find(Pt.base);
                if (it == sub.end()) sub[Pt.base] = Tt.base;
                else if (it->second != Tt.base) return false;
            }
            auto it = sub.find(pat);
            if (it == sub.end()) { sub[pat] = tgt; return true; }
            return it->second == tgt;
        }
        if (Pt.op != Tt.op || Pt.k != Tt.k) return false;
        if (Pt.op == Op::RConst) return pat == tgt;
        for (auto [x, y] : {std::pair{Pt.a, Tt.a}, std::pair{Pt.b, Tt.b}, std::pair{Pt.c, Tt.c}}) {
            if ((x < 0) != (y < 0)) return false;
            if (x >= 0 && !match(x, y, sub)) return false;
        }
        return true;
    };
    auto instance_of = [&](const Law& k, const Law& l) {
        for (int swap = 0; swap < (k.le ? 1 : 2); ++swap) {
            std::map<int, int> sub;
            if (match(swap ? k.b : k.a, l.a, sub) && match(swap ? k.a : k.b, l.b, sub)) return true;
        }
        return false;
    };
    struct Closure {
        std::vector<int> up;
        int find(int x) { while (up[x] != x) x = up[x] = up[up[x]]; return x; }
        void join(int x, int y) { x = find(x); y = find(y); if (x != y) up[std::max(x, y)] = std::min(x, y); }
    };
    std::map<Mask, Closure> closures;
    auto close = [&](Closure& c) {
        for (bool changed = true; changed;) {
            changed = false;
            std::map<std::tuple<int, double, int, int, int>, int> seen;
            for (int t = 0; t < n; ++t) {
                const auto& T = w.pool[t];
                if (T.a < 0) continue;
                const auto k = std::tuple{int(T.op), T.k, c.find(T.a), T.b < 0 ? -1 : c.find(T.b), T.c < 0 ? -1 : c.find(T.c)};
                auto [it, fresh] = seen.emplace(k, t);
                if (!fresh && c.find(it->second) != c.find(t)) { c.join(it->second, t); changed = true; }
            }
        }
    };

    std::vector<Law> kept;
    // The merges are laws of every space: their instances are no news either.
    for (auto& [t, o] : w.ghosts) {
        w.pool.push_back(t);
        kept.push_back({int(w.pool.size()) - 1, o, false, ~Mask{0}, ~Mask{0}, t.size + w.pool[o].size});
    }
    const std::size_t seeded = kept.size();
    for (const auto& L : found) {
        bool redundant = false;
        for (const auto& K_ : kept)
            if (K_.size < L.size && (K_.hold & L.hold) == L.hold && K_.le == L.le && instance_of(K_, L)) { redundant = true; break; }
        if (redundant) continue;
        auto [it, fresh] = closures.try_emplace(L.hold);
        auto& c = it->second;
        if (fresh) {
            c.up.resize(n);
            for (int t = 0; t < n; ++t) c.up[t] = t;
            for (const auto& K_ : kept)
                if (!K_.le && K_.a < n && K_.b < n && (K_.hold & L.hold) == L.hold) c.join(K_.a, K_.b);
        }
        close(c);
        if (c.find(L.a) == c.find(L.b)) continue;
        if (L.le) {
            // the same pair up to the equalities, one term on both sides up to a constant factor,
            // or 0 below what is known not to be negative
            auto strip = [&](int t) {
                for (;;) {
                    const auto& T = w.pool[t];
                    if (T.op != Op::RMul || T.a == T.b) break;
                    if (w.pool[T.a].op == Op::RConst) t = T.b;
                    else if (w.pool[T.b].op == Op::RConst) t = T.a;
                    else break;
                }
                return c.find(t);
            };
            std::function<bool(int)> nonneg = [&](int t) -> bool {
                const auto& T = w.pool[t];
                if (T.op == Op::RConst) return T.k >= 0;
                if (T.op == Op::RMul && c.find(T.a) == c.find(T.b)) return true;
                if (T.op == Op::RMul && w.pool[T.a].op == Op::RConst) return w.pool[T.a].k >= 0 && nonneg(T.b);
                if (T.op == Op::RMul && w.pool[T.b].op == Op::RConst) return w.pool[T.b].k >= 0 && nonneg(T.a);
                if (T.op == Op::RAdd) return nonneg(T.a) && nonneg(T.b);
                for (const auto& K_ : kept)
                    if (K_.le && (K_.hold & L.hold) == L.hold && w.pool[K_.a].op == Op::RConst && w.pool[K_.a].k == 0 &&
                        K_.b < n && c.find(K_.b) == c.find(t))
                        return true;
                return false;
            };
            bool follows = strip(L.a) == strip(L.b);
            if (!follows && w.pool[L.a].op == Op::RConst && w.pool[L.a].k == 0) follows = nonneg(L.b);
            for (const auto& K_ : kept)
                if (!follows && K_.le && K_.a < n && K_.b < n && (K_.hold & L.hold) == L.hold &&
                    c.find(K_.a) == c.find(L.a) && c.find(K_.b) == c.find(L.b))
                    follows = true;
            if (follows) continue;
        }
        kept.push_back(L);
        if (!L.le)
            for (auto& [m, cl] : closures) if ((L.hold & m) == m) cl.join(L.a, L.b);
    }

    auto bits = [&](Mask m) {
        std::string s;
        for (int k = 0; k < ns; ++k) s += (m >> k & 1) ? '1' : '0';
        return s;
    };
    std::printf("{\"tuples\": %d, \"max_size\": %d, \"terms\": %d, \"merged\": %zu,\n \"spaces\": [", K, max_size, n,
                w.ghosts.size());
    for (int k = 0; k < ns; ++k) std::printf("%s%s", k ? ", " : "", json_string(w.spaces[k]->name).c_str());
    std::printf("],\n \"missing\": [");
    bool first_missing = true;
    for (auto& s : w.spaces)
        for (auto& m : s->missing()) {
            std::printf("%s\n  {\"space\": %s, \"operation\": %s}", first_missing ? "" : ",", json_string(s->name).c_str(),
                        json_string(m).c_str());
            first_missing = false;
        }
    std::printf("],\n \"laws\": [");
    bool printed = false;
    for (std::size_t i = 0; i < kept.size(); ++i) {
        const auto& L = kept[i];
        const bool everywhere = i < seeded;
        if (everywhere && !space_ops(w.pool[L.a].text) && !space_ops(w.pool[L.b].text)) continue;
        std::printf("%s\n  {\"lhs\": %s, \"rel\": \"%s\", \"rhs\": %s, \"size\": %d, \"hold\": \"%s\", \"defined\": \"%s\",\n"
                    "   \"lhs_form\": %s, \"rhs_form\": %s}",
                    printed ? "," : "", json_string(w.pool[L.a].text).c_str(), L.le ? "<=" : "=", json_string(w.pool[L.b].text).c_str(),
                    L.size, everywhere ? "all" : bits(L.hold).c_str(), everywhere ? "all" : bits(L.def).c_str(),
                    json_string(w.pool[L.a].sexp).c_str(), json_string(w.pool[L.b].sexp).c_str());
        printed = true;
    }
    std::printf("\n]}\n");
    return 0;
}

// Laws read back, one a line on stdin -- an id, the left side, the relation, the right side, in
// prefix form, tab-separated -- and evaluated on the named spaces alone: for each, an id and on
// which of those spaces it holds and is defined. What `discover` found once is checked again column
// by column, so a change pays for the spaces it touches, not for the whole enumeration.
int check(const std::vector<std::string>& only) {
    World w;
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

    std::string line;
    while (std::getline(std::cin, line)) {
        std::vector<std::string> f;
        std::size_t from = 0;
        for (std::size_t tab; (tab = line.find('\t', from)) != std::string::npos; from = tab + 1) f.push_back(line.substr(from, tab - from));
        f.push_back(line.substr(from));
        if (f.size() != 4) continue;
        std::size_t at = 0;
        const int a = parse(f[1], at);
        at = 0;
        const int b = parse(f[3], at);
        std::string hold, def;
        for (auto& s : w.spaces) {
            const Rows &x = s->pool[a], &y = s->pool[b];
            const bool d = x.defined() && y.defined();
            const bool h = d && (f[2] == "=" ? (a == b || rows_eq(x, y)) : rows_le(x, y));
            hold += h ? '1' : '0';
            def += d ? '1' : '0';
        }
        std::printf("%s\t%s\t%s\n", f[0].c_str(), hold.c_str(), def.c_str());
    }
    return 0;
}

}  // namespace laws

// rsc_laws [max size]            discover: enumerate, print the laws as JSON
// rsc_laws check space...        check: laws on stdin, evaluated on those spaces alone
int main(int argc, char** argv) {
    if (argc > 1 && std::string(argv[1]) == "check") return laws::check(std::vector<std::string>(argv + 2, argv + argc));
    return laws::run(argc > 1 ? std::atoi(argv[1]) : 6);
}
