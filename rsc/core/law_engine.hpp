// Laws found by enumeration: terms over a project's operations, evaluated on every space of a matrix, equal ones
// merged, and what is left read off as laws -- each with the spaces it holds on and is defined on. What a project
// supplies is a signature (its operations and sorts), its spaces (how a batch of terms is evaluated on one), and the
// enumeration of terms by size; the rest is here and knows nothing of what a term means.
//
// A signature `Sig` gives the term's two enumerations and says what the engine may assume of them:
//
//   using Op = ...; using Sort = ...;             enumerations of the operations and the sorts
//   static bool is_var(Op);                       a variable: matching binds it (with its base, when it has one)
//   static bool is_const(Op);                     a constant, its value in Term::k
//   static bool is_mul(Op), is_add(Op);           the product and sum of the ordered sort, for the order's pruning
//   static bool ordered(Sort);                    values of this sort are numbers, so <= is asked of them too
//   static bool interesting(const Term<Sig>&);    a law is reported only when one side is
//
// The steps: World::add evaluates a batch of terms on every space at once and keeps each term unless it equals,
// wherever both are defined, one already kept (the merge is a law of every space); report() then takes every pair of
// kept terms equal or ordered on some spaces, drops a law that is an instance of a smaller one or follows from
// smaller ones by equality, congruence or a constant factor, and prints the rest as JSON.
#pragma once

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace rsc::laws {

template<class Sig> struct Term {
    typename Sig::Op op;
    int a = -1, b = -1, c = -1;
    double k = 0;
    typename Sig::Sort sort;
    int base = -1;   // for a sort that lives at a point: the term it is at
    int size = 1;
    std::string text;
    int deg = 0;
    std::string sexp;   // the same term in prefix form, which `check` reads back
};

// A term's values on one space: its coordinates on each tuple, empty where it is undefined there.
struct Rows {
    std::vector<std::vector<double>> at;
    int count = 0;
    // Defined on a space: on at least half of its tuples.
    bool defined() const { return !at.empty() && count * 2 >= int(at.size()); }
};

// Two terms compared on the tuples where both are defined, a quarter of them at least.
template<class Cmp> bool all_common(const Rows& x, const Rows& y, Cmp cmp) {
    if (!x.defined() || !y.defined()) return false;
    const int n = int(x.at.size());
    int common = 0;
    for (int t = 0; t < n; ++t) {
        const auto &a = x.at[t], &b = y.at[t];
        if (a.empty() || b.empty()) continue;
        if (a.size() != b.size() || !cmp(a, b)) return false;
        ++common;
    }
    return common * 4 >= n;
}
inline bool rows_eq(const Rows& x, const Rows& y, double tol) {
    return all_common(x, y, [tol](const auto& a, const auto& b) {
        for (std::size_t m = 0; m < a.size(); ++m)
            if (std::abs(a[m] - b[m]) > tol * (1 + std::abs(a[m]))) return false;
        return true;
    });
}
inline bool rows_le(const Rows& x, const Rows& y, double tol) {
    return all_common(x, y, [tol](const auto& a, const auto& b) { return a[0] <= b[0] + tol * (1 + std::abs(b[0])); });
}
// What merge candidates share, read on the first two tuples.
inline double rows_key(const Rows& x) {
    if (!x.defined() || x.at[0].empty() || x.at[1].empty()) return std::numeric_limits<double>::quiet_NaN();
    double h = 0;
    for (int t = 0; t < 2; ++t)
        for (double y : x.at[t]) h = h * 1.618 + y;
    return h;
}

template<class Sig> struct Space {
    std::string name;
    std::vector<Rows> pool, staged;
    virtual ~Space() = default;
    // Evaluates a batch of terms whose arguments are in the pool, into `staged`.
    virtual void stage(const std::vector<Term<Sig>>& batch) = 0;
    // Moves one staged term to the end of the pool.
    virtual void accept(std::size_t i) = 0;
    virtual void clear_stage() = 0;
    // The operations this space lacks, reported beside the laws.
    virtual std::vector<std::string> missing() const { return {}; }
};

template<class Sig> struct World {
    using T = Term<Sig>;
    double tol = 1e-6;
    std::vector<std::unique_ptr<Space<Sig>>> spaces;
    std::vector<T> pool;
    std::vector<std::pair<T, int>> ghosts;   // each merged term and the kept term it equals
    std::unordered_map<std::string, std::vector<int>> buckets;

    // Merge candidates share a key read on spaces whose operations are closed forms, where a
    // term is defined on every tuple or on none.
    std::vector<int> key_spaces;
    std::string key(const std::vector<const Rows*>& rows, const T& t) const {
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
    void add(std::vector<T> batch, bool merge = true) {
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
                    if (!rows_eq(a, b, tol)) differ = true;
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

inline std::string json_string(const std::string& s) {
    std::string out = "\"";
    for (char ch : s) {
        if (ch == '"' || ch == '\\') out += '\\';
        out += ch;
    }
    return out + "\"";
}

// The laws of an enumerated world, printed as JSON on stdout: the spaces, the operations they lack, and every law
// kept with the spaces it holds and is defined on. `tuples` and `max_size` are what the enumeration used. With
// `derivations`, each law the closure already held is written there as a question with its chain.
template<class Sig> int report(World<Sig>& w, int tuples, int max_size, const char* derivations = nullptr) {
    using Mask = std::uint64_t;
    const int ns = int(w.spaces.size());
    if (ns > 64) throw std::runtime_error("rsc::laws::report: more than 64 spaces");
    const int n = int(w.pool.size());
    const double tol = w.tol;
    auto is_zero = [&](int i) { return Sig::is_const(w.pool[i].op) && w.pool[i].k == 0; };

    struct Law { int a, b; bool le; Mask hold, def; int size; };
    std::vector<Law> found;
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) {
            const auto &A = w.pool[i], &B = w.pool[j];
            if (A.sort != B.sort || A.base != B.base) continue;
            if (!Sig::interesting(A) && !Sig::interesting(B)) continue;
            const bool zero = is_zero(i) || is_zero(j);
            if (Sig::ordered(A.sort) && A.deg != B.deg && !zero) continue;
            Mask def = 0, eq = 0, l = 0, g = 0;
            for (int k = 0; k < ns; ++k) {
                const Rows &x = w.spaces[k]->pool[i], &y = w.spaces[k]->pool[j];
                if (!(x.defined() && y.defined())) continue;
                def |= Mask{1} << k;
                if (rows_eq(x, y, tol)) eq |= Mask{1} << k;
                if (Sig::ordered(A.sort)) {
                    if (rows_le(x, y, tol)) l |= Mask{1} << k;
                    if (rows_le(y, x, tol)) g |= Mask{1} << k;
                }
            }
            const int sz = A.size + B.size;
            if (eq) found.push_back({i, j, false, eq, def, sz});
            if (l & ~eq) found.push_back({i, j, true, l, def, sz});
            if (g & ~eq) found.push_back({j, i, true, g, def, sz});
        }
    std::stable_sort(found.begin(), found.end(), [](const Law& x, const Law& y) { return x.size < y.size; });

    // A kept law matched onto a candidate: its variables bound to terms, a variable's base with it.
    std::function<bool(int, int, std::map<int, int>&)> match = [&](int pat, int tgt, std::map<int, int>& sub) -> bool {
        const auto &Pt = w.pool[pat], &Tt = w.pool[tgt];
        if (Pt.sort != Tt.sort) return false;
        if (Sig::is_var(Pt.op)) {
            if (Pt.base >= 0) {
                auto it = sub.find(Pt.base);
                if (it == sub.end()) sub[Pt.base] = Tt.base;
                else if (it->second != Tt.base) return false;
            }
            auto it = sub.find(pat);
            if (it == sub.end()) { sub[pat] = tgt; return true; }
            return it->second == tgt;
        }
        if (Pt.op != Tt.op || Pt.k != Tt.k) return false;
        if (Sig::is_const(Pt.op)) return pat == tgt;
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
    // why two terms were joined: a kept law (its index), or congruence (-1): the same operation on
    // arguments already equal
    struct Edge { int x, y, why; int ax = -1, ay = -1; };   // congruence: the arguments already equal
    struct Closure {
        std::vector<int> up;
        std::vector<Edge> edges;
        int find(int x) { while (up[x] != x) x = up[x] = up[up[x]]; return x; }
        void join(int x, int y, int why = -1, int ax = -1, int ay = -1) {
            const int rx = find(x), ry = find(y);
            if (rx == ry) return;
            edges.push_back({x, y, why, ax, ay});
            up[std::max(rx, ry)] = std::min(rx, ry);
        }
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
                if (!fresh && c.find(it->second) != c.find(t)) {
                    // the first argument the two differ in, already equal: what makes them equal
                    const auto& S = w.pool[it->second];
                    int ax = -1, ay = -1;
                    for (auto [x, y] : {std::pair{S.a, T.a}, std::pair{S.b, T.b}, std::pair{S.c, T.c}})
                        if (x >= 0 && x != y) { ax = x, ay = y; break; }
                    c.join(it->second, t, -1, ax, ay);
                    changed = true;
                }
            }
        }
    };

    // A law the closure already holds, with the chain that shows it: each hop an equality, by a kept
    // law or by congruence. `instance` is the kept law it is an instance of, in one hop.
    struct Derived { Law law; std::vector<Edge> chain; int instance; };
    std::vector<Derived> derived;
    // the hops from a to b through the joins that made them one: a breadth-first search over the
    // recorded edges, each taken in the direction walked
    auto path = [&](const Closure& c, int a, int b) {
        std::map<int, std::vector<std::pair<int, int>>> adj;   // term -> (other term, edge index)
        for (int e = 0; e < int(c.edges.size()); ++e) {
            adj[c.edges[e].x].push_back({c.edges[e].y, e});
            adj[c.edges[e].y].push_back({c.edges[e].x, e});
        }
        std::map<int, std::pair<int, int>> from;   // term -> (previous term, edge)
        std::vector<int> todo{a};
        from[a] = {-1, -1};
        for (std::size_t i = 0; i < todo.size() && !from.count(b); ++i)
            for (auto [y, e] : adj[todo[i]])
                if (!from.count(y)) from[y] = {todo[i], e}, todo.push_back(y);
        std::vector<Edge> hops;
        if (!from.count(b)) return hops;
        for (int t = b; t != a; t = from[t].first) {
            const Edge& e = c.edges[from[t].second];
            const bool forward = e.x == from[t].first;
            hops.push_back({from[t].first, t, e.why, forward ? e.ax : e.ay, forward ? e.ay : e.ax});
        }
        std::reverse(hops.begin(), hops.end());
        return hops;
    };
    std::vector<Law> kept;
    // The merges are laws of every space: their instances are no news either.
    for (auto& [t, o] : w.ghosts) {
        w.pool.push_back(t);
        kept.push_back({int(w.pool.size()) - 1, o, false, ~Mask{0}, ~Mask{0}, t.size + w.pool[o].size});
    }
    const std::size_t seeded = kept.size();
    for (const auto& L : found) {
        int instance = -1;
        for (std::size_t q = 0; q < kept.size(); ++q) {
            const auto& K_ = kept[q];
            if (K_.size < L.size && (K_.hold & L.hold) == L.hold && K_.le == L.le && instance_of(K_, L)) { instance = int(q); break; }
        }
        if (instance >= 0) {
            if (!L.le) derived.push_back({L, {{L.a, L.b, instance}}, instance});
            continue;
        }
        auto [it, fresh] = closures.try_emplace(L.hold);
        auto& c = it->second;
        if (fresh) {
            c.up.resize(n);
            for (int t = 0; t < n; ++t) c.up[t] = t;
            for (std::size_t q = 0; q < kept.size(); ++q) {
                const auto& K_ = kept[q];
                if (!K_.le && K_.a < n && K_.b < n && (K_.hold & L.hold) == L.hold) c.join(K_.a, K_.b, int(q));
            }
        }
        close(c);
        if (c.find(L.a) == c.find(L.b)) {
            if (!L.le && derivations) derived.push_back({L, path(c, L.a, L.b), -1});
            continue;
        }
        if (L.le) {
            // the same pair up to the equalities, one term on both sides up to a constant factor,
            // or 0 below what is known not to be negative
            auto strip = [&](int t) {
                for (;;) {
                    const auto& T = w.pool[t];
                    if (!Sig::is_mul(T.op) || T.a == T.b) break;
                    if (Sig::is_const(w.pool[T.a].op)) t = T.b;
                    else if (Sig::is_const(w.pool[T.b].op)) t = T.a;
                    else break;
                }
                return c.find(t);
            };
            std::function<bool(int)> nonneg = [&](int t) -> bool {
                const auto& T = w.pool[t];
                if (Sig::is_const(T.op)) return T.k >= 0;
                if (Sig::is_mul(T.op) && c.find(T.a) == c.find(T.b)) return true;
                if (Sig::is_mul(T.op) && Sig::is_const(w.pool[T.a].op)) return w.pool[T.a].k >= 0 && nonneg(T.b);
                if (Sig::is_mul(T.op) && Sig::is_const(w.pool[T.b].op)) return w.pool[T.b].k >= 0 && nonneg(T.a);
                if (Sig::is_add(T.op)) return nonneg(T.a) && nonneg(T.b);
                for (const auto& K_ : kept)
                    if (K_.le && (K_.hold & L.hold) == L.hold && is_zero(K_.a) && K_.b < n && c.find(K_.b) == c.find(t))
                        return true;
                return false;
            };
            bool follows = strip(L.a) == strip(L.b);
            if (!follows && is_zero(L.a)) follows = nonneg(L.b);
            for (const auto& K_ : kept)
                if (!follows && K_.le && K_.a < n && K_.b < n && (K_.hold & L.hold) == L.hold &&
                    c.find(K_.a) == c.find(L.a) && c.find(K_.b) == c.find(L.b))
                    follows = true;
            if (follows) continue;
        }
        kept.push_back(L);
        if (!L.le)
            for (auto& [m, cl] : closures) if ((L.hold & m) == m) cl.join(L.a, L.b, int(kept.size()) - 1);
    }

    auto bits = [&](Mask m) {
        std::string s;
        for (int k = 0; k < ns; ++k) s += (m >> k & 1) ? '1' : '0';
        return s;
    };
    // Each law the closure already held, as a question with its derivation: the premises its chain
    // uses and a few others true on the same spaces, the equality asked, the chain an equality a hop,
    // its depth, what a one-hop rule would answer, and near misses -- terms of the right side's form
    // that a tuple of these spaces shows unequal to the left side.
    if (derivations) {
        std::FILE* out = std::fopen(derivations, "w");
        if (!out) throw std::runtime_error(std::string("cannot write ") + derivations);
        auto eq_text = [&](int x, int y) { return w.pool[x].text + " = " + w.pool[y].text; };
        auto unequal_somewhere = [&](int x, int y, Mask m) {
            for (int k = 0; k < ns; ++k) {
                if (!(m >> k & 1)) continue;
                const Rows &rx = w.spaces[k]->pool[x], &ry = w.spaces[k]->pool[y];
                if (rx.defined() && ry.defined() && !rows_eq(rx, ry, tol)) return true;
            }
            return false;
        };
        std::mt19937 rng(7);
        std::map<int, int> histogram;
        int id = 0;
        for (const auto& D : derived) {
            const auto& L = D.law;
            const int depth = int(D.chain.size());
            if (depth == 0 || L.a >= n || L.b >= n) continue;
            if (!Sig::interesting(w.pool[L.a]) && !Sig::interesting(w.pool[L.b])) continue;
            ++histogram[depth];
            std::vector<int> used;
            for (const auto& h : D.chain)
                if (h.why >= 0 && std::find(used.begin(), used.end(), h.why) == used.end()) used.push_back(h.why);
            std::vector<int> others;
            for (int q = 0; q < int(kept.size()); ++q)
                if (!kept[q].le && kept[q].a < n && kept[q].b < n && (kept[q].hold & L.hold) == L.hold &&
                    std::find(used.begin(), used.end(), q) == used.end() &&
                    (Sig::interesting(w.pool[kept[q].a]) || Sig::interesting(w.pool[kept[q].b])))
                    others.push_back(q);
            std::shuffle(others.begin(), others.end(), rng);
            others.resize(std::min<std::size_t>(others.size(), 4));
            std::vector<int> premises = used;
            premises.insert(premises.end(), others.begin(), others.end());
            std::shuffle(premises.begin(), premises.end(), rng);
            std::vector<int> misses;
            for (int t = 0; t < n && misses.size() < 2; ++t) {
                const auto &T = w.pool[t], &B = w.pool[L.b];
                if (t == L.b || t == L.a || T.sort != B.sort || T.base != B.base || T.op != B.op || std::abs(T.size - B.size) > 1) continue;
                if (unequal_somewhere(L.a, t, L.hold)) misses.push_back(t);
            }
            std::string spaces, prem, chain, miss;
            for (int k = 0; k < ns; ++k)
                if (L.hold >> k & 1) spaces += (spaces.empty() ? "" : ", ") + json_string(w.spaces[k]->name);
            for (int q : premises) prem += (prem.empty() ? "" : ", ") + json_string(eq_text(kept[q].a, kept[q].b));
            for (const auto& h : D.chain)
                chain += std::string(chain.empty() ? "" : ", ") + "{\"from\": " + json_string(w.pool[h.x].text) + ", \"to\": " +
                         json_string(w.pool[h.y].text) + ", \"by\": " +
                         json_string(h.why >= 0 ? eq_text(kept[h.why].a, kept[h.why].b) : std::string("congruence")) +
                         (h.why < 0 && h.ax >= 0 ? ", \"because\": " + json_string(eq_text(h.ax, h.ay)) : std::string()) + "}";
            for (int t : misses) miss += (miss.empty() ? "" : ", ") + json_string(eq_text(L.a, t));
            std::fprintf(out,
                         "{\"id\": %d, \"spaces\": [%s], \"premises\": [%s], \"question\": %s, \"answer\": true, \"depth\": %d, "
                         "\"chain\": [%s], \"shallow\": %s, \"false_candidates\": [%s]}\n",
                         id++, spaces.c_str(), prem.c_str(), json_string(eq_text(L.a, L.b)).c_str(), depth, chain.c_str(),
                         depth == 1 ? "true" : "false", miss.c_str());
        }
        std::fclose(out);
        std::fprintf(stderr, "derivations: %d;", id);
        for (auto [d, c] : histogram) std::fprintf(stderr, " depth %d: %d", d, c);
        std::fprintf(stderr, "\n");
    }
    std::printf("{\"tuples\": %d, \"max_size\": %d, \"terms\": %d, \"merged\": %zu,\n \"spaces\": [", tuples, max_size, n,
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
        if (everywhere && !Sig::interesting(w.pool[L.a]) && !Sig::interesting(w.pool[L.b])) continue;
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

// Laws read back, one a line on stdin -- an id, the left side, the relation, the right side, in prefix form,
// tab-separated -- and evaluated on the world's spaces: for each, an id and on which spaces it holds and is defined.
// `parse` reads one side into the world's pool (the project's own syntax) and returns its index.
template<class Sig> int check(World<Sig>& w, const std::function<int(const std::string&)>& parse) {
    std::string line;
    while (std::getline(std::cin, line)) {
        std::vector<std::string> f;
        std::size_t from = 0;
        for (std::size_t tab; (tab = line.find('\t', from)) != std::string::npos; from = tab + 1) f.push_back(line.substr(from, tab - from));
        f.push_back(line.substr(from));
        if (f.size() != 4) continue;
        const int a = parse(f[1]);
        const int b = parse(f[3]);
        std::string hold, def;
        for (auto& s : w.spaces) {
            const Rows &x = s->pool[a], &y = s->pool[b];
            const bool d = x.defined() && y.defined();
            const bool h = d && (f[2] == "=" ? (a == b || rows_eq(x, y, w.tol)) : rows_le(x, y, w.tol));
            hold += h ? '1' : '0';
            def += d ? '1' : '0';
        }
        std::printf("%s\t%s\t%s\n", f[0].c_str(), hold.c_str(), def.c_str());
    }
    return 0;
}

}  // namespace rsc::laws
