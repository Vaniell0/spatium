#pragma once
// A typed graph of a library, and the question "by what routes can I get this?" asked of it.
//
// Nodes are concepts (the library's interfaces), types (its concrete spaces, shapes, groups) and
// scalars. Edges are
//
//   refines    concept -> concept   `from` is built from `to` (MetricSpace refines TopologicalSpace)
//   satisfies  type -> concept      the compiler said `Concept<Type>` holds
//   operation  one way to get a quantity out of a type: `inputs[0]` is the type, the rest are
//              the concepts that make the way available; `needs` are the quantities the way is
//              built from; `output` is the quantity it yields; `level` says how far the library's
//              own tests have verified it (0 declared, 1 runs and is finite, 2 an axiom or a
//              symmetry holds, 3 its numbers agree across scalars).
//
// A route to a quantity on a type is a set of operation edges that closes: every edge's needs are
// yielded by another edge of the route, and the last yields the goal. Two routes to one quantity
// (a closed form and a derivation from the basis) are two different sets.
//
// This header is the generic core: it knows nothing of the library it describes, includes nothing
// of it, and reads the graph from JSON (rsc/spatium/graph.json is the library's, written by
// rsc/spatium/gen_graph.py). CI greps rsc/core for any include of the library.

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <charconv>
#include <deque>
#include <expected>
#include <filesystem>
#include <fstream>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

namespace rsc {

// ── A JSON reader ──────────────────────────────────────────────
// Only what the graph needs: no comments, no NaN, \uXXXX read as UTF-8, an object keeps its keys in order.

namespace json {

struct Value;
using Array = std::vector<Value>;
using Object = std::vector<std::pair<std::string, Value>>;

struct Value {
    std::variant<std::nullptr_t, bool, double, std::string, Array, Object> v{nullptr};

    bool is_object() const { return std::holds_alternative<Object>(v); }
    bool is_array() const { return std::holds_alternative<Array>(v); }
    bool is_string() const { return std::holds_alternative<std::string>(v); }
    bool is_number() const { return std::holds_alternative<double>(v); }

    const Object& object() const { return std::get<Object>(v); }
    const Array& array() const { return std::get<Array>(v); }
    const std::string& string() const { return std::get<std::string>(v); }
    double number() const { return std::get<double>(v); }

    const Value* get(std::string_view key) const {
        if (!is_object()) return nullptr;
        for (const auto& [k, x] : object())
            if (k == key) return &x;
        return nullptr;
    }
};

class Parser {
public:
    explicit Parser(std::string_view s) : s_(s) {}

    std::expected<Value, std::string> run() {
        auto v = value();
        if (!v) return v;
        skip();
        if (i_ != s_.size()) return fail("trailing characters");
        return v;
    }

private:
    std::string_view s_;
    std::size_t i_ = 0;
    int depth_ = 0;

    std::unexpected<std::string> fail(const std::string& what) const {
        return std::unexpected("json: " + what + " at byte " + std::to_string(i_));
    }
    void skip() {
        while (i_ < s_.size() && (s_[i_] == ' ' || s_[i_] == '\n' || s_[i_] == '\t' || s_[i_] == '\r')) ++i_;
    }
    bool eat(char c) {
        skip();
        if (i_ < s_.size() && s_[i_] == c) { ++i_; return true; }
        return false;
    }
    bool word(std::string_view w) {
        if (s_.substr(i_, w.size()) != w) return false;
        i_ += w.size();
        return true;
    }

    std::expected<Value, std::string> value() {
        skip();
        if (i_ >= s_.size()) return fail("unexpected end");
        if (++depth_ > 64) return fail("nested too deep");
        struct Leave { int& d; ~Leave() { --d; } } leave{depth_};
        const char c = s_[i_];
        if (c == '{') return object();
        if (c == '[') return array();
        if (c == '"') {
            auto str = string();
            if (!str) return std::unexpected(str.error());
            return Value{std::move(*str)};
        }
        if (word("true")) return Value{true};
        if (word("false")) return Value{false};
        if (word("null")) return Value{nullptr};
        return number();
    }

    std::expected<Value, std::string> number() {
        const std::size_t start = i_;
        while (i_ < s_.size() && (std::string_view("+-0123456789.eE").find(s_[i_]) != std::string_view::npos)) ++i_;
        double x = 0;
        const auto [p, ec] = std::from_chars(s_.data() + start, s_.data() + i_, x);
        if (start == i_ || ec != std::errc{} || p != s_.data() + i_) { i_ = start; return fail("bad value"); }
        return Value{x};
    }

    std::expected<std::string, std::string> string() {
        ++i_;   // the opening quote
        std::string out;
        while (i_ < s_.size() && s_[i_] != '"') {
            char c = s_[i_++];
            if (c != '\\') { out += c; continue; }
            if (i_ >= s_.size()) break;
            c = s_[i_++];
            switch (c) {
                case '"': case '\\': case '/': out += c; break;
                case 'b': out += '\b'; break;
                case 'f': out += '\f'; break;
                case 'n': out += '\n'; break;
                case 'r': out += '\r'; break;
                case 't': out += '\t'; break;
                case 'u': {
                    if (i_ + 4 > s_.size()) return fail("short \\u escape");
                    unsigned cp = 0;
                    const auto [p, ec] = std::from_chars(s_.data() + i_, s_.data() + i_ + 4, cp, 16);
                    if (ec != std::errc{} || p != s_.data() + i_ + 4) return fail("bad \\u escape");
                    i_ += 4;
                    if (cp < 0x80) out += static_cast<char>(cp);
                    else if (cp < 0x800) { out += static_cast<char>(0xC0 | (cp >> 6)); out += static_cast<char>(0x80 | (cp & 0x3F)); }
                    else { out += static_cast<char>(0xE0 | (cp >> 12)); out += static_cast<char>(0x80 | ((cp >> 6) & 0x3F)); out += static_cast<char>(0x80 | (cp & 0x3F)); }
                    break;
                }
                default: return fail("bad escape");
            }
        }
        if (i_ >= s_.size()) return fail("unterminated string");
        ++i_;   // the closing quote
        return out;
    }

    std::expected<Value, std::string> array() {
        ++i_;
        Array a;
        if (eat(']')) return Value{std::move(a)};
        for (;;) {
            auto v = value();
            if (!v) return v;
            a.push_back(std::move(*v));
            if (eat(',')) continue;
            if (eat(']')) return Value{std::move(a)};
            return fail("expected , or ]");
        }
    }

    std::expected<Value, std::string> object() {
        ++i_;
        Object o;
        if (eat('}')) return Value{std::move(o)};
        for (;;) {
            skip();
            if (i_ >= s_.size() || s_[i_] != '"') return fail("expected a key");
            auto k = string();
            if (!k) return std::unexpected(k.error());
            if (!eat(':')) return fail("expected :");
            auto v = value();
            if (!v) return v;
            o.emplace_back(std::move(*k), std::move(*v));
            if (eat(',')) continue;
            if (eat('}')) return Value{std::move(o)};
            return fail("expected , or }");
        }
    }
};

inline std::expected<Value, std::string> parse(std::string_view text) { return Parser(text).run(); }

}  // namespace json

// ── The graph ──────────────────────────────────────────────────

enum class NodeKind { Concept, Type, Scalar };
enum class EdgeKind { Refines, Satisfies, Operation };

struct Node {
    std::string id;
    NodeKind kind = NodeKind::Concept;
    std::string cpp;                       // how the library spells it
    std::vector<std::string> aliases;      // other names a query may use
};

struct Edge {
    EdgeKind kind = EdgeKind::Satisfies;
    std::string from, to;                  // refines: `from` is built from `to`; satisfies: type, concept
    std::string op, variant, output;       // operation: what it is, which way, the quantity it yields
    std::vector<std::string> inputs;       // operation: the type, then the concepts that make the way available
    std::vector<std::string> needs;        // operation: quantities the way is built from
    int level = 0;                         // operation: 0 declared .. 3 agrees across scalars (on double)
    std::map<std::string, int> levels;     // operation: the level per scalar the matrix ran
};

// The edges of a route, in an order where every edge comes after the ones that yield its needs.
using Route = std::vector<std::size_t>;

class Graph {
public:
    static std::expected<Graph, std::string> parse(std::string_view text) {
        auto doc = json::parse(text);
        if (!doc) return std::unexpected(doc.error());
        return from_value(*doc);
    }

    static std::expected<Graph, std::string> load(const std::filesystem::path& path) {
        std::ifstream in(path, std::ios::binary);
        if (!in) return std::unexpected("cannot open " + path.string());
        std::ostringstream ss;
        ss << in.rdbuf();
        return parse(ss.str());
    }

    const std::vector<Node>& nodes() const { return nodes_; }
    const std::vector<Edge>& edges() const { return edges_; }

    std::size_t count(NodeKind k) const {
        return static_cast<std::size_t>(std::ranges::count_if(nodes_, [k](const Node& n) { return n.kind == k; }));
    }
    std::size_t count(EdgeKind k) const {
        return static_cast<std::size_t>(std::ranges::count_if(edges_, [k](const Edge& e) { return e.kind == k; }));
    }

    // By id or by alias; null when there is no such node.
    const Node* find(std::string_view name) const {
        const auto it = index_.find(std::string(name));
        return it == index_.end() ? nullptr : &nodes_[it->second];
    }

    // Declared directly: the compiler said so for this type and this concept.
    bool satisfies(std::string_view type, std::string_view concept_name) const {
        const Node* t = find(type);
        const Node* c = find(concept_name);
        return t && c && satisfied_.contains({t->id, c->id});
    }

    // All routes from `from` (a type) to a quantity `goal`, each at most `max_depth` edges long,
    // found breadth first over partial routes: a state is the edges chosen so far and the
    // quantities still open; an open quantity is expanded by every edge that yields it. Distinct as
    // sets of edges, shortest first.
    std::vector<Route> routes(std::string_view from, std::string_view goal, std::size_t max_depth) const {
        const Node* start = find(from);
        if (!start || start->kind == NodeKind::Concept || max_depth == 0) return {};

        struct State {
            Route chosen;
            std::vector<std::string> open;
        };
        const auto yielded = [this](const Route& chosen, const std::string& q) {
            return std::ranges::any_of(chosen, [&](std::size_t i) { return edges_[i].output == q; });
        };

        std::set<Route> seen;
        std::vector<Route> found;
        std::deque<State> queue;
        queue.push_back({{}, {std::string(goal)}});
        while (!queue.empty()) {
            State s = std::move(queue.front());
            queue.pop_front();
            while (!s.open.empty() && yielded(s.chosen, s.open.back())) s.open.pop_back();
            if (s.open.empty()) {
                Route key = s.chosen;
                std::ranges::sort(key);
                if (seen.insert(key).second)
                    if (auto ordered = in_order(s.chosen)) found.push_back(std::move(*ordered));
                continue;
            }
            const std::string next = s.open.back();
            s.open.pop_back();
            if (s.chosen.size() >= max_depth) continue;
            const auto by = by_subject_.find(start->id);
            if (by == by_subject_.end()) continue;
            for (const std::size_t e : by->second) {
                if (edges_[e].output != next) continue;
                State t = s;
                t.chosen.push_back(e);
                for (auto need = edges_[e].needs.rbegin(); need != edges_[e].needs.rend(); ++need) t.open.push_back(*need);   // the first need is expanded first
                queue.push_back(std::move(t));
            }
        }
        std::ranges::stable_sort(found, [](const Route& a, const Route& b) {
            return a.size() != b.size() ? a.size() < b.size() : a < b;
        });
        return found;
    }

    // A route is as verified as its least verified edge.
    int route_level(const Route& r) const {
        if (r.empty()) return 0;
        int lo = 3;
        for (const std::size_t i : r) lo = std::min(lo, edges_[i].level);
        return lo;
    }

private:
    std::vector<Node> nodes_;
    std::vector<Edge> edges_;
    std::unordered_map<std::string, std::size_t> index_;                       // id and alias -> node
    std::set<std::pair<std::string, std::string>> satisfied_;                  // (type id, concept id)
    std::unordered_map<std::string, std::vector<std::size_t>> by_subject_;     // type id -> operation edges

    // Chosen edges sorted so each follows the edges yielding its needs; none if they cannot be.
    std::optional<Route> in_order(const Route& chosen) const {
        Route out;
        std::vector<bool> placed(chosen.size(), false);
        while (out.size() < chosen.size()) {
            bool moved = false;
            for (std::size_t k = 0; k < chosen.size(); ++k) {
                if (placed[k]) continue;
                const bool ready = std::ranges::all_of(edges_[chosen[k]].needs, [&](const std::string& q) {
                    return std::ranges::any_of(out, [&](std::size_t i) { return edges_[i].output == q; });
                });
                if (!ready) continue;
                out.push_back(chosen[k]);
                placed[k] = true;
                moved = true;
            }
            if (!moved) return std::nullopt;
        }
        return out;
    }

    static std::expected<std::string, std::string> text(const json::Value& o, std::string_view key) {
        const json::Value* v = o.get(key);
        if (!v || !v->is_string()) return std::unexpected("missing string \"" + std::string(key) + "\"");
        return v->string();
    }
    static std::expected<std::vector<std::string>, std::string> texts(const json::Value& o, std::string_view key,
                                                                      bool required) {
        const json::Value* v = o.get(key);
        if (!v) {
            if (required) return std::unexpected("missing list \"" + std::string(key) + "\"");
            return std::vector<std::string>{};
        }
        if (!v->is_array()) return std::unexpected("\"" + std::string(key) + "\" is not a list");
        std::vector<std::string> out;
        for (const json::Value& x : v->array()) {
            if (!x.is_string()) return std::unexpected("\"" + std::string(key) + "\" holds a non-string");
            out.push_back(x.string());
        }
        return out;
    }

#define RSC_TRY(var, expr) \
    auto var = (expr);     \
    if (!var) return std::unexpected(var.error())

    static std::expected<Graph, std::string> from_value(const json::Value& doc) {
        if (!doc.is_object()) return std::unexpected("graph: the document is not an object");
        const json::Value* version = doc.get("version");
        if (!version || !version->is_number() || version->number() != 1) return std::unexpected("graph: version must be 1");
        const json::Value* ns = doc.get("nodes");
        const json::Value* es = doc.get("edges");
        if (!ns || !ns->is_array() || !es || !es->is_array()) return std::unexpected("graph: nodes and edges must be lists");

        Graph g;
        for (const json::Value& n : ns->array()) {
            RSC_TRY(id, text(n, "id"));
            RSC_TRY(kind, text(n, "kind"));
            RSC_TRY(aliases, texts(n, "aliases", false));
            Node node;
            node.id = *id;
            if (*kind == "concept") node.kind = NodeKind::Concept;
            else if (*kind == "type") node.kind = NodeKind::Type;
            else if (*kind == "scalar") node.kind = NodeKind::Scalar;
            else return std::unexpected("graph: node " + *id + " has kind " + *kind);
            if (const json::Value* c = n.get("cpp"); c && c->is_string()) node.cpp = c->string();
            node.aliases = std::move(*aliases);
            const std::size_t at = g.nodes_.size();
            for (const std::string& name : [&] { auto v = node.aliases; v.push_back(node.id); return v; }())
                if (!g.index_.emplace(name, at).second) return std::unexpected("graph: the name " + name + " is used twice");
            g.nodes_.push_back(std::move(node));
        }

        const auto node_of = [&g](const std::string& name) { return g.find(name); };
        for (const json::Value& e : es->array()) {
            RSC_TRY(kind, text(e, "kind"));
            Edge edge;
            if (*kind == "refines" || *kind == "satisfies") {
                edge.kind = *kind == "refines" ? EdgeKind::Refines : EdgeKind::Satisfies;
                RSC_TRY(from, text(e, "from"));
                RSC_TRY(to, text(e, "to"));
                const Node* a = node_of(*from);
                const Node* b = node_of(*to);
                if (!a || !b) return std::unexpected("graph: " + *kind + " edge " + *from + " -> " + *to + " names no node");
                if (b->kind != NodeKind::Concept) return std::unexpected("graph: " + *to + " is not a concept");
                if (edge.kind == EdgeKind::Refines && a->kind != NodeKind::Concept)
                    return std::unexpected("graph: " + *from + " does not refine, it is not a concept");
                if (edge.kind == EdgeKind::Satisfies && a->kind == NodeKind::Concept)
                    return std::unexpected("graph: " + *from + " satisfies, it is a concept");
                edge.from = a->id;
                edge.to = b->id;
                if (edge.kind == EdgeKind::Satisfies) g.satisfied_.insert({edge.from, edge.to});
            } else if (*kind == "operation") {
                edge.kind = EdgeKind::Operation;
                RSC_TRY(op, text(e, "op"));
                RSC_TRY(variant, text(e, "variant"));
                RSC_TRY(output, text(e, "output"));
                RSC_TRY(inputs, texts(e, "inputs", true));
                RSC_TRY(needs, texts(e, "needs", true));
                edge.op = *op;
                edge.variant = *variant;
                edge.output = *output;
                edge.inputs = std::move(*inputs);
                edge.needs = std::move(*needs);
                if (edge.inputs.empty()) return std::unexpected("graph: operation " + edge.op + " has no inputs");
                for (std::string& in : edge.inputs) {
                    const Node* n = node_of(in);
                    if (!n) return std::unexpected("graph: operation " + edge.op + " names no node " + in);
                    in = n->id;
                }
                const Node* subject = node_of(edge.inputs[0]);
                if (subject->kind == NodeKind::Concept)
                    return std::unexpected("graph: operation " + edge.op + " starts from a concept");
                for (std::size_t k = 1; k < edge.inputs.size(); ++k)
                    if (node_of(edge.inputs[k])->kind != NodeKind::Concept)
                        return std::unexpected("graph: operation " + edge.op + " lists " + edge.inputs[k] + " after the type");
                if (std::ranges::count(edge.needs, edge.output))
                    return std::unexpected("graph: operation " + edge.op + " needs its own output");
                const json::Value* lv = e.get("level");
                if (!lv || !lv->is_number() || lv->number() < 0 || lv->number() > 3)
                    return std::unexpected("graph: operation " + edge.op + " has no level 0..3");
                edge.level = static_cast<int>(lv->number());
                if (const json::Value* ls = e.get("levels"); ls && ls->is_object())
                    for (const auto& [scalar, x] : ls->object())
                        if (x.is_number()) edge.levels[scalar] = static_cast<int>(x.number());
                g.by_subject_[subject->id].push_back(g.edges_.size());
            } else {
                return std::unexpected("graph: an edge of kind " + *kind);
            }
            g.edges_.push_back(std::move(edge));
        }

        // An operation lists the concepts that make it available; the type must satisfy each one.
        for (const Edge& e : g.edges_) {
            if (e.kind != EdgeKind::Operation) continue;
            for (std::size_t k = 1; k < e.inputs.size(); ++k)
                if (!g.satisfied_.contains({e.inputs[0], e.inputs[k]}))
                    return std::unexpected("graph: " + e.inputs[0] + " " + e.op + " (" + e.variant + ") needs the concept " +
                                           e.inputs[k] + ", which the type does not satisfy");
        }
        return g;
    }
#undef RSC_TRY
};

}  // namespace rsc
