// ccd_policy -- which SurfaceCcdPolicy a query of physics/mechanics/
// surface_ccd.hpp wants, learned from measurement and distilled to a tree
// over SurfaceCcdFeatures.
//
// Every policy is right on every query -- none changes whether contact is
// found, only what finding it costs -- so the oracle is free and nothing
// is labelled by hand: each query is run under all 32 policies (four slab
// directions on or off, two split rules) and the cost of each is recorded.
// Cost is the pairs of cells bounded plus the slab directions tried across
// them, deterministic, so a label does not move with the machine's load;
// time is measured separately at the end, on queries the tree never saw.
//
// The tree is cost-sensitive rather than a classifier over "the best
// policy": a leaf takes the policy of least total cost over its queries,
// and a split is the one that lowers the sum of those totals most. With
// thirty-two policies and many near ties, a classifier would spend its
// splits separating policies that cost the same. Each query's costs are
// taken relative to its cheapest policy, so every query weighs the same:
// in absolute terms a cubic patch costs a thousand vertex-triangle
// queries, and the first measurement of this tree spent every split on
// patches and made vertex-face 2.15x slower than the default.
//
// Queries: tensor Bezier patches of order 1, 2 and 3 with moving control
// points, as Chen et al. test on; vertex-triangle and edge-edge pairs, the
// queries of cloth; a torus from its formula against a patch. The class a
// query came from is never shown to the tree -- only the features a solver
// could compute -- and is used only to report per class.
//
// Build: part of the rsc tools target.

#include <ccd_queries.hpp>
#include <spatium/physics/mechanics/surface_ccd.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <format>
#include <functional>
#include <memory>
#include <numbers>
#include <print>
#include <random>
#include <string>
#include <vector>

namespace mech = spatium::physics::mechanics;
using V3 = spatium::Vec<double, 3>;
using Chart = mech::MovingChart<double>;
using Policy = mech::SurfaceCcdPolicy;

using namespace rsc::ccd;

namespace {

constexpr int kPolicies = 32;
Policy policy_of(int k) {
    Policy p;
    p.axes = static_cast<std::uint32_t>(k & 15);
    p.split_by_reach = (k & 16) == 0;
    return p;
}
std::string name_of(int k) {
    std::string s;
    const auto p = policy_of(k);
    if (p.axes & Policy::kNormals) s += "N";
    if (p.axes & Policy::kCentres) s += "C";
    if (p.axes & Policy::kTangentCross) s += "X";
    if (p.axes & Policy::kWorld) s += "W";
    if (s.empty()) s = "ball";
    return s + (p.split_by_reach ? "/reach" : "/width");
}

constexpr double kWidth = 1e-5;
constexpr std::size_t kBudget = std::size_t{1} << 22;

struct Row {
    std::array<double, 5> x;
    std::array<double, kPolicies> cost;
    std::array<double, kPolicies> rel;   // cost over the query's cheapest
    int kind;
};

std::array<double, 5> features_of(const SurfaceQuery& q) {
    const auto f = mech::surface_ccd_features(q.a, q.b);
    return {double(f.dim_a), double(f.dim_b), f.bend_a, f.bend_b, f.speed};
}

std::vector<Row> measure(const std::vector<SurfaceQuery>& qs) {
    std::vector<Row> rows(qs.size());
    for (std::size_t i = 0; i < qs.size(); ++i) {
        rows[i].x = features_of(qs[i]);
        rows[i].kind = qs[i].kind;
        int hit = -1;
        for (int k = 0; k < kPolicies; ++k) {
            const auto r = mech::first_contact(qs[i].a, qs[i].b, kWidth, kBudget, policy_of(k));
            rows[i].cost[k] = double(r.pairs + r.slabs);
            if (hit < 0) hit = r.hit;
            else if (hit != int(r.hit) && r.pairs < kBudget)
                std::println(stderr, "query {}: policy {} answers {} where the first answered {}", i, name_of(k), r.hit, hit);
        }
        const double least = *std::min_element(rows[i].cost.begin(), rows[i].cost.end());
        for (int k = 0; k < kPolicies; ++k) rows[i].rel[k] = rows[i].cost[k] / std::max(least, 1.0);
    }
    return rows;
}

// ── A cost-sensitive tree ───────────────────────────────────────

struct Node {
    int policy = 0;
    int feature = -1;
    double threshold = 0;
    std::unique_ptr<Node> lo, hi;
};

std::pair<int, double> best_policy(const std::vector<const Row*>& rows) {
    std::array<double, kPolicies> sum{};
    for (const auto* r : rows)
        for (int k = 0; k < kPolicies; ++k) sum[k] += r->rel[k];
    const int k = int(std::min_element(sum.begin(), sum.end()) - sum.begin());
    return {k, sum[k]};
}

std::unique_ptr<Node> fit(const std::vector<const Row*>& rows, int depth, int max_depth, std::size_t min_rows) {
    auto node = std::make_unique<Node>();
    const auto [k, total] = best_policy(rows);
    node->policy = k;
    if (depth >= max_depth || rows.size() < 2 * min_rows) return node;
    double best = total;
    int bf = -1;
    double bt = 0;
    for (int f = 0; f < 5; ++f) {
        std::vector<double> vals;
        for (const auto* r : rows) vals.push_back(r->x[f]);
        std::sort(vals.begin(), vals.end());
        vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
        const std::size_t steps = std::min<std::size_t>(vals.size() - 1, 16);
        for (std::size_t q = 1; q <= steps; ++q) {
            const double t = vals[(vals.size() - 1) * q / (steps + 1)];
            std::vector<const Row*> lo, hi;
            for (const auto* r : rows) (r->x[f] <= t ? lo : hi).push_back(r);
            if (lo.size() < min_rows || hi.size() < min_rows) continue;
            const double c = best_policy(lo).second + best_policy(hi).second;
            if (c < best * (1 - 1e-3)) { best = c; bf = f; bt = t; }
        }
    }
    if (bf < 0) return node;
    std::vector<const Row*> lo, hi;
    for (const auto* r : rows) (r->x[bf] <= bt ? lo : hi).push_back(r);
    node->feature = bf;
    node->threshold = bt;
    node->lo = fit(lo, depth + 1, max_depth, min_rows);
    node->hi = fit(hi, depth + 1, max_depth, min_rows);
    return node;
}

int predict(const Node& n, const std::array<double, 5>& x) {
    if (n.feature < 0) return n.policy;
    return predict(x[n.feature] <= n.threshold ? *n.lo : *n.hi, x);
}

const char* kFeature[] = {"f.dim_a", "f.dim_b", "f.bend_a", "f.bend_b", "f.speed"};
void emit(std::FILE* out, const Node& n, int indent) {
    const std::string pad(indent, ' ');
    if (n.feature < 0) {
        const auto p = policy_of(n.policy);
        std::println(out, "{}return policy({}u, {});   // {}", pad, p.axes, p.split_by_reach ? "true" : "false",
                     name_of(n.policy));
        return;
    }
    std::println(out, "{}if ({} <= {:.17g}) {{", pad, kFeature[n.feature], n.threshold);
    emit(out, *n.lo, indent + 4);
    std::println(out, "{}}}", pad);
    emit(out, *n.hi, indent);
}

}  // namespace

int main(int argc, char** argv) {
    const int per_kind = argc > 1 ? std::atoi(argv[1]) : 60;
    const int max_depth = argc > 2 ? std::atoi(argv[2]) : 3;
    const auto train_q = generate(1, per_kind), test_q = generate(2, per_kind);
    const auto train = measure(train_q), test = measure(test_q);

    std::vector<const Row*> tr;
    for (const auto& r : train) tr.push_back(&r);
    const auto tree = fit(tr, 0, max_depth, std::max<std::size_t>(8, train.size() / 40));
    const int single = best_policy(tr).first;
    const int dflt = [] {
        const Policy d{};
        for (int k = 0; k < kPolicies; ++k)
            if (policy_of(k).axes == d.axes && policy_of(k).split_by_reach == d.split_by_reach) return k;
        return 0;
    }();

    std::println("// ccd_policy: {} queries a kind, depth {}; cost = pairs + slab directions", per_kind, max_depth);
    emit(stdout, *tree, 4);

    // Held-out cost and time per kind: default, best single, tree, oracle.
    std::vector<std::string> report;
    std::println("\n{:<12} | {:>12} {:>12} {:>12} {:>12} | {:>9} {:>9} {:>9}", "kind", "default", "best single", "tree",
                 "oracle", "us dflt", "us tree", "tree/dflt");
    for (int kind = 0; kind < 6; ++kind) {
        double c_d = 0, c_s = 0, c_t = 0, c_o = 0, us_d = 0, us_t = 0;
        int n = 0;
        for (std::size_t i = 0; i < test.size(); ++i) {
            if (test[i].kind != kind) continue;
            ++n;
            const auto& r = test[i];
            const int k = predict(*tree, r.x);
            c_d += r.cost[dflt];
            c_s += r.cost[single];
            c_t += r.cost[k];
            c_o += *std::min_element(r.cost.begin(), r.cost.end());
            using clk = std::chrono::steady_clock;
            for (int rep = 0; rep < 3; ++rep) {
                auto t0 = clk::now();
                (void)mech::first_contact(test_q[i].a, test_q[i].b, kWidth, kBudget, policy_of(dflt));
                auto t1 = clk::now();
                (void)mech::first_contact(test_q[i].a, test_q[i].b, kWidth, kBudget, policy_of(k));
                auto t2 = clk::now();
                us_d += std::chrono::duration<double, std::micro>(t1 - t0).count() / 3;
                us_t += std::chrono::duration<double, std::micro>(t2 - t1).count() / 3;
            }
        }
        std::println("{:<12} | {:>12.0f} {:>12.0f} {:>12.0f} {:>12.0f} | {:>9.1f} {:>9.1f} {:>9.2f}", kKinds[kind], c_d / n,
                     c_s / n, c_t / n, c_o / n, us_d / n, us_t / n, us_t / us_d);
        report.push_back(std::format("{:<12} {:>9.0f} / {:>9.0f} / {:>9.0f}", kKinds[kind], c_d / n, c_t / n, c_o / n));
    }
    std::println("best single policy: {}; default: {}", name_of(single), name_of(dflt));

    // The tree as a header, so what was found is kept, compiled and
    // checked (tests/test_rsc_ccd.cpp) rather than left in a log.
    if (argc > 3) {
        std::FILE* out = std::fopen(argv[3], "w");
        if (!out) { std::println(stderr, "cannot write {}", argv[3]); return 1; }
        std::println(out, "#pragma once");
        std::println(out, "// Generated by rsc/tools/ccd_policy {} {} -- do not edit; rerun it.", per_kind, max_depth);
        std::println(out, "// A SurfaceCcdPolicy per query, from a cost-sensitive tree over SurfaceCcdFeatures,");
        std::println(out, "// fitted on seed 1 and held out on seed 2. Held-out cost per kind (pairs + slab");
        std::println(out, "// directions), default / tree / per-query best:");
        for (const auto& line : report) std::println(out, "//   {}", line);
        std::println(out, "#include <spatium/physics/mechanics/surface_ccd.hpp>\n");
        std::println(out, "namespace rsc::ccd {{\n");
        std::println(out, "struct LearnedSurfaceChooser {{");
        std::println(out, "    static spatium::physics::mechanics::SurfaceCcdPolicy policy(std::uint32_t axes, bool by_reach) {{");
        std::println(out, "        spatium::physics::mechanics::SurfaceCcdPolicy p;");
        std::println(out, "        p.axes = axes;");
        std::println(out, "        p.split_by_reach = by_reach;");
        std::println(out, "        return p;");
        std::println(out, "    }}");
        std::println(out, "    spatium::physics::mechanics::SurfaceCcdPolicy choose(");
        std::println(out, "        const spatium::physics::mechanics::SurfaceCcdFeatures<double>& f) const {{");
        emit(out, *tree, 8);
        std::println(out, "    }}");
        std::println(out, "}};\n");
        std::println(out, "}}  // namespace rsc::ccd");
        std::fclose(out);
    }
}
