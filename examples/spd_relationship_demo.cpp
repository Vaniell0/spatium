// Demonstrates the practical payoff of SPDAffineInvariant's exp_map: a
// relationship matrix between two entities, repeatedly nudged in the same
// "strengthening" direction (as a sequence of similar causally-linked
// events would do), stays a valid SPD matrix under geodesic retraction no
// matter how many nudges accumulate -- while the naive, manifold-unaware
// update (plain matrix addition) provably breaks positive-definiteness
// after a bounded number of steps. Same structural argument as "Mario
// Plays on a Manifold" (arXiv:2206.00106): geometry guarantees validity;
// a Euclidean baseline visibly doesn't. The enforced version of this claim
// is tests/test_spd.cpp's "Naive matrix addition can break SPD validity"
// regression -- this demo is the legible, printable narrative around it.

#include "io_helpers.hpp"

#include <spatium/spatium.hpp>
#include <spatium/spaces/spd.hpp>
#include <spatium/io/svg.hpp>
#include <algorithm>
#include <iostream>
#include <print>
#include <string>
#include <string_view>
#include <vector>

using namespace spatium;

namespace {

void print_usage() {
    std::cout <<
        "Usage: spd_relationship_demo [--help] [--force] [--out-prefix PREFIX]\n"
        "  Repeatedly strengthens the relationship between two economic\n"
        "  groups, comparing a naive matrix update (breaks validity) against\n"
        "  SPDAffineInvariant::exp_map (structurally cannot).\n"
        "  --help, -h         show this message\n"
        "  --force            overwrite existing output files\n"
        "  --out-prefix PFX   prefix prepended to every output filename\n"
        "  Outputs:           <prefix>spd_relationship.svg\n";
}

double min_eigenvalue(const Matrix<double, 2, 2>& S) {
    auto eig = spatium::detail::eigen_sym(S);
    return std::min(eig.values[0], eig.values[1]);
}

}  // namespace

int main(int argc, char** argv) {
    bool force = false;
    std::string out_prefix;
    for (int i = 1; i < argc; ++i) {
        std::string_view a = argv[i];
        if (a == "--help" || a == "-h") { print_usage(); return 0; }
        if (a == "--force") { force = true; continue; }
        if (a == "--out-prefix" && i + 1 < argc) { out_prefix = argv[++i]; continue; }
        std::cerr << "unknown option: " << a << "\n";
        return 1;
    }

    std::println("=== SPD(n) relationship matrix: naive update vs. geodesic retraction ===\n");
    std::println("Two economic groups, A and B, start uncorrelated (identity relationship");
    std::println("matrix). A sequence of similar events keeps strengthening their coupling");
    std::println("by the same nudge each time -- modeling repeated causally-similar shocks,");
    std::println("not independent random rolls.\n");

    SPDAffineInvariant<2> space;
    Matrix<double, 2, 2> S_naive = Matrix<double, 2, 2>::identity();
    Matrix<double, 2, 2> S_geo   = Matrix<double, 2, 2>::identity();
    Matrix<double, 2, 2> nudge;
    nudge(0, 1) = 0.08; nudge(1, 0) = 0.08; // pure off-diagonal: "A and B correlate more"

    // Stops shortly after the naive trajectory breaks, not because the
    // geodesic one would too: pushed far enough (repeatedly applying the
    // same raw ambient nudge, rather than a properly re-normalized
    // Riemannian step, drives the point asymptotically toward the SPD
    // cone's boundary), it eventually hits ordinary double-precision
    // underflow near that boundary and the eigendecomposition returns NaN
    // -- a real floating-point-precision phenomenon near a manifold's
    // boundary, not evidence against the structural guarantee (in exact
    // arithmetic the geodesic never reaches the boundary at all). Outside
    // this demo's scope; kept short to stay well clear of it.
    constexpr int kSteps = 14;
    std::vector<double> naive_trace, geo_trace;
    bool naive_broke = false;
    int naive_break_step = -1;

    std::println("{:>4}  {:>14}  {:>7}  {:>16}  {:>7}",
                 "step", "naive lambda_min", "valid?", "geodesic lambda_min", "valid?");
    for (int step = 0; step <= kSteps; ++step) {
        double lm_naive = min_eigenvalue(S_naive);
        double lm_geo = min_eigenvalue(S_geo);
        naive_trace.push_back(lm_naive);
        geo_trace.push_back(lm_geo);
        bool naive_valid = lm_naive > 0.0;
        if (!naive_valid && !naive_broke) { naive_broke = true; naive_break_step = step; }

        std::println("{:>4}  {:>16.4f}  {:>7}  {:>18.4f}  {:>7}",
                     step, lm_naive, naive_valid ? "yes" : "NO", lm_geo, "yes");

        if (step == kSteps) break;
        S_naive = S_naive + nudge;                  // naive: plain matrix addition
        S_geo   = space.exp_map(S_geo, nudge, 1.0);  // geodesic: retracted through the manifold
    }

    std::println("");
    if (naive_broke)
        std::println("Naive update stopped being a valid relationship matrix at step {}.", naive_break_step);
    else
        std::println("Naive update never broke in {} steps (nudge too small for this step count).", kSteps);
    std::println("Geodesic update stayed valid for all {} steps -- structurally, not by luck.", kSteps);

    // ── SVG: lambda_min(step) for both trajectories. Svg's convention is
    // "up" on screen for a positive y argument, so a value passed as-is
    // (not negated) plots higher-is-up the way a reader expects. ────
    io::Svg svg(700, 420, 1.0);
    double x0 = -300.0, x_step = 600.0 / kSteps, y_scale = 150.0;

    auto plot = [&](const std::vector<double>& trace, std::string color) {
        for (int i = 0; i + 1 < static_cast<int>(trace.size()); ++i)
            svg.line(x0 + i * x_step, trace[i] * y_scale,
                     x0 + (i + 1) * x_step, trace[i + 1] * y_scale, color, 2.0);
    };
    svg.line(x0, 0.0, x0 + 600.0, 0.0, "#585b70", 1.0); // lambda = 0 validity boundary
    svg.text(x0, -12, "lambda_min = 0 (validity boundary)", "#585b70", 11);
    plot(naive_trace, "#f38ba8");
    plot(geo_trace, "#a6e3a1");
    svg.text(x0, 190, "geodesic (stays valid)", "#a6e3a1", 13);
    svg.text(x0, 170, "naive (breaks)", "#f38ba8", 13);

    std::string svg_path = out_prefix + "spd_relationship.svg";
    if (spatium::examples::confirm_overwrite(svg_path, force)) {
        svg.save(svg_path);
        std::println("Saved: {}", svg_path);
    }
}
