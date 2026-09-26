// napkins -- several napkins dropped onto a sphere and onto each other
// (napkin_scene.hpp), measured: what each part of a substep costs, how
// many pairs the continuous pass asks and stops, and the independent
// count of edges piercing triangles, with the continuous pass on and off.
//
//   napkins [--seconds S] [--napkins K] [--grid N] [--throw SPEED] [--rate SUBSTEPS_A_SECOND] [--trace]

#include "napkin_scene.hpp"

#include <cstdlib>
#include <algorithm>
#include <print>
#include <string>

int main(int argc, char** argv) {
    napkins::Config c;
    double seconds = 2.0;
    bool trace = false;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&] { return std::atof(i + 1 < argc ? argv[++i] : "0"); };
        if (a == "--seconds") seconds = next();
        else if (a == "--napkins") c.napkins = int(next());
        else if (a == "--grid") c.grid = int(next());
        else if (a == "--trace") trace = true;
        else if (a == "--throw") c.throw_speed = next();
        else if (a == "--rate") c.substep_dt = 1.0 / next();
    }
    std::println("{} napkins of {}x{}, {} s at {} substeps a second", c.napkins, c.grid, c.grid, seconds,
                 1.0 / c.substep_dt);
    std::println("{:<8} | {:>8} {:>8} {:>8} {:>8} | {:>9} {:>9} {:>7} {:>8} | {:>10} {:>8}", "ccd", "predict", "broad",
                 "solve", "ccd", "vf/step", "ee/step", "stops", "reverts", "piercings", "exploded");
    std::println("(asked: pairs past the coplanarity filter that reached surface_ccd)");
    for (bool ccd : {false, true}) {
        c.ccd = ccd;
        auto cloth = napkins::make_cloth(c);
        // The claim is that nothing starts to cross; a start that already
        // crosses would be counted every substep and prove nothing.
        if (const long p0 = napkins::count_piercings(cloth, 2.5 * c.spacing); p0) {
            std::println(stderr, "the starting layout already has {} piercings", p0);
            return 1;
        }
        napkins::Stats st;
        const int steps = int(seconds / c.substep_dt);
        const int every = int(0.5 / c.substep_dt);
        for (int k = 0; k < steps && !st.exploded; ++k) {
            napkins::substep(cloth, c, st);
            st.piercings += napkins::count_piercings(cloth, 2.5 * c.spacing);
            if (trace && ccd && (k + 1) % every == 0) {
                // Kinetic energy and the cloth's lowest and highest points:
                // a stack that has settled is still and keeps its height.
                double ke = 0, lo = 1e300, hi = -1e300;
                for (const auto& p : cloth.parts) {
                    const auto v = napkins::V3{(p.x - p.x_prev) * (1.0 / c.substep_dt)};
                    ke += 0.5 * v.dot(v);
                    lo = std::min(lo, p.x[1]);
                    hi = std::max(hi, p.x[1]);
                }
                std::println("  t {:5.2f}  kinetic {:10.4g}  y [{:.3f} {:.3f}]  stops {}  reverts {}  piercings {}",
                             (k + 1) * c.substep_dt, ke, lo, hi, st.hits, st.reverted, st.piercings);
            }
        }
        const double n = double(st.substeps);
        std::println("{:<8} | {:>7.3f}m {:>7.3f}m {:>7.3f}m {:>7.3f}m | {:>9.0f} {:>9.0f} {:>7} {:>8} | {:>10} {:>8}",
                     ccd ? "on" : "off", st.ms_xpbd / n, st.ms_broad / n, st.ms_push / n, st.ms_ccd / n, double(st.vf) / n,
                     double(st.ee) / n, st.hits, st.reverted, st.piercings, st.exploded ? "yes" : "no");
        std::println("         near {:.0f}, asked {:.0f} a substep", double(st.near) / n, double(st.asked) / n);
    }
}
