// N bodies pulling each other on a torus.
//
// The consumer of physics/mechanics/gravity_field.hpp and gravity_tree.hpp:
// a few thousand softened point masses in a periodic box (2D or 3D), two
// kinds of body with different gravitational charge per unit of inertia,
// advanced by kick-drift-kick leapfrog (velocity Verlet) with the force from
// the Barnes-Hut tree and, for comparison, from the direct sum over all pairs
// and all images (Ewald). It prints the drift of total energy and momentum
// and the milliseconds per step of each, and writes an SVG of the final
// positions.
//
//   nbody_demo [--n 2000] [--dim 3] [--steps 40] [--dt 0.004] [--theta 0.4]
//              [--exponent P] [--free] [--seed 1] [--alpha-L 6.0]
//              [--out nbody_density.svg] [--force] [--no-direct] [--no-energy]
//   nbody_demo --bench [--free] [--n N]   time one force evaluation, tree vs direct
//   nbody_demo --theta-scan         tree error against the direct sum vs theta
//
// Units: box side 1, total inertia 1, G = 1.

#include "io_helpers.hpp"
#include <spatium/io/svg.hpp>
#include <spatium/physics/mechanics/body.hpp>
#include <spatium/physics/mechanics/gravity_field.hpp>
#include <spatium/physics/mechanics/gravity_tree.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>

using namespace spatium;
using namespace spatium::physics::mechanics;

namespace {

using Clock = std::chrono::steady_clock;
double ms_since(Clock::time_point t) {
    return std::chrono::duration<double, std::milli>(Clock::now() - t).count();
}

struct Options {
    std::size_t n = 2000;
    bool n_set = false;
    int dim = 3;
    int steps = 40;
    double dt = 0.004;
    double theta = GravityTree<3>::default_theta;
    int exponent = 0;           // 0: Newtonian for the dimension
    bool periodic = true;
    unsigned seed = 1;
    double alpha_L = 6.0;
    std::string out = "nbody_density.svg";
    bool force = false;
    bool direct = true;
    bool energy = true;         // the potential energy is a direct O(N^2) sum
    bool bench = false;
    bool scan = false;
};

template<std::size_t N>
struct World {
    std::vector<PointMass<N>> bodies;     // inertia and phase state
    std::vector<double> coupling;         // charge / inertia, per body
    GravityField<N> field;
    GravityTree<N> tree;

    World(GravityField<N> f) : field(std::move(f)), tree(field) {}

    void load() {
        std::vector<GravitySource<N>> s(bodies.size());
        for (std::size_t i = 0; i < bodies.size(); ++i)
            s[i] = {bodies[i].state.position, coupling[i] * bodies[i].mass};
        field.set_sources(std::move(s));
    }
};

template<std::size_t N>
GravityField<N> make_field(const Options& o) {
    const unsigned p = o.exponent > 0 ? static_cast<unsigned>(o.exponent) : static_cast<unsigned>(N - 1);
    // Softening: a fraction of the mean spacing, so close pairs stay resolved by the step.
    const double eps = 0.3 * std::pow(1.0 / static_cast<double>(o.n), 1.0 / N);
    GravityLaw<N> law{p, 1.0, eps, 0.0};
    if (!o.periodic) return GravityField<N>(law);
    Vec<double, N> box;
    for (std::size_t i = 0; i < N; ++i) box[i] = 1.0;
    return GravityField<N>(law, box, EwaldOptions{o.alpha_L, 4.4});
}

// Two kinds of body: half pulled with ratio 1, half with 0.5. A source's mass is
// coupling * inertia, so the third law holds in inertia: sum(M v) is conserved.
template<std::size_t N>
void fill(World<N>& w, const Options& o, std::mt19937& rng) {
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::normal_distribution<double> g(0.0, 1.0);
    const double v0 = 0.15;
    w.bodies.resize(o.n);
    w.coupling.resize(o.n);
    for (std::size_t i = 0; i < o.n; ++i) {
        Vec<double, N> x, v;
        for (std::size_t d = 0; d < N; ++d) { x[d] = u(rng); v[d] = v0 * g(rng); }
        w.bodies[i] = PointMass<N>(1.0 / o.n, x, v);
        w.coupling[i] = (i % 2 == 0) ? 1.0 : 0.5;
    }
    // Remove the net drift.
    Vec<double, N> p = total_momentum(w.bodies);
    for (auto& b : w.bodies) b.state.velocity = Vec<double, N>{b.state.velocity - p};
    w.load();
    w.tree.rebuild();
}

template<std::size_t N>
double energy(const World<N>& w) {
    return total_kinetic_energy(w.bodies) + w.field.potential_energy();
}

enum class Method { Tree, Direct };

template<std::size_t N>
std::vector<Vec<double, N>> accelerations(World<N>& w, Method m, double theta) {
    std::vector<Vec<double, N>> a(w.bodies.size());
    if (m == Method::Tree) w.tree.rebuild();
    for (std::size_t i = 0; i < a.size(); ++i) {
        const auto& x = w.field.sources()[i].position;
        const auto g = m == Method::Tree ? w.tree.acceleration(x, theta, i) : w.field.acceleration(x, i);
        a[i] = Vec<double, N>{g * w.coupling[i]};
    }
    return a;
}

// Kick-drift-kick; `a` holds the accelerations at the current positions on entry and exit.
template<std::size_t N>
void leapfrog(World<N>& w, std::vector<Vec<double, N>>& a, Method m, double theta, double dt) {
    for (std::size_t i = 0; i < a.size(); ++i) {
        auto& s = w.bodies[i].state;
        s.velocity = Vec<double, N>{s.velocity + a[i] * (0.5 * dt)};
        s.position = w.field.wrap(Vec<double, N>{s.position + s.velocity * dt});
    }
    w.load();
    a = accelerations(w, m, theta);
    for (std::size_t i = 0; i < a.size(); ++i) {
        auto& s = w.bodies[i].state;
        s.velocity = Vec<double, N>{s.velocity + a[i] * (0.5 * dt)};
    }
}

template<std::size_t N>
void write_svg(const World<N>& w, const Options& o) {
    if (!examples::confirm_overwrite(o.out, o.force)) return;
    const double size = 800;
    io::Svg svg(size, size, size);
    svg.line(-0.5, -0.5, 0.5, -0.5, "#585b70", 1).line(0.5, -0.5, 0.5, 0.5, "#585b70", 1)
       .line(0.5, 0.5, -0.5, 0.5, "#585b70", 1).line(-0.5, 0.5, -0.5, -0.5, "#585b70", 1);
    for (std::size_t i = 0; i < w.bodies.size(); ++i) {
        const auto& x = w.bodies[i].state.position;
        svg.point(x[0] - 0.5, x[1] - 0.5, 1.6, w.coupling[i] == 1.0 ? "#89b4fa" : "#f9e2af");
    }
    svg.text(-0.49, 0.47, std::to_string(w.bodies.size()) + " bodies, " + std::to_string(N) +
             "D torus, blue: coupling 1, yellow: coupling 0.5", "#cdd6f4", 12);
    if (svg.save(o.out)) std::printf("wrote %s\n", o.out.c_str());
}

template<std::size_t N>
int simulate(const Options& o) {
    std::mt19937 rng(o.seed);
    World<N> tw(make_field<N>(o));
    fill(tw, o, rng);
    World<N> dw = [&] {
        World<N> d(make_field<N>(o));
        d.bodies = tw.bodies;
        d.coupling = tw.coupling;
        d.load();
        return d;
    }();

    const auto& law = tw.field.law();
    std::printf("%zu bodies, %zuD, %s, law 1/r^%u, softening %.4f, dt %.4g, theta %.2f\n",
                o.n, N, o.periodic ? "torus (Ewald)" : "free space", law.exponent, law.softening, o.dt, o.theta);
    if (o.periodic)
        std::printf("Ewald: alpha %.2f, real cutoff %.3f, %zu reciprocal vectors\n",
                    tw.field.alpha(), tw.field.cutoff(), tw.field.reciprocal_vectors());

    const double e0 = o.energy ? energy(tw) : 0.0;
    const auto p0 = total_momentum(tw.bodies);
    const double pscale = total_kinetic_energy(tw.bodies) > 0
        ? std::sqrt(2.0 * total_kinetic_energy(tw.bodies)) : 1.0;   // |P| against sqrt(2 M K)
    if (o.energy) std::printf("E0 = %.6g (K = %.6g)\n", e0, total_kinetic_energy(tw.bodies));

    auto at = accelerations(tw, Method::Tree, o.theta);
    std::vector<Vec<double, N>> ad;
    if (o.direct) ad = accelerations(dw, Method::Direct, o.theta);

    double t_tree = 0, t_direct = 0;
    const int every = std::max(1, o.steps / 5);
    std::printf("%6s %14s %14s %14s %14s\n", "step", "dE/E0 tree", "|dP| tree", "dE/E0 direct", "|dP| direct");
    for (int s = 1; s <= o.steps; ++s) {
        auto t0 = Clock::now();
        leapfrog(tw, at, Method::Tree, o.theta, o.dt);
        t_tree += ms_since(t0);
        if (o.direct) {
            t0 = Clock::now();
            leapfrog(dw, ad, Method::Direct, o.theta, o.dt);
            t_direct += ms_since(t0);
        }
        if (o.energy && (s % every == 0 || s == o.steps)) {
            const double de = (energy(tw) - e0) / std::abs(e0);
            const double dp = Vec<double, N>{total_momentum(tw.bodies) - p0}.norm() / pscale;
            if (o.direct) {
                const double de2 = (energy(dw) - e0) / std::abs(e0);
                const double dp2 = Vec<double, N>{total_momentum(dw.bodies) - p0}.norm() / pscale;
                std::printf("%6d %14.3e %14.3e %14.3e %14.3e\n", s, de, dp, de2, dp2);
            } else {
                std::printf("%6d %14.3e %14.3e %14s %14s\n", s, de, dp, "-", "-");
            }
        }
    }
    std::printf("ms/step: tree %.2f", t_tree / o.steps);
    if (o.direct) {
        std::printf(", direct %.2f (x%.1f)\n", t_direct / o.steps, t_direct / std::max(t_tree, 1e-9));
        double worst = 0;
        for (std::size_t i = 0; i < o.n; ++i) {
            Vec<double, N> d{tw.bodies[i].state.position - dw.bodies[i].state.position};
            for (std::size_t k = 0; k < N; ++k) d[k] -= std::round(d[k]) * (o.periodic ? 1.0 : 0.0);
            worst = std::max(worst, d.norm());
        }
        std::printf("largest position difference, tree run vs direct run: %.3e\n", worst);
    } else {
        std::printf("\n");
    }
    if constexpr (N >= 2) write_svg(tw, o);
    return 0;
}

// ── Measurements ──────────────────────────────────────────────

template<std::size_t N>
World<N> random_world(std::size_t n, const Options& o) {
    Options q = o;
    q.n = n;
    World<N> w(make_field<N>(q));
    std::mt19937 rng(o.seed);
    fill(w, q, rng);
    return w;
}

double median(std::vector<double> v) {
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

template<std::size_t N>
int theta_scan(const Options& o) {
    World<N> w = random_world<N>(o.n, o);
    std::vector<Vec<double, N>> ref(o.n);
    double rms = 0;
    for (std::size_t i = 0; i < o.n; ++i) {
        ref[i] = w.field.acceleration(w.field.sources()[i].position, i);
        rms += ref[i].dot(ref[i]);
    }
    rms = std::sqrt(rms / o.n);
    std::printf("%zuD, %s, n = %zu: tree force against the direct sum\n", N,
                o.periodic ? "torus" : "free", o.n);
    std::printf("%6s %12s %12s %14s %14s\n", "theta", "mean rel", "max rel", "mean (/rms)", "max (/rms)");
    for (double th : {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0}) {
        GravityTree<N> t(w.field, th);
        double mean = 0, mx = 0, mean_r = 0, mx_r = 0;
        for (std::size_t i = 0; i < o.n; ++i) {
            const auto a = t.acceleration(w.field.sources()[i].position, th, i);
            const double e = Vec<double, N>{a - ref[i]}.norm();
            mean += e / ref[i].norm() / o.n;
            mx = std::max(mx, e / ref[i].norm());
            mean_r += e / rms / o.n;
            mx_r = std::max(mx_r, e / rms);
        }
        std::printf("%6.2f %12.2e %12.2e %14.2e %14.2e\n", th, mean, mx, mean_r, mx_r);
    }
    return 0;
}

// One evaluation of the force on every body: set_sources (the structure factor
// on a torus) + tree build + N queries, against set_sources + N direct sums.
// Past `sample` bodies the direct sum is timed on that many targets and scaled.
template<std::size_t N>
int bench(const Options& o) {
    std::printf("%zuD, %s, theta %.2f, medians of 5; direct timed on at most 256 targets and scaled\n",
                N, o.periodic ? "torus" : "free", o.theta);
    std::printf("%8s %12s %12s %9s\n", "N", "tree ms", "direct ms", "speedup");
    const std::size_t sample = 256;
    std::vector<std::size_t> sizes = {1000, 3000, 10000, 30000, 100000};
    if (o.n_set) sizes = {o.n};
    for (std::size_t n : sizes) {
        World<N> w = random_world<N>(n, o);
        std::vector<double> tt, td;
        for (int rep = 0; rep < 5; ++rep) {
            auto t0 = Clock::now();
            w.load();
            w.tree.rebuild();
            volatile double sink = 0;
            for (std::size_t i = 0; i < n; ++i) sink = sink + w.tree.acceleration(w.field.sources()[i].position, o.theta, i)[0];
            tt.push_back(ms_since(t0));

            t0 = Clock::now();
            w.load();
            const std::size_t m = std::min(n, sample);
            const std::size_t stride = n / m;
            for (std::size_t k = 0; k < m; ++k) {
                const std::size_t i = k * stride;
                sink = sink + w.field.acceleration(w.field.sources()[i].position, i)[0];
            }
            td.push_back(ms_since(t0) * static_cast<double>(n) / static_cast<double>(m));
        }
        const double a = median(tt), b = median(td);
        std::printf("%8zu %12.1f %12.1f %8.1fx\n", n, a, b, b / a);
    }
    return 0;
}

} // namespace

int main(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        auto is = [&](const char* f) { return std::strcmp(argv[i], f) == 0; };
        auto val = [&]() -> const char* {
            if (i + 1 >= argc) { std::fprintf(stderr, "missing value for %s\n", argv[i]); std::exit(2); }
            return argv[++i];
        };
        if (is("--n")) { o.n = std::strtoull(val(), nullptr, 10); o.n_set = true; }
        else if (is("--dim")) o.dim = std::atoi(val());
        else if (is("--steps")) o.steps = std::atoi(val());
        else if (is("--dt")) o.dt = std::atof(val());
        else if (is("--theta")) o.theta = std::atof(val());
        else if (is("--exponent")) o.exponent = std::atoi(val());
        else if (is("--seed")) o.seed = static_cast<unsigned>(std::atoi(val()));
        else if (is("--alpha-L")) o.alpha_L = std::atof(val());
        else if (is("--out")) o.out = val();
        else if (is("--free")) o.periodic = false;
        else if (is("--force")) o.force = true;
        else if (is("--no-direct")) o.direct = false;
        else if (is("--no-energy")) o.energy = false;
        else if (is("--bench")) o.bench = true;
        else if (is("--theta-scan")) o.scan = true;
        else { std::fprintf(stderr, "unknown option %s\n", argv[i]); return 2; }
    }
    if (o.dim != 2 && o.dim != 3) { std::fprintf(stderr, "--dim must be 2 or 3\n"); return 2; }
    try {
        if (o.bench) return o.dim == 2 ? bench<2>(o) : bench<3>(o);
        if (o.scan) return o.dim == 2 ? theta_scan<2>(o) : theta_scan<3>(o);
        return o.dim == 2 ? simulate<2>(o) : simulate<3>(o);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
