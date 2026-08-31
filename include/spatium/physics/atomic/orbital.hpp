#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <spatium/spaces/implicit.hpp>
#  include <array>
#  include <atomic>
#  include <cmath>
#  include <mutex>
#  include <numbers>
#  include <random>
#  include <thread>
#  include <tuple>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::physics::atomic {

// ── Coordinate conversion ─────────────────────────────────────

template<Scalar T>
std::tuple<T, T, T> cartesian_to_spherical(T x, T y, T z) {
    T r = std::sqrt(x * x + y * y + z * z);
    T theta = (r > T{1e-30}) ? std::acos(std::clamp(z / r, T{-1}, T{1})) : T{0};
    T phi = std::atan2(y, x);
    return {r, theta, phi};
}

// ── Factorial (small values only, n ≤ 170) ────────────────────
// 170! is the last factorial that fits in IEEE 754 double; 21!
// already overflows int64. Legendre / spherical-harmonic recursions
// reach at most n = 2·l_max for high-Z atoms (≈ 14 for the
// lanthanides), so the practical input range stays well under that
// ceiling. Above n = 170 the result is +Inf for double — callers
// should not feed values that high.

template<Scalar T>
T factorial(int n) {
    if (n < 0) return T{0};
    T result{1};
    for (int i = 2; i <= n; ++i)
        result *= static_cast<T>(i);
    return result;
}

// ── Associated Legendre polynomial P_l^m(x) ──────────────────
// Condon-Shortley phase included.

template<Scalar T>
T assoc_legendre(int l, int m, T x) {
    if (m < 0) {
        m = -m;
        T sign = (m % 2 == 0) ? T{1} : T{-1};
        return sign * factorial<T>(l - m) / factorial<T>(l + m) * assoc_legendre<T>(l, m, x);
    }

    T pmm{1};
    if (m > 0) {
        T somx2 = std::sqrt(std::abs(T{1} - x * x));
        T fact = T{1};
        for (int i = 1; i <= m; ++i) {
            pmm *= -fact * somx2;
            fact += T{2};
        }
    }

    if (l == m) return pmm;

    T pmmp1 = x * static_cast<T>(2 * m + 1) * pmm;
    if (l == m + 1) return pmmp1;

    T pll{};
    for (int ll = m + 2; ll <= l; ++ll) {
        pll = (static_cast<T>(2 * ll - 1) * x * pmmp1
             - static_cast<T>(ll + m - 1) * pmm)
            / static_cast<T>(ll - m);
        pmm = pmmp1;
        pmmp1 = pll;
    }
    return pll;
}

// ── Real spherical harmonic Y_l^m(θ, φ) ──────────────────────

template<Scalar T>
T real_spherical_harmonic(int l, int m, T theta, T phi) {
    constexpr T pi = std::numbers::pi_v<T>;
    int am = std::abs(m);

    T K = std::sqrt(
        static_cast<T>(2 * l + 1) / (T{4} * pi)
        * factorial<T>(l - am) / factorial<T>(l + am)
    );

    T P = assoc_legendre<T>(l, am, std::cos(theta));

    if (m > 0)      return K * std::sqrt(T{2}) * P * std::cos(static_cast<T>(m) * phi);
    else if (m < 0) return K * std::sqrt(T{2}) * P * std::sin(static_cast<T>(am) * phi);
    else            return K * P;
}

// ── Generalized Laguerre polynomial L_n^α(x) ─────────────────

template<Scalar T>
T gen_laguerre(int n, T alpha, T x) {
    if (n == 0) return T{1};
    T L0{1};
    T L1 = T{1} + alpha - x;
    if (n == 1) return L1;

    for (int k = 2; k <= n; ++k) {
        T Lk = ((static_cast<T>(2 * k - 1) + alpha - x) * L1
              - (static_cast<T>(k - 1) + alpha) * L0)
             / static_cast<T>(k);
        L0 = L1;
        L1 = Lk;
    }
    return L1;
}

// ── Hydrogenic radial wavefunction R_nl(r) ────────────────────
// Atomic units: a₀ = 1. Exact for one-electron atoms (H, He⁺, Li²⁺, ...).
// For many-electron atoms see `slater_radial_wavefunction` below.

template<Scalar T>
T hydrogenic_radial_wavefunction(int n, int l, T r, T Z = T{1}) {
    T rho = T{2} * Z * r / static_cast<T>(n);
    T na0 = static_cast<T>(n);
    T norm = std::sqrt(
        std::pow(T{2} * Z / na0, 3)
        * factorial<T>(n - l - 1)
        / (T{2} * na0 * std::pow(factorial<T>(n + l), 3))
    );
    T L = gen_laguerre<T>(n - l - 1, T(2 * l + 1), rho);
    return norm * std::exp(-rho / T{2}) * std::pow(rho, l) * L;
}

// Backwards-compatible alias (Z=1, hydrogen).
template<Scalar T>
T radial_wavefunction(int n, int l, T r) {
    return hydrogenic_radial_wavefunction<T>(n, l, r, T{1});
}

// ── Slater's rules: effective nuclear charge Z_eff ────────────
// Shielding σ depends on electron group. For electron in (n, l):
//   - Electrons with larger n  → σ += 0
//   - Electrons in same group  → σ += 0.35 (0.30 for 1s)
//   - Electrons one shell below (n-1) with s/p target → σ += 0.85
//   - Electrons two or more shells below → σ += 1.00
//   - For d/f targets: all electrons with smaller n → σ += 1.00
// Group boundaries: (1s), (2s,2p), (3s,3p), (3d), (4s,4p), (4d), (4f), ...

inline int slater_group_id(int n, int l) {
    // Groups: (1s)=0, (2s,2p)=1, (3s,3p)=2, (3d)=3, (4s,4p)=4, (4d)=5, (4f)=6,
    //         (5s,5p)=7, (5d)=8, (5f)=9, (6s,6p)=10, (6d)=11, (7s,7p)=12, ...
    // s/p group per shell, then each d/f/... gets its own group.
    if (n == 1) return 0;
    // Explicit table keyed by (n, l), following common Slater's rules
    // convention where f group is counted with lower shells.
    //   1:                                 1s → 0
    //   2:                          2s, 2p → 1
    //   3:                  3s, 3p → 2;  3d → 3
    //   4:    4s, 4p → 4;  4d → 5;  4f → 6
    //   5:    5s, 5p → 7;  5d → 8;  5f → 9
    //   6:   6s, 6p → 10;          6d, 6f → 11
    //   7:                          7s, 7p → 12
    if (n == 2) return 1;
    if (n == 3) return (l < 2) ? 2 : 3;
    if (n == 4) {
        if (l < 2) return 4;   // 4s, 4p
        if (l == 2) return 5;  // 4d
        return 6;              // 4f
    }
    if (n == 5) {
        if (l < 2) return 7;   // 5s, 5p
        if (l == 2) return 8;  // 5d
        return 9;              // 5f
    }
    if (n == 6) {
        if (l < 2) return 10;  // 6s, 6p
        return 11;             // 6d, 6f bundled for Slater's rules
    }
    if (n == 7) return 12;
    return 13;  // fallback for exotic configs
}

inline double slater_n_star(int n) {
    // Effective principal quantum number.
    switch (n) {
        case 1: case 2: case 3: return static_cast<double>(n);
        case 4: return 3.7;
        case 5: return 4.0;
        case 6: return 4.2;
        case 7: return 4.3;
        default: return static_cast<double>(n);
    }
}

// Compute Z_eff for an electron in subshell (n_target, l_target) given the
// full electron configuration (list of subshells, each with its occupancy).
template<typename Subshells>
double slater_z_eff(int Z, int n_target, int l_target, const Subshells& cfg) {
    int my_group = slater_group_id(n_target, l_target);
    bool target_is_sp = (l_target < 2);
    double sigma = 0.0;

    for (const auto& sub : cfg) {
        int n  = sub.n;
        int l  = sub.l;
        int ne = sub.electrons;
        int g  = slater_group_id(n, l);

        if (g > my_group) continue;                 // outer shells don't shield
        if (g == my_group) {
            // Same group — subtract one ONLY from the target's own (n,l) subshell.
            int shield_count = ne;
            if (n == n_target && l == l_target) shield_count -= 1;
            if (shield_count <= 0) continue;
            double contrib = (n_target == 1 && n == 1) ? 0.30 : 0.35;
            sigma += shield_count * contrib;
            continue;
        }
        // g < my_group — lower group.
        if (target_is_sp) {
            if (n == n_target - 1) sigma += 0.85 * ne;  // shell just below
            else                   sigma += 1.00 * ne;  // two or more below
        } else {
            // d or f target — every electron in a lower group screens fully.
            sigma += 1.00 * ne;
        }
    }
    double z_eff = static_cast<double>(Z) - sigma;
    return (z_eff > 0.1) ? z_eff : 0.1;  // floor for numerical safety
}

// ── Slater-type orbital (STO) radial ──────────────────────────
// R(r) = N · r^(n*-1) · exp(-ζr),  ζ = Z_eff / n*,  N = (2ζ)^(n*+1/2) / sqrt(Γ(2n*+1)).
// Far more accurate than hydrogenic for multi-electron atoms — it reproduces
// correct atomic radii and first ionization energies to 10-20% without any
// self-consistent iteration. This is the simplest upgrade that matters.

template<Scalar T>
T slater_radial_wavefunction(double n_star, T zeta, T r) {
    // Normalization: (2ζ)^(n*+1/2) / sqrt(Γ(2n*+1))
    // For half-integer n*, use std::tgamma.
    using std::pow; using std::tgamma; using std::sqrt;
    T two_zeta = T{2} * zeta;
    T exp_np = static_cast<T>(n_star) + T{0.5};
    T norm = pow(two_zeta, exp_np) / sqrt(tgamma(T{2} * static_cast<T>(n_star) + T{1}));
    return norm * pow(r, static_cast<T>(n_star) - T{1}) * std::exp(-zeta * r);
}

// Convenience: STO radial for electron (n_target, l_target) in given atom.
template<Scalar T, typename Subshells>
T slater_radial(int Z, int n_target, int l_target, const Subshells& cfg, T r) {
    double zeff    = slater_z_eff(Z, n_target, l_target, cfg);
    double n_star  = slater_n_star(n_target);
    T      zeta    = static_cast<T>(zeff / n_star);
    return slater_radial_wavefunction<T>(n_star, zeta, r);
}

// ── Per-orbital spin occupancy (Hund's rule) ──────────────────
// Subshell (n, l) has 2l+1 m-values, each holds ↑ or ↓ (spin-1/2).
// Hund: fill all m's with ↑ first (parallel spins), then pair with ↓.
// Returns spin count {up, down} for each m in [-l..l].

struct SpinOccupancy {
    int m;            // magnetic quantum number
    int up;           // 0 or 1
    int down;         // 0 or 1
    [[nodiscard]] int total() const { return up + down; }
    [[nodiscard]] bool paired() const { return up == 1 && down == 1; }
    [[nodiscard]] bool unpaired() const { return total() == 1; }
};

// Fill occupancy for a subshell with `electrons` electrons, l ∈ [0, ...].
template<typename Out>
void fill_subshell_spins(int l, int electrons, Out& out) {
    int orbitals = 2 * l + 1;
    out.clear();
    out.reserve(orbitals);
    for (int m = -l; m <= l; ++m) out.push_back({m, 0, 0});

    // Step 1: single ↑ in each m (up to `orbitals` or `electrons`, whichever smaller).
    int remaining = electrons;
    for (int i = 0; i < orbitals && remaining > 0; ++i) {
        out[i].up = 1;
        --remaining;
    }
    // Step 2: pair with ↓.
    for (int i = 0; i < orbitals && remaining > 0; ++i) {
        out[i].down = 1;
        --remaining;
    }
}

// Total unpaired electrons in an atom (paramagnetic signature).
template<typename Subshells>
int unpaired_electrons(const Subshells& cfg) {
    int total = 0;
    for (const auto& sub : cfg) {
        int orbitals = 2 * sub.l + 1;
        int e = sub.electrons;
        // Hund: unpaired count = min(e, orbitals) - max(0, e - orbitals)
        int filled_up  = (e < orbitals) ? e : orbitals;
        int paired     = (e > orbitals) ? (e - orbitals) : 0;
        total += filled_up - paired;
    }
    return total;
}

// ── Probability density |ψ_nlm|² ─────────────────────────────

template<Scalar T>
T orbital_density(int n, int l, int m, T r, T theta, T phi) {
    T R = radial_wavefunction<T>(n, l, r);
    T Y = real_spherical_harmonic<T>(l, m, theta, phi);
    return R * R * Y * Y;
}

// Cartesian convenience
template<Scalar T>
T orbital_density_xyz(int n, int l, int m, T x, T y, T z) {
    auto [r, theta, phi] = cartesian_to_spherical(x, y, z);
    return orbital_density(n, l, m, r, theta, phi);
}

// ── Signed wavefunction (for +/- lobes) ───────────────────────

template<Scalar T>
T orbital_wavefunction_xyz(int n, int l, int m, T x, T y, T z) {
    auto [r, theta, phi] = cartesian_to_spherical(x, y, z);
    return radial_wavefunction<T>(n, l, r) * real_spherical_harmonic<T>(l, m, theta, phi);
}

// ── ImplicitSurface factory for orbital isosurface ────────────

template<Scalar T = double>
ImplicitSurface<T> make_orbital(int n, int l, int m, T iso_level = T{0.0001}) {
    // Orbital spatial extent scales as n²
    T bound = static_cast<T>(n * n) * T{4} + T{6};

    return ImplicitSurface<T>(
        [=](T x, T y, T z) -> T {
            return orbital_density_xyz(n, l, m, x, y, z) - iso_level;
        },
        {-bound, bound, -bound, bound, -bound, bound}
    );
}

// Signed version: isosurface of ψ = ±iso for +/- lobes
template<Scalar T = double>
std::pair<ImplicitSurface<T>, ImplicitSurface<T>> make_orbital_lobes(
    int n, int l, int m, T iso_level = T{0.01})
{
    T bound = static_cast<T>(n * n) * T{4} + T{6};
    typename ImplicitSurface<T>::Bounds b{-bound, bound, -bound, bound, -bound, bound};

    auto positive = ImplicitSurface<T>(
        [=](T x, T y, T z) -> T {
            return orbital_wavefunction_xyz(n, l, m, x, y, z) - iso_level;
        }, b);

    auto negative = ImplicitSurface<T>(
        [=](T x, T y, T z) -> T {
            return -orbital_wavefunction_xyz(n, l, m, x, y, z) - iso_level;
        }, b);

    return {std::move(positive), std::move(negative)};
}

// ── Point cloud from rejection sampling |ψ|² ─────────────────

template<Scalar T = double>
struct OrbitalPointCloud {
    std::vector<Vec<T, 3>> positions;
    std::vector<bool> positive_lobe;  // true if ψ > 0 at this point
};

template<Scalar T = double>
OrbitalPointCloud<T> sample_orbital_points(int n, int l, int m,
    std::size_t num_points, unsigned seed = 42,
    unsigned threads = 0 /* 0 = hardware_concurrency */)
{
    // Tight bound: density negligible beyond this radius
    T bound = static_cast<T>(n * (n + 1)) + T{3};

    // Find max density: 1D radial scan x 2D angular scan (much cheaper than 3D grid)
    constexpr std::size_t r_steps = 150;
    constexpr std::size_t ang_steps = 16;
    constexpr T pi = std::numbers::pi_v<T>;
    T max_density{0};

    for (std::size_t ir = 0; ir <= r_steps; ++ir) {
        T r = bound * static_cast<T>(ir) / static_cast<T>(r_steps);
        for (std::size_t it = 0; it <= ang_steps; ++it) {
            T theta = pi * static_cast<T>(it) / static_cast<T>(ang_steps);
            for (std::size_t ip = 0; ip < ang_steps; ++ip) {
                T phi = T{2} * pi * static_cast<T>(ip) / static_cast<T>(ang_steps);
                T d = orbital_density(n, l, m, r, theta, phi);
                if (d > max_density) max_density = d;
            }
        }
    }
    max_density *= T{1.2};  // safety margin

    if (max_density <= T{0}) return {};

    // Rejection sampling in sphere — parallel.
    // Each worker keeps a local cloud; stop as soon as the total count crosses
    // num_points. A single mutex gates the final merge so no hot shared path.
    if (threads == 0) threads = std::thread::hardware_concurrency();
    if (threads < 1) threads = 1;
    if (num_points < 256) threads = 1;  // overhead not worth it for tiny atoms

    OrbitalPointCloud<T> cloud;
    cloud.positions.reserve(num_points);
    cloud.positive_lobe.reserve(num_points);

    if (threads == 1) {
        std::mt19937 rng(seed);
        std::uniform_real_distribution<T> pos_dist(-bound, bound);
        std::uniform_real_distribution<T> unit_dist(T{0}, T{1});
        while (cloud.positions.size() < num_points) {
            T x = pos_dist(rng), y = pos_dist(rng), z = pos_dist(rng);
            if (x * x + y * y + z * z > bound * bound) continue;
            T density = orbital_density_xyz(n, l, m, x, y, z);
            if (unit_dist(rng) * max_density < density) {
                cloud.positions.push_back(Vec<T, 3>{x, y, z});
                cloud.positive_lobe.push_back(
                    orbital_wavefunction_xyz(n, l, m, x, y, z) >= T{0});
            }
        }
        return cloud;
    }

    std::atomic<std::size_t> accepted{0};
    std::vector<std::vector<Vec<T, 3>>> worker_pos(threads);
    std::vector<std::vector<char>>      worker_sign(threads);  // bool compressed
    for (auto& v : worker_pos)  v.reserve(num_points / threads + 64);
    for (auto& v : worker_sign) v.reserve(num_points / threads + 64);

    std::vector<std::thread> pool;
    pool.reserve(threads);
    for (unsigned tid = 0; tid < threads; ++tid) {
        pool.emplace_back([&, tid, num_points, bound, max_density]() {
            std::mt19937 rng(seed + tid * 104729u);  // distinct stream per worker
            std::uniform_real_distribution<T> pos_dist(-bound, bound);
            std::uniform_real_distribution<T> unit_dist(T{0}, T{1});
            while (accepted.load(std::memory_order_relaxed) < num_points) {
                T x = pos_dist(rng), y = pos_dist(rng), z = pos_dist(rng);
                if (x * x + y * y + z * z > bound * bound) continue;
                T density = orbital_density_xyz(n, l, m, x, y, z);
                if (unit_dist(rng) * max_density < density) {
                    auto prev = accepted.fetch_add(1, std::memory_order_relaxed);
                    if (prev >= num_points) return;
                    worker_pos[tid].push_back(Vec<T, 3>{x, y, z});
                    worker_sign[tid].push_back(
                        orbital_wavefunction_xyz(n, l, m, x, y, z) >= T{0} ? 1 : 0);
                }
            }
        });
    }
    for (auto& t : pool) t.join();

    // Merge — accept may have over-shot num_points by up to (threads-1) samples
    // because workers check the counter before pushing; cap at num_points.
    for (unsigned tid = 0; tid < threads; ++tid) {
        for (std::size_t i = 0; i < worker_pos[tid].size(); ++i) {
            if (cloud.positions.size() >= num_points) break;
            cloud.positions.push_back(worker_pos[tid][i]);
            cloud.positive_lobe.push_back(worker_sign[tid][i] != 0);
        }
    }
    return cloud;
}

} // namespace spatium::physics::atomic
