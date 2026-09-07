#pragma once

// Explicit finite-difference time-stepping for the 2D wave equation on a
// rectangular membrane with fixed (Dirichlet) boundary on all four edges
// -- a drumhead struck at a point and released. Same numerical family as
// wave_string.hpp (see that header's comment for why an explicit
// central-difference "leapfrog" scheme is used directly instead of
// routing through algebra/ode.hpp), generalized from a 1D stencil to a
// 2D 5-point one:
//
//   u_tt(x,y,t) = c^2 (u_xx + u_yy),      u = 0 on the rectangle boundary
//
// c = sqrt(surface_tension / areal_density) is the 2D wave speed
// (surface_tension in N/m, areal_density in kg/m^2 -- the 2D analogues
// of a string's tension/linear_density).
//
// Stability: the 2D CFL condition is stricter than the 1D case --
// c*dt*sqrt(1/dx^2 + 1/dy^2) <= 1 (courant_number_2d() below reports it).
//
// Closed-form reference used to verify this against: a RECTANGULAR
// membrane (deliberately chosen over the more familiar circular
// drumhead) has mode shapes sin(m*pi*x/lx)*sin(n*pi*y/ly) with frequency
//
//   f_mn = (c/2) * sqrt((m/lx)^2 + (n/ly)^2),   m,n = 1,2,3,...
//
// A circular membrane's mode frequencies are zeros of Bessel functions
// (no closed form as simple as the rectangular case) -- the rectangle
// keeps this header verifiable against an exact formula the same way
// wave_string.hpp is, at the cost of the visually-iconic circular shape.
//
// Initial condition contrast with wave_string.hpp's pluck: a *struck*
// membrane starts at zero displacement with a localized initial
// VELOCITY impulse (a mallet imparts momentum, not a "held then released"
// displacement) -- see make_struck_membrane()'s own comment for the
// phantom-previous-state derivation this implies, which is different
// from (and, for this initial condition, simpler than) the pluck case's.
//
// Damping: `damping` (gamma, units 1/s) adds the same velocity-
// proportional loss term as wave_string.hpp, u_tt = c^2(u_xx+u_yy) -
// 2*gamma*u_t -- see that header's own comment for why this form (rather
// than a flat per-step amplitude multiplier) is the one that doesn't
// bias the closed-form frequency check below.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <cmath>
#  include <cstddef>
#  include <numbers>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::physics::mechanics {

// ── Closed-form reference ────────────────────────────────────────
// Frequency of mode (m,n) of an ideal rectangular membrane clamped on
// all four edges. (m,n) = (1,1) is the fundamental.
template<Scalar T>
T membrane_mode_frequency(T surface_tension, T areal_density, T lx, T ly, int m, int n) {
    using std::sqrt;
    T c  = sqrt(surface_tension / areal_density);
    T mx = T(m) / lx;
    T ny = T(n) / ly;
    return (c / T{2}) * sqrt(mx * mx + ny * ny);
}

// c*dt*sqrt(1/dx^2 + 1/dy^2) -- must stay <= 1 for VibratingMembrane::step's
// explicit scheme to be numerically stable.
template<Scalar T>
T courant_number_2d(T wave_speed, T dt, T dx, T dy) {
    using std::sqrt;
    return wave_speed * dt * sqrt(T{1} / (dx * dx) + T{1} / (dy * dy));
}

// ── Simulation state ─────────────────────────────────────────────
// A rectangular membrane discretized on an nx*ny grid (row-major, index
// (i,j) at offset j*nx+i), i in [0,nx-1] spanning [0,lx], j in [0,ny-1]
// spanning [0,ly]; the four edges (i=0, i=nx-1, j=0, j=ny-1) are the
// fixed Dirichlet boundary. Build one via make_struck_membrane() below.
template<Scalar T = double>
struct VibratingMembrane {
    std::size_t nx, ny;
    T lx, ly;
    T wave_speed;        // c = sqrt(surface_tension / areal_density)
    T damping = T{0};    // gamma in u_tt = c^2(u_xx+u_yy) - 2*gamma*u_t, units 1/s

    std::vector<T> current;   // u^n,   size nx*ny
    std::vector<T> previous;  // u^{n-1}

    [[nodiscard]] T dx() const { return lx / T(nx - 1); }
    [[nodiscard]] T dy() const { return ly / T(ny - 1); }

    T&       at(std::vector<T>& f, std::size_t i, std::size_t j) const { return f[j * nx + i]; }
    const T& at(const std::vector<T>& f, std::size_t i, std::size_t j) const { return f[j * nx + i]; }

    // Advance one time step of size dt. dt should satisfy
    // courant_number_2d(wave_speed, dt, dx(), dy()) <= 1.
    //
    // Same central-difference-in-time damping derivation as
    // wave_string.hpp's step() (u_t(t_n) approximated centrally, solved
    // for the single unknown u^{n+1}), generalized to the 2D Laplacian.
    void step(T dt) {
        const T hx    = dx();
        const T hy    = dy();
        const T rx2   = (wave_speed * dt / hx) * (wave_speed * dt / hx);
        const T ry2   = (wave_speed * dt / hy) * (wave_speed * dt / hy);
        const T gdt   = damping * dt;
        const T denom = T{1} + gdt;

        std::vector<T> next(nx * ny, T{0});
        for (std::size_t j = 1; j + 1 < ny; ++j) {
            for (std::size_t i = 1; i + 1 < nx; ++i) {
                T lap_x = at(current, i + 1, j) - T{2} * at(current, i, j) + at(current, i - 1, j);
                T lap_y = at(current, i, j + 1) - T{2} * at(current, i, j) + at(current, i, j - 1);
                T v = (T{2} * at(current, i, j) - (T{1} - gdt) * at(previous, i, j)
                       + rx2 * lap_x + ry2 * lap_y) / denom;
                at(next, i, j) = v;
            }
        }

        previous = std::move(current);
        current  = std::move(next);
    }
};

// Builds a struck membrane: zero initial displacement, a localized
// initial-velocity impulse (a raised-cosine bump, radius `strike_radius`
// around (strike_x*lx, strike_y*ly), amplitude `strike_height`) standing
// in for a mallet strike. The raised-cosine profile (rather than a hard
// cutoff disc) keeps the velocity field's spatial derivative continuous
// at the bump's edge, which avoids seeding spurious high-frequency
// content the sharp edge of a flat disc would inject into the stencil.
//
// Phantom-previous-state derivation for this (displacement=0, velocity=
// v0) initial condition, at damping=0 (undamped case; see below for
// damping != 0): Taylor-expanding u(dt) = u(0) + dt*u_t(0) +
// (dt^2/2)*u_tt(0), with u(0) = 0 and u_tt(0) = c^2*Laplacian(u(0)) = 0
// (the Laplacian of an identically-zero field is zero), collapses to
// u(dt) = dt*v0 -- no diffusion term survives at the very first step,
// unlike wave_string.hpp's pluck case where u(0) itself already has
// curvature. Plugging that into VibratingMembrane::step's own recurrence
// at n=0 (2*u^0 - u^{-1} + r*Laplacian(u^0) = -u^{-1} since u^0 = 0) and
// matching against the true u^1 = dt*v0 gives the phantom state directly:
// previous = -dt*v0. With damping != 0, step()'s own (1-gdt)/(1+gdt)
// factors apply to this same phantom state, giving a first step of
// dt*v0*(1-gdt)/(1+gdt) -- a small, physically real reduction (damping
// already acts on the very first step, since it opposes velocity and
// velocity is exactly what the strike sets up) rather than an artifact;
// it collapses back to the plain dt*v0 above at damping=0.
template<Scalar T = double>
VibratingMembrane<T> make_struck_membrane(std::size_t nx, std::size_t ny, T lx, T ly,
                                           T wave_speed, T dt,
                                           T strike_x = T{0.5}, T strike_y = T{0.5},
                                           T strike_radius = T{0.1},
                                           T strike_height = T{1},
                                           T damping = T{0}) {
    using std::cos;
    using std::sqrt;

    VibratingMembrane<T> m;
    m.nx = nx; m.ny = ny;
    m.lx = lx; m.ly = ly;
    m.wave_speed = wave_speed;
    m.damping = damping;
    m.current.assign(nx * ny, T{0});
    m.previous.assign(nx * ny, T{0});

    const T hx = m.dx();
    const T hy = m.dy();
    const T cx = strike_x * lx;
    const T cy = strike_y * ly;
    constexpr T half_pi = std::numbers::pi_v<T> / T{2};

    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            bool boundary = (i == 0 || i + 1 == nx || j == 0 || j + 1 == ny);
            T v0 = T{0};
            if (!boundary) {
                T x = T(i) * hx - cx;
                T y = T(j) * hy - cy;
                T r = sqrt(x * x + y * y);
                if (r < strike_radius) {
                    T c = cos(half_pi * r / strike_radius);
                    v0 = strike_height * c * c;
                }
            }
            m.at(m.previous, i, j) = -dt * v0;
        }
    }

    return m;
}

} // namespace spatium::physics::mechanics
