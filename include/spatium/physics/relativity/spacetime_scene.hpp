#pragma once
// A black-hole scene described rather than coded: which holes, how many,
// how they move, and what is around them. The spacetime is one metric
// built from the holes -- Minkowski plus one Kerr-Schild term per hole,
// each centred where its hole is at the event's time -- so a hole's path
// is a motion field (in the scene's ordinary vocabulary, on the motion
// clock) that `Field::in_spacetime()` puts on the metric's clock.
//
// What is and is not physics here, since a picture will not say:
//
//   - One hole at rest is Kerr, exactly (the Kerr-Schild form of
//     metric_field.hpp).
//   - Several holes are a superposition, g = eta + sum f_i l_i (x) l_i,
//     which is not a solution of the field equations: the Ricci residual,
//     measured in tests/test_metric_field.cpp's terms, is half the tidal
//     curvature between two holes 40 M apart and grows as they close.
//     `residual_at()` reports it at an event so a caller can see where the
//     approximation holds.
//   - A moving hole's term is the static one at the moving centre, not
//     boosted, so its field is carried along rather than Lorentz
//     contracted -- a further approximation of order v, a few percent at
//     the separations a binary is shown at.
//   - A binary's orbit is Newtonian in shape and shrinks by the
//     quadrupole formula (Peters 1964): a^4 = a0^4 - (256/5) m1 m2 M t.
//     The shrinking stops at `stop` separation, after which the pair
//     keeps circling there; the merger itself is not modelled.
//
// Disk, dust, sky and camera are settings the renderer reads; they carry
// no physics of their own here.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/concepts.hpp>
#  include <spatium/io/field.hpp>
#  include <spatium/physics/relativity/metric_field.hpp>
#  include <array>
#  include <cmath>
#  include <cstdint>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::physics::relativity {

template<Scalar T = double>
struct HoleSpec {
    using F = io::build::Field<T>;
    T mass{1};
    T spin{0};    // about the z axis, |spin| < mass
    // Where the hole's centre is, as motion fields of time (Field::t()).
    F x{T{0}}, y{T{0}}, z{T{0}};
    // Its velocity, as fields of time too: what the metric's boost reads
    // (SpacetimeScene::boost). The field language has no derivative, so a
    // path's velocity is written beside it; binary() writes both.
    F vx{T{0}}, vy{T{0}}, vz{T{0}};
};

template<Scalar T = double>
struct DiskSettings {
    bool on = true;
    T inner = T{6}, outer = T{30};   // in units of the total mass
    T temperature = T{1};            // scales the blackbody colour
};

template<Scalar T = double>
struct DustSettings {
    std::uint32_t count = 0;
    T inner = T{8}, outer = T{60};
    std::uint32_t seed = 7;
};

template<Scalar T = double>
struct SkySettings {
    std::uint32_t stars = 20000;
    bool galaxy_band = true;
    std::uint32_t seed = 11;
};

template<Scalar T = double>
struct CameraSettings {
    T distance = T{60}, azimuth_deg = T{30}, elevation_deg = T{10}, fov_deg = T{50};
};

template<Scalar T = double>
class SpacetimeScene {
public:
    using F = io::build::Field<T>;

    // A hole at rest at the origin; give it a path with the returned spec.
    HoleSpec<T>& hole(T mass, T spin) {
        holes_.push_back(HoleSpec<T>{mass, spin, F(T{0}), F(T{0}), F(T{0})});
        return holes_.back();
    }

    // Two holes circling their centre of mass in the z = 0 plane,
    // separation a0, shrinking by the quadrupole formula down to `stop`
    // (or never, with `inspiral` false).
    void binary(T m1, T m2, T a0, bool inspiral = true, T stop = T{10}, T spin1 = T{0}, T spin2 = T{0}) {
        const T M = m1 + m2;
        const T beta = T{64} / T{5} * m1 * m2 * M;
        const F t = F::t();
        F a, phase;
        if (!inspiral || stop >= a0) {
            a = F(a0);
            phase = F(std::sqrt(M / (a0 * a0 * a0))) * t;
        } else {
            // Until t_stop the separation follows Peters; past it, the pair
            // circles at `stop`. Both pieces are fields of t, joined by min
            // and max, so nothing here is a closure.
            const T t_stop = (a0 * a0 * a0 * a0 - stop * stop * stop * stop) / (T{4} * beta);
            const F ts = min(t, F(t_stop));
            const F a4 = F(a0 * a0 * a0 * a0) - F(T{4} * beta) * ts;
            a = sqrt(sqrt(a4));
            // dphi/dt = sqrt(M) a^(-3/2) and da/dt = -beta / a^3 give
            // phi = (2/5) sqrt(M) / beta (a0^(5/2) - a^(5/2)).
            const auto pow52 = [](const F& v) { return v * v * sqrt(v); };
            const T a052 = a0 * a0 * std::sqrt(a0);
            const F inspiral_phase = F(T{2} / T{5} * std::sqrt(M) / beta) * (F(a052) - pow52(a));
            const T omega_stop = std::sqrt(M / (stop * stop * stop));
            phase = inspiral_phase + F(omega_stop) * max(t - F(t_stop), F(T{0}));
        }
        const F c = cos(phase), s = sin(phase);
        // Velocities: phi' = sqrt(M / a^3) on both pieces (a is `stop` on
        // the second), a' = -beta / a^3 while the pair is shrinking and 0
        // after -- the switch a ramp of width 1e-6 at t_stop, since the
        // field language has no step.
        const F omega = sqrt(F(M) / (a * a * a));
        F adot(T{0});
        if (inspiral && stop < a0) {
            const T t_stop = (a0 * a0 * a0 * a0 - stop * stop * stop * stop) / (T{4} * beta);
            const F on = min(F(T{1}), max(F(T{0}), (F(t_stop) - t) * F(T{1e6})));
            adot = F(T{0}) - F(beta) / (a * a * a) * on;
        }
        const F dx = adot * c - a * omega * s, dy = adot * s + a * omega * c;
        auto& h1 = hole(m1, spin1);
        h1.x = F(m2 / M) * a * c;
        h1.y = F(m2 / M) * a * s;
        h1.vx = F(m2 / M) * dx;
        h1.vy = F(m2 / M) * dy;
        auto& h2 = hole(m2, spin2);
        h2.x = F(T{0}) - F(m1 / M) * a * c;
        h2.y = F(T{0}) - F(m1 / M) * a * s;
        h2.vx = F(T{0}) - F(m1 / M) * dx;
        h2.vy = F(T{0}) - F(m1 / M) * dy;
    }

    // How much of each hole's velocity its Kerr-Schild term is boosted by:
    // 0 superposes terms at rest wherever the holes are (as this scene did
    // first), 1 boosts each into the frame it moves in -- exact for a hole
    // in uniform motion, and what rsc/tools/metric_calibrate finds by
    // minimising the vacuum residual alone.
    void boost(T b) { boost_ = b; }
    T boost() const { return boost_; }

    DiskSettings<T>& disk() { return disk_; }
    DustSettings<T>& dust() { return dust_; }
    SkySettings<T>& sky() { return sky_; }
    CameraSettings<T>& camera() { return camera_; }
    const DiskSettings<T>& disk() const { return disk_; }
    const DustSettings<T>& dust() const { return dust_; }
    const SkySettings<T>& sky() const { return sky_; }
    const CameraSettings<T>& camera() const { return camera_; }
    const std::vector<HoleSpec<T>>& holes() const { return holes_; }
    T total_mass() const {
        T m{0};
        for (const auto& h : holes_) m += h.mass;
        return m;
    }

    // Minkowski plus one Kerr-Schild term per hole, each at its centre at
    // the event's time.
    //
    // Which of Kerr-Schild's two forms: `ingoing` is regular across the
    // future horizon, where a body or a photon falls in; `outgoing` across
    // the past horizon, which is what a ray traced backwards from a camera
    // approaches -- in the ingoing form it meets a coordinate singularity
    // there, its tangent grows without bound and the step jumps the
    // horizon. The outgoing form is the ingoing one with time reversed,
    // and time reversal turns a hole's rotation around, so it is the
    // ingoing form for spin -a with the sign of l's spatial part flipped.
    enum class Form { ingoing, outgoing };
    Result<MetricField<T>> metric(Form form = Form::ingoing) const {
        std::array<F, 10> g;
        for (std::size_t e = 0; e < 10; ++e) {
            const auto [i, j] = MetricField<T>::kEntries[e];
            g[e] = F(i != j ? T{0} : (i == 0 ? T{-1} : T{1}));
        }
        for (const auto& h : holes_) {
            auto cx = h.x.in_spacetime(), cy = h.y.in_spacetime(), cz = h.z.in_spacetime();
            if (!cx) return std::unexpected(cx.error());
            if (!cy) return std::unexpected(cy.error());
            if (!cz) return std::unexpected(cz.error());
            auto ux = h.vx.in_spacetime(), uy = h.vy.in_spacetime(), uz = h.vz.in_spacetime();
            if (!ux) return std::unexpected(ux.error());
            if (!uy) return std::unexpected(uy.error());
            if (!uz) return std::unexpected(uz.error());
            // The event in the hole's rest frame. With b = boost * velocity
            // and gamma = 1 / sqrt(1 - b.b), an offset d taken at one lab
            // time is x' = d + k (b.d) b in the frame the hole is at rest
            // in, k = gamma^2 / (gamma + 1) -- which is (gamma - 1) / b.b,
            // written so it stays finite at b = 0. Kerr-Schild is stationary
            // there, so x' is all it reads.
            // At boost 0 none of this is written, so the pool -- and the
            // shader made from it -- is what it was before the boost existed.
            const bool boosted = boost_ != T{0};
            const F bs(boost_);
            const F bx = bs * *ux, by = bs * *uy, bz = bs * *uz;
            const F gam = boosted ? F(T{1}) / sqrt(F(T{1}) - (bx * bx + by * by + bz * bz)) : F(T{1});
            const F k = gam * gam / (gam + F(T{1}));
            const F dx0 = F::coord(1) - *cx, dy0 = F::coord(2) - *cy, dz0 = F::coord(3) - *cz;
            const F bd = bx * dx0 + by * dy0 + bz * dz0;
            const F x = boosted ? dx0 + k * bd * bx : dx0, y = boosted ? dy0 + k * bd * by : dy0,
                    z = boosted ? dz0 + k * bd * bz : dz0;
            const F M(h.mass), a(h.spin), two(T{2}), half(T{0.5}), quarter(T{0.25});
            const F w = x * x + y * y + z * z - a * a;
            const F r2 = half * w + sqrt(quarter * w * w + a * a * z * z);
            const F r = sqrt(r2);
            const F den = r2 + a * a;
            const bool out = form == Form::outgoing;
            const F sa = out ? F(T{0}) - a : a;          // the spin the reversed form is written with
            const F sign(out ? T{-1} : T{1});
            const std::array<F, 4> lr{F(T{1}), sign * (r * x + sa * y) / den, sign * (r * y - sa * x) / den,
                                      sign * z / r};
            // Back to the lab: a covector goes by the transpose of the boost,
            // l_t = gamma (l'_t - b.l'), l_j = -gamma b_j l'_t + l'_j + k b_j (b.l').
            const F bl = bx * lr[1] + by * lr[2] + bz * lr[3];
            const std::array<F, 4> l = boosted
                ? std::array<F, 4>{gam * (lr[0] - bl), lr[1] - gam * bx * lr[0] + k * bx * bl,
                                   lr[2] - gam * by * lr[0] + k * by * bl, lr[3] - gam * bz * lr[0] + k * bz * bl}
                : lr;
            const F f = two * M * r2 * r / (r2 * r2 + a * a * z * z);
            for (std::size_t e = 0; e < 10; ++e) {
                const auto [i, j] = MetricField<T>::kEntries[e];
                g[e] = g[e] + f * l[static_cast<std::size_t>(i)] * l[static_cast<std::size_t>(j)];
            }
        }
        return MetricField<T>::make(g);
    }

    // How far the scene's metric is from a vacuum solution at an event;
    // see vacuum_residual().
    //
    // Exact (nested Dual, no step size). An event where the superposition
    // is not a spacetime at all -- det g >= 0, between two close holes --
    // is an error rather than a number, since no residual describes it.
    Result<T> residual_at(const Vec<T, 4>& event) const {
        auto g = metric();
        if (!g) return std::unexpected(g.error());
        if (!lorentzian_at(*g, event))
            return std::unexpected(Error{ErrorCode::DegenerateInput,
                                         "metric is not Lorentzian at this event"});
        return vacuum_residual(*g, event);
    }

private:
    std::vector<HoleSpec<T>> holes_;
    T boost_{0};
    DiskSettings<T> disk_;
    DustSettings<T> dust_;
    SkySettings<T> sky_;
    CameraSettings<T> camera_;
};

} // namespace spatium::physics::relativity
