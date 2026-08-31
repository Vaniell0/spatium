// Narrow-phase point ↔ analytical surface ContactQuery driving the IPC
// barrier.
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <spatium/physics/mechanics/narrow_phase.hpp>
#include <spatium/spaces/parametric.hpp>
#include <cmath>
#include <numbers>

using namespace spatium;
using namespace spatium::physics::mechanics;
using Catch::Approx;

// ── point ↔ sphere ─────────────────────────────────────────────

TEST_CASE("point_to_sphere: outside point, exact distance and normal",
          "[narrow_phase][sphere]") {
    Vec<double, 3> centre{1.0, 0.0, 0.0};
    auto q = point_to_sphere(Vec<double, 3>{4.0, 0.0, 0.0}, centre, 1.0);

    REQUIRE_FALSE(q.inside);
    REQUIRE(q.distance == Approx(2.0).epsilon(1e-12));
    // Closest point lies on +x face of the unit sphere centred at (1,0,0).
    REQUIRE(q.closest_point[0] == Approx(2.0).epsilon(1e-12));
    REQUIRE(q.closest_point[1] == Approx(0.0).margin(1e-12));
    REQUIRE(q.closest_point[2] == Approx(0.0).margin(1e-12));
    // Outward normal points along +x.
    REQUIRE(q.normal[0] == Approx(1.0).epsilon(1e-12));
}

TEST_CASE("point_to_sphere: inside flag for points within radius",
          "[narrow_phase][sphere]") {
    auto q = point_to_sphere(Vec<double, 3>{0.3, 0.0, 0.0},
                             Vec<double, 3>{}, 1.0);
    REQUIRE(q.inside);
    REQUIRE(q.signed_distance() < 0.0);
    REQUIRE(q.distance == Approx(0.7).epsilon(1e-12));
}

// ── point ↔ torus ──────────────────────────────────────────────

TEST_CASE("point_to_torus: point on a meridional circle just outside the tube",
          "[narrow_phase][torus]") {
    geometry::Torus<double> t;
    t.major_radius = 1.0;
    t.minor_radius = 0.25;

    // Point sits on +x axis at distance 1 + 0.25 + 0.05 from origin.
    // Closest point should be the outer rim of the tube, distance 0.05.
    Vec<double, 3> p{1.0 + 0.25 + 0.05, 0.0, 0.0};
    auto q = point_to_torus(p, t);

    REQUIRE_FALSE(q.inside);
    REQUIRE(q.distance == Approx(0.05).epsilon(1e-10));
    // Normal points outward in +x.
    REQUIRE(q.normal[0] == Approx(1.0).epsilon(1e-10));
}

TEST_CASE("point_to_torus: point inside the tube reports inside=true",
          "[narrow_phase][torus]") {
    geometry::Torus<double> t;
    t.major_radius = 1.0;
    t.minor_radius = 0.25;

    // Place a point just inside the tube on +x.
    Vec<double, 3> p{1.0 + 0.20, 0.0, 0.0};
    auto q = point_to_torus(p, t);

    REQUIRE(q.inside);
    REQUIRE(q.distance == Approx(0.05).epsilon(1e-10));
    REQUIRE(q.signed_distance() < 0.0);
}

TEST_CASE("point_to_torus: distance through the donut hole equals R",
          "[narrow_phase][torus]") {
    geometry::Torus<double> t;
    t.major_radius = 1.0;
    t.minor_radius = 0.25;

    // The torus axis is +z; the centre of the hole is at the origin.
    // Closest tube point to (0,0,0) sits on the centreline circle
    // (any meridional rim point), distance R - r = 0.75.
    auto q = point_to_torus(Vec<double, 3>{}, t);

    REQUIRE_FALSE(q.inside);
    REQUIRE(q.distance == Approx(t.major_radius - t.minor_radius)
                          .epsilon(1e-10));
}

// ── point ↔ parametric ─────────────────────────────────────────

TEST_CASE("point_to: unit sphere via parametric chart matches "
          "analytical sphere",
          "[narrow_phase][parametric]") {
    using std::sin; using std::cos;
    constexpr double pi = std::numbers::pi_v<double>;

    ParametricSurface<double> sphere(
        [](double u, double v) -> Vec<double, 3> {
            return {sin(v) * cos(u), sin(v) * sin(u), cos(v)};
        },
        ParametricSurface<double>::Domain{0.0, 2.0 * pi, 0.0, pi},
        /*periodic_u=*/true, /*periodic_v=*/false);

    Vec<double, 3> p{2.0, 0.0, 0.0};
    auto qp = point_to(p, sphere);
    auto qs = point_to_sphere(p, Vec<double, 3>{}, 1.0);

    REQUIRE(qp.distance == Approx(qs.distance).epsilon(1e-3));
    // Closest point should be near (1,0,0).
    REQUIRE(qp.closest_point[0] == Approx(1.0).epsilon(1e-3));
    REQUIRE(qp.closest_point[1] == Approx(0.0).margin(1e-3));
    REQUIRE(qp.closest_point[2] == Approx(0.0).margin(1e-3));
}

// ── IPC contact energy / force on top of the queries ──────────

TEST_CASE("ipc_contact_energy is zero outside d_hat, positive inside",
          "[narrow_phase][ipc]") {
    geometry::Torus<double> t;
    t.major_radius = 1.0;
    t.minor_radius = 0.25;
    constexpr double d_hat = 0.05;

    // Far from the tube: distance > d_hat, energy = 0.
    auto far = point_to_torus(Vec<double, 3>{2.0, 0.0, 0.0}, t);
    REQUIRE(ipc_contact_energy(far, d_hat) == Approx(0.0).margin(1e-15));

    // Just outside the tube but inside the band: energy > 0.
    auto near = point_to_torus(
        Vec<double, 3>{1.0 + 0.25 + d_hat * 0.5, 0.0, 0.0}, t);
    REQUIRE(near.distance < d_hat);
    REQUIRE(ipc_contact_energy(near, d_hat) > 0.0);
}

TEST_CASE("ipc_contact_force points along the outward normal and grows "
          "as the point approaches the surface",
          "[narrow_phase][ipc]") {
    Vec<double, 3> centre{};
    constexpr double r = 1.0;
    constexpr double d_hat = 0.1;

    // Three points along +x at distances 0.5·d_hat, 0.1·d_hat, 0.01·d_hat
    // outside the unit sphere.
    auto qa = point_to_sphere(
        Vec<double, 3>{r + d_hat * 0.5, 0.0, 0.0}, centre, r);
    auto qb = point_to_sphere(
        Vec<double, 3>{r + d_hat * 0.1, 0.0, 0.0}, centre, r);
    auto qc = point_to_sphere(
        Vec<double, 3>{r + d_hat * 0.01, 0.0, 0.0}, centre, r);

    auto fa = ipc_contact_force(qa, d_hat);
    auto fb = ipc_contact_force(qb, d_hat);
    auto fc = ipc_contact_force(qc, d_hat);

    // All forces point outward in +x.
    REQUIRE(fa[0] > 0.0); REQUIRE(fb[0] > 0.0); REQUIRE(fc[0] > 0.0);
    REQUIRE(fa[1] == Approx(0.0).margin(1e-15));
    REQUIRE(fb[1] == Approx(0.0).margin(1e-15));
    REQUIRE(fc[1] == Approx(0.0).margin(1e-15));

    // Force magnitude grows as we approach the surface (barrier gets steep).
    REQUIRE(fb[0] > fa[0]);
    REQUIRE(fc[0] > fb[0]);
}

TEST_CASE("ipc_contact_force on inside point pushes outward (repair step)",
          "[narrow_phase][ipc]") {
    auto q = point_to_sphere(Vec<double, 3>{0.5, 0.0, 0.0},
                             Vec<double, 3>{}, 1.0);
    REQUIRE(q.inside);
    auto f = ipc_contact_force(q, 0.1);
    // Outward normal is +x, force must push along it.
    REQUIRE(f[0] > 0.0);
}

// ── property tests: closest_point, normal, penetration ─────────

TEST_CASE("point_to_sphere: closest_point lies exactly on the sphere",
          "[narrow_phase][sphere][property]") {
    Vec<double, 3> c{0.3, -0.7, 1.2};
    constexpr double r = 1.5;
    // Cover both inside and outside query points.
    for (auto p : {Vec<double, 3>{2.0, 1.0, 3.0},
                   Vec<double, 3>{0.4, -0.5, 1.0},
                   Vec<double, 3>{-3.0, -2.0, 0.0},
                   Vec<double, 3>{0.3, -0.7, 1.2 + 0.1}}) {
        auto q = point_to_sphere(p, c, r);
        double d_to_centre = (q.closest_point - c).norm();
        REQUIRE(d_to_centre == Approx(r).epsilon(1e-12));
    }
}

TEST_CASE("point_to_sphere: outward normal is the radial direction",
          "[narrow_phase][sphere][property]") {
    Vec<double, 3> c{};
    constexpr double r = 1.0;
    Vec<double, 3> p{2.0, 1.0, 0.5};
    auto q = point_to_sphere(p, c, r);
    Vec<double, 3> radial = (q.closest_point - c) / r;
    REQUIRE(q.normal[0] == Approx(radial[0]).epsilon(1e-12));
    REQUIRE(q.normal[1] == Approx(radial[1]).epsilon(1e-12));
    REQUIRE(q.normal[2] == Approx(radial[2]).epsilon(1e-12));
}

TEST_CASE("point_to_torus: closest_point sits exactly on the tube surface",
          "[narrow_phase][torus][property]") {
    geometry::Torus<double> t;
    t.major_radius = 1.0;
    t.minor_radius = 0.3;
    // Build a few tilt-axis tori too.
    geometry::Torus<double> tt;
    tt.center = Vec<double, 3>{0.5, 0.0, 0.0};
    tt.axis   = Vec<double, 3>{1.0, 0.0, 0.0};
    tt.major_radius = 0.7; tt.minor_radius = 0.1;

    for (auto* torus : {&t, &tt}) {
        for (auto p : {Vec<double, 3>{2.0, 1.0, 0.5},
                       Vec<double, 3>{0.0, 0.0, 0.0},
                       Vec<double, 3>{1.1, 0.0, 0.05}}) {
            auto q = point_to_torus(p, *torus);
            // Reproduce the torus implicit form on the closest_point:
            //   (|loc|² + R² − r²)² − 4R²(x² + y²) = 0  (in local frame)
            Vec<double, 3> u, v, w = torus->axis;
            geometry::torus_basis(w, u, v);
            Vec<double, 3> d = q.closest_point - torus->center;
            double x = d.dot(u), y = d.dot(v), z = d.dot(w);
            double R = torus->major_radius, rr = torus->minor_radius;
            double lhs = (x*x + y*y + z*z + R*R - rr*rr);
            double impl = lhs * lhs - 4.0 * R*R * (x*x + y*y);
            REQUIRE(impl == Approx(0.0).margin(1e-8));
        }
    }
}

TEST_CASE("point_to_sphere: penetration depth equals r − |p − c| for inside",
          "[narrow_phase][sphere][property]") {
    Vec<double, 3> c{};
    constexpr double r = 2.0;
    Vec<double, 3> p{0.5, 0.7, -0.3};
    auto q = point_to_sphere(p, c, r);
    REQUIRE(q.inside);
    double expected = r - (p - c).norm();
    REQUIRE(q.distance == Approx(expected).epsilon(1e-12));
}

TEST_CASE("ContactQuery::signed_distance flips with inside flag",
          "[narrow_phase][api]") {
    Vec<double, 3> c{};
    auto qo = point_to_sphere(Vec<double, 3>{2.0, 0.0, 0.0}, c, 1.0);
    auto qi = point_to_sphere(Vec<double, 3>{0.5, 0.0, 0.0}, c, 1.0);
    REQUIRE(qo.signed_distance() == Approx( 1.0).epsilon(1e-12));
    REQUIRE(qi.signed_distance() == Approx(-0.5).epsilon(1e-12));
}

// ── ipc_contact_force matches finite difference of energy ──────
// The IPC barrier is a real potential, so the force on the query
// point must equal −∂E/∂x. Verifies sign + magnitude in one shot —
// catches the most likely class of bugs in this header.

TEST_CASE("ipc_contact_force = −∂(ipc_contact_energy)/∂x  (sphere)",
          "[narrow_phase][ipc][property]") {
    Vec<double, 3> c{};
    constexpr double r = 1.0;
    constexpr double d_hat = 0.1;
    constexpr double h = 1e-5;

    auto E = [&](const Vec<double, 3>& x) {
        return ipc_contact_energy(point_to_sphere(x, c, r), d_hat);
    };

    Vec<double, 3> x0{1.0 + d_hat * 0.3, 0.05, -0.02};
    auto q = point_to_sphere(x0, c, r);
    Vec<double, 3> f_ana = ipc_contact_force(q, d_hat);

    Vec<double, 3> f_fd{};
    for (int i = 0; i < 3; ++i) {
        Vec<double, 3> xp = x0, xm = x0;
        xp[i] += h; xm[i] -= h;
        f_fd[i] = -(E(xp) - E(xm)) / (2.0 * h);
    }

    // Finite-diff tolerance — barrier is steep but well-conditioned at
    // d ≈ 0.3·d_hat; 1e-5 is comfortable.
    for (int i = 0; i < 3; ++i)
        REQUIRE(f_ana[i] == Approx(f_fd[i]).margin(5e-5));
}

TEST_CASE("ipc_contact_force = −∂(ipc_contact_energy)/∂x  (torus)",
          "[narrow_phase][ipc][property]") {
    geometry::Torus<double> t;
    t.major_radius = 1.0;
    t.minor_radius = 0.25;
    constexpr double d_hat = 0.05;
    constexpr double h = 1e-5;

    auto E = [&](const Vec<double, 3>& x) {
        return ipc_contact_energy(point_to_torus(x, t), d_hat);
    };

    Vec<double, 3> x0{1.0 + 0.25 + d_hat * 0.4, 0.02, 0.01};
    auto q = point_to_torus(x0, t);
    Vec<double, 3> f_ana = ipc_contact_force(q, d_hat);

    Vec<double, 3> f_fd{};
    for (int i = 0; i < 3; ++i) {
        Vec<double, 3> xp = x0, xm = x0;
        xp[i] += h; xm[i] -= h;
        f_fd[i] = -(E(xp) - E(xm)) / (2.0 * h);
    }
    for (int i = 0; i < 3; ++i)
        REQUIRE(f_ana[i] == Approx(f_fd[i]).margin(5e-4));
}

// ── ContactSurface concept dispatch ─────────────────────────────
// `point_to(p, surface)` resolves at compile time to the right
// analytical / parametric / fallback path. No runtime dispatch.

TEST_CASE("ContactSurface concept satisfied by Sphere<2>, Torus, "
          "ParametricSurface", "[narrow_phase][concept]") {
    using std::sin; using std::cos;
    constexpr double pi = std::numbers::pi_v<double>;

    static_assert(ContactSurface<spatium::Sphere<2, double>>);
    static_assert(ContactSurface<geometry::Torus<double>>);
    static_assert(ContactSurface<ParametricSurface<double>>);
    // Type without a `point_to` overload must NOT satisfy.
    static_assert(!ContactSurface<int>);

    // Functional check — generic API resolves to specialised paths.
    spatium::Sphere<2, double> sph{1.0};
    geometry::Torus<double> tor;
    tor.major_radius = 1.0; tor.minor_radius = 0.25;
    ParametricSurface<double> psurf(
        [](double u, double v) -> Vec<double, 3> {
            return {sin(v) * cos(u), sin(v) * sin(u), cos(v)};
        },
        ParametricSurface<double>::Domain{0.0, 2.0 * pi, 0.0, pi},
        true, false);

    Vec<double, 3> p{2.0, 0.0, 0.0};
    auto qs = point_to(p, sph);
    auto qt = point_to(p, tor);
    auto qp = point_to(p, psurf);

    REQUIRE(qs.distance == Approx(1.0).epsilon(1e-12));
    REQUIRE(qt.distance == Approx(2.0 - (1.0 + 0.25)).epsilon(1e-12));
    REQUIRE(qp.distance == Approx(1.0).epsilon(1e-3));
}

TEST_CASE("ipc_contact_force_on dispatches via the concept",
          "[narrow_phase][concept][ipc]") {
    spatium::Sphere<2, double> sph{1.0};
    Vec<double, 3> p{1.05, 0.0, 0.0};
    auto F = ipc_contact_force_on(p, sph, 0.1);
    REQUIRE(F[0] > 0.0);                       // outward push
    REQUIRE(F[1] == Approx(0.0).margin(1e-15));
    REQUIRE(F[2] == Approx(0.0).margin(1e-15));
}

// ── parametric torus reproduces analytical torus distance ───────

TEST_CASE("point_to on torus chart matches analytical "
          "point_to_torus",
          "[narrow_phase][parametric]") {
    using std::sin; using std::cos;
    constexpr double pi = std::numbers::pi_v<double>;
    constexpr double R = 1.0, r = 0.25;

    ParametricSurface<double> torus_surf(
        [=](double u, double v) -> Vec<double, 3> {
            double cu = cos(u), su = sin(u);
            double cv = cos(v), sv = sin(v);
            return {(R + r * cv) * cu, (R + r * cv) * su, r * sv};
        },
        ParametricSurface<double>::Domain{0.0, 2.0 * pi, 0.0, 2.0 * pi},
        /*periodic_u=*/true, /*periodic_v=*/true);

    geometry::Torus<double> torus_geom;
    torus_geom.major_radius = R;
    torus_geom.minor_radius = r;

    // Same query points should yield the same unsigned distance, up to
    // the parametric Newton tolerance from ParametricSurface::find_params
    // (8×8 grid seed + 5 Newton iters — gets within ~5 % of analytical
    // for points away from the singular meridional axis).
    for (auto p : {Vec<double, 3>{1.4, 0.0, 0.05},
                   Vec<double, 3>{0.0, 1.5, 0.10},
                   Vec<double, 3>{2.0, 0.5, 0.30}}) {
        auto qa = point_to_torus(p, torus_geom);
        auto qp = point_to(p, torus_surf);
        REQUIRE(qp.distance == Approx(qa.distance).epsilon(5e-2));
    }
}
