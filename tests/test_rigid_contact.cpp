// Pairwise sphere-sphere / sphere-box contact queries, the broad-phase
// AABB sweep, and the XPBD sphere-collision constraint.
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <spatium/physics/mechanics/rigid_contact.hpp>
#include <cmath>

using namespace spatium;
using namespace spatium::physics::mechanics;
using Catch::Approx;

// ── sphere ↔ sphere ──────────────────────────────────────────────

TEST_CASE("sphere_sphere_contact: separated spheres report positive gap",
          "[rigid_contact][sphere_sphere]") {
    Vec<double, 3> a{0, 0, 0}, b{5, 0, 0};
    auto q = sphere_sphere_contact(a, 1.0, b, 1.0);

    REQUIRE_FALSE(q.inside);
    // Centers 5 apart, radii 1 each -> surfaces 3 apart.
    REQUIRE(q.distance == Approx(3.0).epsilon(1e-12));
    REQUIRE(q.normal[0] == Approx(1.0).epsilon(1e-12));   // A -> B is +x
    REQUIRE(q.normal[1] == Approx(0.0).margin(1e-12));
    REQUIRE(q.normal[2] == Approx(0.0).margin(1e-12));
    // Contact point midway between the two surface points: (1,0,0) and (4,0,0).
    REQUIRE(q.closest_point[0] == Approx(2.5).epsilon(1e-12));
}

TEST_CASE("sphere_sphere_contact: exactly touching spheres report zero gap",
          "[rigid_contact][sphere_sphere]") {
    Vec<double, 3> a{0, 0, 0}, b{3, 0, 0};
    auto q = sphere_sphere_contact(a, 1.0, b, 2.0);   // 1 + 2 == 3

    REQUIRE_FALSE(q.inside);
    REQUIRE(q.distance == Approx(0.0).margin(1e-12));
    // Both surface points coincide at (1,0,0).
    REQUIRE(q.closest_point[0] == Approx(1.0).epsilon(1e-12));
}

TEST_CASE("sphere_sphere_contact: overlapping spheres report inside=true "
          "and the correct penetration depth",
          "[rigid_contact][sphere_sphere]") {
    Vec<double, 3> a{0, 0, 0}, b{1.5, 0, 0};
    auto q = sphere_sphere_contact(a, 1.0, b, 1.0);   // rest 2.0, d 1.5

    REQUIRE(q.inside);
    REQUIRE(q.distance == Approx(0.5).epsilon(1e-12));
    REQUIRE(q.signed_distance() == Approx(-0.5).epsilon(1e-12));
    REQUIRE(q.normal[0] == Approx(1.0).epsilon(1e-12));
}

TEST_CASE("sphere_sphere_contact: normal always points from A toward B",
          "[rigid_contact][sphere_sphere][property]") {
    Vec<double, 3> a{1, 2, 3}, b{4, -1, 5};
    auto q = sphere_sphere_contact(a, 0.5, b, 0.7);
    Vec<double, 3> expected = Vec<double, 3>{(b - a).normalized()};
    REQUIRE(q.normal[0] == Approx(expected[0]).epsilon(1e-12));
    REQUIRE(q.normal[1] == Approx(expected[1]).epsilon(1e-12));
    REQUIRE(q.normal[2] == Approx(expected[2]).epsilon(1e-12));
}

// ── swept (continuous) sphere ↔ sphere ──────────────────────────

TEST_CASE("sweep_sphere_sphere catches a fast body tunnelling straight "
          "through a static one that the discrete query misses at both "
          "the start and end of the step",
          "[rigid_contact][sweep][ccd]") {
    // A moves from x=-2 to x=+2 in one step (displacement 4); B sits
    // still at the origin. Both radii 0.1, so they touch when the
    // gap between centers is 0.2 -- which happens at x = ±0.2, well
    // inside the step, even though A starts and ends 2 units away
    // (>> the 0.2 contact distance) on either side.
    Vec<double, 3> a0{-2.0, 0.0, 0.0}, disp_a{4.0, 0.0, 0.0};
    Vec<double, 3> b0{0.0, 0.0, 0.0}, disp_b{0.0, 0.0, 0.0};
    double ra = 0.1, rb = 0.1;

    // The old, discrete-only approach: check start and end positions
    // only (exactly what a per-substep positional/AABB check without
    // a sweep would see). Both report clean separation -- the crossing
    // in between is invisible to it.
    auto q_start = sphere_sphere_contact(a0, ra, b0, rb);
    Vec<double, 3> a1 = Vec<double, 3>{a0 + disp_a};
    auto q_end = sphere_sphere_contact(a1, ra, b0, rb);
    REQUIRE_FALSE(q_start.inside);
    REQUIRE_FALSE(q_end.inside);
    REQUIRE(q_start.distance > 1.0);
    REQUIRE(q_end.distance > 1.0);

    // The swept query catches the crossing the discrete snapshots miss.
    auto sweep = sweep_sphere_sphere(a0, ra, disp_a, b0, rb, disp_b);
    REQUIRE(sweep.hit);
    // Entering root: |2 - 4*toi| = 0.2 -> toi = 0.45 (first crossing).
    REQUIRE(sweep.toi == Approx(0.45).epsilon(1e-9));
    REQUIRE(sweep.contact.distance == Approx(0.0).margin(1e-9));
    // At toi = 0.45, A's center sits at x = -2 + 4*0.45 = -0.2, B stays
    // at the origin; the two surface points at exact touching coincide
    // at x = -0.2 + radius_a = -0.1.
    REQUIRE(sweep.contact.closest_point[0] == Approx(-0.1).margin(1e-6));
}

TEST_CASE("sweep_sphere_sphere: already overlapping at the start reports "
          "immediate contact (toi = 0), not the future separation point",
          "[rigid_contact][sweep][ccd]") {
    Vec<double, 3> a0{0.0, 0.0, 0.0}, disp_a{0.0, 0.0, 0.0};
    Vec<double, 3> b0{0.5, 0.0, 0.0}, disp_b{5.0, 0.0, 0.0};   // separating fast
    double ra = 1.0, rb = 1.0;   // rest = 2.0, overlapping at start (d = 0.5)

    auto sweep = sweep_sphere_sphere(a0, ra, disp_a, b0, rb, disp_b);
    REQUIRE(sweep.hit);
    REQUIRE(sweep.toi == 0.0);
    REQUIRE(sweep.contact.inside);
}

TEST_CASE("sweep_sphere_sphere: converging but falling short reports no "
          "hit this step",
          "[rigid_contact][sweep][ccd]") {
    // A approaches B but stops 0.5 short of contact distance.
    Vec<double, 3> a0{-3.0, 0.0, 0.0}, disp_a{1.5, 0.0, 0.0};   // ends at -1.5
    Vec<double, 3> b0{0.0, 0.0, 0.0}, disp_b{0.0, 0.0, 0.0};
    double ra = 1.0, rb = 0.0;   // contact distance 1.0; gap at end = 0.5

    auto sweep = sweep_sphere_sphere(a0, ra, disp_a, b0, rb, disp_b);
    REQUIRE_FALSE(sweep.hit);
    REQUIRE(sweep.toi == 1.0);
    REQUIRE(sweep.contact.distance == Approx(0.5).epsilon(1e-9));
}

TEST_CASE("sweep_sphere_sphere: zero relative motion while separated "
          "never reports a hit (and does not crash on the degenerate "
          "zero-direction sweep)",
          "[rigid_contact][sweep][ccd]") {
    Vec<double, 3> a0{0.0, 0.0, 0.0}, disp_a{2.0, 3.0, -1.0};
    Vec<double, 3> b0{10.0, 0.0, 0.0}, disp_b = disp_a;   // identical motion -> zero relative
    double ra = 0.5, rb = 0.5;

    auto sweep = sweep_sphere_sphere(a0, ra, disp_a, b0, rb, disp_b);
    REQUIRE_FALSE(sweep.hit);
    REQUIRE(std::isfinite(sweep.contact.distance));
}

TEST_CASE("sweep_sphere_sphere reduces to sphere_sphere_contact at the "
          "impact positions",
          "[rigid_contact][sweep][ccd][property]") {
    Vec<double, 3> a0{-1.0, 0.2, 0.0}, disp_a{2.0, 0.0, 0.0};
    Vec<double, 3> b0{1.0, 0.2, 0.0}, disp_b{-2.0, 0.0, 0.0};   // head-on approach
    double ra = 0.3, rb = 0.4;

    auto sweep = sweep_sphere_sphere(a0, ra, disp_a, b0, rb, disp_b);
    REQUIRE(sweep.hit);
    Vec<double, 3> pa = Vec<double, 3>{a0 + disp_a * sweep.toi};
    Vec<double, 3> pb = Vec<double, 3>{b0 + disp_b * sweep.toi};
    auto direct = sphere_sphere_contact(pa, ra, pb, rb);
    REQUIRE(sweep.contact.distance == Approx(direct.distance).margin(1e-9));
    REQUIRE(sweep.contact.normal[0] == Approx(direct.normal[0]).epsilon(1e-9));
}

// ── sphere ↔ AABB ─────────────────────────────────────────────────

TEST_CASE("sphere_aabb_contact: sphere well outside the box",
          "[rigid_contact][sphere_aabb]") {
    Vec<double, 3> box_min{-1, -1, -1}, box_max{1, 1, 1};
    Vec<double, 3> c{4, 0, 0};
    auto q = sphere_aabb_contact(c, 1.0, box_min, box_max);

    REQUIRE_FALSE(q.inside);
    // Closest box point is (1,0,0); center-to-point distance 3, minus radius 1.
    REQUIRE(q.distance == Approx(2.0).epsilon(1e-12));
    REQUIRE(q.normal[0] == Approx(1.0).epsilon(1e-12));
    REQUIRE(q.closest_point[0] == Approx(1.0).epsilon(1e-12));
    REQUIRE(q.closest_point[1] == Approx(0.0).margin(1e-12));
    REQUIRE(q.closest_point[2] == Approx(0.0).margin(1e-12));
}

TEST_CASE("sphere_aabb_contact: sphere exactly touching a face",
          "[rigid_contact][sphere_aabb]") {
    Vec<double, 3> box_min{-1, -1, -1}, box_max{1, 1, 1};
    Vec<double, 3> c{2, 0, 0};   // distance to face = 1 == radius
    auto q = sphere_aabb_contact(c, 1.0, box_min, box_max);

    REQUIRE_FALSE(q.inside);
    REQUIRE(q.distance == Approx(0.0).margin(1e-12));
}

TEST_CASE("sphere_aabb_contact: sphere overlapping a face (shallow)",
          "[rigid_contact][sphere_aabb]") {
    Vec<double, 3> box_min{-1, -1, -1}, box_max{1, 1, 1};
    Vec<double, 3> c{1.4, 0, 0};   // distance to face = 0.4 < radius 1
    auto q = sphere_aabb_contact(c, 1.0, box_min, box_max);

    REQUIRE(q.inside);
    REQUIRE(q.distance == Approx(0.6).epsilon(1e-12));
    REQUIRE(q.normal[0] == Approx(1.0).epsilon(1e-12));
    REQUIRE(q.closest_point[0] == Approx(1.0).epsilon(1e-12));
}

TEST_CASE("sphere_aabb_contact: corner case, closest point is a box corner",
          "[rigid_contact][sphere_aabb]") {
    Vec<double, 3> box_min{-1, -1, -1}, box_max{1, 1, 1};
    Vec<double, 3> c{2, 2, 2};
    auto q = sphere_aabb_contact(c, 1.0, box_min, box_max);

    REQUIRE_FALSE(q.inside);
    Vec<double, 3> corner{1, 1, 1};
    double expected = (c - corner).norm() - 1.0;
    REQUIRE(q.distance == Approx(expected).epsilon(1e-12));
    REQUIRE(q.closest_point[0] == Approx(1.0).epsilon(1e-12));
    REQUIRE(q.closest_point[1] == Approx(1.0).epsilon(1e-12));
    REQUIRE(q.closest_point[2] == Approx(1.0).epsilon(1e-12));
}

TEST_CASE("sphere_aabb_contact: deep penetration falls back to the nearest "
          "face's minimum-translation axis",
          "[rigid_contact][sphere_aabb]") {
    // A small sphere fully inside a big box, off-center toward +x: the
    // nearest face is +x (distance 1) rather than -x (distance 9).
    Vec<double, 3> box_min{-10, -10, -10}, box_max{10, 10, 10};
    Vec<double, 3> c{9, 0, 0};
    auto q = sphere_aabb_contact(c, 0.5, box_min, box_max);

    REQUIRE(q.inside);
    REQUIRE(q.normal[0] == Approx(1.0).epsilon(1e-12));
    REQUIRE(q.normal[1] == Approx(0.0).margin(1e-12));
    // Penetration depth = distance-to-face (1) + radius (0.5).
    REQUIRE(q.distance == Approx(1.5).epsilon(1e-12));
    REQUIRE(q.closest_point[0] == Approx(10.0).epsilon(1e-12));
}

// ── broad phase ──────────────────────────────────────────────────

TEST_CASE("broad_phase_aabb_pairs: reports exactly the overlapping AABBs",
          "[rigid_contact][broad_phase]") {
    std::vector<geometry::Box<3, double>> boxes = {
        sphere_aabb(Vec<double, 3>{0, 0, 0}, 1.0),      // 0
        sphere_aabb(Vec<double, 3>{1.5, 0, 0}, 1.0),    // 1: overlaps 0
        sphere_aabb(Vec<double, 3>{100, 0, 0}, 1.0),    // 2: isolated
        sphere_aabb(Vec<double, 3>{1.5, 1.5, 0}, 1.0),  // 3: overlaps 0 and 1
                                                         //    (AABB overlap is
                                                         //    a per-axis interval
                                                         //    test, not a
                                                         //    center-distance
                                                         //    one -- diagonal
                                                         //    offset alone
                                                         //    doesn't separate
                                                         //    two unit boxes).
    };
    auto pairs = broad_phase_aabb_pairs(boxes);

    auto has_pair = [&](std::size_t i, std::size_t j) {
        for (auto& p : pairs)
            if (p.first == i && p.second == j) return true;
        return false;
    };
    REQUIRE(has_pair(0, 1));
    REQUIRE(has_pair(0, 3));
    REQUIRE(has_pair(1, 3));
    REQUIRE_FALSE(has_pair(0, 2));
    REQUIRE_FALSE(has_pair(1, 2));
    REQUIRE_FALSE(has_pair(2, 3));
}

TEST_CASE("swept_sphere_aabb: catches a fast mover's crossing that the "
          "end-of-step-only AABB would miss",
          "[rigid_contact][broad_phase][sweep][ccd]") {
    // Mover A crosses right past stationary B, ending up far away --
    // same scenario as the tunnelling narrow-phase test above, checked
    // one level up at the broad phase.
    Vec<double, 3> a0{-2.0, 0.0, 0.0}, disp_a{4.0, 0.0, 0.0};
    Vec<double, 3> b0{0.0, 0.0, 0.0};
    double ra = 0.1, rb = 0.1;

    auto end_only_box = sphere_aabb(Vec<double, 3>{a0 + disp_a}, ra);
    auto b_box = sphere_aabb(b0, rb);
    REQUIRE_FALSE(end_only_box.intersects(b_box));   // end-position-only box misses it

    auto swept_box = swept_sphere_aabb(a0, ra, disp_a);
    REQUIRE(swept_box.intersects(b_box));            // swept box catches it
}

// ── XPBD sphere collision constraint ────────────────────────────

TEST_CASE("XPBD sphere collision: inactive while not penetrating",
          "[rigid_contact][xpbd]") {
    using Vec3d = Vec<double, 3>;
    std::vector<XpbdParticle<3>> parts(2);
    parts[0].x = parts[0].x_prev = Vec3d{0, 0, 0};
    parts[0].w = 1.0;
    parts[1].x = parts[1].x_prev = Vec3d{3, 0, 0};   // radii sum 2, well separated
    parts[1].w = 1.0;

    XpbdSphereCollisionConstraint<3> c{0, 1, 1.0, 1.0, 0.0, 0.0};
    xpbd_solve_sphere_collision(c, parts, 0.01);

    REQUIRE(parts[0].x[0] == 0.0);
    REQUIRE(parts[1].x[0] == 3.0);
    REQUIRE(c.lambda == 0.0);
}

TEST_CASE("XPBD sphere collision converges to non-overlapping within a few "
          "Gauss-Seidel iterations",
          "[rigid_contact][xpbd]") {
    using Vec3d = Vec<double, 3>;
    std::vector<XpbdParticle<3>> parts(2);
    parts[0].x = parts[0].x_prev = Vec3d{0, 0, 0};
    parts[0].w = 1.0;
    parts[1].x = parts[1].x_prev = Vec3d{0.5, 0, 0};   // deeply overlapping
    parts[1].w = 1.0;

    XpbdSphereCollisionConstraint<3> c{0, 1, 0.6, 0.6, 0.0, 0.0};
    for (int it = 0; it < 8; ++it) {
        c.reset();
        xpbd_solve_sphere_collision(c, parts, 1.0);
    }

    double dist = (parts[0].x - parts[1].x).norm();
    REQUIRE(dist == Approx(1.2).epsilon(1e-6));   // radius_i + radius_j
}

TEST_CASE("XPBD sphere collision: heavier (lower w) body moves less",
          "[rigid_contact][xpbd]") {
    using Vec3d = Vec<double, 3>;
    std::vector<XpbdParticle<3>> parts(2);
    parts[0].x = parts[0].x_prev = Vec3d{0, 0, 0};
    parts[0].w = 0.0;       // pinned / infinite mass
    parts[1].x = parts[1].x_prev = Vec3d{0.5, 0, 0};
    parts[1].w = 1.0;

    XpbdSphereCollisionConstraint<3> c{0, 1, 0.6, 0.6, 0.0, 0.0};
    for (int it = 0; it < 8; ++it) {
        c.reset();
        xpbd_solve_sphere_collision(c, parts, 1.0);
    }

    // Pinned sphere never moves; the free one absorbs the whole correction.
    REQUIRE(parts[0].x[0] == 0.0);
    double dist = (parts[0].x - parts[1].x).norm();
    REQUIRE(dist == Approx(1.2).epsilon(1e-6));
}

TEST_CASE("XPBD sphere collision conserves the center of mass for equal "
          "inverse masses",
          "[rigid_contact][xpbd][conservation]") {
    // The Gauss-Seidel correction is +correction*w_i on one particle and
    // -correction*w_j on the other; with w_i == w_j the two displacements
    // are equal and opposite, so (x_i + x_j)/2 is an exact invariant of
    // the projection — the discrete analogue of the conservation checks
    // test_symplectic_integrators.cpp runs for its own integrators.
    using Vec3d = Vec<double, 3>;
    std::vector<XpbdParticle<3>> parts(2);
    parts[0].x = parts[0].x_prev = Vec3d{-0.3, 0.2, 0.0};
    parts[0].w = 1.0;
    parts[1].x = parts[1].x_prev = Vec3d{0.3, 0.2, 0.1};
    parts[1].w = 1.0;

    Vec3d com_before = Vec3d{(parts[0].x + parts[1].x) * 0.5};

    XpbdSphereCollisionConstraint<3> c{0, 1, 0.5, 0.5, 0.0, 0.0};
    for (int it = 0; it < 10; ++it) {
        c.reset();
        xpbd_solve_sphere_collision(c, parts, 1.0);
    }

    Vec3d com_after = Vec3d{(parts[0].x + parts[1].x) * 0.5};
    REQUIRE(com_after[0] == Approx(com_before[0]).margin(1e-12));
    REQUIRE(com_after[1] == Approx(com_before[1]).margin(1e-12));
    REQUIRE(com_after[2] == Approx(com_before[2]).margin(1e-12));

    // And the pair is genuinely separated now.
    double dist = (parts[0].x - parts[1].x).norm();
    REQUIRE(dist == Approx(1.0).epsilon(1e-6));
}

TEST_CASE("XPBD sphere collision: compliance > 0 leaves residual overlap "
          "compared to the rigid (compliance = 0) limit",
          "[rigid_contact][xpbd][compliance]") {
    using Vec3d = Vec<double, 3>;
    auto run = [](double compliance) {
        std::vector<XpbdParticle<3>> parts(2);
        parts[0].x = parts[0].x_prev = Vec3d{0, 0, 0};
        parts[0].w = 1.0;
        parts[1].x = parts[1].x_prev = Vec3d{0.5, 0, 0};
        parts[1].w = 1.0;

        XpbdSphereCollisionConstraint<3> c{0, 1, 0.6, 0.6, compliance, 0.0};
        double dt = 0.01;
        for (int it = 0; it < 20; ++it) {
            c.reset();
            xpbd_solve_sphere_collision(c, parts, dt);
        }
        return (parts[0].x - parts[1].x).norm();
    };

    double dist_rigid = run(0.0);
    double dist_soft  = run(1e-2);
    REQUIRE(dist_rigid == Approx(1.2).epsilon(1e-6));
    REQUIRE(dist_soft < dist_rigid);   // soft contact still overlapping some
}
