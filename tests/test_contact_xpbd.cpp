// IPC barrier kernel + XPBD distance solver.
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <spatium/physics/mechanics/contact.hpp>
#include <spatium/physics/mechanics/xpbd.hpp>
#include <cmath>

using namespace spatium;
using namespace spatium::physics::mechanics;
using Catch::Approx;

// ── IPC barrier ────────────────────────────────────────────────

TEST_CASE("IPC barrier vanishes at and beyond d_hat", "[ipc][barrier]") {
    constexpr double d_hat = 0.1;
    REQUIRE(ipc_barrier(d_hat,            d_hat) == Approx(0.0).margin(1e-15));
    REQUIRE(ipc_barrier(d_hat * 2.0,      d_hat) == Approx(0.0).margin(1e-15));
    REQUIRE(ipc_barrier(d_hat * 100.0,    d_hat) == Approx(0.0).margin(1e-15));
}

TEST_CASE("IPC barrier grows monotonically as d → 0+", "[ipc][barrier]") {
    constexpr double d_hat = 0.1;
    // B ~ d_hat² · |log(d/d_hat)| as d → 0 — grows logarithmically.
    double b_far    = ipc_barrier(d_hat * 0.5,  d_hat);
    double b_mid    = ipc_barrier(d_hat * 0.1,  d_hat);
    double b_near   = ipc_barrier(d_hat * 1e-6, d_hat);
    REQUIRE(b_far  > 0.0);
    REQUIRE(b_mid  > b_far);
    REQUIRE(b_near > b_mid);
    REQUIRE(std::isinf(ipc_barrier(0.0, d_hat)));
}

TEST_CASE("IPC barrier monotonic decrease with d", "[ipc][barrier]") {
    constexpr double d_hat = 0.1;
    double prev = std::numeric_limits<double>::infinity();
    for (double d = 0.01; d < d_hat; d += 0.005) {
        double b = ipc_barrier(d, d_hat);
        REQUIRE(b >= 0.0);
        REQUIRE(b < prev);
        prev = b;
    }
}

TEST_CASE("IPC barrier gradient direction: pushes away from contact",
          "[ipc][gradient]") {
    constexpr double d_hat = 0.1;
    // dB/dd < 0 in the active region — a step in +d direction reduces energy.
    double g = ipc_barrier_grad(d_hat * 0.5, d_hat);
    REQUIRE(g < 0.0);
}

TEST_CASE("IPC default stiffness is positive and finite", "[ipc][stiffness]") {
    constexpr double d_hat = 0.1;
    double k = ipc_default_stiffness(d_hat);
    REQUIRE(k > 0.0);
    REQUIRE(std::isfinite(k));
}

// ── XPBD ───────────────────────────────────────────────────────

TEST_CASE("XPBD distance constraint pulls particles to rest length",
          "[xpbd][distance]") {
    using Vec3d = Vec<double, 3>;
    std::vector<XpbdParticle<3>> parts(2);
    parts[0].x = parts[0].x_prev = Vec3d{0, 0, 0};
    parts[0].w = 1.0;
    parts[1].x = parts[1].x_prev = Vec3d{2.0, 0, 0};        // stretched 2× rest
    parts[1].w = 1.0;

    std::vector<XpbdDistanceConstraint<3>> cons{{0, 1, 1.0, 0.0, 0.0}};

    // No external accel: gravity-free. Just project the constraint.
    auto no_force = [](const XpbdParticle<3>&, double) {
        return Vec3d{};
    };
    xpbd_step(parts, cons.begin(), cons.end(), no_force, 0.01, 8);

    Vec3d diff = Vec3d{parts[0].x - parts[1].x};
    double dist = std::sqrt(diff.dot(diff));
    REQUIRE(dist == Approx(1.0).margin(1e-6));   // converges to rest after 8 iter
}

TEST_CASE("XPBD pinned particle stays in place", "[xpbd][pin]") {
    using Vec3d = Vec<double, 3>;
    std::vector<XpbdParticle<3>> parts(2);
    parts[0].x = parts[0].x_prev = Vec3d{0, 0, 0};
    parts[0].w = 0.0;       // pinned
    parts[1].x = parts[1].x_prev = Vec3d{1.0, 0, 0};
    parts[1].w = 1.0;

    std::vector<XpbdDistanceConstraint<3>> cons{{0, 1, 1.0, 0.0, 0.0}};
    auto gravity = [](const XpbdParticle<3>&, double) {
        return Vec3d{0, -9.81, 0};
    };

    for (int i = 0; i < 100; ++i)
        xpbd_step(parts, cons.begin(), cons.end(), gravity, 0.01, 4);

    // Pinned particle never moves.
    REQUIRE(parts[0].x[0] == 0.0);
    REQUIRE(parts[0].x[1] == 0.0);
    REQUIRE(parts[0].x[2] == 0.0);

    // Free particle settles below pin at distance ≈ rest length.
    Vec3d diff = Vec3d{parts[0].x - parts[1].x};
    double dist = std::sqrt(diff.dot(diff));
    REQUIRE(std::abs(dist - 1.0) < 0.05);
    REQUIRE(parts[1].x[1] < 0.0);                // hanging down under gravity
}

TEST_CASE("XPBD compliance > 0 makes the spring soft", "[xpbd][compliance]") {
    using Vec3d = Vec<double, 3>;
    auto run = [](double compliance) {
        std::vector<XpbdParticle<3>> parts(2);
        parts[0].x = parts[0].x_prev = Vec3d{0, 0, 0};
        parts[0].w = 0.0;       // pinned
        parts[1].x = parts[1].x_prev = Vec3d{1.0, 0, 0};
        parts[1].w = 1.0;

        std::vector<XpbdDistanceConstraint<3>> cons{
            {0, 1, 1.0, compliance, 0.0}
        };
        auto gravity = [](const XpbdParticle<3>&, double) {
            return Vec3d{0, -9.81, 0};
        };
        for (int i = 0; i < 200; ++i)
            xpbd_step(parts, cons.begin(), cons.end(), gravity, 0.01, 4);
        Vec3d diff = Vec3d{parts[0].x - parts[1].x};
        return std::sqrt(diff.dot(diff));
    };

    double dist_stiff  = run(0.0);     // exactly enforces rest length
    double dist_soft   = run(1e-3);    // some give

    REQUIRE(dist_stiff < dist_soft);   // soft spring stretches further
}

TEST_CASE("XPBD build_distance_constraints from triangle list",
          "[xpbd][cloth-setup]") {
    using Vec3d = Vec<double, 3>;
    std::vector<XpbdParticle<3>> parts(4);
    parts[0].x = parts[0].x_prev = Vec3d{0, 0, 0};
    parts[1].x = parts[1].x_prev = Vec3d{1, 0, 0};
    parts[2].x = parts[2].x_prev = Vec3d{0, 1, 0};
    parts[3].x = parts[3].x_prev = Vec3d{1, 1, 0};
    for (auto& p : parts) p.w = 1.0;

    std::vector<std::array<std::uint32_t, 3>> faces = {
        {0, 1, 2}, {1, 3, 2}        // 2-triangle quad sharing edge 1-2
    };

    auto cons = build_distance_constraints<3>(parts, faces);
    // Quad has 5 unique edges: 0-1, 1-2, 0-2, 1-3, 2-3.
    REQUIRE(cons.size() == 5);
}

// ── XPBD distance-based bending ───────────────────────────────

TEST_CASE("build_bending_distance_constraints picks interior edges only",
          "[xpbd][bending]") {
    using Vec3d = Vec<double, 3>;
    std::vector<XpbdParticle<3>> parts(4);
    parts[0].x = parts[0].x_prev = Vec3d{0, 0, 0};
    parts[1].x = parts[1].x_prev = Vec3d{1, 0, 0};
    parts[2].x = parts[2].x_prev = Vec3d{0, 1, 0};
    parts[3].x = parts[3].x_prev = Vec3d{1, 1, 0};
    for (auto& p : parts) p.w = 1.0;

    // The quad has one interior edge (1-2) shared by both triangles.
    std::vector<std::array<std::uint32_t, 3>> faces = {
        {0, 1, 2}, {1, 3, 2}
    };
    auto bend = build_bending_distance_constraints<3, double>(
        parts, faces);
    REQUIRE(bend.size() == 1);
    // Constraint links the two opposite vertices: 0 and 3.
    REQUIRE(((bend[0].i == 0 && bend[0].j == 3) ||
             (bend[0].i == 3 && bend[0].j == 0)));
    // Rest distance equals √2 (diagonal of the unit square).
    REQUIRE(bend[0].rest == Approx(std::sqrt(2.0)).epsilon(1e-12));
}

TEST_CASE("XPBD distance-bending drives the opposite-vertex distance "
          "back to its rest value when the cloth is folded",
          "[xpbd][bending]") {
    using Vec3d = Vec<double, 3>;
    // Flat butterfly: two triangles sharing edge 0-1, opposite
    // vertices 2 and 3 sit at distance 2 along ±y.
    std::vector<XpbdParticle<3>> parts(4);
    parts[0].x = parts[0].x_prev = Vec3d{0,  0, 0};
    parts[1].x = parts[1].x_prev = Vec3d{1,  0, 0};
    parts[2].x = parts[2].x_prev = Vec3d{0.5,  1, 0};
    parts[3].x = parts[3].x_prev = Vec3d{0.5, -1, 0};
    for (auto& p : parts) p.w = 1.0;

    std::vector<std::array<std::uint32_t, 3>> faces = {
        {0, 1, 2}, {1, 0, 3}
    };
    auto bend = build_bending_distance_constraints<3, double>(
        parts, faces);
    REQUIRE(bend.size() == 1);
    REQUIRE(bend[0].rest == Approx(2.0).epsilon(1e-12));

    // Fold the sheet by rotating the upper triangle 90° about the
    // shared edge (the x-axis): p2 → (0.5, 0, 1). The opposite-
    // vertex distance |p2 − p3| collapses from 2 to √2.
    parts[2].x = Vec3d{0.5, 0.0, 1.0};
    parts[2].x_prev = parts[2].x;
    double d_before = (parts[2].x - parts[3].x).norm();
    REQUIRE(d_before == Approx(std::sqrt(2.0)).epsilon(1e-12));

    for (int it = 0; it < 30; ++it) {
        for (auto& c : bend) c.reset();
        for (auto& c : bend) xpbd_solve_distance(c, parts, 1.0);
    }

    double d_after = (parts[2].x - parts[3].x).norm();
    REQUIRE(d_after == Approx(2.0).epsilon(1e-6));     // restored to rest
    // Distance bending alone does not pull the cloth back to flat
    // — it only enforces |p_k − p_l| = rest. The OTHER (structural)
    // distance constraints in a real cloth handle the local
    // rigidity that, combined with bending, locks in the flat
    // shape. Standalone test guarantees only the distance invariant.
}
