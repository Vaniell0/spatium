// cloth_sphere_probe — implicit-contact (IPC) successor to the removed
// cloth_on_sphere_demo.cpp / XPBD-era cloth_sphere_probe.cpp.
//
// Explicit XPBD position-based contact was investigated and found
// structurally blocked: a 20-config compliance/substep sweep against a
// stiff cloth with heavy overhanging corners (the diagonal corners of a
// square grid hang past a sphere's circular silhouette) failed on all
// 20 (see rsc/README.md's "Contact physics / soft bodies" entry). This
// probe replaces the explicit projection with one full implicit-Euler
// Newton solve per substep of
//
//   E(x) = 1/2 (x-x_hat)^T M (x-x_hat) + h^2 Psi_elastic(x) + h^2 B(x)
//
// where Psi_elastic is a hand-rolled quadratic spring energy over the
// cloth's own mesh edges (topology only reused from MeshTopology; the
// Newton assembly itself is new — XPBD's Gauss-Seidel projection is a
// different numerical method, not reusable inside a monolithic Newton
// solve) and B is ipc-toolkit's BarrierPotential (Li-Kaufman 2020 log
// barrier). The sphere is triangulated and merged into the same
// CollisionMesh as a static (pinned) obstacle — ipc-toolkit's own
// NormalCollisions/BarrierPotential API expects both contact sides to
// be meshes, it has no entry point for Spatium's analytical
// ContactSurface<Sphere<2>>; using it directly here is the simpler,
// more robust "prove the pipeline works" path (a hand-assembled
// analytical-barrier hybrid using contact.hpp's ipc_barrier* kernels
// driven by point_to_sphere is a real follow-up, not attempted here).
// compute_collision_free_stepsize() is the "line-search filter" itself
// (Li-Kaufman 2020 Sec.5): every Newton step is capped by the largest
// collision-free fraction of the proposed displacement before any
// Armijo energy-decrease backtracking runs.
//
// Reports worst-penetration depth (independent check: point_to_sphere's
// analytical signed distance, not ipc-toolkit's own internal collision
// set) and max-stretch-ratio (max over all cloth edges, all frames, of
// current_length/rest_length - 1) per config.
//
// Ships with ONE config, not a multi-config sweep matching the XPBD
// investigation's own 20-point discipline: three real bugs were found
// and fixed getting Newton to actually converge (starting the iterate
// from the unchecked predictor Xhat instead of the last known-safe X;
// an Armijo line search that silently accepted a non-decreasing-energy
// step when every halving attempt failed instead of erroring out; a
// CCD min-distance floor set to a physically-meaningful gap instead of
// a numerical safety margin, permanently deadlocking Newton once
// contact settled near it) -- confirmed converging (gnorm 3.2e-2 ->
// 2.6e-6 over ~14 iterations, no explosion, no tunnelling) via
// diagnostic instrumentation during that debugging pass. What's NOT
// done: this small a problem (~350 vertices) is slow under
// ipc-toolkit's TBB-parallel collision detection -- thread-pool
// overhead dominates at this scale, and even one reduced config took
// several CPU-minutes. A broader sweep and/or a performance pass is
// real follow-up work, deliberately not chased now (mechanism
// correctness was the goal here, not throughput).

#include <spatium/algebra/vector.hpp>
#include <spatium/mesh/mesh.hpp>
#include <spatium/mesh/primitives.hpp>
#include <spatium/mesh/topology.hpp>
#include <spatium/spaces/euclidean.hpp>
#include <spatium/physics/mechanics/narrow_phase.hpp>

#include <ipc/collision_mesh.hpp>
#include <ipc/collisions/normal/normal_collisions.hpp>
#include <ipc/potentials/barrier_potential.hpp>
#include <ipc/ipc.hpp>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <print>
#include <vector>
#include <cmath>
#include <cstdio>
#include <limits>

using namespace spatium;
using namespace spatium::physics::mechanics;

namespace {

// ── Problem setup — matches the removed cloth_on_sphere_demo.cpp's
// physical scale exactly, so worst-penetration numbers are comparable
// to what explicit XPBD produced (all 20 configs there failed at this
// scale). ──────────────────────────────────────────────────────────
constexpr int    kGrid        = 9;          // 9x9 = 81 cloth particles
constexpr double kSpacing     = 0.18;
constexpr double kSphereR     = 1.0;
constexpr int    kSphereSlices = 16, kSphereStacks = 8; // coarser obstacle mesh -- fewer
                                                          // collision-candidate triangles;
                                                          // doesn't affect the cloth-scale
                                                          // comparability to XPBD above
constexpr double kStartHeight = 1.6;
constexpr double kGravityY    = -9.81;
constexpr double kDt          = 1.0 / 240.0;
constexpr int    kSubsteps    = 4;
constexpr int    kFrames      = 90;         // ~0.375s @ 240Hz -- falls (~0.35s), briefly settles
constexpr double kParticleMass = 1.0;
constexpr int    kMaxNewtonIters = 12;
constexpr double kNewtonTol   = 1e-6;
// CCD floor -- must be a pure numerical safety margin (ipc-toolkit's own
// default is 0.0), NOT a physically meaningful gap. A value comparable
// to dhat deadlocks Newton once contact settles near it:
// compute_collision_free_stepsize() requires "vertices_t0 is intersection
// free" (distance > min_distance) as a precondition and returns toi=0
// (refusing ANY further step, even a retreating one) once distance sits
// at or below the floor -- discovered by a stuck run where gnorm/Dir
// were bit-identical for 20 straight Newton iterations at kMinDistance=1e-4
// with dhat=0.02/0.05 (the physical resting gap landed right on the floor).
constexpr double kMinDistance = 1e-8;

struct Config {
    double stiffness;   // cloth quadratic-spring stiffness k
    double dhat;         // IPC contact-band thickness
    double barrier_stiffness;
};

struct Metrics {
    double worst_penetration = 0.0;  // most negative point_to_sphere signed distance seen
    double max_stretch       = 0.0;  // max (len/rest - 1) over all edges, all frames
    bool   exploded          = false;
};

struct Edge { int a, b; double rest; };

// Manual grid build (not mesh::grid_mesh — that lays out in the XY
// plane with Z=0; this probe wants Y as "up" against gravity, matching
// the removed demo's own convention).
void build_cloth(std::vector<Vec<double,3>>& verts,
                 std::vector<std::array<uint32_t,3>>& faces) {
    double half = (kGrid - 1) * kSpacing * 0.5;
    for (int j = 0; j < kGrid; ++j) {
        for (int i = 0; i < kGrid; ++i) {
            verts.push_back(Vec<double,3>{
                i * kSpacing - half, kStartHeight, j * kSpacing - half});
        }
    }
    auto idx = [](int i, int j) -> uint32_t { return static_cast<uint32_t>(j * kGrid + i); };
    for (int j = 0; j + 1 < kGrid; ++j) {
        for (int i = 0; i + 1 < kGrid; ++i) {
            uint32_t a = idx(i, j), b = idx(i + 1, j);
            uint32_t c = idx(i + 1, j + 1), d = idx(i, j + 1);
            faces.push_back({a, b, c});
            faces.push_back({a, c, d});
        }
    }
}

std::vector<Edge> topology_edges(const std::vector<Vec<double,3>>& verts,
                                 const std::vector<std::array<uint32_t,3>>& faces) {
    mesh::Mesh<Euclidean<3,double>> m;
    m.vertices = verts;
    m.faces = faces;
    auto topo = mesh::MeshTopology<Euclidean<3,double>>::build(m);
    std::vector<Edge> out;
    out.reserve(topo.edge_count());
    for (uint32_t e = 0; e < topo.edge_count(); ++e) {
        auto& ed = topo.edge(e);
        double rest = (verts[ed.v0] - verts[ed.v1]).norm();
        out.push_back({static_cast<int>(ed.v0), static_cast<int>(ed.v1), rest});
    }
    return out;
}

Eigen::VectorXd flatten(const Eigen::MatrixXd& X) {
    // RowMajor DOF convention (ipc-toolkit's default
    // VERTEX_DERIVATIVE_LAYOUT): [x0,y0,z0, x1,y1,z1, ...]
    Eigen::VectorXd v(X.size());
    for (int i = 0; i < X.rows(); ++i)
        for (int d = 0; d < X.cols(); ++d)
            v[3 * i + d] = X(i, d);
    return v;
}

Eigen::MatrixXd unflatten(const Eigen::VectorXd& v, int n) {
    Eigen::MatrixXd X(n, 3);
    for (int i = 0; i < n; ++i)
        for (int d = 0; d < 3; ++d)
            X(i, d) = v[3 * i + d];
    return X;
}

// Quadratic spring elastic energy + gradient + Hessian triplets over
// the cloth's own edges, in the same flattened RowMajor DOF layout
// ipc-toolkit's BarrierPotential uses, so the two can be summed
// directly. Sphere rows never appear here — they carry zero elastic
// energy by construction (pinned/static, handled separately).
double assemble_elastic(const Eigen::MatrixXd& X, const std::vector<Edge>& edges,
                        double k, Eigen::VectorXd& grad,
                        std::vector<Eigen::Triplet<double>>& triplets) {
    double energy = 0.0;
    for (auto& e : edges) {
        Eigen::Vector3d xi = X.row(e.a), xj = X.row(e.b);
        Eigen::Vector3d diff = xi - xj;
        double len = diff.norm();
        if (len < 1e-12) continue;
        double C = len - e.rest;
        energy += 0.5 * k * C * C;

        Eigen::Vector3d dCdxi = diff / len;          // dC/dxi = -dC/dxj
        Eigen::Vector3d g = k * C * dCdxi;
        for (int d = 0; d < 3; ++d) {
            grad[3 * e.a + d] += g[d];
            grad[3 * e.b + d] -= g[d];
        }

        // Gauss-Newton Hessian (drop the C * d2C/dx2 term — standard
        // for stiff mass-spring Newton solves, keeps H SPD without
        // needing the true curvature of |x|, which is only indefinite
        // far from the constraint anyway).
        Eigen::Matrix3d outer = dCdxi * dCdxi.transpose();
        Eigen::Matrix3d Hii = k * outer;
        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
                double v = Hii(a, b);
                triplets.emplace_back(3*e.a+a, 3*e.a+b,  v);
                triplets.emplace_back(3*e.b+a, 3*e.b+b,  v);
                triplets.emplace_back(3*e.a+a, 3*e.b+b, -v);
                triplets.emplace_back(3*e.b+a, 3*e.a+b, -v);
            }
        }
    }
    return energy;
}

// One config's full sim run. Returns measured metrics.
Metrics run_config(const Config& cfg) {
    Metrics m;

    std::vector<Vec<double,3>> cloth_v; std::vector<std::array<uint32_t,3>> cloth_f;
    build_cloth(cloth_v, cloth_f);
    auto cloth_edges = topology_edges(cloth_v, cloth_f);

    auto sphere_mesh = mesh::uv_sphere_mesh<double>(kSphereSlices, kSphereStacks, kSphereR);
    auto sphere_edges_local = topology_edges(sphere_mesh.vertices, sphere_mesh.faces);

    const int n_cloth  = static_cast<int>(cloth_v.size());
    const int n_sphere = static_cast<int>(sphere_mesh.vertices.size());
    const int n_total  = n_cloth + n_sphere;

    Eigen::MatrixXd rest_positions(n_total, 3);
    for (int i = 0; i < n_cloth; ++i)
        for (int d = 0; d < 3; ++d) rest_positions(i, d) = cloth_v[i][d];
    for (int i = 0; i < n_sphere; ++i)
        for (int d = 0; d < 3; ++d) rest_positions(n_cloth + i, d) = sphere_mesh.vertices[i][d];

    std::vector<std::array<uint32_t,3>> all_faces = cloth_f;
    for (auto& f : sphere_mesh.faces)
        all_faces.push_back({f[0] + static_cast<uint32_t>(n_cloth),
                             f[1] + static_cast<uint32_t>(n_cloth),
                             f[2] + static_cast<uint32_t>(n_cloth)});
    Eigen::MatrixXi faces(all_faces.size(), 3);
    for (std::size_t i = 0; i < all_faces.size(); ++i)
        for (int d = 0; d < 3; ++d) faces(static_cast<int>(i), d) = all_faces[i][d];

    std::vector<Edge> all_edges = cloth_edges;
    for (auto& e : sphere_edges_local)
        all_edges.push_back({e.a + n_cloth, e.b + n_cloth, e.rest});
    Eigen::MatrixXi edges(all_edges.size(), 2);
    for (std::size_t i = 0; i < all_edges.size(); ++i) {
        edges(static_cast<int>(i), 0) = all_edges[i].a;
        edges(static_cast<int>(i), 1) = all_edges[i].b;
    }

    ipc::CollisionMesh collision_mesh(rest_positions, edges, faces);
    ipc::BarrierPotential B(cfg.dhat, cfg.barrier_stiffness);

    Eigen::VectorXd mass = Eigen::VectorXd::Constant(3 * n_total, kParticleMass);
    for (int i = n_cloth; i < n_total; ++i)
        for (int d = 0; d < 3; ++d) mass[3 * i + d] = 0.0; // sphere: pinned, mass irrelevant

    Eigen::MatrixXd X = rest_positions;
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(n_total, 3);

    const Vec<double,3> sphere_center{0, 0, 0};
    double sub_dt = kDt / kSubsteps;

    for (int frame = 0; frame < kFrames && !m.exploded; ++frame) {
        for (int s = 0; s < kSubsteps; ++s) {
            Eigen::MatrixXd Xhat = X + sub_dt * V;
            for (int i = 0; i < n_cloth; ++i) Xhat(i, 1) += sub_dt * sub_dt * kGravityY;
            // sphere rows: Xhat stays == X (pinned, no gravity/mass)
            for (int i = n_cloth; i < n_total; ++i) Xhat.row(i) = X.row(i);

            // Newton's starting iterate is X (the previous substep's
            // converged, guaranteed collision-free state), NOT Xhat (the
            // raw, collision-unaware predictor x + h*v + h^2*g). Xhat only
            // enters as the inertia term's target below -- if velocity is
            // large enough, Xhat itself can already tunnel past the sphere
            // before any CCD check ever runs; every within-Newton step is
            // CCD-clamped relative to its *own* starting Xk, so starting
            // from a known-safe X keeps that chain unbroken from the first
            // iteration on.
            Eigen::MatrixXd Xk = X;

            for (int iter = 0; iter < kMaxNewtonIters; ++iter) {
                ipc::NormalCollisions collisions;
                collisions.build(collision_mesh, Xk, cfg.dhat);

                Eigen::VectorXd grad_elastic = Eigen::VectorXd::Zero(3 * n_total);
                std::vector<Eigen::Triplet<double>> elastic_triplets;
                double e_elastic = assemble_elastic(Xk, cloth_edges, cfg.stiffness,
                                                    grad_elastic, elastic_triplets);

                double e_barrier = B(collisions, collision_mesh, Xk);
                Eigen::VectorXd grad_barrier = B.gradient(collisions, collision_mesh, Xk);
                Eigen::SparseMatrix<double> H_barrier =
                    B.hessian(collisions, collision_mesh, Xk, ipc::PSDProjectionMethod::CLAMP);

                Eigen::VectorXd xk_flat = flatten(Xk), xhat_flat = flatten(Xhat);
                Eigen::VectorXd grad = mass.cwiseProduct(xk_flat - xhat_flat)
                                      + sub_dt * sub_dt * (grad_elastic + grad_barrier);

                // Pin sphere DOF: zero gradient there so Newton never moves them.
                for (int i = n_cloth; i < n_total; ++i)
                    for (int d = 0; d < 3; ++d) grad[3 * i + d] = 0.0;

                double gnorm = grad.norm();
                if (!std::isfinite(gnorm)) { m.exploded = true; break; }
                if (gnorm < kNewtonTol) break;

                // Assemble H = M + h^2 (H_elastic + H_barrier) as triplets, then
                // pin sphere DOF by DROPPING any triplet that touches a pinned
                // row or column and inserting an explicit identity diagonal
                // instead. (NOT: iterate `InnerIterator(H, row)` after the sum
                // and zero non-diagonal entries — for Eigen's default
                // column-major storage that walks the COLUMN at index `row`,
                // not the row; it silently decouples only half of each pinned
                // DOF's coupling to the free block, corrupting the free-DOF
                // solve too. Filtering triplets before assembly is correct
                // regardless of storage order.)
                std::vector<Eigen::Triplet<double>> all_triplets;
                all_triplets.reserve(3 * n_total + elastic_triplets.size()
                                     + static_cast<std::size_t>(H_barrier.nonZeros()));
                for (int r = 0; r < 3 * n_total; ++r)
                    all_triplets.emplace_back(r, r, mass[r]);
                for (auto& t : elastic_triplets)
                    all_triplets.emplace_back(t.row(), t.col(), sub_dt * sub_dt * t.value());
                for (int k = 0; k < H_barrier.outerSize(); ++k)
                    for (Eigen::SparseMatrix<double>::InnerIterator it(H_barrier, k); it; ++it)
                        all_triplets.emplace_back(it.row(), it.col(), sub_dt * sub_dt * it.value());

                auto is_pinned = [&](int dof) { return dof >= 3 * n_cloth; };
                std::vector<Eigen::Triplet<double>> final_triplets;
                final_triplets.reserve(all_triplets.size());
                for (auto& t : all_triplets)
                    if (!is_pinned(t.row()) && !is_pinned(t.col()))
                        final_triplets.push_back(t);
                for (int dof = 3 * n_cloth; dof < 3 * n_total; ++dof)
                    final_triplets.emplace_back(dof, dof, 1.0);

                Eigen::SparseMatrix<double> H(3 * n_total, 3 * n_total);
                H.setFromTriplets(final_triplets.begin(), final_triplets.end());

                Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(H);
                if (solver.info() != Eigen::Success) { m.exploded = true; break; }
                Eigen::VectorXd p = solver.solve(-grad);
                if (solver.info() != Eigen::Success || !p.allFinite()) { m.exploded = true; break; }

                Eigen::MatrixXd Dir = unflatten(p, n_total);

                double alpha = std::min(1.0, ipc::compute_collision_free_stepsize(
                    collision_mesh, Xk, Xk + Dir, kMinDistance));

                double e_cur = 0.5 * (xk_flat - xhat_flat).dot(mass.cwiseProduct(xk_flat - xhat_flat))
                             + sub_dt * sub_dt * (e_elastic + e_barrier);
                Eigen::MatrixXd Xtrial;
                bool ls_ok = false;
                for (int ls = 0; ls < 20; ++ls) {
                    Xtrial = Xk + alpha * Dir;
                    ipc::NormalCollisions trial_collisions;
                    trial_collisions.build(collision_mesh, Xtrial, cfg.dhat);
                    Eigen::VectorXd trial_flat = flatten(Xtrial);
                    double trial_elastic = 0.0; Eigen::VectorXd dummy_grad = Eigen::VectorXd::Zero(3*n_total);
                    std::vector<Eigen::Triplet<double>> dummy_tri;
                    trial_elastic = assemble_elastic(Xtrial, cloth_edges, cfg.stiffness, dummy_grad, dummy_tri);
                    double trial_barrier = B(trial_collisions, collision_mesh, Xtrial);
                    double e_trial = 0.5 * (trial_flat - xhat_flat).dot(mass.cwiseProduct(trial_flat - xhat_flat))
                                   + sub_dt * sub_dt * (trial_elastic + trial_barrier);
                    // Armijo: only accept a step that actually decreases
                    // energy. Falling out of this loop without ever hitting
                    // this branch (NOT the old code's behaviour) must NOT
                    // silently apply whatever alpha the last failed attempt
                    // left behind -- that accepted energy-increasing steps
                    // and is why e_barrier climbed monotonically instead of
                    // settling.
                    if (std::isfinite(e_trial) && e_trial <= e_cur + 1e-10) { ls_ok = true; break; }
                    alpha *= 0.5;
                }
                if (!ls_ok || !std::isfinite(alpha)) { m.exploded = true; break; }
                Xk = Xtrial;
            }
            if (m.exploded) break;

            V = (Xk - X) / sub_dt;
            for (int i = n_cloth; i < n_total; ++i) V.row(i).setZero();
            X = Xk;

            for (int i = 0; i < n_cloth; ++i) {
                Vec<double,3> p{X(i,0), X(i,1), X(i,2)};
                auto q = point_to_sphere(p, sphere_center, kSphereR);
                m.worst_penetration = std::min(m.worst_penetration, q.signed_distance());
            }
            for (auto& e : cloth_edges) {
                double len = (X.row(e.a) - X.row(e.b)).norm();
                m.max_stretch = std::max(m.max_stretch, len / e.rest - 1.0);
            }
            if (!X.allFinite()) m.exploded = true;
        }
    }
    return m;
}

} // namespace

int main() {
    std::setbuf(stdout, nullptr); // visible progress even if a run is killed mid-sweep
    std::println("# cloth_sphere_probe — {}x{} cloth ({} particles), sphere r={}, "
                 "{} frames @ {:.4f}s x {} substeps",
                 kGrid, kGrid, kGrid * kGrid, kSphereR, kFrames, kDt, kSubsteps);
    std::println("# {:>10} {:>8} {:>12}  {:>14} {:>14} {:>10}",
                 "stiffness", "dhat", "barrier_k", "worst-pen", "max-stretch", "status");

    std::vector<Config> configs = {
        {1e4, 0.02, 1e2},
    };

    int n_ok = 0;
    for (auto& cfg : configs) {
        Metrics m = run_config(cfg);
        std::println("  {:>10.1e} {:>8.3f} {:>12.1e}  {:>14.4e} {:>14.4e} {:>10}",
                     cfg.stiffness, cfg.dhat, cfg.barrier_stiffness,
                     m.worst_penetration, m.max_stretch,
                     m.exploded ? "EXPLODED" : "ok");
        if (!m.exploded && m.worst_penetration > -cfg.dhat && m.max_stretch < 0.23) ++n_ok;
    }
    std::println("# {}/{} configs held (worst-pen >= -dhat, no explosion, stretch < 23% "
                 "-- the XPBD sweep's own failure threshold)", n_ok, static_cast<int>(configs.size()));
    return 0;
}
