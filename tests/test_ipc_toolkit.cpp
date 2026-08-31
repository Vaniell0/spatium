#include <catch2/catch_test_macros.hpp>

#if SPATIUM_HAS_IPC_TOOLKIT

#include <ipc/collision_mesh.hpp>
#include <Eigen/Dense>

TEST_CASE("ipc-toolkit CollisionMesh builds from a single triangle", "[ipc_toolkit]") {
    Eigen::MatrixXd vertices(3, 3);
    vertices << 0.0, 0.0, 0.0,
                1.0, 0.0, 0.0,
                0.0, 1.0, 0.0;
    Eigen::MatrixXi edges(3, 2);
    edges << 0, 1,
              1, 2,
              0, 2;
    Eigen::MatrixXi faces(1, 3);
    faces << 0, 1, 2;

    ipc::CollisionMesh mesh(vertices, edges, faces);

    CHECK(mesh.num_vertices() == 3);
    CHECK(mesh.num_faces() == 1);
    CHECK(mesh.dim() == 3);
}

#endif
