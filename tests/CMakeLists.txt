add_executable(spatium_tests
    test_concepts.cpp
    test_vector.cpp
    test_euclidean.cpp
    test_box.cpp
    test_line.cpp
    test_hyperplane.cpp
    test_triangle.cpp
    test_simplex.cpp
    test_circle.cpp
    test_polygon.cpp
    test_intersection.cpp
    test_sphere.cpp
    test_hyperbolic.cpp
    test_morphism.cpp
    test_format.cpp
    test_mesh.cpp
    test_precision.cpp
    test_verify.cpp
    test_hardening.cpp
    test_matrix.cpp
    test_transform.cpp
    test_product.cpp
    test_convex_hull.cpp
    test_algebra.cpp
    test_surface_adapter.cpp
    test_clip.cpp
    test_boolean.cpp
    test_bvh.cpp
    test_topology.cpp
    test_geodesic.cpp
    test_transport.cpp
    test_voronoi.cpp
    test_parametric.cpp
    test_mesh_ops.cpp
    test_quaternion.cpp
    test_implicit.cpp
    test_adaptive.cpp
    test_svg_projection.cpp
    test_orbital.cpp
    test_elements.cpp
    test_literals.cpp
    test_distance.cpp
    test_complex.cpp
    test_dual.cpp
    test_calculus.cpp
    test_ode.cpp
    test_linear_solve.cpp
    test_render.cpp
    test_polynomial.cpp
    test_ray_surface.cpp
    test_ray_parametric.cpp
    test_eigen_interop.cpp
    test_heat_geodesic.cpp
    test_atom_upgrade.cpp
    test_mechanics.cpp
    test_dec.cpp
    test_lie_integrator.cpp
    test_manifold_body.cpp
    test_variational.cpp
    test_symplectic_integrators.cpp
    test_geometric_mechanics.cpp
    test_contact_xpbd.cpp
    test_ipc_toolkit.cpp
    test_narrow_phase.cpp
    test_relativity.cpp
    test_rsc_registry.cpp
    test_rsc_task.cpp
    test_rsc_dispatcher.cpp
    test_rsc_describe.cpp
    test_rsc_train.cpp
    test_rsc_chain.cpp
    test_rsc_precision.cpp
    test_rsc_geodesic.cpp
    test_rsc_rootfind.cpp
    test_rsc_ode.cpp
    test_rsc_integrator.cpp
    test_rsc_mesh.cpp
    test_rsc_linear.cpp
    test_rsc_base.cpp
    test_rsc_checkpoint.cpp
    test_rsc_embedded_base.cpp
    test_rsc_custom_layer.cpp
)

target_include_directories(spatium_tests PRIVATE
    ${CMAKE_SOURCE_DIR}/rsc/include
    ${CMAKE_BINARY_DIR}/generated
    ${CMAKE_SOURCE_DIR}/gpu
)

target_link_libraries(spatium_tests PRIVATE
    spatium_core
    Catch2::Catch2WithMain
)

if(SPATIUM_IPC_TOOLKIT)
    target_link_libraries(spatium_tests PRIVATE ipc::toolkit)
endif()

target_compile_options(spatium_tests PRIVATE
    -Wall -Wextra -Wpedantic -Wno-unused-parameter
)

# ── Modules tests ──────────────────────────────────────────────
# Compiled only when modules are enabled. Each phase consumes the matching
# module via `import`; the smoke check (test_modules_smoke.cpp) imports
# spatium.core and proves the helper + std.cppm chain is wired correctly.
if(SPATIUM_USE_MODULES AND TARGET spatium_core_module)
    target_sources(spatium_tests PRIVATE test_modules_smoke.cpp)
endif()
if(SPATIUM_USE_MODULES AND TARGET spatium_algebra_module)
    target_sources(spatium_tests PRIVATE test_modules_phase1.cpp)
    target_link_libraries(spatium_tests PRIVATE spatium_algebra_module)
endif()
if(SPATIUM_USE_MODULES AND TARGET spatium_spaces_module)
    target_sources(spatium_tests PRIVATE test_modules_phase2.cpp)
    target_link_libraries(spatium_tests PRIVATE
        spatium_mesh_module
        spatium_spaces_module
        spatium_point_module)
endif()
if(SPATIUM_USE_MODULES AND TARGET spatium_geometry_module)
    target_sources(spatium_tests PRIVATE test_modules_phase3.cpp)
    target_link_libraries(spatium_tests PRIVATE spatium_geometry_module)
endif()
if(SPATIUM_USE_MODULES AND TARGET spatium_mesh_module)
    target_sources(spatium_tests PRIVATE test_modules_phase4.cpp)
    # spatium_mesh_module already linked via phase2 deps; listed for clarity.
    target_link_libraries(spatium_tests PRIVATE spatium_mesh_module)
endif()
if(SPATIUM_USE_MODULES AND TARGET spatium_io_module)
    target_sources(spatium_tests PRIVATE test_modules_phase5.cpp)
    target_link_libraries(spatium_tests PRIVATE
        spatium_spatial_module
        spatium_discrete_module
        spatium_io_module)
endif()
if(SPATIUM_USE_MODULES AND TARGET spatium_physics_module)
    target_sources(spatium_tests PRIVATE test_modules_phase6.cpp)
    target_link_libraries(spatium_tests PRIVATE spatium_physics_module)
endif()
if(SPATIUM_USE_MODULES AND TARGET spatium_render_module)
    target_sources(spatium_tests PRIVATE test_modules_render.cpp)
    target_link_libraries(spatium_tests PRIVATE spatium_render_module)
endif()

include(CTest)
include(Catch)
catch_discover_tests(spatium_tests)
