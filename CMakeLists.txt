cmake_minimum_required(VERSION 3.28)
project(spatium VERSION 1.0.0 LANGUAGES CXX
    DESCRIPTION "C++23 mathematical spaces, geometry and mesh library"
    HOMEPAGE_URL "https://github.com/Vaniell0/spatium")

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ── RSC reproducibility: commit SHA baked in at configure time ────
# rsc/checkpoint.hpp pins a saved base checkpoint to the exact commit it
# was trained against (ground truth IS Spatium's own code, so a Spatium
# change alters both the training environment and the definition of
# "correct" at once -- see rsc/README.md's Reproducibility section).
find_package(Git QUIET)
if(GIT_FOUND)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE SPATIUM_GIT_SHA
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
endif()
if(NOT SPATIUM_GIT_SHA)
    set(SPATIUM_GIT_SHA "unknown")
endif()
configure_file(
    ${CMAKE_SOURCE_DIR}/rsc/include/build_info.hpp.in
    ${CMAKE_BINARY_DIR}/generated/build_info.hpp
    @ONLY
)

# ── RSC base weights: embedded into the binary, not loaded at runtime ──
# Base is compiled in (this generated header), the per-deployment custom
# layer (rsc/include/custom_layer.hpp) is a real, separate, mutable file
# instead -- it has to update in the field without recompiling, the base
# doesn't. Re-runs at every configure so the embedded weights always
# match whatever's checked into rsc/checkpoints/.
set(SPATIUM_RSC_BASE_CHECKPOINT
    ${CMAKE_SOURCE_DIR}/rsc/checkpoints/base_v1.checkpoint
    CACHE FILEPATH "Base checkpoint to embed into the binary at configure time")
set(SPATIUM_EMBEDDED_BASE_HEADER ${CMAKE_BINARY_DIR}/generated/embedded_base_data.hpp)
if(EXISTS ${SPATIUM_RSC_BASE_CHECKPOINT})
    find_package(Python3 COMPONENTS Interpreter QUIET)
    if(Python3_FOUND)
        execute_process(
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_SOURCE_DIR}/rsc/tools/embed_checkpoint.py
                ${SPATIUM_RSC_BASE_CHECKPOINT}
                ${SPATIUM_EMBEDDED_BASE_HEADER}
            RESULT_VARIABLE SPATIUM_EMBED_RESULT
        )
        if(NOT SPATIUM_EMBED_RESULT EQUAL 0)
            message(WARNING "rsc: embed_checkpoint.py failed, embedded base will be unavailable")
        endif()
    else()
        message(WARNING "rsc: python3 not found, cannot embed base checkpoint")
    endif()
else()
    message(STATUS "rsc: no checkpoint at ${SPATIUM_RSC_BASE_CHECKPOINT} yet "
                    "(run rsc/tools/train_base, then re-run cmake) -- embedded base unavailable")
    file(WRITE ${SPATIUM_EMBEDDED_BASE_HEADER}
        "#pragma once\n#define SPATIUM_HAS_EMBEDDED_BASE 0\n")
endif()

# ── C++23 modules infrastructure ───────────────────────────────
# Caught up 2026-08-28 to cover the full header tree: core/algebra/spaces/
# mesh/geometry/point/spatial/discrete/io/physics/render partitions in
# modules/. algebra/{dual,ode,linear_solve,calculus}.hpp and
# physics/relativity/* (added since the 2026-04-24 migration, plus
# physics/atomic/ after the 2026-08-28 split) now have partitions.
# Two headers deliberately stay header-only, same as before, for
# documented reasons rather than neglect: algebra/eigen_interop.hpp (SSE
# intrinsics pulled by Eigen/Core conflict with vec_simd's BMI — see
# algebra.cppm) and render/write_image.hpp (its
# STB_IMAGE_WRITE_IMPLEMENTATION contract is a link-time, not module,
# concern — see render.cppm). viewer/* also stays out: Vulkan/GLFW/ImGui
# pull C-linkage that would hit the same kind of BMI conflict.
# Off by default not because of open compiler bugs, but because it still
# isn't wired into CI (no workflow passes this ON yet).
option(SPATIUM_USE_MODULES "Enable C++23 named modules build (behind current headers, see comment above)" OFF)

if(SPATIUM_USE_MODULES)
    if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
       OR CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14.0)
        message(FATAL_ERROR
            "SPATIUM_USE_MODULES requires GCC >= 14 (got ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION})")
    endif()
    set(CMAKE_CXX_SCAN_FOR_MODULES ON)
    list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
    include(SpatiumModule)

    # `import std;` BMI producer. GCC has no CMAKE_CXX_MODULE_STD=ON support yet,
    # so we ship a local copy of bits/std.cc as modules/std.cppm.
    add_library(spatium_std STATIC)
    target_sources(spatium_std
        PUBLIC
            FILE_SET CXX_MODULES
            BASE_DIRS ${CMAKE_SOURCE_DIR}/modules
            FILES
                ${CMAKE_SOURCE_DIR}/modules/std.cppm
                ${CMAKE_SOURCE_DIR}/modules/std.compat.cppm
    )
    target_compile_features(spatium_std PUBLIC cxx_std_23)
    target_compile_options(spatium_std PUBLIC -fmodules-ts)
    set_target_properties(spatium_std PROPERTIES POSITION_INDEPENDENT_CODE ON)

    # ── spatium.core (Phase 1) ─────────────────────────────────
    # Primary module first, then partitions in dep order. Every .cppm wraps
    # the matching .hpp via SPATIUM_BUILDING_MODULE; the helper supplies the
    # `import std;` link.
    spatium_module(NAME spatium_core_module
        PARTITIONS
            ${CMAKE_SOURCE_DIR}/modules/core.cppm
            ${CMAKE_SOURCE_DIR}/modules/core/concepts.cppm
            ${CMAKE_SOURCE_DIR}/modules/core/error.cppm
            ${CMAKE_SOURCE_DIR}/modules/core/epsilon.cppm
            ${CMAKE_SOURCE_DIR}/modules/core/precision.cppm
            ${CMAKE_SOURCE_DIR}/modules/core/verify.cppm
        DEPS spatium_sdk)
    add_library(Spatium::core_module ALIAS spatium_core_module)

    # ── spatium.algebra (Phase 1) ──────────────────────────────
    spatium_module(NAME spatium_algebra_module
        PARTITIONS
            ${CMAKE_SOURCE_DIR}/modules/algebra.cppm
            ${CMAKE_SOURCE_DIR}/modules/algebra/concepts.cppm
            ${CMAKE_SOURCE_DIR}/modules/algebra/vec_simd.cppm
            ${CMAKE_SOURCE_DIR}/modules/algebra/vector.cppm
            ${CMAKE_SOURCE_DIR}/modules/algebra/matrix.cppm
            ${CMAKE_SOURCE_DIR}/modules/algebra/format.cppm
            ${CMAKE_SOURCE_DIR}/modules/algebra/literals.cppm
            ${CMAKE_SOURCE_DIR}/modules/algebra/complex.cppm
            ${CMAKE_SOURCE_DIR}/modules/algebra/polynomial.cppm
            ${CMAKE_SOURCE_DIR}/modules/algebra/functions.cppm
            ${CMAKE_SOURCE_DIR}/modules/algebra/quaternion.cppm
            ${CMAKE_SOURCE_DIR}/modules/algebra/verify.cppm
            ${CMAKE_SOURCE_DIR}/modules/algebra/groups_so3.cppm
            ${CMAKE_SOURCE_DIR}/modules/algebra/groups_se3.cppm
            ${CMAKE_SOURCE_DIR}/modules/algebra/dual.cppm
            ${CMAKE_SOURCE_DIR}/modules/algebra/ode.cppm
            ${CMAKE_SOURCE_DIR}/modules/algebra/linear_solve.cppm
            ${CMAKE_SOURCE_DIR}/modules/algebra/calculus.cppm
        DEPS spatium_sdk spatium_core_module)
    add_library(Spatium::algebra_module ALIAS spatium_algebra_module)

    # ── spatium.mesh (Phase 4 — single partition, upstream of spaces/geometry)
    # Excludes primitives.hpp (depends on spatium.spaces) and operations.hpp
    # (depends on spatium.geometry) to avoid cycles; both stay header-only.
    spatium_module(NAME spatium_mesh_module
        PARTITIONS ${CMAKE_SOURCE_DIR}/modules/mesh.cppm
        DEPS spatium_sdk spatium_core_module spatium_algebra_module)
    add_library(Spatium::mesh_module ALIAS spatium_mesh_module)

    # ── spatium.render: single partition, see render.cppm for the
    # write_image.hpp exclusion rationale (STB implementation macro is a
    # link-time contract, not a module one).
    spatium_module(NAME spatium_render_module
        PARTITIONS ${CMAKE_SOURCE_DIR}/modules/render.cppm
        DEPS spatium_sdk spatium_core_module spatium_algebra_module)
    add_library(Spatium::render_module ALIAS spatium_render_module)

    # ── spatium.spaces (Phase 2): single partition, see spaces.cppm comment.
    spatium_module(NAME spatium_spaces_module
        PARTITIONS ${CMAKE_SOURCE_DIR}/modules/spaces.cppm
        DEPS spatium_sdk spatium_core_module spatium_algebra_module
             spatium_mesh_module)
    add_library(Spatium::spaces_module ALIAS spatium_spaces_module)

    # ── spatium.point (Phase 2): point + morphism
    spatium_module(NAME spatium_point_module
        PARTITIONS
            ${CMAKE_SOURCE_DIR}/modules/point_primary.cppm
            ${CMAKE_SOURCE_DIR}/modules/point.cppm
            ${CMAKE_SOURCE_DIR}/modules/morphism.cppm
        DEPS spatium_sdk spatium_core_module spatium_algebra_module)
    add_library(Spatium::point_module ALIAS spatium_point_module)

    # ── spatium.geometry (Phase 3): single partition, see geometry.cppm
    spatium_module(NAME spatium_geometry_module
        PARTITIONS ${CMAKE_SOURCE_DIR}/modules/geometry.cppm
        DEPS spatium_sdk spatium_core_module spatium_algebra_module
             spatium_mesh_module spatium_spaces_module spatium_point_module)
    add_library(Spatium::geometry_module ALIAS spatium_geometry_module)

    # ── spatium.spatial / discrete / io (Phase 5)
    spatium_module(NAME spatium_spatial_module
        PARTITIONS ${CMAKE_SOURCE_DIR}/modules/spatial.cppm
        DEPS spatium_sdk spatium_core_module spatium_algebra_module
             spatium_geometry_module)
    add_library(Spatium::spatial_module ALIAS spatium_spatial_module)

    spatium_module(NAME spatium_discrete_module
        PARTITIONS ${CMAKE_SOURCE_DIR}/modules/discrete.cppm
        DEPS spatium_sdk spatium_core_module spatium_algebra_module
             spatium_geometry_module)
    add_library(Spatium::discrete_module ALIAS spatium_discrete_module)

    spatium_module(NAME spatium_io_module
        PARTITIONS ${CMAKE_SOURCE_DIR}/modules/io.cppm
        DEPS spatium_sdk spatium_core_module spatium_algebra_module
             spatium_spaces_module spatium_mesh_module)
    add_library(Spatium::io_module ALIAS spatium_io_module)

    # ── spatium.physics (Phase 6): interface + impl unit
    # viewer/* stays out of modules: Vulkan/GLFW/ImGui pull C-linkage and
    # would conflict with algebra :vec_simd BMI (same pitfall as :eigen).
    spatium_module(NAME spatium_physics_module
        PARTITIONS ${CMAKE_SOURCE_DIR}/modules/physics.cppm
        IMPL       ${CMAKE_SOURCE_DIR}/modules/physics_impl.cpp
        DEPS spatium_sdk spatium_core_module spatium_algebra_module
             spatium_spaces_module spatium_mesh_module spatium_io_module
             spatium_geometry_module)
    add_library(Spatium::physics_module ALIAS spatium_physics_module)

    message(STATUS "Spatium: C++23 modules ENABLED (GCC ${CMAKE_CXX_COMPILER_VERSION})")
else()
    message(STATUS "Spatium: C++23 modules DISABLED (set -DSPATIUM_USE_MODULES=ON to enable)")
endif()

# ── SDK (INTERFACE — header-only) ──────────────────────────────
add_library(spatium_sdk INTERFACE)
add_library(Spatium::sdk ALIAS spatium_sdk)
set_target_properties(spatium_sdk PROPERTIES EXPORT_NAME sdk)

target_include_directories(spatium_sdk INTERFACE
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_compile_features(spatium_sdk INTERFACE cxx_std_23)

# ── SIMD (opt-in native arch) ─────────────────────────────────
option(SPATIUM_NATIVE_ARCH "Use -march=native for SIMD (NOT portable across CPUs)" OFF)
if(SPATIUM_NATIVE_ARCH)
    target_compile_options(spatium_sdk INTERFACE -march=native)
    message(WARNING
        "SPATIUM_NATIVE_ARCH is ON: binaries are tuned for the build host's CPU "
        "(-march=native) and may fault with SIGILL on different machines. "
        "Disable for distributable builds, packaging, or CI.")
endif()

# ── Eigen3 interop (optional) ────────────────────────────────
option(SPATIUM_EIGEN "Enable Eigen3 interop (required for heat method)" OFF)
if(SPATIUM_EIGEN)
    find_package(Eigen3 3.4 REQUIRED)
    target_link_libraries(spatium_sdk INTERFACE Eigen3::Eigen)
    target_compile_definitions(spatium_sdk INTERFACE SPATIUM_HAS_EIGEN=1)
    message(STATUS "Spatium: Eigen3 ENABLED")
else()
    target_compile_definitions(spatium_sdk INTERFACE SPATIUM_HAS_EIGEN=0)
endif()

# ── ipc-toolkit (optional) ────────────────────────────────────
# Implicit contact backend for contact physics/soft bodies (Newton +
# log-barrier + line-search filter, Li-Kaufman 2020 §5) -- replaces
# explicit XPBD position-based contact, which can't hold stiff structural
# constraints against contact under a heavy overhang (20/20 swept configs
# failed, see contact.hpp).
# When SPATIUM_EIGEN is also ON, ipc-toolkit reuses Spatium's
# Eigen3::Eigen target instead of fetching its own (its eigen.cmake
# recipe checks `if(TARGET Eigen3::Eigen) return()` first).
#
# NOT linked into spatium_sdk: unlike Eigen3::Eigen (an IMPORTED target
# from find_package, safe to expose from an EXPORT-installed interface
# target), ipc::toolkit is built locally via FetchContent and isn't part
# of any install export set -- installing spatium_sdk INTERFACE-linked to
# it breaks `install(EXPORT SpatiumTargets ...)`. Same reasoning as
# Vulkan/GLFW living on spatium_viewer instead of spatium_sdk. Consumers
# (tests, the future contact-physics demo) link ipc::toolkit directly.
option(SPATIUM_IPC_TOOLKIT "Enable ipc-toolkit for implicit contact physics" OFF)
if(SPATIUM_IPC_TOOLKIT)
    include(FetchContent)
    FetchContent_Declare(
        ipc_toolkit
        GIT_REPOSITORY https://github.com/ipc-sim/ipc-toolkit.git
        GIT_TAG v1.6.0
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(ipc_toolkit)
    target_compile_definitions(spatium_sdk INTERFACE SPATIUM_HAS_IPC_TOOLKIT=1)
    message(STATUS "Spatium: ipc-toolkit ENABLED (v1.6.0)")
else()
    target_compile_definitions(spatium_sdk INTERFACE SPATIUM_HAS_IPC_TOOLKIT=0)
endif()

# ── Core library (SHARED) ──────────────────────────────────────
add_library(spatium_core SHARED
    src/physics/elements.cpp
)
add_library(Spatium::core ALIAS spatium_core)
set_target_properties(spatium_core PROPERTIES
    EXPORT_NAME core
    POSITION_INDEPENDENT_CODE ON
)
target_link_libraries(spatium_core PUBLIC spatium_sdk)
target_compile_options(spatium_core PRIVATE
    -Wall -Wextra -Wpedantic -Wno-unused-parameter
)

# ── Viewer library (Vulkan + GLFW) ────────────────────────────
option(SPATIUM_BUILD_VIEWER "Build Vulkan viewer" ON)
if(SPATIUM_BUILD_VIEWER)
    find_package(Vulkan QUIET)
    find_package(glfw3 3.3 QUIET)
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
        pkg_check_modules(SHADERC QUIET shaderc)
    endif()

    if(NOT Vulkan_FOUND OR NOT glfw3_FOUND OR NOT SHADERC_FOUND)
        set(_missing "")
        if(NOT Vulkan_FOUND)
            list(APPEND _missing "Vulkan")
        endif()
        if(NOT glfw3_FOUND)
            list(APPEND _missing "glfw3")
        endif()
        if(NOT SHADERC_FOUND)
            list(APPEND _missing "shaderc")
        endif()
        message(WARNING
            "SPATIUM_BUILD_VIEWER=ON but viewer dependencies are missing: ${_missing}. "
            "The viewer target will NOT be built; pass -DSPATIUM_BUILD_VIEWER=OFF "
            "to silence this warning, or install the missing packages.")
    endif()

    if(Vulkan_FOUND AND glfw3_FOUND AND SHADERC_FOUND)
        add_library(spatium_viewer SHARED
            src/viewer/app.cpp
            src/viewer/camera.cpp
        )
        add_library(Spatium::viewer ALIAS spatium_viewer)
        set_target_properties(spatium_viewer PROPERTIES
            EXPORT_NAME viewer
            POSITION_INDEPENDENT_CODE ON
        )
        target_link_libraries(spatium_viewer PUBLIC spatium_core Vulkan::Vulkan glfw)
        target_link_libraries(spatium_viewer PRIVATE ${SHADERC_LIBRARIES})
        target_include_directories(spatium_viewer PRIVATE ${SHADERC_INCLUDE_DIRS})
        target_compile_options(spatium_viewer PRIVATE
            -Wall -Wextra -Wpedantic -Wno-unused-parameter
        )
        # ImGui (optional, from source directory)
        if(NOT IMGUI_DIR AND DEFINED ENV{IMGUI_DIR})
            set(IMGUI_DIR $ENV{IMGUI_DIR})
        endif()
        if(IMGUI_DIR)
            set(IMGUI_SOURCES
                ${IMGUI_DIR}/imgui.cpp
                ${IMGUI_DIR}/imgui_draw.cpp
                ${IMGUI_DIR}/imgui_tables.cpp
                ${IMGUI_DIR}/imgui_widgets.cpp
                ${IMGUI_DIR}/backends/imgui_impl_glfw.cpp
                ${IMGUI_DIR}/backends/imgui_impl_vulkan.cpp
            )
            target_sources(spatium_viewer PRIVATE ${IMGUI_SOURCES})
            target_include_directories(spatium_viewer PUBLIC ${IMGUI_DIR} ${IMGUI_DIR}/backends)
            target_compile_definitions(spatium_viewer PUBLIC SPATIUM_HAS_IMGUI=1)
            set_source_files_properties(${IMGUI_SOURCES} PROPERTIES COMPILE_FLAGS "-w")
            set(SPATIUM_HAS_IMGUI TRUE)
            message(STATUS "Spatium: ImGui ENABLED (${IMGUI_DIR})")
        else()
            target_compile_definitions(spatium_viewer PUBLIC SPATIUM_HAS_IMGUI=0)
            set(SPATIUM_HAS_IMGUI FALSE)
            message(STATUS "Spatium: ImGui DISABLED (set IMGUI_DIR to enable)")
        endif()

        set(SPATIUM_HAS_VIEWER TRUE)
        message(STATUS "Spatium: Vulkan viewer ENABLED")
    else()
        set(SPATIUM_HAS_VIEWER FALSE)
        message(STATUS "Spatium: Vulkan viewer DISABLED (missing deps)")
    endif()
endif()

# ── Backward compat alias ─────────────────────────────────────
add_library(spatium ALIAS spatium_sdk)

# ── Tests ──────────────────────────────────────────────────────
option(SPATIUM_BUILD_TESTS "Build tests" ON)
if(SPATIUM_BUILD_TESTS)
    enable_testing()
    find_package(Catch2 3 REQUIRED)
    add_subdirectory(tests)
endif()

# ── Benchmarks ────────────────────────────────────────────────
option(SPATIUM_BUILD_BENCHMARKS "Build benchmarks" OFF)
if(SPATIUM_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

# ── Examples ───────────────────────────────────────────────────
option(SPATIUM_BUILD_EXAMPLES "Build examples" ON)
if(SPATIUM_BUILD_EXAMPLES)
    include(${CMAKE_SOURCE_DIR}/cmake/SpatiumTarget.cmake)
    add_subdirectory(examples)
endif()

# ── RSC tools ──────────────────────────────────────────────────
option(SPATIUM_BUILD_RSC_TOOLS "Build RSC tools (train_base, ...)" ON)
if(SPATIUM_BUILD_RSC_TOOLS)
    add_subdirectory(rsc/tools)
endif()

# ── CUDA geodesic kernel (optional) ───────────────────────────
# Closed-form (hand/symbolically derived, not the generic Dual<T> +
# general 4x4 invert() machinery geodesic.hpp uses on the CPU) RK4
# geodesic integration for Schwarzschild/Kerr, one thread per ray. Off
# by default: needs an actual CUDA toolchain (nvcc) and, to run
# anything, an NVIDIA GPU -- neither assumed present on a default build
# machine.
option(SPATIUM_CUDA "Build CUDA geodesic kernel (gpu/)" OFF)
if(SPATIUM_CUDA)
    enable_language(CUDA)
    add_subdirectory(gpu)
endif()

# ── Install ────────────────────────────────────────────────────
include(CMakePackageConfigHelpers)

set(_install_targets spatium_sdk spatium_core)
if(SPATIUM_HAS_VIEWER)
    list(APPEND _install_targets spatium_viewer)
endif()

install(TARGETS ${_install_targets} EXPORT SpatiumTargets
    LIBRARY     DESTINATION lib
    ARCHIVE     DESTINATION lib
    RUNTIME     DESTINATION bin
    INCLUDES    DESTINATION include
)
install(DIRECTORY include/spatium DESTINATION include
    FILES_MATCHING PATTERN "*.hpp"
)
install(EXPORT SpatiumTargets
    FILE SpatiumTargets.cmake
    NAMESPACE Spatium::
    DESTINATION lib/cmake/Spatium
)

# Track which optional pieces shipped so downstream consumers can react.
set(SPATIUM_INSTALLED_HAS_VIEWER ${SPATIUM_HAS_VIEWER})
set(SPATIUM_INSTALLED_HAS_EIGEN  ${SPATIUM_EIGEN})
set(SPATIUM_INSTALLED_HAS_IMGUI  ${SPATIUM_HAS_IMGUI})

configure_package_config_file(
    ${CMAKE_SOURCE_DIR}/cmake/SpatiumConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/SpatiumConfig.cmake
    INSTALL_DESTINATION lib/cmake/Spatium
)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/SpatiumConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/SpatiumConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/SpatiumConfigVersion.cmake
    DESTINATION lib/cmake/Spatium
)
