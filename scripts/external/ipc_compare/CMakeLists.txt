# Out-of-tree comparison of Spatium against ipc-toolkit (Li-Kaufman 2020
# implicit contact: Newton + log barrier + CCD line search) on the
# cloth-on-sphere scene. ipc-toolkit is NOT a Spatium dependency; this
# project fetches it itself and is never part of the Spatium build or CI.
#
#   cmake -S scripts/external/ipc_compare -B /tmp/ipc_compare -G Ninja \
#         -DCMAKE_BUILD_TYPE=Release
#   ninja -C /tmp/ipc_compare && /tmp/ipc_compare/cloth_sphere_probe
#
# Last in-tree version: git tag ipc-toolkit-final.
cmake_minimum_required(VERSION 3.25)
project(spatium_ipc_compare CXX)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
    ipc_toolkit
    GIT_REPOSITORY https://github.com/ipc-sim/ipc-toolkit.git
    GIT_TAG v1.6.0
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(ipc_toolkit)

get_filename_component(SPATIUM_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../.." ABSOLUTE)
add_executable(cloth_sphere_probe cloth_sphere_probe.cpp)
target_include_directories(cloth_sphere_probe PRIVATE "${SPATIUM_ROOT}/include")
target_link_libraries(cloth_sphere_probe PRIVATE ipc::toolkit)
