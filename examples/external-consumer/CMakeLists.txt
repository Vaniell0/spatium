# Standalone downstream project consuming Spatium via FetchContent.
#
# This is a REAL, independent CMake project -- it is not wired into
# Spatium's own examples/CMakeLists.txt or spatium_add_example() machinery.
# It exists to prove (and document) the one-command path for pulling
# Spatium into an unrelated project without Nix and without a manual
# `cmake --install` + `find_package(... PATHS ...)` step. See README.md in
# this directory for build instructions.
cmake_minimum_required(VERSION 3.28)
project(spatium_external_consumer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FetchContent)
FetchContent_Declare(spatium
    GIT_REPOSITORY https://github.com/Vaniell0/spatium.git
    GIT_TAG main
)
FetchContent_MakeAvailable(spatium)

add_executable(sphere_distance_demo main.cpp)
target_link_libraries(sphere_distance_demo PRIVATE Spatium::sdk)
