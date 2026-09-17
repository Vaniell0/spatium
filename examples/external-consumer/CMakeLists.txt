# Standalone downstream project consuming Spatium two different ways.
#
# This is a REAL, independent CMake project -- it is not wired into
# Spatium's own examples/CMakeLists.txt or spatium_add_example() machinery.
# It exists to prove (and document) how an unrelated project picks Spatium
# up, and there are exactly two such ways:
#
#   -DSPATIUM_CONSUME=fetchcontent  (default) pull the sources at configure
#                                   time; no install step, nothing on the
#                                   system, the one-command path.
#   -DSPATIUM_CONSUME=find_package  consume an *installed* Spatium, the way
#                                   a system or distro package is consumed.
#
# One project and one main.cpp for both, deliberately. They are the same
# program asking the same question of the same library; a second directory
# would be a second thing to keep in step, and the first time it drifted
# nobody would notice, because a build that is not run is not a check.
#
# The find_package path is here because it was the untested one. Spatium
# has shipped install(EXPORT SpatiumTargets), a generated SpatiumConfig and
# a version file since before its first release, and CI installed none of
# it -- so the package config was, in the exact sense conventions.md means,
# a mechanism with nothing running it. See README.md for both commands.
cmake_minimum_required(VERSION 3.28)
project(spatium_external_consumer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(SPATIUM_CONSUME "fetchcontent" CACHE STRING
    "How to obtain Spatium: fetchcontent | find_package")
set_property(CACHE SPATIUM_CONSUME PROPERTY STRINGS fetchcontent find_package)

if(SPATIUM_CONSUME STREQUAL "find_package")
    # Point at an install tree with -DCMAKE_PREFIX_PATH=<prefix>.
    find_package(Spatium REQUIRED)
    message(STATUS "spatium_external_consumer: found installed Spatium ${Spatium_VERSION}")
elseif(SPATIUM_CONSUME STREQUAL "fetchcontent")
    include(FetchContent)
    FetchContent_Declare(spatium
        GIT_REPOSITORY https://github.com/Vaniell0/spatium.git
        GIT_TAG main
    )
    FetchContent_MakeAvailable(spatium)
else()
    message(FATAL_ERROR
        "SPATIUM_CONSUME must be 'fetchcontent' or 'find_package', got '${SPATIUM_CONSUME}'")
endif()

add_executable(sphere_distance_demo main.cpp)
target_link_libraries(sphere_distance_demo PRIVATE Spatium::sdk)
