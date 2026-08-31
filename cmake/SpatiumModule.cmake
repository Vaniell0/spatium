# SpatiumModule.cmake — helper for named C++23 module targets.
#
# Usage:
#   spatium_module(NAME <target>
#                  PARTITIONS <file.cppm> [file.cppm ...]
#                  [DEPS <target> [target ...]]
#                  [IMPL <file.cpp> [file.cpp ...]]
#                  [KIND STATIC|SHARED|INTERFACE])
#
# - Primary module file must be listed first in PARTITIONS.
# - The helper applies -fmodules-ts, scans for modules via FILE_SET CXX_MODULES,
#   links std module (via the global `spatium_std` target populated in the
#   top-level CMakeLists), and links any DEPS.
# - IMPL sources are normal .cpp translation units that may `import` the module.

function(spatium_module)
    set(_opts)
    set(_single NAME KIND)
    set(_multi PARTITIONS DEPS IMPL)
    cmake_parse_arguments(SM "${_opts}" "${_single}" "${_multi}" ${ARGN})

    if(NOT SM_NAME)
        message(FATAL_ERROR "spatium_module: NAME required")
    endif()
    if(NOT SM_PARTITIONS)
        message(FATAL_ERROR "spatium_module(${SM_NAME}): PARTITIONS required")
    endif()
    if(NOT SM_KIND)
        set(SM_KIND STATIC)
    endif()

    add_library(${SM_NAME} ${SM_KIND} ${SM_IMPL})

    target_sources(${SM_NAME}
        PUBLIC
            FILE_SET CXX_MODULES
            BASE_DIRS ${CMAKE_SOURCE_DIR}/modules
            FILES ${SM_PARTITIONS}
    )

    target_compile_features(${SM_NAME} PUBLIC cxx_std_23)
    target_compile_options(${SM_NAME} PUBLIC -fmodules-ts)

    # Every module consumes `import std;`, so depend on the std BMI target.
    if(TARGET spatium_std)
        target_link_libraries(${SM_NAME} PUBLIC spatium_std)
    endif()

    if(SM_DEPS)
        target_link_libraries(${SM_NAME} PUBLIC ${SM_DEPS})
    endif()

    set_target_properties(${SM_NAME} PROPERTIES
        POSITION_INDEPENDENT_CODE ON
    )
endfunction()
