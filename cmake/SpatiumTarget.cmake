# SpatiumTarget.cmake — helper for example/demo executables.
#
# Usage:
#   spatium_add_example(NAME <target>
#                        [SOURCES <file.cpp> ...]     # defaults to <target>.cpp
#                        [LIBS <target> ...]           # defaults to spatium_core
#                        [OPTIONS <flag> ...]          # extra compile options, e.g. -O3 -march=native
#                        [INCLUDE_DIRS <dir> ...])
#
# Every Spatium example target repeats the same three lines
# (add_executable + target_link_libraries + -Wall -Wextra -Wpedantic);
# this collects that pattern in one place instead of copy-pasting it once
# per demo. Per-demo specifics (extra libs, -O3 -march=native for the
# offline raytracers, Vulkan include dirs) stay as explicit arguments,
# not hidden defaults.

function(spatium_add_example)
    set(_opts)
    set(_single NAME)
    set(_multi SOURCES LIBS OPTIONS INCLUDE_DIRS)
    cmake_parse_arguments(EX "${_opts}" "${_single}" "${_multi}" ${ARGN})

    if(NOT EX_NAME)
        message(FATAL_ERROR "spatium_add_example: NAME required")
    endif()
    if(NOT EX_SOURCES)
        set(EX_SOURCES ${EX_NAME}.cpp)
    endif()
    if(NOT EX_LIBS)
        set(EX_LIBS spatium_core)
    endif()

    add_executable(${EX_NAME} ${EX_SOURCES})
    target_link_libraries(${EX_NAME} PRIVATE ${EX_LIBS})
    # -Wno-missing-field-initializers: stb_image_write.h's own
    # `= { 0 }` zero-init idiom trips this under -Wextra in every example
    # that writes PNGs (render/write_image.hpp's one shared include site)
    # -- vendor code we don't own and wouldn't want to edit to appease it.
    target_compile_options(${EX_NAME} PRIVATE
        -Wall -Wextra -Wpedantic -Wno-missing-field-initializers ${EX_OPTIONS})

    if(EX_INCLUDE_DIRS)
        target_include_directories(${EX_NAME} PRIVATE ${EX_INCLUDE_DIRS})
    endif()
endfunction()
