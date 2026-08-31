# Only reached when -DSPATIUM_CUDA=ON (see top-level CMakeLists.txt),
# which also does enable_language(CUDA) first. Flat-directory layout
# matching benchmarks/'s own convention (one CMakeLists.txt here, no
# further nesting).

add_library(spatium_gpu_kernel STATIC geodesic_kernel.cu)
set_target_properties(spatium_gpu_kernel PROPERTIES
    CUDA_ARCHITECTURES "75"  # T4 = sm_75; widen once other GPUs matter
    CUDA_SEPARABLE_COMPILATION ON
    POSITION_INDEPENDENT_CODE ON
    # Without an explicit CUDA_STANDARD, CMake propagates the top-level
    # CMAKE_CXX_STANDARD (23) onto CUDA_STANDARD too -- nvcc (up through
    # at least 12.9) has no c++23 device dialect, so that fails at
    # configure time ("requires the language dialect CUDA23"). The .cu
    # files here only use c++20 features anyway (matches every manual
    # `nvcc -std=c++20 ...` invocation this project's VM builds have
    # used from the start).
    CUDA_STANDARD 20
    CUDA_STANDARD_REQUIRED ON
)
target_include_directories(spatium_gpu_kernel PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(verify_cuda_geodesic verify_cuda_geodesic.cpp)
target_link_libraries(verify_cuda_geodesic PRIVATE spatium_gpu_kernel spatium_sdk)

# Phase 3/4: the full per-ray render kernels (geodesic + volumetric disk
# emission-absorption) -- previously only ever built by hand-typed
# `g++`/`nvcc` command sequences run over SSH on the render VM each
# session. Wired into CMake here so the ONLY thing that changes between
# a local build (compile-check, no GPU needed) and a run on the VM is
# which machine `cmake --build` runs on -- rsync the repo, configure
# once with -DSPATIUM_CUDA=ON, build. No manual per-file compiler
# invocations, no post-rsync header patching (see io_helpers.hpp's and
# write_image.hpp's own <print> fallback guards for why that's gone too).
add_library(spatium_gpu_render_kernel STATIC
    schwarzschild_render_kernel.cu
    kerr_render_kernel.cu
)
set_target_properties(spatium_gpu_render_kernel PROPERTIES
    CUDA_ARCHITECTURES "75"
    CUDA_SEPARABLE_COMPILATION ON
    POSITION_INDEPENDENT_CODE ON
    CUDA_STANDARD 20
    CUDA_STANDARD_REQUIRED ON
)
target_include_directories(spatium_gpu_render_kernel PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(verify_cuda_render verify_cuda_render.cpp)
target_link_libraries(verify_cuda_render PRIVATE spatium_gpu_render_kernel spatium_sdk)

add_executable(verify_cuda_kerr_render verify_cuda_kerr_render.cpp)
target_link_libraries(verify_cuda_kerr_render PRIVATE spatium_gpu_render_kernel spatium_sdk)

# The actual final-video render driver (examples/blackhole_gr_demo.cpp's
# choreography, executed per-ray on GPU). Needs -O2: unoptimized host-side
# per-pixel compositing/sky-lookup is the dominant cost at 1920x1080, not
# the GPU kernel launch itself.
add_executable(render_4k_frame render_4k_frame.cpp)
target_link_libraries(render_4k_frame PRIVATE spatium_gpu_render_kernel spatium_sdk)
target_compile_options(render_4k_frame PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-O2>)
