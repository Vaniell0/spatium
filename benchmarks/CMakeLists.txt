find_package(benchmark REQUIRED)

add_executable(spatium_bench
    bench_vec.cpp
    bench_intersection.cpp
    bench_mesh.cpp
    bench_bvh.cpp
    bench_orbital.cpp
    bench_raycast.cpp
    bench_narrow_phase.cpp
)

target_link_libraries(spatium_bench PRIVATE
    spatium_core
    benchmark::benchmark_main
)

target_compile_options(spatium_bench PRIVATE
    -Wall -Wextra -Wpedantic -Wno-unused-parameter -O3
)
