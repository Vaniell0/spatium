#pragma once

// Work-stealing row parallelism for CPU raytracers.
//
// Every raytracer example (blackhole_demo.cpp, wormhole_demo.cpp,
// parametric_analytical_demo.cpp's three render functions) independently
// reimplemented the same "compute nthreads, split rows into contiguous
// per-thread blocks, launch jthreads, join" boilerplate -- found
// duplicated during the 2026-08-26 architecture audit. Worse than just
// duplicated: static contiguous blocks assume every row costs the same,
// which is false for these renderers specifically -- rows crossing a
// black hole's photon ring, a wormhole's throat, or a Newton-UV
// surface's silhouette cost far more per-ray than an empty-sky row (see
// `ray_parametric.hpp`'s own "~568,000 ns/ray for miss pixels" note), so
// a static block can leave some threads idle while others are still
// grinding through the expensive rows.
//
// `parallel_for_rows` fixes both at once: one shared function here
// instead of four copies, and every thread pulls the next unclaimed row
// off a single atomic counter instead of owning a fixed range -- a
// thread that finishes its (cheap) row immediately steals the next one,
// so uneven per-row cost gets load-balanced automatically. No new
// dependency (no TBB, no std::execution::par -- libstdc++'s parallel
// algorithms need TBB as a backend, deliberately not pulled in here);
// this is plain std::jthread + std::atomic, same building blocks the
// four duplicated versions already used.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <algorithm>
#  include <atomic>
#  include <cstddef>
#  include <thread>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::render {

// Calls `row_fn(y)` once for every y in [0, height), across
// hardware_concurrency() worker threads, via atomic work-stealing.
// `row_fn` must be safe to call concurrently from multiple threads with
// different `y` values -- the usual contract for a raytracer whose rows
// write disjoint slices of one output buffer. Any per-row accumulation
// the caller needs (hit counts, captured-pixel counts, ...) should use
// its own `std::atomic` captured by reference inside `row_fn`; the
// per-row increment cost is negligible next to a row's actual ray-trace
// work.
template<class RowFn>
void parallel_for_rows(int height, RowFn&& row_fn) {
    std::atomic<int> next_row{0};
    unsigned nthreads = std::max(1u, std::thread::hardware_concurrency());
    nthreads = std::min(nthreads, static_cast<unsigned>(std::max(1, height)));

    std::vector<std::jthread> workers;
    workers.reserve(nthreads);
    for (unsigned t = 0; t < nthreads; ++t) {
        workers.emplace_back([&next_row, height, &row_fn] {
            int y;
            while ((y = next_row.fetch_add(1, std::memory_order_relaxed)) < height) {
                row_fn(y);
            }
        });
    }
}  // jthreads join here

}  // namespace spatium::render
