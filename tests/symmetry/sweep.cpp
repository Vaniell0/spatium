// The sweep: every registered pair on many more inputs than the test takes.
//
//     ninja -C build symmetry_sweep && ./build/tests/symmetry_sweep 300000 [name-fragment]
//
// A rare failure (one ray in twenty thousand was off by 4e-4) is found by volume, and this is the
// loop a scheduled run repeats. Prints one line per pair and exits 1 if any failed.

#include "symmetry/pair.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <string>

int main(int argc, char** argv) {
    const std::size_t n = argc > 1 ? std::strtoull(argv[1], nullptr, 10) : 100000;
    const std::string only = argc > 2 ? argv[2] : "";
    int failed = 0;
    auto all = symmetry::registry();
    std::sort(all.begin(), all.end(), [](const auto& a, const auto& b) { return a.name < b.name; });
    for (const auto& registered : all) {
        if (!only.empty() && registered.name.find(only) == std::string::npos) continue;
        const auto r = registered.build(n).run();
        std::printf("%-44s n=%zu failures=%zu first=%zu worst/bound=%.3g\n", registered.name.c_str(), r.samples,
                    r.failures, r.first_failure, r.worst_ratio);
        std::fflush(stdout);
        if (r.failures) ++failed;
    }
    return failed ? 1 : 0;
}
