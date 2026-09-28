// The connectivity matrix, held (docs/connectivity.md,
// scripts/gen_connectivity.py). Every cell the generator found compiling is
// instantiated here -- so a body that stops compiling stops the build -- and
// run, and must reach the level recorded for it: L1 and L2 as the cell
// reports them, L3 by the same comparison the generator makes, against the
// reference scalar's cell. A cell that got better is not a failure; the
// generator is rerun to record it.

#include <catch2/catch_test_macros.hpp>

#include "connectivity/probes.hpp"
#include "connectivity/expected.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <limits>
#include <map>
#include <string>
#include <tuple>
#include <vector>

namespace {

struct Ran {
    int level;
    std::vector<double> sig;
    double eps;
};

using Key = std::tuple<std::string, std::string, std::string>;

std::map<Key, Ran>& runs() {
    static std::map<Key, Ran> r;
    return r;
}

std::map<Key, int>& expected() {
    static std::map<Key, int> e;
    return e;
}

template<class T, class Space, class Probe>
void run_one(const char* scalar, const char* space, const char* probe, int level) {
    std::vector<double> sig;
    const int got = connectivity::run_cell<T, Space, Probe>(sig);
    runs()[{scalar, space, probe}] = {got, sig, connectivity::base_epsilon<T>()};
    expected()[{scalar, space, probe}] = level;
}

void run_all() {
    if (!runs().empty()) return;
#define CONNECTIVITY_RUN(T, NAME, SPACE, PROBE, LEVEL) \
    run_one<T, connectivity::SPACE, connectivity::PROBE>(NAME, #SPACE, #PROBE, LEVEL);
    CONNECTIVITY_CELLS(CONNECTIVITY_RUN)
#undef CONNECTIVITY_RUN
}

}  // namespace

TEST_CASE("Every cell of the connectivity matrix reaches its recorded level", "[connectivity]") {
    run_all();
    REQUIRE_FALSE(runs().empty());
    for (const auto& [key, want] : expected()) {
        const auto& [scalar, space, probe] = key;
        const auto& got = runs().at(key);
        int level = got.level;
        // L3: the same numbers as the reference scalar's cell.
        const std::string ref = scalar == "double" ? "Real50" : "double";
        const auto it = runs().find({ref, space, probe});
        if (level == 2 && it != runs().end() && it->second.level >= 1 && it->second.sig.size() == got.sig.size() &&
            !got.sig.empty()) {
            const double eps = std::max(got.eps, it->second.eps);
            bool agree = true;
            for (std::size_t i = 0; i < got.sig.size(); ++i)
                agree = agree && std::abs(got.sig[i] - it->second.sig[i]) <=
                                     std::sqrt(eps) * 64 * (std::abs(it->second.sig[i]) + 1);
            if (agree) level = 3;
        }
        // Without the reference scalar in this build (double's is Real50,
        // which needs Boost), L3 cannot be confirmed here; L2 can.
        const int reachable = it == runs().end() ? std::min(want, 2) : want;
        INFO(std::format("{} / {} / {}: recorded L{}, now L{}", scalar, space, probe, want, level));
        CHECK(level >= reachable);
    }
}
