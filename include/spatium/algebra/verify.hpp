#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/concepts.hpp>
#  include <spatium/core/verify.hpp>
#  include <span>
#endif

SPATIUM_EXPORT namespace spatium::inline algebra {

// Verify for matrix groups (elements support operator-): uses Frobenius norm for comparison

template<Group G>
    requires requires(typename G::ElementType a, typename G::ElementType b) {
        { a - b };
    }
VerifyResult verify_matrix_group(const G& g,
                                 std::span<const typename G::ElementType> samples,
                                 double tolerance = 1e-8) {
    auto e = g.identity();

    for (const auto& a : samples) {
        // Inverse check: a * inv(a) should be close to identity
        auto inv = g.inverse(a);
        auto prod = g.compose(a, inv);
        auto diff = prod - e;
        // Check all elements near zero
        for (std::size_t i = 0; i < diff.data.size(); ++i) {
            using std::abs;
            if (static_cast<double>(abs(diff.data[i])) > tolerance)
                return VerifyResult::fail("a * inv(a) != identity");
        }
    }

    // Associativity
    for (std::size_t i = 0; i + 2 < samples.size(); ++i) {
        auto ab = g.compose(samples[i], samples[i + 1]);
        auto bc = g.compose(samples[i + 1], samples[i + 2]);
        auto ab_c = g.compose(ab, samples[i + 2]);
        auto a_bc = g.compose(samples[i], bc);
        auto diff = ab_c - a_bc;
        for (std::size_t j = 0; j < diff.data.size(); ++j) {
            using std::abs;
            if (static_cast<double>(abs(diff.data[j])) > tolerance)
                return VerifyResult::fail("associativity violated");
        }
    }

    return VerifyResult::ok();
}

} // namespace spatium::algebra
