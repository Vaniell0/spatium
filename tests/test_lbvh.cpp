#include <catch2/catch_test_macros.hpp>
#include <spatium/algebra/groups/so3.hpp>
#include <spatium/geometry/ray_surface.hpp>
#include <spatium/render/lbvh.hpp>

#include <cstring>
#include <random>
#include <vector>

using namespace spatium;
namespace g = spatium::render::gpu;

namespace {

g::Quadric quadric_of(const geometry::BoundedQuadric<double>& bq) {
    g::Quadric q{};
    for (std::size_t c = 0; c < 4; ++c)
        for (std::size_t r = 0; r < 4; ++r) q.q[c * 4 + r] = static_cast<float>(bq.surface.Q(r, c));
    for (std::size_t k = 0; k < 3; ++k) {
        q.lo[k] = static_cast<float>(bq.clip.min_corner[k]);
        q.hi[k] = static_cast<float>(bq.clip.max_corner[k]);
    }
    q.lo[3] = bq.closed ? 1.0f : 0.0f;
    return q;
}

std::uint32_t bits(float f) { std::uint32_t u; std::memcpy(&u, &f, 4); return u; }

}  // namespace

// The tree is judged by the one thing a traversal promises: the nearest hit
// it finds is the nearest hit there is. Brute force over every instance is
// the reference, with the same ray test, so agreement is to the bit -- a
// culled branch that held the true nearest hit would show as a different
// index, not as a small error.
TEST_CASE("An LBVH finds the same nearest instance as brute force", "[lbvh]") {
    std::vector<g::Quadric> quadrics{
        quadric_of(geometry::BoundedQuadric<double>::sphere(1.0)),
        quadric_of([] { auto c = geometry::BoundedQuadric<double>::cylinder_z(0.3, 0.0, 2.0);
                        c.closed = true; return c; }())};
    std::mt19937 rng(7);
    std::uniform_real_distribution<float> pos(-10.0f, 10.0f), sc(0.05f, 0.4f), ang(-3.0f, 3.0f);

    for (std::size_t count : {1u, 2u, 3u, 17u, 5000u}) {
        std::vector<g::Instance> inst(count);
        for (std::size_t i = 0; i < count; ++i) {
            auto R = algebra::SO3<double>{}.exp(Vec<double, 3>{ang(rng), ang(rng), ang(rng)});
            auto& x = inst[i];
            float* rows[3] = {x.r0, x.r1, x.r2};
            for (int r = 0; r < 3; ++r) {
                for (int c = 0; c < 3; ++c) rows[r][c] = static_cast<float>(R(r, c));
                rows[r][3] = pos(rng);
            }
            // Every seventh one scaled to nothing, which the tree must leave out.
            x.scale_quadric[0] = (i % 7 == 6) ? 0.0f : sc(rng);
            const std::uint32_t qi = static_cast<std::uint32_t>(i % 2);
            std::memcpy(&x.scale_quadric[1], &qi, 4);
        }
        const auto tree = g::build_lbvh(inst, quadrics);
        CHECK(tree.leaves + tree.culled == count);
        if (tree.leaves > 0) CHECK(tree.nodes.size() == 2 * tree.leaves - 1);

        std::size_t hits = 0;
        for (int k = 0; k < 400; ++k) {
            g::V3 o{pos(rng) * 1.5f, pos(rng) * 1.5f, pos(rng) * 1.5f};
            g::V3 tgt{pos(rng) * 0.5f, pos(rng) * 0.5f, pos(rng) * 0.5f};
            g::V3 d = tgt - o;
            d = d * (1.0f / std::sqrt(g::dot(d, d)));
            const g::Ray32 ray{o, d};

            g::Hit brute;
            for (std::uint32_t i = 0; i < count; ++i)
                if (g::hit_instance(ray, inst[i], quadrics[i % 2], brute)) {
                    brute.kind = g::HitKind::Instance;
                    brute.index = i;
                }
            g::Hit fast;
            g::trace_lbvh(tree, inst, quadrics, ray, fast);
            INFO("count " << count << " ray " << k);
            REQUIRE(fast.kind == brute.kind);
            if (brute.kind == g::HitKind::None) continue;
            ++hits;
            CHECK(fast.index == brute.index);
            CHECK(bits(fast.t) == bits(brute.t));
        }
        if (count == 5000) CHECK(hits > 50);   // the rays do meet things
    }
}
