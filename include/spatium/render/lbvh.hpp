#pragma once
// A tree over instances that move every frame, built in milliseconds
// rather than rebuilt with care: a linear BVH (Karras 2012).
//
// `spatial::BVH` builds by binned SAH, which is the right tree for geometry
// that holds still -- and seconds for two million instances, which is the
// wrong cost for a cloud that moves every frame. An LBVH gives up some tree
// quality for a build that is a sort and two parallel passes: every
// instance gets a Morton code from where its box sits, the codes are
// sorted, and the hierarchy falls out of where adjacent codes first differ.
// Each internal node is decided independently of the others, which is what
// makes the build parallel, on the host's cores now and on a device later
// by the same steps.
//
// Written for the host first, and that is not a stopgap. On a device with
// memory shared with the host -- this one -- the instances a kernel moved
// are readable here without a copy, and the nodes written here are
// readable by the trace kernel without one either; and on a machine with no
// usable GPU this is the tree the CPU path traces.
//
// **The layout differs from `gpu::Node`.** An LBVH's children are wherever
// the split puts them, not "left is next", so an internal node names both.
// Leaves hold one instance each.
//
// Not part of the `spatium.render` module, for the reason
// `render/cooked_scene.hpp` is not.
#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/render/gpu_types.hpp>
#  include <algorithm>
#  include <array>
#  include <atomic>
#  include <bit>
#  include <chrono>
#  include <cmath>
#  include <cstdint>
#  include <limits>
#  include <thread>
#  include <vector>
#endif

namespace spatium::render::gpu {

struct Lbvh {
    std::vector<LNode> nodes;           // root first; empty when nothing is visible
    std::size_t leaves = 0;             // instances in the tree, after culling
    std::size_t culled = 0;             // instances left out: scaled to nothing
    // Where the build's time went, in ms: boxes, keys, sort, hierarchy, bounds.
    double ms[5] = {0, 0, 0, 0, 0};
};

namespace detail {

// Run `fn(begin, end)` over [0, n) in one contiguous chunk per thread. The
// chunks are contiguous, not interleaved, because the radix sort below
// needs each thread's share of the input to be a range.
template<typename Fn>
void parallel_chunks(std::size_t n, Fn&& fn) {
    const std::size_t threads =
        std::max<std::size_t>(1, std::min<std::size_t>(std::thread::hardware_concurrency(),
                                                        (n + 16383) / 16384));
    if (threads == 1) { fn(std::size_t{0}, n, std::size_t{0}); return; }
    std::vector<std::jthread> pool;
    pool.reserve(threads);
    for (std::size_t t = 0; t < threads; ++t)
        pool.emplace_back([&, t] { fn(n * t / threads, n * (t + 1) / threads, t); });
}

inline std::size_t chunk_threads(std::size_t n) {
    return std::max<std::size_t>(1, std::min<std::size_t>(std::thread::hardware_concurrency(),
                                                          (n + 16383) / 16384));
}

// Spread the low 10 bits of `v` to every third bit.
inline std::uint32_t spread10(std::uint32_t v) {
    v &= 0x3ffu;
    v = (v | (v << 16)) & 0x030000FFu;
    v = (v | (v << 8)) & 0x0300F00Fu;
    v = (v | (v << 4)) & 0x030C30C3u;
    v = (v | (v << 2)) & 0x09249249u;
    return v;
}

// An instance's box in world space: the quadric's clip box through the
// uniform scale and the rotation -- whose absolute value takes a local
// half-extent to a world one -- then the translation. Widened by a few
// ulps, for the reason `convert_nodes` rounds outward.
inline void instance_box(const Instance& g, const Quadric& q, float lo[3], float hi[3]) {
    const float s = std::abs(g.scale_quadric[0]);
    const float* rows[3] = {g.r0, g.r1, g.r2};
    float c[3], h[3];
    for (int k = 0; k < 3; ++k) {
        c[k] = 0.5f * (q.lo[k] + q.hi[k]) * g.scale_quadric[0];
        h[k] = 0.5f * (q.hi[k] - q.lo[k]) * s;
    }
    for (int r = 0; r < 3; ++r) {
        float wc = rows[r][3], wh = 0.0f;
        for (int k = 0; k < 3; ++k) {
            wc += rows[r][k] * c[k];
            wh += std::abs(rows[r][k]) * h[k];
        }
        const float pad = wh * 1e-6f + std::abs(wc) * 1e-7f;
        lo[r] = wc - wh - pad;
        hi[r] = wc + wh + pad;
    }
}

// LSD radix sort of 64-bit keys, eight bits a pass, `passes` passes, with
// every thread counting and scattering its own contiguous share -- stable,
// which is what lets the low passes survive the high ones.
inline void radix_sort(std::vector<std::uint64_t>& keys, int passes) {
    const std::size_t n = keys.size();
    std::vector<std::uint64_t> tmp(n);
    const std::size_t T = chunk_threads(n);
    std::vector<std::array<std::size_t, 256>> count(T);
    for (int pass = 0; pass < passes; ++pass) {
        const int shift = pass * 8;
        for (auto& c : count) c.fill(0);
        parallel_chunks(n, [&](std::size_t b, std::size_t e, std::size_t t) {
            for (std::size_t i = b; i < e; ++i) ++count[t][(keys[i] >> shift) & 0xffu];
        });
        // Offsets: bucket-major, thread-minor, so each thread's share of a
        // bucket lands after the shares of the threads before it.
        std::size_t sum = 0;
        for (std::size_t d = 0; d < 256; ++d)
            for (std::size_t t = 0; t < T; ++t) {
                const std::size_t c = count[t][d];
                count[t][d] = sum;
                sum += c;
            }
        parallel_chunks(n, [&](std::size_t b, std::size_t e, std::size_t t) {
            auto& off = count[t];
            for (std::size_t i = b; i < e; ++i) tmp[off[(keys[i] >> shift) & 0xffu]++] = keys[i];
        });
        keys.swap(tmp);
    }
}

// The same sort over keys held elsewhere -- a mapped device buffer -- read
// once into host memory, where the passes run, and written back once.
inline void radix_sort(std::uint64_t* data, std::size_t n, int passes) {
    std::vector<std::uint64_t> keys(data, data + n);
    radix_sort(keys, passes);
    std::copy(keys.begin(), keys.end(), data);
}

}  // namespace detail

// The tree over `instances`, whose shapes are `quadrics`. Instances scaled
// to zero are left out: a ray test refuses them, so a box for them is only
// a node every ray has to walk past.
inline Lbvh build_lbvh(const std::vector<Instance>& instances, const std::vector<Quadric>& quadrics) {
    Lbvh out;
    auto clock = std::chrono::steady_clock::now();
    auto lap = [&](int phase) {
        const auto now = std::chrono::steady_clock::now();
        out.ms[phase] = std::chrono::duration<double, std::milli>(now - clock).count();
        clock = now;
    };
    const std::size_t total = instances.size();
    std::vector<float> box(total * 6);
    std::vector<char> live(total, 0);

    // Boxes, and the bounds of their centres, one partial bound per thread.
    const std::size_t T = detail::chunk_threads(total);
    std::vector<std::array<float, 6>> part(T, {std::numeric_limits<float>::max(),
                                               std::numeric_limits<float>::max(),
                                               std::numeric_limits<float>::max(),
                                               -std::numeric_limits<float>::max(),
                                               -std::numeric_limits<float>::max(),
                                               -std::numeric_limits<float>::max()});
    detail::parallel_chunks(total, [&](std::size_t b, std::size_t e, std::size_t t) {
        for (std::size_t i = b; i < e; ++i) {
            const auto& g = instances[i];
            if (g.scale_quadric[0] == 0.0f) continue;
            const auto qi = detail::bits_of(g.scale_quadric[1]);
            float* bx = &box[i * 6];
            detail::instance_box(g, quadrics[qi], bx, bx + 3);
            live[i] = 1;
            for (int k = 0; k < 3; ++k) {
                const float c = 0.5f * (bx[k] + bx[3 + k]);
                part[t][k] = std::min(part[t][k], c);
                part[t][3 + k] = std::max(part[t][3 + k], c);
            }
        }
    });
    float lo[3], ext[3];
    for (int k = 0; k < 3; ++k) {
        float a = std::numeric_limits<float>::max(), z = -std::numeric_limits<float>::max();
        for (const auto& p : part) { a = std::min(a, p[k]); z = std::max(z, p[3 + k]); }
        lo[k] = a;
        ext[k] = z > a ? z - a : 1.0f;
    }

    lap(0);
    // Keys: the Morton code of the centre in the high bits, the instance
    // index in the low 32, so no two keys are equal -- the split search
    // needs a strict order, and an index is the cheapest tiebreak.
    std::vector<std::uint64_t> keys;
    keys.reserve(total);
    for (std::size_t i = 0; i < total; ++i) {
        if (!live[i]) { ++out.culled; continue; }
        std::uint32_t code = 0;
        for (int k = 0; k < 3; ++k) {
            const float c = 0.5f * (box[i * 6 + k] + box[i * 6 + 3 + k]);
            const float u = std::clamp((c - lo[k]) / ext[k], 0.0f, 1.0f);
            code |= detail::spread10(static_cast<std::uint32_t>(u * 1023.0f)) << (2 - k);
        }
        keys.push_back((static_cast<std::uint64_t>(code) << 32) | static_cast<std::uint32_t>(i));
    }
    const std::size_t n = keys.size();
    out.leaves = n;
    lap(1);
    if (n == 0) return out;
    detail::radix_sort(keys, 8);
    lap(2);   // all 64 bits: simple before it is fast

    auto leaf = [&](std::size_t i) { return static_cast<std::uint32_t>(n - 1 + i); };
    out.nodes.resize(2 * n - 1);
    std::vector<std::uint32_t> parent(2 * n - 1, 0xffffffffu);

    // Leaves.
    detail::parallel_chunks(n, [&](std::size_t b, std::size_t e, std::size_t) {
        for (std::size_t i = b; i < e; ++i) {
            const auto inst = static_cast<std::uint32_t>(keys[i] & 0xffffffffu);
            auto& L = out.nodes[leaf(i)];
            for (int k = 0; k < 3; ++k) {
                L.lo[k] = box[inst * 6 + k];
                L.hi[k] = box[inst * 6 + 3 + k];
            }
            L.left = inst | kLeafBit;
            L.right = 0;
        }
    });

    if (n == 1) return out;   // 2n - 1 = 1: the root is the one leaf

    // Internal nodes, each on its own: Karras 2012, section 4.
    auto delta = [&](std::int64_t i, std::int64_t j) -> int {
        if (j < 0 || j >= static_cast<std::int64_t>(n)) return -1;
        return std::countl_zero(keys[static_cast<std::size_t>(i)] ^ keys[static_cast<std::size_t>(j)]);
    };
    detail::parallel_chunks(n - 1, [&](std::size_t b, std::size_t e, std::size_t) {
        for (std::size_t ui = b; ui < e; ++ui) {
            const auto i = static_cast<std::int64_t>(ui);
            const int d = delta(i, i + 1) - delta(i, i - 1) > 0 ? 1 : -1;
            const int dmin = delta(i, i - d);
            std::int64_t lmax = 2;
            while (delta(i, i + lmax * d) > dmin) lmax *= 2;
            std::int64_t l = 0;
            for (std::int64_t t = lmax / 2; t >= 1; t /= 2)
                if (delta(i, i + (l + t) * d) > dmin) l += t;
            const std::int64_t j = i + l * d;
            const int dnode = delta(i, j);
            std::int64_t s = 0;
            for (std::int64_t div = 2;; div *= 2) {
                const std::int64_t t = (l + div - 1) / div;
                if (delta(i, i + (s + t) * d) > dnode) s += t;
                if (t <= 1) break;
            }
            const std::int64_t g = i + s * d + std::min(d, 0);
            const auto lo_i = std::min(i, j), hi_i = std::max(i, j);
            const auto left = lo_i == g ? leaf(static_cast<std::size_t>(g)) : static_cast<std::uint32_t>(g);
            const auto right = hi_i == g + 1 ? leaf(static_cast<std::size_t>(g + 1))
                                             : static_cast<std::uint32_t>(g + 1);
            out.nodes[ui].left = left;
            out.nodes[ui].right = right;
            parent[left] = static_cast<std::uint32_t>(ui);
            parent[right] = static_cast<std::uint32_t>(ui);
        }
    });

    lap(3);
    // Bounds, bottom up: a leaf climbs, and at each parent the first child
    // to arrive stops while the second, which knows both are written,
    // takes the union and climbs on.
    std::vector<std::atomic<std::uint32_t>> arrived(n - 1);
    for (auto& a : arrived) a.store(0, std::memory_order_relaxed);
    detail::parallel_chunks(n, [&](std::size_t b, std::size_t e, std::size_t) {
        for (std::size_t i = b; i < e; ++i) {
            std::uint32_t at = parent[leaf(i)];
            while (at != 0xffffffffu) {
                if (arrived[at].fetch_add(1, std::memory_order_acq_rel) == 0) break;
                auto& N = out.nodes[at];
                const auto& A = out.nodes[N.left];
                const auto& B = out.nodes[N.right];
                for (int k = 0; k < 3; ++k) {
                    N.lo[k] = std::min(A.lo[k], B.lo[k]);
                    N.hi[k] = std::max(A.hi[k], B.hi[k]);
                }
                at = parent[at];
            }
        }
    });
    lap(4);
    return out;
}

}  // namespace spatium::render::gpu
