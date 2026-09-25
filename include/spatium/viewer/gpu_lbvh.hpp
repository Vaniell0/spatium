#pragma once
// A linear BVH built on the device over instances that already live there:
// the kernels of `render/gpu_lbvh_glsl.hpp`, with the sort between them on
// the host.
//
// Two submissions a build, because the sort sits between them: boxes and
// keys, then -- after the host has sorted the keys in place in shared
// memory -- leaves, internal nodes and bounds. The node buffer it leaves
// behind is the `LNode` array the trace kernel reads, so a frame is: move
// the instances, build this, trace.

#include <spatium/viewer/compute.hpp>

#include <cstdint>
#include <memory>

namespace spatium::viewer::compute {

class DeviceLbvh {
public:
    explicit DeviceLbvh(Context& ctx);
    ~DeviceLbvh();
    DeviceLbvh(const DeviceLbvh&) = delete;
    DeviceLbvh& operator=(const DeviceLbvh&) = delete;

    struct Timing {
        double boxes_keys_ms = 0;   // device
        double read_sort_ms = 0;    // host, including the read of the keys
        double tree_ms = 0;         // device: leaves, nodes, bounds
    };

    // Over `count` instances in `instances` (gpu::Instance[]), whose shapes
    // are in `quadrics` (gpu::Quadric[]). Instances scaled to zero are left
    // out, as build_lbvh() leaves them out.
    Timing build(Buffer& instances, Buffer& quadrics, std::uint32_t count);

    // The same, leaving out instances a projection pass draws instead:
    // `cam` and `fwd` as render/gpu_splat_glsl.hpp's splat_small() reads
    // them -- position and pixels-per-unit, forward and the threshold.
    Timing build(Buffer& instances, Buffer& quadrics, std::uint32_t count, const float cam[4],
                 const float fwd[4]);

    // Each instance's box, lo (w = live) and hi, as the build left them:
    // what a projection pass reads.
    Buffer& boxes() { return *box_buf_; }

    // The tree: `2 * leaves() - 1` nodes, root first; none when no instance
    // is visible.
    Buffer& nodes();
    std::uint32_t leaves() const { return leaves_; }

private:
    void reserve_(std::uint32_t count);

    Context& ctx_;
    std::unique_ptr<Kernel> boxes_, keys_, tree_, bounds_;
    std::unique_ptr<Buffer> box_buf_, bounds_buf_, key_buf_, node_buf_, parent_buf_, arrived_buf_;
    std::uint32_t capacity_ = 0, leaves_ = 0;
    Buffer* bound_instances_ = nullptr;
    Buffer* bound_quadrics_ = nullptr;
};

}  // namespace spatium::viewer::compute
