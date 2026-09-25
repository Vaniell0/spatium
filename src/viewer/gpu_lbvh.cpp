#include <spatium/viewer/gpu_lbvh.hpp>

#include <spatium/render/gpu_lbvh_glsl.hpp>
#include <spatium/render/gpu_types.hpp>
#include <spatium/render/lbvh.hpp>

#include <chrono>
#include <string>

namespace spatium::viewer::compute {

namespace {

namespace gpu = spatium::render::gpu;

// One dispatch's writes visible to the next dispatch's reads, and a fill's
// to a shader's.
void barrier(VkCommandBuffer cmd, VkPipelineStageFlags from, VkAccessFlags src) {
    VkMemoryBarrier mb{VK_STRUCTURE_TYPE_MEMORY_BARRIER};
    mb.srcAccessMask = src;
    mb.dstAccessMask = VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT;
    vkCmdPipelineBarrier(cmd, from, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &mb, 0, nullptr, 0,
                         nullptr);
}

std::uint32_t groups(std::uint32_t n) { return (n + 255) / 256; }

struct Push { std::uint32_t info[4]; };

}  // namespace

DeviceLbvh::DeviceLbvh(Context& ctx) : ctx_(ctx) {
    const std::string common = gpu::kLbvhCommonGlsl;
    boxes_ = std::make_unique<Kernel>(ctx_, (common + gpu::kBoxesGlsl).c_str(), "lbvh_boxes.comp", 4,
                                      sizeof(Push));
    keys_ = std::make_unique<Kernel>(ctx_, (common + gpu::kKeysGlsl).c_str(), "lbvh_keys.comp", 3,
                                     sizeof(Push));
    tree_ = std::make_unique<Kernel>(ctx_, (common + gpu::kTreeGlsl).c_str(), "lbvh_tree.comp", 4,
                                     sizeof(Push));
    bounds_ = std::make_unique<Kernel>(ctx_, (common + gpu::kBoundsGlsl).c_str(), "lbvh_bounds.comp",
                                       3, sizeof(Push));
    bounds_buf_ = std::make_unique<Buffer>(ctx_, 8 * sizeof(std::uint32_t));
}

DeviceLbvh::~DeviceLbvh() = default;

Buffer& DeviceLbvh::nodes() { return *node_buf_; }

void DeviceLbvh::reserve_(std::uint32_t count) {
    if (count <= capacity_ && node_buf_) return;
    capacity_ = std::max<std::uint32_t>(count, 1);
    box_buf_ = std::make_unique<Buffer>(ctx_, std::size_t{capacity_} * 2 * 16);
    key_buf_ = std::make_unique<Buffer>(ctx_, std::size_t{capacity_} * 8);
    node_buf_ = std::make_unique<Buffer>(ctx_, (2 * std::size_t{capacity_} - 1) * sizeof(gpu::LNode));
    parent_buf_ = std::make_unique<Buffer>(ctx_, (2 * std::size_t{capacity_} - 1) * 4);
    arrived_buf_ = std::make_unique<Buffer>(ctx_, std::size_t{capacity_} * 4);
    bound_instances_ = bound_quadrics_ = nullptr;   // descriptor sets name the old buffers
}

DeviceLbvh::Timing DeviceLbvh::build(Buffer& instances, Buffer& quadrics, std::uint32_t count) {
    Timing t;
    reserve_(count);
    if (bound_instances_ != &instances || bound_quadrics_ != &quadrics) {
        Buffer* b0[] = {&instances, &quadrics, box_buf_.get(), bounds_buf_.get()};
        boxes_->bind(b0);
        Buffer* b1[] = {box_buf_.get(), bounds_buf_.get(), key_buf_.get()};
        keys_->bind(b1);
        Buffer* b2[] = {key_buf_.get(), box_buf_.get(), node_buf_.get(), parent_buf_.get()};
        tree_->bind(b2);
        Buffer* b3[] = {node_buf_.get(), parent_buf_.get(), arrived_buf_.get()};
        bounds_->bind(b3);
        bound_instances_ = &instances;
        bound_quadrics_ = &quadrics;
    }

    // Boxes and keys. The centre bounds start at (max, min) in the ordered
    // encoding: all ones for the minima, zero for the maxima.
    Push p{{count, 0, 0, 0}};
    t.boxes_keys_ms = ctx_.run([&](VkCommandBuffer cmd) {
        vkCmdFillBuffer(cmd, bounds_buf_->handle(), 0, 12, 0xffffffffu);
        vkCmdFillBuffer(cmd, bounds_buf_->handle(), 12, 12, 0u);
        barrier(cmd, VK_PIPELINE_STAGE_TRANSFER_BIT, VK_ACCESS_TRANSFER_WRITE_BIT);
        boxes_->dispatch(cmd, &p, groups(count), 1);
        barrier(cmd, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_ACCESS_SHADER_WRITE_BIT);
        keys_->dispatch(cmd, &p, groups(count), 1);
    });

    // The sort, in place in the shared buffer. Dead instances carry an
    // all-ones code and land at the end, so the live count is where they
    // start.
    const auto t0 = std::chrono::steady_clock::now();
    auto* keys = static_cast<std::uint64_t*>(key_buf_->data());
    render::gpu::detail::radix_sort(keys, count, 8);
    std::uint32_t n = count;
    while (n > 0 && (keys[n - 1] >> 32) == 0xffffffffu) --n;
    leaves_ = n;
    t.read_sort_ms =
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    if (n == 0) return t;

    // Leaves and nodes, then bounds. The root has no parent, and every
    // counter starts at zero.
    Push q{{n, 0, 0, 0}};
    t.tree_ms = ctx_.run([&](VkCommandBuffer cmd) {
        vkCmdFillBuffer(cmd, parent_buf_->handle(), 0, (2 * VkDeviceSize{n} - 1) * 4, 0xffffffffu);
        vkCmdFillBuffer(cmd, arrived_buf_->handle(), 0, VkDeviceSize{n} * 4, 0u);
        barrier(cmd, VK_PIPELINE_STAGE_TRANSFER_BIT, VK_ACCESS_TRANSFER_WRITE_BIT);
        tree_->dispatch(cmd, &q, groups(n), 1);
        barrier(cmd, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_ACCESS_SHADER_WRITE_BIT);
        bounds_->dispatch(cmd, &q, groups(n), 1);
    });
    return t;
}

}  // namespace spatium::viewer::compute
