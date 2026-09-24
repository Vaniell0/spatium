#pragma once
// Headless Vulkan compute: a device, host-visible storage buffers, and a
// kernel compiled from GLSL at run time.
//
// Deliberately not `viewer::App`. App is a rasterizer -- render pass,
// depth buffer, mesh and point pipelines -- and a ray tracer needs none of
// it: an image written by a compute shader and, later, copied to a
// swapchain. Keeping the two apart costs some repeated setup and keeps the
// compute path free of a window, which is what lets it be checked on a
// machine with no display at all.
//
// Buffers are host-visible and coherent. On the integrated GPUs this is
// first written for, device memory *is* host memory, so a staging copy
// would be a copy of the same bytes into the same RAM.

#include <vulkan/vulkan.h>

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <span>
#include <string>

namespace spatium::viewer::compute {

class Context {
public:
    // Picks a discrete GPU, then an integrated one, and refuses a CPU
    // implementation (llvmpipe) unless SPATIUM_VK_ALLOW_CPU is set -- a
    // software device would pass every check here and measure nothing.
    explicit Context(const char* app_name);
    ~Context();
    Context(const Context&) = delete;
    Context& operator=(const Context&) = delete;

    VkDevice device() const { return device_; }
    VkPhysicalDevice physical() const { return physical_; }
    const std::string& device_name() const { return name_; }

    // Records through `record`, submits, and waits. Returns the GPU time
    // between the two timestamps written around the recording, in ms.
    template<typename F>
    double run(F&& record) {
        begin_();
        record(cmd_);
        return end_and_wait_();
    }

private:
    friend class Buffer;
    friend class Kernel;
    std::uint32_t memory_type_(std::uint32_t bits) const;
    void begin_();
    double end_and_wait_();

    VkInstance instance_{};
    VkPhysicalDevice physical_{};
    VkDevice device_{};
    VkQueue queue_{};
    std::uint32_t family_ = 0;
    VkCommandPool pool_{};
    VkCommandBuffer cmd_{};
    VkFence fence_{};
    VkQueryPool timestamps_{};
    double ns_per_tick_ = 1.0;
    std::string name_;
};

class Buffer {
public:
    Buffer(Context& ctx, std::size_t bytes);
    ~Buffer();
    Buffer(Buffer&& o) noexcept;
    Buffer& operator=(Buffer&&) = delete;
    Buffer(const Buffer&) = delete;
    Buffer& operator=(const Buffer&) = delete;

    void* data() { return mapped_; }
    std::size_t size() const { return bytes_; }
    VkBuffer handle() const { return buffer_; }

    template<typename T>
    static Buffer from(Context& ctx, std::span<const T> src) {
        // A zero-sized storage buffer is not allowed; an empty array still
        // needs a binding for the shader to compile against.
        Buffer b(ctx, src.empty() ? 16 : src.size_bytes());
        if (!src.empty()) std::memcpy(b.data(), src.data(), src.size_bytes());
        return b;
    }

private:
    VkDevice device_{};
    VkBuffer buffer_{};
    VkDeviceMemory memory_{};
    void* mapped_ = nullptr;
    std::size_t bytes_ = 0;
};

// A compute pipeline whose bindings are `buffers` storage buffers at
// set 0, bindings 0..n-1, plus `push_bytes` of push constants.
class Kernel {
public:
    Kernel(Context& ctx, const char* glsl, const char* tag, std::uint32_t buffers,
           std::uint32_t push_bytes);
    ~Kernel();
    Kernel(const Kernel&) = delete;
    Kernel& operator=(const Kernel&) = delete;

    void bind(std::span<Buffer* const> buffers);
    void dispatch(VkCommandBuffer cmd, const void* push, std::uint32_t gx, std::uint32_t gy,
                  std::uint32_t gz = 1) const;

private:
    VkDevice device_{};
    VkDescriptorSetLayout set_layout_{};
    VkPipelineLayout layout_{};
    VkPipeline pipeline_{};
    VkDescriptorPool pool_{};
    VkDescriptorSet set_{};
    std::uint32_t buffers_ = 0, push_bytes_ = 0;
};

}  // namespace spatium::viewer::compute
