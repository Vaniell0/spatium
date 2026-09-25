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

struct GLFWwindow;

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <span>
#include <string>
#include <vector>

namespace spatium::viewer::compute {

class Context {
public:
    // Picks a discrete GPU, then an integrated one, and refuses a CPU
    // implementation (llvmpipe) unless SPATIUM_VK_ALLOW_CPU is set -- a
    // software device would pass every check here and measure nothing.
    explicit Context(const char* app_name);

    // The same, able to present into `window`: the instance gets the
    // surface extensions, the device the swapchain one, and the queue is
    // one that can compute, draw (ImGui) and present.
    Context(const char* app_name, GLFWwindow* window);
    ~Context();
    Context(const Context&) = delete;
    Context& operator=(const Context&) = delete;

    VkDevice device() const { return device_; }
    VkPhysicalDevice physical() const { return physical_; }
    const std::string& device_name() const { return name_; }
    VkInstance instance() const { return instance_; }
    VkQueue queue() const { return queue_; }
    std::uint32_t queue_family() const { return family_; }
    VkSurfaceKHR surface() const { return surface_; }

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
    friend class Presenter;
    void init_(const char* app_name, GLFWwindow* window);
    std::uint32_t memory_type_(std::uint32_t bits) const;
    void begin_();
    double end_and_wait_();

    VkInstance instance_{};
    VkSurfaceKHR surface_{};
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

// A window's swapchain, fed by a compute kernel's pixel buffer.
//
// Each frame: the kernel writes RGBA8 into a storage buffer, the buffer is
// copied into the acquired swapchain image, and ImGui draws on top in a
// render pass that loads rather than clears -- so the traced picture is
// the background and the controls are the only thing rasterized. One frame
// in flight: the trace dominates the frame, and a second in flight would
// only queue a second trace behind the first.
class Presenter {
public:
    Presenter(Context& ctx, GLFWwindow* window);
    ~Presenter();
    Presenter(const Presenter&) = delete;
    Presenter& operator=(const Presenter&) = delete;

    std::uint32_t width() const { return extent_.width; }
    std::uint32_t height() const { return extent_.height; }
    // True when the swapchain stores B,G,R,A, so the kernel can write
    // pixels in the order the image expects and the copy stays a copy.
    bool bgra() const { return bgra_; }

    // Records `trace` (which must fill `pixels`, width*height uint32),
    // copies, draws the ImGui frame the caller has already built, presents.
    // Returns false when the swapchain had to be rebuilt -- the caller
    // then resizes its pixel buffer to the new extent and tries again.
    // `gpu_ms` receives the time of the trace alone.
    template<typename F>
    bool frame(Buffer& pixels, F&& trace, double& gpu_ms) {
        if (!begin_frame_()) return false;
        trace(cmd_);
        return end_frame_(pixels, gpu_ms);
    }

private:
    bool begin_frame_();
    bool end_frame_(Buffer& pixels, double& gpu_ms);
    void create_swapchain_();
    void create_framebuffers_();
    void destroy_swapchain_();

    Context& ctx_;
    GLFWwindow* window_ = nullptr;
    VkSwapchainKHR swapchain_{};
    VkFormat format_{};
    bool bgra_ = false;
    VkExtent2D extent_{};
    std::vector<VkImage> images_;
    std::vector<VkImageView> views_;
    std::vector<VkFramebuffer> framebuffers_;
    std::vector<VkSemaphore> finished_;   // one per image, signalled by the submit
    VkRenderPass render_pass_{};
    VkSemaphore acquired_{};
    VkFence fence_{};
    VkCommandBuffer cmd_{};
    VkQueryPool timestamps_{};
    VkDescriptorPool imgui_pool_{};
    std::uint32_t image_ = 0;
};

}  // namespace spatium::viewer::compute
