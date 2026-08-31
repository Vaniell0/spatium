#include <spatium/viewer/app.hpp>

#define GLFW_INCLUDE_VULKAN
#include <GLFW/glfw3.h>

#include <vulkan/vulkan.h>
#include <shaderc/shaderc.hpp>

#if SPATIUM_HAS_IMGUI
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_vulkan.h>
#endif

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <format>
#include <memory>
#include <stdexcept>
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <spatium/vendor/stb_image_write.h>

#define VK_CHECK(expr) do { \
    VkResult _r = (expr); \
    if (_r != VK_SUCCESS) \
        throw std::runtime_error(std::format("Vulkan error {} at {}:{}", (int)_r, __FILE__, __LINE__)); \
} while(0)

namespace spatium::viewer {

// ── Shader sources ─────────────────────────────────────────────

static const char* VERT_SRC = R"(
#version 450
layout(push_constant) uniform PushConstants {
    mat4 mvp;
    vec4 mesh_color;
} pc;
layout(location = 0) in vec3 inPos;
layout(location = 1) in vec3 inNormal;
layout(location = 0) out vec3 fragNormal;
layout(location = 1) out vec3 fragPos;
layout(location = 2) out vec4 fragColor;
void main() {
    gl_Position = pc.mvp * vec4(inPos, 1.0);
    fragNormal = inNormal;
    fragPos = inPos;
    fragColor = pc.mesh_color;
}
)";

static const char* FRAG_SRC = R"(
#version 450
layout(location = 0) in vec3 fragNormal;
layout(location = 1) in vec3 fragPos;
layout(location = 2) in vec4 fragColor;
layout(location = 0) out vec4 outColor;

vec3 hsv2rgb(float h, float s, float v) {
    vec3 p = abs(fract(vec3(h) + vec3(1.0, 2.0/3.0, 1.0/3.0)) * 6.0 - 3.0);
    return v * mix(vec3(1.0), clamp(p - 1.0, 0.0, 1.0), s);
}

void main() {
    vec3 n = normalize(fragNormal);
    if (!gl_FrontFacing) n = -n;
    vec3 lightDir = normalize(vec3(1.0, 1.0, 1.0));
    float diff = max(dot(n, lightDir), 0.2);

    // Per-face color from triangle ID (golden ratio hash)
    uint id = gl_PrimitiveID;
    float hue = fract(float(id) * 0.618034);
    float sat = 0.55 + 0.15 * fract(float(id) * 0.3819);
    float val = 0.75 + 0.15 * fract(float(id) * 0.2473);

    vec3 faceColor = hsv2rgb(hue, sat, val);
    // Tint with mesh color
    vec3 color = mix(faceColor, fragColor.rgb, 0.3) * diff;
    outColor = vec4(color, fragColor.a);
}
)";

static const char* FRAG_EDGE_SRC = R"(
#version 450
layout(location = 0) in vec3 fragNormal;
layout(location = 1) in vec3 fragPos;
layout(location = 0) out vec4 outColor;
void main() {
    outColor = vec4(0.08, 0.08, 0.12, 1.0);
}
)";

static const char* POINT_VERT_SRC = R"(
#version 450
layout(push_constant) uniform PushConstants {
    mat4 mvp;        // offset 0, 64 bytes
    vec4 mesh_color; // offset 64, 16 bytes
    vec4 clip_plane; // offset 80, 16 bytes (xyz = normal, w = offset)
    float point_size;// offset 96, 4 bytes
} pc;
layout(location = 0) in vec3 inPos;
layout(location = 0) out vec4 fragColor;
void main() {
    // Quarter cutaway: discard points in the world +X/+Z quadrant (fixed, camera-independent)
    if (pc.clip_plane.w > -9998.0 && inPos.x > 0.0 && inPos.z > 0.0) {
        gl_Position = vec4(0.0);
        gl_PointSize = 0.0;
        fragColor = vec4(0.0);
        return;
    }
    gl_Position = pc.mvp * vec4(inPos, 1.0);
    gl_PointSize = pc.point_size;
    fragColor = pc.mesh_color;
}
)";

static const char* POINT_FRAG_SRC = R"(
#version 450
layout(location = 0) in vec4 fragColor;
layout(location = 0) out vec4 outColor;
void main() {
    outColor = fragColor;
}
)";

// ── Compile GLSL to SPIR-V ─────────────────────────────────────

static std::vector<uint32_t> compile_glsl(const char* src, shaderc_shader_kind kind) {
    shaderc::Compiler compiler;
    shaderc::CompileOptions options;
    options.SetTargetEnvironment(shaderc_target_env_vulkan, shaderc_env_version_vulkan_1_3);
    auto result = compiler.CompileGlslToSpv(src, kind, "shader", options);
    if (result.GetCompilationStatus() != shaderc_compilation_status_success)
        throw std::runtime_error(std::string("Shader compile error: ") + result.GetErrorMessage());
    return {result.cbegin(), result.cend()};
}

// ── Vulkan implementation ──────────────────────────────────────

struct App::Impl {
    GLFWwindow* window = nullptr;
    VkInstance instance = VK_NULL_HANDLE;
    VkSurfaceKHR surface = VK_NULL_HANDLE;
    VkPhysicalDevice physical_device = VK_NULL_HANDLE;
    VkDevice device = VK_NULL_HANDLE;
    uint32_t graphics_family = 0;
    VkQueue graphics_queue = VK_NULL_HANDLE;
    VkQueue present_queue = VK_NULL_HANDLE;
    VkSwapchainKHR swapchain = VK_NULL_HANDLE;
    VkFormat swapchain_format{};
    VkExtent2D swapchain_extent{};
    std::vector<VkImage> swapchain_images;
    std::vector<VkImageView> swapchain_views;
    VkRenderPass render_pass = VK_NULL_HANDLE;
    std::vector<VkFramebuffer> framebuffers;
    VkPipelineLayout pipeline_layout = VK_NULL_HANDLE;
    VkPipeline pipeline_solid = VK_NULL_HANDLE;
    VkPipeline pipeline_wire = VK_NULL_HANDLE;
    VkCommandPool cmd_pool = VK_NULL_HANDLE;
    std::vector<VkCommandBuffer> cmd_buffers;
    std::vector<VkSemaphore> sem_available;
    std::vector<VkSemaphore> sem_finished;
    std::vector<VkFence> fences_in_flight;
    uint32_t current_frame = 0;
    struct MeshBuffer {
        VkBuffer vertex_buf = VK_NULL_HANDLE;
        VkDeviceMemory vertex_mem = VK_NULL_HANDLE;
        VkBuffer index_buf = VK_NULL_HANDLE;
        VkDeviceMemory index_mem = VK_NULL_HANDLE;
    };
    std::vector<MeshBuffer> mesh_buffers;

    VkPipeline pipeline_points = VK_NULL_HANDLE;
    struct PointBuffer {
        VkBuffer vertex_buf = VK_NULL_HANDLE;
        VkDeviceMemory vertex_mem = VK_NULL_HANDLE;
    };
    std::vector<PointBuffer> point_buffers;

    VkImage depth_image = VK_NULL_HANDLE;
    VkDeviceMemory depth_memory = VK_NULL_HANDLE;
    VkImageView depth_view = VK_NULL_HANDLE;

    VkDescriptorPool imgui_pool = VK_NULL_HANDLE;

    // Mouse state
    bool dragging = false;
    double last_x = 0, last_y = 0;
    bool framebuffer_resized = false;
};

// ── Helper: find memory type ───────────────────────────────────

static uint32_t find_memory_type(VkPhysicalDevice pd, uint32_t filter, VkMemoryPropertyFlags props) {
    VkPhysicalDeviceMemoryProperties mem;
    vkGetPhysicalDeviceMemoryProperties(pd, &mem);
    for (uint32_t i = 0; i < mem.memoryTypeCount; ++i)
        if ((filter & (1u << i)) && (mem.memoryTypes[i].propertyFlags & props) == props)
            return i;
    throw std::runtime_error("No suitable memory type");
}

// ── Helper: create buffer ──────────────────────────────────────

static void create_buffer(VkDevice dev, VkPhysicalDevice pd,
                          VkDeviceSize size, VkBufferUsageFlags usage,
                          VkMemoryPropertyFlags props,
                          VkBuffer& buffer, VkDeviceMemory& memory) {
    VkBufferCreateInfo bi{VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO};
    bi.size = size;
    bi.usage = usage;
    bi.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
    VK_CHECK(vkCreateBuffer(dev, &bi, nullptr, &buffer));

    VkMemoryRequirements req;
    vkGetBufferMemoryRequirements(dev, buffer, &req);

    VkMemoryAllocateInfo ai{VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
    ai.allocationSize = req.size;
    ai.memoryTypeIndex = find_memory_type(pd, req.memoryTypeBits, props);
    VK_CHECK(vkAllocateMemory(dev, &ai, nullptr, &memory));
    VK_CHECK(vkBindBufferMemory(dev, buffer, memory, 0));
}

// ── App implementation ─────────────────────────────────────────

App::App(const std::string& title, int width, int height)
    : width_(width), height_(height) {
    impl_ = std::make_unique<Impl>();

    glfwInit();
    glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API);
    glfwWindowHint(GLFW_RESIZABLE, GLFW_TRUE);
    impl_->window = glfwCreateWindow(width, height, title.c_str(), nullptr, nullptr);
    glfwSetWindowUserPointer(impl_->window, this);
    glfwSetKeyCallback(impl_->window, key_callback);
    glfwSetScrollCallback(impl_->window, scroll_callback);
    glfwSetCursorPosCallback(impl_->window, cursor_callback);
    glfwSetMouseButtonCallback(impl_->window, mouse_button_callback);
    glfwSetFramebufferSizeCallback(impl_->window, framebuffer_resize_callback);

    init_vulkan();
}

App::~App() {
    cleanup();
}

void App::enable_imgui() {
    imgui_enabled_ = true;
}

void App::add_mesh(MeshData data, Vec4f color) {
    meshes_.push_back({std::move(data), color, true});
    buffers_dirty_ = true;
}

void App::clear_meshes() {
    meshes_.clear();
    buffers_dirty_ = true;
}

void App::set_mesh_visible(std::size_t index, bool visible) {
    if (index < meshes_.size())
        meshes_[index].visible = visible;
}

void App::update_mesh(std::size_t index, MeshData data) {
    if (index >= meshes_.size()) return;
    meshes_[index].data = std::move(data);

    auto* d = impl_.get();
    if (index >= d->mesh_buffers.size()) return;
    auto& mb = d->mesh_buffers[index];
    auto& entry = meshes_[index];

    if (mb.vertex_buf && !entry.data.vertices.empty()) {
        void* mapped;
        VkDeviceSize vsize = entry.data.vertices.size() * sizeof(float);
        vkMapMemory(d->device, mb.vertex_mem, 0, vsize, 0, &mapped);
        std::memcpy(mapped, entry.data.vertices.data(), vsize);
        vkUnmapMemory(d->device, mb.vertex_mem);
    }
}

void App::add_point_cloud(PointCloudData data, Vec4f color) {
    points_.push_back({std::move(data), color, true});
    buffers_dirty_ = true;
}

void App::clear_point_clouds() {
    points_.clear();
    buffers_dirty_ = true;
}

void App::set_point_visible(std::size_t index, bool visible) {
    if (index < points_.size())
        points_[index].visible = visible;
}

void App::set_point_color(std::size_t index, Vec4f color) {
    if (index < points_.size())
        points_[index].color = color;
}

void App::fit_camera(float radius) {
    camera.fit(radius);
}

void App::save_screenshot(const std::string& path) {
    auto* d = impl_.get();
    vkDeviceWaitIdle(d->device);

    uint32_t w = d->swapchain_extent.width;
    uint32_t h = d->swapchain_extent.height;

    // Determine bytes per pixel from swapchain format
    uint32_t bpp = 4;
    bool is_float16 = false;
    bool is_bgra = false;
    switch (d->swapchain_format) {
        case VK_FORMAT_R16G16B16A16_SFLOAT: bpp = 8; is_float16 = true; break;
        case VK_FORMAT_B8G8R8A8_UNORM:
        case VK_FORMAT_B8G8R8A8_SRGB:       bpp = 4; is_bgra = true; break;
        default:                             bpp = 4; break;
    }
    VkDeviceSize size = static_cast<VkDeviceSize>(w) * h * bpp;

    // Create staging buffer
    VkBuffer staging_buf;
    VkDeviceMemory staging_mem;
    create_buffer(d->device, d->physical_device, size,
                  VK_BUFFER_USAGE_TRANSFER_DST_BIT,
                  VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                  staging_buf, staging_mem);

    // Allocate one-shot command buffer
    VkCommandBufferAllocateInfo cbai{VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO};
    cbai.commandPool = d->cmd_pool;
    cbai.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
    cbai.commandBufferCount = 1;
    VkCommandBuffer cmd;
    vkAllocateCommandBuffers(d->device, &cbai, &cmd);

    VkCommandBufferBeginInfo begin{VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO};
    begin.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
    vkBeginCommandBuffer(cmd, &begin);

    // Transition swapchain image: PRESENT_SRC → TRANSFER_SRC
    VkImageMemoryBarrier barrier{VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER};
    barrier.oldLayout = VK_IMAGE_LAYOUT_PRESENT_SRC_KHR;
    barrier.newLayout = VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL;
    barrier.srcAccessMask = VK_ACCESS_MEMORY_READ_BIT;
    barrier.dstAccessMask = VK_ACCESS_TRANSFER_READ_BIT;
    barrier.image = d->swapchain_images[d->current_frame % d->swapchain_images.size()];
    barrier.subresourceRange = {VK_IMAGE_ASPECT_COLOR_BIT, 0, 1, 0, 1};
    vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_TRANSFER_BIT, VK_PIPELINE_STAGE_TRANSFER_BIT,
                         0, 0, nullptr, 0, nullptr, 1, &barrier);

    // Copy image to buffer
    VkBufferImageCopy region{};
    region.imageSubresource = {VK_IMAGE_ASPECT_COLOR_BIT, 0, 0, 1};
    region.imageExtent = {w, h, 1};
    vkCmdCopyImageToBuffer(cmd, barrier.image, VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL,
                           staging_buf, 1, &region);

    // Transition back: TRANSFER_SRC → PRESENT_SRC
    barrier.oldLayout = VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL;
    barrier.newLayout = VK_IMAGE_LAYOUT_PRESENT_SRC_KHR;
    barrier.srcAccessMask = VK_ACCESS_TRANSFER_READ_BIT;
    barrier.dstAccessMask = VK_ACCESS_MEMORY_READ_BIT;
    vkCmdPipelineBarrier(cmd, VK_PIPELINE_STAGE_TRANSFER_BIT, VK_PIPELINE_STAGE_BOTTOM_OF_PIPE_BIT,
                         0, 0, nullptr, 0, nullptr, 1, &barrier);

    vkEndCommandBuffer(cmd);

    VkSubmitInfo si{VK_STRUCTURE_TYPE_SUBMIT_INFO};
    si.commandBufferCount = 1;
    si.pCommandBuffers = &cmd;
    vkQueueSubmit(d->graphics_queue, 1, &si, VK_NULL_HANDLE);
    vkQueueWaitIdle(d->graphics_queue);

    // Read pixels and convert to RGBA8 for PNG
    void* data;
    vkMapMemory(d->device, staging_mem, 0, size, 0, &data);

    std::vector<uint8_t> rgba(w * h * 4);

    if (is_float16) {
        // R16G16B16A16_SFLOAT → RGBA8: convert half-float to byte
        auto* src = static_cast<const uint16_t*>(data);
        auto half_to_float = [](uint16_t h) -> float {
            uint32_t sign = (h >> 15) & 1;
            uint32_t exp = (h >> 10) & 0x1F;
            uint32_t mant = h & 0x3FF;
            if (exp == 0) return sign ? -0.0f : 0.0f;
            if (exp == 31) return sign ? -1e30f : 1e30f;
            float f = std::ldexp(static_cast<float>(mant | 0x400), static_cast<int>(exp) - 25);
            return sign ? -f : f;
        };
        for (uint32_t i = 0; i < w * h; ++i) {
            float r = half_to_float(src[i * 4 + 0]);
            float g = half_to_float(src[i * 4 + 1]);
            float b = half_to_float(src[i * 4 + 2]);
            rgba[i * 4 + 0] = static_cast<uint8_t>(std::clamp(r, 0.0f, 1.0f) * 255.0f);
            rgba[i * 4 + 1] = static_cast<uint8_t>(std::clamp(g, 0.0f, 1.0f) * 255.0f);
            rgba[i * 4 + 2] = static_cast<uint8_t>(std::clamp(b, 0.0f, 1.0f) * 255.0f);
            rgba[i * 4 + 3] = 255;
        }
    } else {
        std::memcpy(rgba.data(), data, w * h * 4);
        if (is_bgra) {
            for (uint32_t i = 0; i < w * h; ++i)
                std::swap(rgba[i * 4 + 0], rgba[i * 4 + 2]);
        }
        for (uint32_t i = 0; i < w * h; ++i)
            rgba[i * 4 + 3] = 255;
    }

    vkUnmapMemory(d->device, staging_mem);

    stbi_write_png(path.c_str(), static_cast<int>(w), static_cast<int>(h), 4, rgba.data(), static_cast<int>(w * 4));

    // Cleanup
    vkFreeCommandBuffers(d->device, d->cmd_pool, 1, &cmd);
    vkDestroyBuffer(d->device, staging_buf, nullptr);
    vkFreeMemory(d->device, staging_mem, nullptr);

    std::fprintf(stderr, "Screenshot saved: %s (%ux%u, %s)\n", path.c_str(), w, h,
                 is_float16 ? "f16->u8" : is_bgra ? "bgra->rgba" : "rgba");
}

void App::init_vulkan() {
    auto* d = impl_.get();

    // Instance
    VkApplicationInfo app_info{VK_STRUCTURE_TYPE_APPLICATION_INFO};
    app_info.pApplicationName = "Spatium Viewer";
    app_info.apiVersion = VK_API_VERSION_1_3;

    uint32_t glfw_ext_count;
    auto** glfw_exts = glfwGetRequiredInstanceExtensions(&glfw_ext_count);

    VkInstanceCreateInfo ci{VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO};
    ci.pApplicationInfo = &app_info;
    ci.enabledExtensionCount = glfw_ext_count;
    ci.ppEnabledExtensionNames = glfw_exts;

    // Try with validation layers first, fall back without
    const char* layers[] = {"VK_LAYER_KHRONOS_validation"};
    ci.enabledLayerCount = 1;
    ci.ppEnabledLayerNames = layers;

    if (vkCreateInstance(&ci, nullptr, &d->instance) != VK_SUCCESS) {
        ci.enabledLayerCount = 0;
        ci.ppEnabledLayerNames = nullptr;
        if (vkCreateInstance(&ci, nullptr, &d->instance) != VK_SUCCESS)
            throw std::runtime_error("Failed to create Vulkan instance");
    }

    // Surface
    glfwCreateWindowSurface(d->instance, d->window, nullptr, &d->surface);

    // Physical device
    uint32_t dev_count = 0;
    vkEnumeratePhysicalDevices(d->instance, &dev_count, nullptr);
    std::vector<VkPhysicalDevice> devs(dev_count);
    vkEnumeratePhysicalDevices(d->instance, &dev_count, devs.data());
    d->physical_device = devs[0]; // pick first

    // Queue family
    uint32_t qf_count = 0;
    vkGetPhysicalDeviceQueueFamilyProperties(d->physical_device, &qf_count, nullptr);
    std::vector<VkQueueFamilyProperties> qfs(qf_count);
    vkGetPhysicalDeviceQueueFamilyProperties(d->physical_device, &qf_count, qfs.data());

    for (uint32_t i = 0; i < qf_count; ++i) {
        VkBool32 present = VK_FALSE;
        vkGetPhysicalDeviceSurfaceSupportKHR(d->physical_device, i, d->surface, &present);
        if ((qfs[i].queueFlags & VK_QUEUE_GRAPHICS_BIT) && present) {
            d->graphics_family = i;
            break;
        }
    }

    // Logical device
    float priority = 1.0f;
    VkDeviceQueueCreateInfo qci{VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO};
    qci.queueFamilyIndex = d->graphics_family;
    qci.queueCount = 1;
    qci.pQueuePriorities = &priority;

    const char* dev_exts[] = {VK_KHR_SWAPCHAIN_EXTENSION_NAME};
    VkPhysicalDeviceFeatures features{};
    features.fillModeNonSolid = VK_TRUE;  // wireframe
    features.geometryShader = VK_TRUE;    // gl_PrimitiveID in fragment

    VkDeviceCreateInfo dci{VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO};
    dci.queueCreateInfoCount = 1;
    dci.pQueueCreateInfos = &qci;
    dci.enabledExtensionCount = 1;
    dci.ppEnabledExtensionNames = dev_exts;
    dci.pEnabledFeatures = &features;

    VK_CHECK(vkCreateDevice(d->physical_device, &dci, nullptr, &d->device));
    vkGetDeviceQueue(d->device, d->graphics_family, 0, &d->graphics_queue);
    d->present_queue = d->graphics_queue;

    // Swapchain
    VkSurfaceCapabilitiesKHR caps;
    vkGetPhysicalDeviceSurfaceCapabilitiesKHR(d->physical_device, d->surface, &caps);

    uint32_t fmt_count;
    vkGetPhysicalDeviceSurfaceFormatsKHR(d->physical_device, d->surface, &fmt_count, nullptr);
    std::vector<VkSurfaceFormatKHR> fmts(fmt_count);
    vkGetPhysicalDeviceSurfaceFormatsKHR(d->physical_device, d->surface, &fmt_count, fmts.data());

    // Prefer B8G8R8A8 or R8G8B8A8 (8-bit) format for compatibility
    d->swapchain_format = fmts[0].format;
    for (auto& f : fmts) {
        if (f.format == VK_FORMAT_B8G8R8A8_SRGB || f.format == VK_FORMAT_B8G8R8A8_UNORM ||
            f.format == VK_FORMAT_R8G8B8A8_SRGB || f.format == VK_FORMAT_R8G8B8A8_UNORM) {
            d->swapchain_format = f.format;
            break;
        }
    }
    d->swapchain_extent = caps.currentExtent;
    if (d->swapchain_extent.width == UINT32_MAX) {
        d->swapchain_extent = {static_cast<uint32_t>(width_), static_cast<uint32_t>(height_)};
    }

    VkSwapchainCreateInfoKHR sci{VK_STRUCTURE_TYPE_SWAPCHAIN_CREATE_INFO_KHR};
    sci.surface = d->surface;
    sci.minImageCount = std::min(caps.minImageCount + 1,
                                 caps.maxImageCount > 0 ? caps.maxImageCount : UINT32_MAX);
    sci.imageFormat = d->swapchain_format;
    sci.imageColorSpace = fmts[0].colorSpace;
    sci.imageExtent = d->swapchain_extent;
    sci.imageArrayLayers = 1;
    sci.imageUsage = VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_TRANSFER_SRC_BIT;
    sci.imageSharingMode = VK_SHARING_MODE_EXCLUSIVE;
    sci.preTransform = caps.currentTransform;
    sci.compositeAlpha = VK_COMPOSITE_ALPHA_OPAQUE_BIT_KHR;
    sci.presentMode = VK_PRESENT_MODE_FIFO_KHR;
    sci.clipped = VK_TRUE;

    VK_CHECK(vkCreateSwapchainKHR(d->device, &sci, nullptr, &d->swapchain));

    uint32_t img_count;
    vkGetSwapchainImagesKHR(d->device, d->swapchain, &img_count, nullptr);
    d->swapchain_images.resize(img_count);
    vkGetSwapchainImagesKHR(d->device, d->swapchain, &img_count, d->swapchain_images.data());

    // Image views
    d->swapchain_views.resize(img_count);
    for (uint32_t i = 0; i < img_count; ++i) {
        VkImageViewCreateInfo ivci{VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO};
        ivci.image = d->swapchain_images[i];
        ivci.viewType = VK_IMAGE_VIEW_TYPE_2D;
        ivci.format = d->swapchain_format;
        ivci.subresourceRange = {VK_IMAGE_ASPECT_COLOR_BIT, 0, 1, 0, 1};
        VK_CHECK(vkCreateImageView(d->device, &ivci, nullptr, &d->swapchain_views[i]));
    }

    // Depth buffer
    VkImageCreateInfo depth_ci{VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO};
    depth_ci.imageType = VK_IMAGE_TYPE_2D;
    depth_ci.format = VK_FORMAT_D32_SFLOAT;
    depth_ci.extent = {d->swapchain_extent.width, d->swapchain_extent.height, 1};
    depth_ci.mipLevels = 1;
    depth_ci.arrayLayers = 1;
    depth_ci.samples = VK_SAMPLE_COUNT_1_BIT;
    depth_ci.tiling = VK_IMAGE_TILING_OPTIMAL;
    depth_ci.usage = VK_IMAGE_USAGE_DEPTH_STENCIL_ATTACHMENT_BIT;
    VK_CHECK(vkCreateImage(d->device, &depth_ci, nullptr, &d->depth_image));

    VkMemoryRequirements depth_req;
    vkGetImageMemoryRequirements(d->device, d->depth_image, &depth_req);
    VkMemoryAllocateInfo depth_ai{VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
    depth_ai.allocationSize = depth_req.size;
    depth_ai.memoryTypeIndex = find_memory_type(d->physical_device, depth_req.memoryTypeBits,
                                                 VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
    VK_CHECK(vkAllocateMemory(d->device, &depth_ai, nullptr, &d->depth_memory));
    VK_CHECK(vkBindImageMemory(d->device, d->depth_image, d->depth_memory, 0));

    VkImageViewCreateInfo depth_vci{VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO};
    depth_vci.image = d->depth_image;
    depth_vci.viewType = VK_IMAGE_VIEW_TYPE_2D;
    depth_vci.format = VK_FORMAT_D32_SFLOAT;
    depth_vci.subresourceRange = {VK_IMAGE_ASPECT_DEPTH_BIT, 0, 1, 0, 1};
    VK_CHECK(vkCreateImageView(d->device, &depth_vci, nullptr, &d->depth_view));

    // Render pass
    VkAttachmentDescription attachments[2]{};
    attachments[0].format = d->swapchain_format;
    attachments[0].samples = VK_SAMPLE_COUNT_1_BIT;
    attachments[0].loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
    attachments[0].storeOp = VK_ATTACHMENT_STORE_OP_STORE;
    attachments[0].initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
    attachments[0].finalLayout = VK_IMAGE_LAYOUT_PRESENT_SRC_KHR;

    attachments[1].format = VK_FORMAT_D32_SFLOAT;
    attachments[1].samples = VK_SAMPLE_COUNT_1_BIT;
    attachments[1].loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
    attachments[1].storeOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
    attachments[1].initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
    attachments[1].finalLayout = VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL;

    VkAttachmentReference color_ref{0, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL};
    VkAttachmentReference depth_ref{1, VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL};

    VkSubpassDescription subpass{};
    subpass.pipelineBindPoint = VK_PIPELINE_BIND_POINT_GRAPHICS;
    subpass.colorAttachmentCount = 1;
    subpass.pColorAttachments = &color_ref;
    subpass.pDepthStencilAttachment = &depth_ref;

    VkRenderPassCreateInfo rpci{VK_STRUCTURE_TYPE_RENDER_PASS_CREATE_INFO};
    rpci.attachmentCount = 2;
    rpci.pAttachments = attachments;
    rpci.subpassCount = 1;
    rpci.pSubpasses = &subpass;

    VK_CHECK(vkCreateRenderPass(d->device, &rpci, nullptr, &d->render_pass));

    // Framebuffers
    d->framebuffers.resize(img_count);
    for (uint32_t i = 0; i < img_count; ++i) {
        VkImageView views[] = {d->swapchain_views[i], d->depth_view};
        VkFramebufferCreateInfo fci{VK_STRUCTURE_TYPE_FRAMEBUFFER_CREATE_INFO};
        fci.renderPass = d->render_pass;
        fci.attachmentCount = 2;
        fci.pAttachments = views;
        fci.width = d->swapchain_extent.width;
        fci.height = d->swapchain_extent.height;
        fci.layers = 1;
        VK_CHECK(vkCreateFramebuffer(d->device, &fci, nullptr, &d->framebuffers[i]));
    }

    // Command pool + buffers
    VkCommandPoolCreateInfo cpci{VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO};
    cpci.flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;
    cpci.queueFamilyIndex = d->graphics_family;
    VK_CHECK(vkCreateCommandPool(d->device, &cpci, nullptr, &d->cmd_pool));

    d->cmd_buffers.resize(img_count);
    VkCommandBufferAllocateInfo cbai{VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO};
    cbai.commandPool = d->cmd_pool;
    cbai.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
    cbai.commandBufferCount = img_count;
    VK_CHECK(vkAllocateCommandBuffers(d->device, &cbai, d->cmd_buffers.data()));

    // Sync — per-frame to avoid semaphore reuse across swapchain images
    d->sem_available.resize(img_count);
    d->sem_finished.resize(img_count);
    d->fences_in_flight.resize(img_count);
    VkSemaphoreCreateInfo semi{VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO};
    VkFenceCreateInfo fi{VK_STRUCTURE_TYPE_FENCE_CREATE_INFO};
    fi.flags = VK_FENCE_CREATE_SIGNALED_BIT;
    for (uint32_t i = 0; i < img_count; ++i) {
        VK_CHECK(vkCreateSemaphore(d->device, &semi, nullptr, &d->sem_available[i]));
        VK_CHECK(vkCreateSemaphore(d->device, &semi, nullptr, &d->sem_finished[i]));
        VK_CHECK(vkCreateFence(d->device, &fi, nullptr, &d->fences_in_flight[i]));
    }

    create_pipeline();
}

void App::create_pipeline() {
    auto* d = impl_.get();

    // Push constants: mvp(64) + color(16) + point_size(4) + clip_plane(16) = 100 bytes
    VkPushConstantRange pc_range{VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, 0, 100};
    VkPipelineLayoutCreateInfo plci{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
    plci.pushConstantRangeCount = 1;
    plci.pPushConstantRanges = &pc_range;
    VK_CHECK(vkCreatePipelineLayout(d->device, &plci, nullptr, &d->pipeline_layout));

    // Shared state for all pipelines
    VkPipelineViewportStateCreateInfo viewport{VK_STRUCTURE_TYPE_PIPELINE_VIEWPORT_STATE_CREATE_INFO};
    viewport.viewportCount = 1;
    viewport.scissorCount = 1;

    VkPipelineMultisampleStateCreateInfo msaa{VK_STRUCTURE_TYPE_PIPELINE_MULTISAMPLE_STATE_CREATE_INFO};
    msaa.rasterizationSamples = VK_SAMPLE_COUNT_1_BIT;

    VkPipelineDepthStencilStateCreateInfo depth{VK_STRUCTURE_TYPE_PIPELINE_DEPTH_STENCIL_STATE_CREATE_INFO};
    depth.depthTestEnable = VK_TRUE;
    depth.depthWriteEnable = VK_TRUE;
    depth.depthCompareOp = VK_COMPARE_OP_LESS;

    VkPipelineColorBlendAttachmentState blend_att{};
    blend_att.colorWriteMask = VK_COLOR_COMPONENT_R_BIT | VK_COLOR_COMPONENT_G_BIT |
                               VK_COLOR_COMPONENT_B_BIT | VK_COLOR_COMPONENT_A_BIT;
    blend_att.blendEnable = VK_TRUE;
    blend_att.srcColorBlendFactor = VK_BLEND_FACTOR_SRC_ALPHA;
    blend_att.dstColorBlendFactor = VK_BLEND_FACTOR_ONE_MINUS_SRC_ALPHA;
    blend_att.colorBlendOp = VK_BLEND_OP_ADD;
    blend_att.srcAlphaBlendFactor = VK_BLEND_FACTOR_ONE;
    blend_att.dstAlphaBlendFactor = VK_BLEND_FACTOR_ZERO;
    blend_att.alphaBlendOp = VK_BLEND_OP_ADD;
    VkPipelineColorBlendStateCreateInfo blend{VK_STRUCTURE_TYPE_PIPELINE_COLOR_BLEND_STATE_CREATE_INFO};
    blend.attachmentCount = 1;
    blend.pAttachments = &blend_att;

    VkDynamicState dyn_states[] = {VK_DYNAMIC_STATE_VIEWPORT, VK_DYNAMIC_STATE_SCISSOR};
    VkPipelineDynamicStateCreateInfo dyn{VK_STRUCTURE_TYPE_PIPELINE_DYNAMIC_STATE_CREATE_INFO};
    dyn.dynamicStateCount = 2;
    dyn.pDynamicStates = dyn_states;

    // Helper: build a pipeline with specific topology, shaders, vertex layout, raster mode
    auto build = [&](VkPrimitiveTopology topo,
                     VkShaderModule vert_mod, VkShaderModule frag_mod,
                     uint32_t stride, const VkVertexInputAttributeDescription* attrs, uint32_t attr_count,
                     VkPolygonMode poly_mode, VkCullModeFlags cull, bool depth_bias) -> VkPipeline
    {
        VkPipelineShaderStageCreateInfo stages[2]{};
        stages[0] = {VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO};
        stages[0].stage = VK_SHADER_STAGE_VERTEX_BIT;
        stages[0].module = vert_mod;
        stages[0].pName = "main";
        stages[1] = {VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO};
        stages[1].stage = VK_SHADER_STAGE_FRAGMENT_BIT;
        stages[1].module = frag_mod;
        stages[1].pName = "main";

        VkVertexInputBindingDescription binding{0, stride, VK_VERTEX_INPUT_RATE_VERTEX};
        VkPipelineVertexInputStateCreateInfo vertex_input{VK_STRUCTURE_TYPE_PIPELINE_VERTEX_INPUT_STATE_CREATE_INFO};
        vertex_input.vertexBindingDescriptionCount = 1;
        vertex_input.pVertexBindingDescriptions = &binding;
        vertex_input.vertexAttributeDescriptionCount = attr_count;
        vertex_input.pVertexAttributeDescriptions = attrs;

        VkPipelineInputAssemblyStateCreateInfo assembly{VK_STRUCTURE_TYPE_PIPELINE_INPUT_ASSEMBLY_STATE_CREATE_INFO};
        assembly.topology = topo;

        VkPipelineRasterizationStateCreateInfo raster{VK_STRUCTURE_TYPE_PIPELINE_RASTERIZATION_STATE_CREATE_INFO};
        raster.polygonMode = poly_mode;
        raster.cullMode = cull;
        raster.frontFace = VK_FRONT_FACE_COUNTER_CLOCKWISE;
        raster.lineWidth = 1.0f;
        if (depth_bias) {
            raster.depthBiasEnable = VK_TRUE;
            raster.depthBiasConstantFactor = -1.0f;
            raster.depthBiasSlopeFactor = -1.0f;
        }

        VkGraphicsPipelineCreateInfo pci{VK_STRUCTURE_TYPE_GRAPHICS_PIPELINE_CREATE_INFO};
        pci.stageCount = 2;
        pci.pStages = stages;
        pci.pVertexInputState = &vertex_input;
        pci.pInputAssemblyState = &assembly;
        pci.pViewportState = &viewport;
        pci.pRasterizationState = &raster;
        pci.pMultisampleState = &msaa;
        pci.pDepthStencilState = &depth;
        pci.pColorBlendState = &blend;
        pci.pDynamicState = &dyn;
        pci.layout = d->pipeline_layout;
        pci.renderPass = d->render_pass;

        VkPipeline pipeline;
        VK_CHECK(vkCreateGraphicsPipelines(d->device, VK_NULL_HANDLE, 1, &pci, nullptr, &pipeline));
        return pipeline;
    };

    // Compile shaders
    VkShaderModuleCreateInfo smci{VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO};
    auto make_module = [&](const char* src, shaderc_shader_kind kind) {
        auto spv = compile_glsl(src, kind);
        smci.codeSize = spv.size() * 4;
        smci.pCode = spv.data();
        VkShaderModule mod;
        VK_CHECK(vkCreateShaderModule(d->device, &smci, nullptr, &mod));
        return mod;
    };

    VkShaderModule mesh_vert  = make_module(VERT_SRC, shaderc_vertex_shader);
    VkShaderModule mesh_frag  = make_module(FRAG_SRC, shaderc_fragment_shader);
    VkShaderModule edge_frag  = make_module(FRAG_EDGE_SRC, shaderc_fragment_shader);
    VkShaderModule point_vert = make_module(POINT_VERT_SRC, shaderc_vertex_shader);
    VkShaderModule point_frag = make_module(POINT_FRAG_SRC, shaderc_fragment_shader);

    // Mesh vertex layout: pos(3f) + normal(3f) = 24 bytes
    VkVertexInputAttributeDescription mesh_attrs[2]{};
    mesh_attrs[0] = {0, 0, VK_FORMAT_R32G32B32_SFLOAT, 0};
    mesh_attrs[1] = {1, 0, VK_FORMAT_R32G32B32_SFLOAT, 12};

    // Point vertex layout: pos(3f) = 12 bytes
    VkVertexInputAttributeDescription point_attr{0, 0, VK_FORMAT_R32G32B32_SFLOAT, 0};

    d->pipeline_solid  = build(VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST, mesh_vert, mesh_frag,
                               24, mesh_attrs, 2, VK_POLYGON_MODE_FILL, VK_CULL_MODE_NONE, false);
    d->pipeline_wire   = build(VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST, mesh_vert, edge_frag,
                               24, mesh_attrs, 2, VK_POLYGON_MODE_LINE, VK_CULL_MODE_NONE, true);
    d->pipeline_points = build(VK_PRIMITIVE_TOPOLOGY_POINT_LIST, point_vert, point_frag,
                               12, &point_attr, 1, VK_POLYGON_MODE_FILL, VK_CULL_MODE_NONE, false);

    vkDestroyShaderModule(d->device, mesh_vert, nullptr);
    vkDestroyShaderModule(d->device, mesh_frag, nullptr);
    vkDestroyShaderModule(d->device, edge_frag, nullptr);
    vkDestroyShaderModule(d->device, point_vert, nullptr);
    vkDestroyShaderModule(d->device, point_frag, nullptr);
}

void App::create_buffers() {
    auto* d = impl_.get();

    // Destroy old mesh buffers
    for (auto& mb : d->mesh_buffers) {
        if (mb.vertex_buf) { vkDestroyBuffer(d->device, mb.vertex_buf, nullptr); vkFreeMemory(d->device, mb.vertex_mem, nullptr); }
        if (mb.index_buf) { vkDestroyBuffer(d->device, mb.index_buf, nullptr); vkFreeMemory(d->device, mb.index_mem, nullptr); }
    }
    d->mesh_buffers.clear();
    d->mesh_buffers.resize(meshes_.size());

    for (std::size_t i = 0; i < meshes_.size(); ++i) {
        auto& entry = meshes_[i];
        auto& mb = d->mesh_buffers[i];
        if (entry.data.vertices.empty()) continue;

        VkDeviceSize vsize = entry.data.vertices.size() * sizeof(float);
        VkDeviceSize isize = entry.data.indices.size() * sizeof(uint32_t);

        create_buffer(d->device, d->physical_device, vsize,
                      VK_BUFFER_USAGE_VERTEX_BUFFER_BIT,
                      VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                      mb.vertex_buf, mb.vertex_mem);

        void* data;
        vkMapMemory(d->device, mb.vertex_mem, 0, vsize, 0, &data);
        std::memcpy(data, entry.data.vertices.data(), vsize);
        vkUnmapMemory(d->device, mb.vertex_mem);

        create_buffer(d->device, d->physical_device, isize,
                      VK_BUFFER_USAGE_INDEX_BUFFER_BIT,
                      VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                      mb.index_buf, mb.index_mem);

        vkMapMemory(d->device, mb.index_mem, 0, isize, 0, &data);
        std::memcpy(data, entry.data.indices.data(), isize);
        vkUnmapMemory(d->device, mb.index_mem);
    }

    // Point cloud buffers
    for (auto& pb : d->point_buffers) {
        if (pb.vertex_buf) { vkDestroyBuffer(d->device, pb.vertex_buf, nullptr); vkFreeMemory(d->device, pb.vertex_mem, nullptr); }
    }
    d->point_buffers.clear();
    d->point_buffers.resize(points_.size());

    for (std::size_t i = 0; i < points_.size(); ++i) {
        auto& entry = points_[i];
        auto& pb = d->point_buffers[i];
        if (entry.data.positions.empty()) continue;

        VkDeviceSize vsize = entry.data.positions.size() * sizeof(float);

        create_buffer(d->device, d->physical_device, vsize,
                      VK_BUFFER_USAGE_VERTEX_BUFFER_BIT,
                      VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                      pb.vertex_buf, pb.vertex_mem);

        void* data;
        vkMapMemory(d->device, pb.vertex_mem, 0, vsize, 0, &data);
        std::memcpy(data, entry.data.positions.data(), vsize);
        vkUnmapMemory(d->device, pb.vertex_mem);
    }
}

void App::draw_frame() {
    auto* d = impl_.get();

    if (buffers_dirty_) {
        vkDeviceWaitIdle(d->device);
        create_buffers();
        buffers_dirty_ = false;
    }

    uint32_t frame = d->current_frame;
    vkWaitForFences(d->device, 1, &d->fences_in_flight[frame], VK_TRUE, UINT64_MAX);

    uint32_t img_idx;
    VkResult acquire_result = vkAcquireNextImageKHR(d->device, d->swapchain, UINT64_MAX, d->sem_available[frame], VK_NULL_HANDLE, &img_idx);
    if (acquire_result == VK_ERROR_OUT_OF_DATE_KHR) {
        recreate_swapchain();
        return;
    }

    vkResetFences(d->device, 1, &d->fences_in_flight[frame]);

    auto cmd = d->cmd_buffers[frame];
    vkResetCommandBuffer(cmd, 0);

    VkCommandBufferBeginInfo begin{VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO};
    vkBeginCommandBuffer(cmd, &begin);

    VkClearValue clears[2]{};
    clears[0].color = {{0.118f, 0.118f, 0.180f, 1.0f}}; // #1e1e2e catppuccin base
    clears[1].depthStencil = {1.0f, 0};

    VkRenderPassBeginInfo rpbi{VK_STRUCTURE_TYPE_RENDER_PASS_BEGIN_INFO};
    rpbi.renderPass = d->render_pass;
    rpbi.framebuffer = d->framebuffers[img_idx];
    rpbi.renderArea = {{0, 0}, d->swapchain_extent};
    rpbi.clearValueCount = 2;
    rpbi.pClearValues = clears;

    vkCmdBeginRenderPass(cmd, &rpbi, VK_SUBPASS_CONTENTS_INLINE);

    {
        int fb_w, fb_h;
        glfwGetFramebufferSize(d->window, &fb_w, &fb_h);
        float aspect = (fb_w > 0 && fb_h > 0) ? (float)fb_w / (float)fb_h : 1.0f;
        auto mvp = camera.view_projection(aspect);

        VkViewport vp{0, 0, (float)d->swapchain_extent.width, (float)d->swapchain_extent.height, 0, 1};
        VkRect2D scissor{{0, 0}, d->swapchain_extent};
        vkCmdSetViewport(cmd, 0, 1, &vp);
        vkCmdSetScissor(cmd, 0, 1, &scissor);

        VkShaderStageFlags pc_stages = VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT;

        // Draw triangle meshes
        for (std::size_t i = 0; i < meshes_.size(); ++i) {
            auto& entry = meshes_[i];
            if (i >= d->mesh_buffers.size()) break;
            auto& mb = d->mesh_buffers[i];
            if (!entry.visible || !mb.vertex_buf || entry.data.index_count == 0) continue;

            VkDeviceSize offset = 0;
            vkCmdBindVertexBuffers(cmd, 0, 1, &mb.vertex_buf, &offset);
            vkCmdBindIndexBuffer(cmd, mb.index_buf, 0, VK_INDEX_TYPE_UINT32);

            struct { float mvp[16]; float color[4]; float point_size; } pc{};
            std::memcpy(pc.mvp, mvp.data(), sizeof(pc.mvp));
            pc.color[0] = entry.color[0];
            pc.color[1] = entry.color[1];
            pc.color[2] = entry.color[2];
            pc.color[3] = entry.color[3];
            pc.point_size = point_size_;

            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, d->pipeline_solid);
            vkCmdPushConstants(cmd, d->pipeline_layout, pc_stages, 0, 84, &pc);
            vkCmdDrawIndexed(cmd, entry.data.index_count, 1, 0, 0, 0);

            if (wireframe_) {
                vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, d->pipeline_wire);
                vkCmdPushConstants(cmd, d->pipeline_layout, pc_stages, 0, 84, &pc);
                vkCmdDrawIndexed(cmd, entry.data.index_count, 1, 0, 0, 0);
            }
        }

        // Draw point clouds
        if (d->pipeline_points && !d->point_buffers.empty()) {
            vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_GRAPHICS, d->pipeline_points);
            for (std::size_t i = 0; i < points_.size(); ++i) {
                auto& entry = points_[i];
                if (i >= d->point_buffers.size()) break;
                auto& pb = d->point_buffers[i];
                if (!entry.visible || !pb.vertex_buf || entry.data.point_count == 0) continue;

                VkDeviceSize offset = 0;
                vkCmdBindVertexBuffers(cmd, 0, 1, &pb.vertex_buf, &offset);

                // Layout: mvp(64) + color(16) + clip_plane(16) + point_size(4) = 100
                struct { float mvp[16]; float color[4]; float clip_plane[4]; float point_size; } pc{};
                std::memcpy(pc.mvp, mvp.data(), sizeof(pc.mvp));
                pc.color[0] = entry.color[0];
                pc.color[1] = entry.color[1];
                pc.color[2] = entry.color[2];
                pc.color[3] = entry.color[3];

                // Quarter cutaway: toggle only; shader uses fixed world +X/+Z quadrant
                pc.clip_plane[3] = cutaway_enabled_ ? 0.0f : -99999.0f;
                pc.point_size = point_size_;

                vkCmdPushConstants(cmd, d->pipeline_layout, pc_stages, 0, 100, &pc);
                vkCmdDraw(cmd, entry.data.point_count, 1, 0, 0);
            }
        }
    }

#if SPATIUM_HAS_IMGUI
    if (imgui_enabled_) {
        ImGui_ImplVulkan_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        if (gui_cb_) gui_cb_();
        ImGui::Render();
        ImGui_ImplVulkan_RenderDrawData(ImGui::GetDrawData(), cmd);
    }
#endif

    vkCmdEndRenderPass(cmd);
    vkEndCommandBuffer(cmd);

    VkPipelineStageFlags wait_stage = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT;
    VkSubmitInfo si{VK_STRUCTURE_TYPE_SUBMIT_INFO};
    si.waitSemaphoreCount = 1;
    si.pWaitSemaphores = &d->sem_available[frame];
    si.pWaitDstStageMask = &wait_stage;
    si.commandBufferCount = 1;
    si.pCommandBuffers = &cmd;
    si.signalSemaphoreCount = 1;
    si.pSignalSemaphores = &d->sem_finished[frame];

    vkQueueSubmit(d->graphics_queue, 1, &si, d->fences_in_flight[frame]);

    VkPresentInfoKHR pi{VK_STRUCTURE_TYPE_PRESENT_INFO_KHR};
    pi.waitSemaphoreCount = 1;
    pi.pWaitSemaphores = &d->sem_finished[frame];
    pi.swapchainCount = 1;
    pi.pSwapchains = &d->swapchain;
    pi.pImageIndices = &img_idx;

    VkResult present_result = vkQueuePresentKHR(d->present_queue, &pi);
    if (present_result == VK_ERROR_OUT_OF_DATE_KHR || present_result == VK_SUBOPTIMAL_KHR || d->framebuffer_resized) {
        d->framebuffer_resized = false;
        recreate_swapchain();
    }
    d->current_frame = (frame + 1) % static_cast<uint32_t>(d->sem_available.size());
}

void App::update_mesh_vertices(std::size_t index, const MeshData& data) {
    if (index >= meshes_.size()) return;
    meshes_[index].data.vertices = data.vertices;
    meshes_[index].data.vertex_count = data.vertex_count;

    auto* d = impl_.get();
    if (index >= d->mesh_buffers.size()) return;
    auto& mb = d->mesh_buffers[index];

    if (mb.vertex_buf && !data.vertices.empty()) {
        void* mapped;
        VkDeviceSize vsize = data.vertices.size() * sizeof(float);
        vkMapMemory(d->device, mb.vertex_mem, 0, vsize, 0, &mapped);
        std::memcpy(mapped, data.vertices.data(), vsize);
        vkUnmapMemory(d->device, mb.vertex_mem);
    }
}

void App::init_imgui() {
#if SPATIUM_HAS_IMGUI
    auto* d = impl_.get();

    // Descriptor pool for ImGui
    VkDescriptorPoolSize pool_size{VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 100};
    VkDescriptorPoolCreateInfo dpci{VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO};
    dpci.flags = VK_DESCRIPTOR_POOL_CREATE_FREE_DESCRIPTOR_SET_BIT;
    dpci.maxSets = 100;
    dpci.poolSizeCount = 1;
    dpci.pPoolSizes = &pool_size;
    VK_CHECK(vkCreateDescriptorPool(d->device, &dpci, nullptr, &d->imgui_pool));

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui::StyleColorsDark();

    // Scale UI for readability
    ImGui::GetStyle().ScaleAllSizes(1.2f);
    ImGui::GetIO().FontGlobalScale = 1.2f;

    ImGui_ImplGlfw_InitForVulkan(d->window, true);

    ImGui_ImplVulkan_InitInfo info{};
    info.Instance = d->instance;
    info.PhysicalDevice = d->physical_device;
    info.Device = d->device;
    info.QueueFamily = d->graphics_family;
    info.Queue = d->graphics_queue;
    info.DescriptorPool = d->imgui_pool;
    info.RenderPass = d->render_pass;
    info.MinImageCount = static_cast<uint32_t>(d->swapchain_images.size());
    info.ImageCount = static_cast<uint32_t>(d->swapchain_images.size());
    info.MSAASamples = VK_SAMPLE_COUNT_1_BIT;
    info.Subpass = 0;
    ImGui_ImplVulkan_Init(&info);

    // Upload font atlas
    ImGui_ImplVulkan_CreateFontsTexture();
#endif
}

void App::shutdown_imgui() {
#if SPATIUM_HAS_IMGUI
    if (!imgui_enabled_) return;
    auto* d = impl_.get();
    vkDeviceWaitIdle(d->device);
    ImGui_ImplVulkan_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    if (d->imgui_pool) {
        vkDestroyDescriptorPool(d->device, d->imgui_pool, nullptr);
        d->imgui_pool = VK_NULL_HANDLE;
    }
#endif
}

void App::run() {
    create_buffers();
    if (imgui_enabled_) init_imgui();

    while (!glfwWindowShouldClose(impl_->window)) {
        glfwPollEvents();
        if (frame_cb_) frame_cb_();
        draw_frame();
    }

    vkDeviceWaitIdle(impl_->device);
}

void App::cleanup() {
    auto* d = impl_.get();
    if (!d || !d->device) return;

    vkDeviceWaitIdle(d->device);
    shutdown_imgui();

    for (auto& mb : d->mesh_buffers) {
        if (mb.vertex_buf) { vkDestroyBuffer(d->device, mb.vertex_buf, nullptr); vkFreeMemory(d->device, mb.vertex_mem, nullptr); }
        if (mb.index_buf) { vkDestroyBuffer(d->device, mb.index_buf, nullptr); vkFreeMemory(d->device, mb.index_mem, nullptr); }
    }
    d->mesh_buffers.clear();

    for (auto& pb : d->point_buffers) {
        if (pb.vertex_buf) { vkDestroyBuffer(d->device, pb.vertex_buf, nullptr); vkFreeMemory(d->device, pb.vertex_mem, nullptr); }
    }
    d->point_buffers.clear();

    for (size_t i = 0; i < d->fences_in_flight.size(); ++i) {
        vkDestroyFence(d->device, d->fences_in_flight[i], nullptr);
        vkDestroySemaphore(d->device, d->sem_finished[i], nullptr);
        vkDestroySemaphore(d->device, d->sem_available[i], nullptr);
    }
    vkDestroyCommandPool(d->device, d->cmd_pool, nullptr);

    if (d->pipeline_points) vkDestroyPipeline(d->device, d->pipeline_points, nullptr);
    if (d->pipeline_wire) vkDestroyPipeline(d->device, d->pipeline_wire, nullptr);
    if (d->pipeline_solid) vkDestroyPipeline(d->device, d->pipeline_solid, nullptr);
    if (d->pipeline_layout) vkDestroyPipelineLayout(d->device, d->pipeline_layout, nullptr);

    for (auto fb : d->framebuffers) vkDestroyFramebuffer(d->device, fb, nullptr);
    vkDestroyRenderPass(d->device, d->render_pass, nullptr);

    vkDestroyImageView(d->device, d->depth_view, nullptr);
    vkDestroyImage(d->device, d->depth_image, nullptr);
    vkFreeMemory(d->device, d->depth_memory, nullptr);

    for (auto iv : d->swapchain_views) vkDestroyImageView(d->device, iv, nullptr);
    vkDestroySwapchainKHR(d->device, d->swapchain, nullptr);
    vkDestroyDevice(d->device, nullptr);
    vkDestroySurfaceKHR(d->instance, d->surface, nullptr);
    vkDestroyInstance(d->instance, nullptr);

    glfwDestroyWindow(d->window);
    glfwTerminate();
}

// ── Input callbacks ────────────────────────────────────────────

void App::key_callback(GLFWwindow* w, int key, int /*scancode*/, int action, int mods) {
    if (action != GLFW_PRESS) return;
#if SPATIUM_HAS_IMGUI
    if (ImGui::GetCurrentContext() && ImGui::GetIO().WantCaptureKeyboard) return;
#endif
    auto* app = static_cast<App*>(glfwGetWindowUserPointer(w));
    if (key == GLFW_KEY_ESCAPE || key == GLFW_KEY_Q) glfwSetWindowShouldClose(w, GLFW_TRUE);
    if (key == GLFW_KEY_W) app->wireframe_ = !app->wireframe_;
    if (key == GLFW_KEY_P) app->point_size_ = std::min(app->point_size_ + 0.5f, 10.0f);
    if (key == GLFW_KEY_O) app->point_size_ = std::max(app->point_size_ - 0.5f, 1.0f);
    if (key == GLFW_KEY_C) app->cutaway_enabled_ = !app->cutaway_enabled_;
    if (key == GLFW_KEY_S) {
        auto now = std::chrono::system_clock::now();
        auto time = std::chrono::system_clock::to_time_t(now);
        char buf[64];
        std::strftime(buf, sizeof(buf), "screenshot_%Y%m%d_%H%M%S.png", std::localtime(&time));
        app->save_screenshot(buf);
    }
    if (app->user_key_cb_) app->user_key_cb_(key, action, mods);
}

void App::scroll_callback(GLFWwindow* w, double /*xoff*/, double yoff) {
#if SPATIUM_HAS_IMGUI
    if (ImGui::GetCurrentContext() && ImGui::GetIO().WantCaptureMouse) return;
#endif
    auto* app = static_cast<App*>(glfwGetWindowUserPointer(w));
    app->camera.zoom(static_cast<float>(yoff));
}

void App::cursor_callback(GLFWwindow* w, double x, double y) {
#if SPATIUM_HAS_IMGUI
    if (ImGui::GetCurrentContext() && ImGui::GetIO().WantCaptureMouse) return;
#endif
    auto* app = static_cast<App*>(glfwGetWindowUserPointer(w));
    auto* d = app->impl_.get();
    if (d->dragging) {
        float dx = static_cast<float>(x - d->last_x);
        float dy = static_cast<float>(y - d->last_y);
        if (glfwGetMouseButton(w, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS)
            app->camera.pan(dx, dy);
        else
            app->camera.orbit(dx, dy);
    }
    d->last_x = x;
    d->last_y = y;
}

void App::mouse_button_callback(GLFWwindow* w, int button, int action, int /*mods*/) {
#if SPATIUM_HAS_IMGUI
    if (ImGui::GetCurrentContext() && ImGui::GetIO().WantCaptureMouse) return;
#endif
    auto* app = static_cast<App*>(glfwGetWindowUserPointer(w));
    if (button == GLFW_MOUSE_BUTTON_LEFT || button == GLFW_MOUSE_BUTTON_RIGHT)
        app->impl_->dragging = (action == GLFW_PRESS);
}

void App::framebuffer_resize_callback(GLFWwindow* w, int /*width*/, int /*height*/) {
    auto* app = static_cast<App*>(glfwGetWindowUserPointer(w));
    app->impl_->framebuffer_resized = true;
}

void App::recreate_swapchain() {
    auto* d = impl_.get();

    int w = 0, h = 0;
    glfwGetFramebufferSize(d->window, &w, &h);
    while (w == 0 || h == 0) {
        glfwGetFramebufferSize(d->window, &w, &h);
        glfwWaitEvents();
    }

    vkDeviceWaitIdle(d->device);

    // Destroy old framebuffers, depth, image views
    for (auto fb : d->framebuffers) vkDestroyFramebuffer(d->device, fb, nullptr);
    vkDestroyImageView(d->device, d->depth_view, nullptr);
    vkDestroyImage(d->device, d->depth_image, nullptr);
    vkFreeMemory(d->device, d->depth_memory, nullptr);
    for (auto iv : d->swapchain_views) vkDestroyImageView(d->device, iv, nullptr);

    // Recreate swapchain
    VkSurfaceCapabilitiesKHR caps;
    vkGetPhysicalDeviceSurfaceCapabilitiesKHR(d->physical_device, d->surface, &caps);

    d->swapchain_extent = caps.currentExtent;
    if (d->swapchain_extent.width == UINT32_MAX) {
        d->swapchain_extent = {static_cast<uint32_t>(w), static_cast<uint32_t>(h)};
    }

    uint32_t fmt_count;
    vkGetPhysicalDeviceSurfaceFormatsKHR(d->physical_device, d->surface, &fmt_count, nullptr);
    std::vector<VkSurfaceFormatKHR> fmts(fmt_count);
    vkGetPhysicalDeviceSurfaceFormatsKHR(d->physical_device, d->surface, &fmt_count, fmts.data());

    VkSwapchainKHR old_swapchain = d->swapchain;
    VkSwapchainCreateInfoKHR sci{VK_STRUCTURE_TYPE_SWAPCHAIN_CREATE_INFO_KHR};
    sci.surface = d->surface;
    sci.minImageCount = std::min(caps.minImageCount + 1,
                                 caps.maxImageCount > 0 ? caps.maxImageCount : UINT32_MAX);
    sci.imageFormat = d->swapchain_format;
    sci.imageColorSpace = fmts[0].colorSpace;
    sci.imageExtent = d->swapchain_extent;
    sci.imageArrayLayers = 1;
    sci.imageUsage = VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_TRANSFER_SRC_BIT;
    sci.imageSharingMode = VK_SHARING_MODE_EXCLUSIVE;
    sci.preTransform = caps.currentTransform;
    sci.compositeAlpha = VK_COMPOSITE_ALPHA_OPAQUE_BIT_KHR;
    sci.presentMode = VK_PRESENT_MODE_FIFO_KHR;
    sci.clipped = VK_TRUE;
    sci.oldSwapchain = old_swapchain;

    VK_CHECK(vkCreateSwapchainKHR(d->device, &sci, nullptr, &d->swapchain));
    vkDestroySwapchainKHR(d->device, old_swapchain, nullptr);

    uint32_t img_count;
    vkGetSwapchainImagesKHR(d->device, d->swapchain, &img_count, nullptr);
    d->swapchain_images.resize(img_count);
    vkGetSwapchainImagesKHR(d->device, d->swapchain, &img_count, d->swapchain_images.data());

    // Image views
    d->swapchain_views.resize(img_count);
    for (uint32_t i = 0; i < img_count; ++i) {
        VkImageViewCreateInfo ivci{VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO};
        ivci.image = d->swapchain_images[i];
        ivci.viewType = VK_IMAGE_VIEW_TYPE_2D;
        ivci.format = d->swapchain_format;
        ivci.subresourceRange = {VK_IMAGE_ASPECT_COLOR_BIT, 0, 1, 0, 1};
        VK_CHECK(vkCreateImageView(d->device, &ivci, nullptr, &d->swapchain_views[i]));
    }

    // Depth buffer
    VkImageCreateInfo depth_ci{VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO};
    depth_ci.imageType = VK_IMAGE_TYPE_2D;
    depth_ci.format = VK_FORMAT_D32_SFLOAT;
    depth_ci.extent = {d->swapchain_extent.width, d->swapchain_extent.height, 1};
    depth_ci.mipLevels = 1;
    depth_ci.arrayLayers = 1;
    depth_ci.samples = VK_SAMPLE_COUNT_1_BIT;
    depth_ci.tiling = VK_IMAGE_TILING_OPTIMAL;
    depth_ci.usage = VK_IMAGE_USAGE_DEPTH_STENCIL_ATTACHMENT_BIT;
    VK_CHECK(vkCreateImage(d->device, &depth_ci, nullptr, &d->depth_image));

    VkMemoryRequirements depth_req;
    vkGetImageMemoryRequirements(d->device, d->depth_image, &depth_req);
    VkMemoryAllocateInfo depth_ai{VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
    depth_ai.allocationSize = depth_req.size;
    depth_ai.memoryTypeIndex = find_memory_type(d->physical_device, depth_req.memoryTypeBits,
                                                 VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
    VK_CHECK(vkAllocateMemory(d->device, &depth_ai, nullptr, &d->depth_memory));
    VK_CHECK(vkBindImageMemory(d->device, d->depth_image, d->depth_memory, 0));

    VkImageViewCreateInfo depth_vci{VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO};
    depth_vci.image = d->depth_image;
    depth_vci.viewType = VK_IMAGE_VIEW_TYPE_2D;
    depth_vci.format = VK_FORMAT_D32_SFLOAT;
    depth_vci.subresourceRange = {VK_IMAGE_ASPECT_DEPTH_BIT, 0, 1, 0, 1};
    VK_CHECK(vkCreateImageView(d->device, &depth_vci, nullptr, &d->depth_view));

    // Framebuffers
    d->framebuffers.resize(img_count);
    for (uint32_t i = 0; i < img_count; ++i) {
        VkImageView views[] = {d->swapchain_views[i], d->depth_view};
        VkFramebufferCreateInfo fci{VK_STRUCTURE_TYPE_FRAMEBUFFER_CREATE_INFO};
        fci.renderPass = d->render_pass;
        fci.attachmentCount = 2;
        fci.pAttachments = views;
        fci.width = d->swapchain_extent.width;
        fci.height = d->swapchain_extent.height;
        fci.layers = 1;
        VK_CHECK(vkCreateFramebuffer(d->device, &fci, nullptr, &d->framebuffers[i]));
    }

    // Resize sync objects if image count changed
    size_t old_count = d->sem_available.size();
    if (img_count != old_count) {
        for (size_t i = 0; i < old_count; ++i) {
            vkDestroyFence(d->device, d->fences_in_flight[i], nullptr);
            vkDestroySemaphore(d->device, d->sem_finished[i], nullptr);
            vkDestroySemaphore(d->device, d->sem_available[i], nullptr);
        }
        d->sem_available.resize(img_count);
        d->sem_finished.resize(img_count);
        d->fences_in_flight.resize(img_count);
        VkSemaphoreCreateInfo semi{VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO};
        VkFenceCreateInfo fi{VK_STRUCTURE_TYPE_FENCE_CREATE_INFO};
        fi.flags = VK_FENCE_CREATE_SIGNALED_BIT;
        for (uint32_t i = 0; i < img_count; ++i) {
            VK_CHECK(vkCreateSemaphore(d->device, &semi, nullptr, &d->sem_available[i]));
            VK_CHECK(vkCreateSemaphore(d->device, &semi, nullptr, &d->sem_finished[i]));
            VK_CHECK(vkCreateFence(d->device, &fi, nullptr, &d->fences_in_flight[i]));
        }
    }

    d->current_frame = 0;
}

} // namespace spatium::viewer
