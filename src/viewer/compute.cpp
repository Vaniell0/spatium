#include <spatium/viewer/compute.hpp>

#include <shaderc/shaderc.hpp>

#define GLFW_INCLUDE_VULKAN
#include <GLFW/glfw3.h>

#if SPATIUM_HAS_IMGUI
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_vulkan.h>
#endif

#include <algorithm>
#include <cstdlib>
#include <stdexcept>
#include <vector>

namespace spatium::viewer::compute {

namespace {

void check(VkResult r, const char* what) {
    if (r != VK_SUCCESS)
        throw std::runtime_error(std::string("Vulkan: ") + what + " failed (" +
                                 std::to_string(static_cast<int>(r)) + ")");
}

std::vector<std::uint32_t> compile(const char* src, const char* tag) {
    shaderc::Compiler compiler;
    shaderc::CompileOptions opts;
    opts.SetTargetEnvironment(shaderc_target_env_vulkan, shaderc_env_version_vulkan_1_3);
    opts.SetOptimizationLevel(shaderc_optimization_level_performance);
    auto res = compiler.CompileGlslToSpv(src, shaderc_compute_shader, tag, opts);
    if (res.GetCompilationStatus() != shaderc_compilation_status_success)
        throw std::runtime_error(std::string("shader ") + tag + ": " + res.GetErrorMessage());
    return {res.cbegin(), res.cend()};
}

int rank(VkPhysicalDeviceType t) {
    switch (t) {
        case VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU:   return 3;
        case VK_PHYSICAL_DEVICE_TYPE_INTEGRATED_GPU: return 2;
        case VK_PHYSICAL_DEVICE_TYPE_VIRTUAL_GPU:    return 1;
        default:                                     return 0;
    }
}

}  // namespace

Context::Context(const char* app_name) { init_(app_name, nullptr); }

Context::Context(const char* app_name, GLFWwindow* window) { init_(app_name, window); }

void Context::init_(const char* app_name, GLFWwindow* window) {
    VkApplicationInfo ai{VK_STRUCTURE_TYPE_APPLICATION_INFO};
    ai.pApplicationName = app_name;
    ai.apiVersion = VK_API_VERSION_1_3;
    VkInstanceCreateInfo ci{VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO};
    ci.pApplicationInfo = &ai;
    std::uint32_t n_ext = 0;
    const char** exts = window ? glfwGetRequiredInstanceExtensions(&n_ext) : nullptr;
    ci.enabledExtensionCount = n_ext;
    ci.ppEnabledExtensionNames = exts;
    check(vkCreateInstance(&ci, nullptr, &instance_), "vkCreateInstance");
    if (window) check(glfwCreateWindowSurface(instance_, window, nullptr, &surface_),
                      "glfwCreateWindowSurface");

    std::uint32_t count = 0;
    vkEnumeratePhysicalDevices(instance_, &count, nullptr);
    std::vector<VkPhysicalDevice> devs(count);
    vkEnumeratePhysicalDevices(instance_, &count, devs.data());

    const bool allow_cpu = std::getenv("SPATIUM_VK_ALLOW_CPU") != nullptr;
    // With a window the one queue also draws ImGui and presents, so it
    // must be a graphics queue the surface accepts, not just a compute one.
    const VkQueueFlags need = window ? (VK_QUEUE_COMPUTE_BIT | VK_QUEUE_GRAPHICS_BIT)
                                     : VK_QUEUE_COMPUTE_BIT;
    int best = -1;
    for (auto pd : devs) {
        VkPhysicalDeviceProperties p;
        vkGetPhysicalDeviceProperties(pd, &p);
        int r = rank(p.deviceType);
        if (r == 0 && !allow_cpu) continue;
        std::uint32_t qn = 0;
        vkGetPhysicalDeviceQueueFamilyProperties(pd, &qn, nullptr);
        std::vector<VkQueueFamilyProperties> qs(qn);
        vkGetPhysicalDeviceQueueFamilyProperties(pd, &qn, qs.data());
        for (std::uint32_t i = 0; i < qn; ++i) {
            if ((qs[i].queueFlags & need) != need || qs[i].timestampValidBits == 0) continue;
            if (surface_) {
                VkBool32 present = VK_FALSE;
                vkGetPhysicalDeviceSurfaceSupportKHR(pd, i, surface_, &present);
                if (!present) continue;
            }
            if (r > best) {
                best = r;
                physical_ = pd;
                family_ = i;
                name_ = p.deviceName;
                ns_per_tick_ = p.limits.timestampPeriod;
            }
            break;
        }
    }
    if (!physical_)
        throw std::runtime_error("no GPU with a suitable queue (a CPU implementation is refused "
                                 "unless SPATIUM_VK_ALLOW_CPU is set)");

    float prio = 1.0f;
    VkDeviceQueueCreateInfo qi{VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO};
    qi.queueFamilyIndex = family_;
    qi.queueCount = 1;
    qi.pQueuePriorities = &prio;
    const char* dev_exts[] = {VK_KHR_SWAPCHAIN_EXTENSION_NAME};
    VkDeviceCreateInfo di{VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO};
    di.queueCreateInfoCount = 1;
    di.pQueueCreateInfos = &qi;
    di.enabledExtensionCount = window ? 1 : 0;
    di.ppEnabledExtensionNames = dev_exts;
    check(vkCreateDevice(physical_, &di, nullptr, &device_), "vkCreateDevice");
    vkGetDeviceQueue(device_, family_, 0, &queue_);

    VkCommandPoolCreateInfo pi{VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO};
    pi.flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;
    pi.queueFamilyIndex = family_;
    check(vkCreateCommandPool(device_, &pi, nullptr, &pool_), "vkCreateCommandPool");
    VkCommandBufferAllocateInfo ca{VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO};
    ca.commandPool = pool_;
    ca.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
    ca.commandBufferCount = 1;
    check(vkAllocateCommandBuffers(device_, &ca, &cmd_), "vkAllocateCommandBuffers");
    VkFenceCreateInfo fi{VK_STRUCTURE_TYPE_FENCE_CREATE_INFO};
    check(vkCreateFence(device_, &fi, nullptr, &fence_), "vkCreateFence");
    VkQueryPoolCreateInfo qp{VK_STRUCTURE_TYPE_QUERY_POOL_CREATE_INFO};
    qp.queryType = VK_QUERY_TYPE_TIMESTAMP;
    qp.queryCount = 2;
    check(vkCreateQueryPool(device_, &qp, nullptr, &timestamps_), "vkCreateQueryPool");
}

Context::~Context() {
    if (device_) {
        vkDeviceWaitIdle(device_);
        vkDestroyQueryPool(device_, timestamps_, nullptr);
        vkDestroyFence(device_, fence_, nullptr);
        vkDestroyCommandPool(device_, pool_, nullptr);
        vkDestroyDevice(device_, nullptr);
    }
    if (surface_) vkDestroySurfaceKHR(instance_, surface_, nullptr);
    if (instance_) vkDestroyInstance(instance_, nullptr);
}

std::uint32_t Context::memory_type_(std::uint32_t bits) const {
    VkPhysicalDeviceMemoryProperties mp;
    vkGetPhysicalDeviceMemoryProperties(physical_, &mp);
    const VkMemoryPropertyFlags need =
        VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT;
    // Device-local as well when it exists, which on an integrated GPU it
    // does for the same heap.
    for (int pass = 0; pass < 2; ++pass) {
        const VkMemoryPropertyFlags want =
            pass == 0 ? need | VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT : need;
        for (std::uint32_t i = 0; i < mp.memoryTypeCount; ++i)
            if ((bits & (1u << i)) && (mp.memoryTypes[i].propertyFlags & want) == want) return i;
    }
    throw std::runtime_error("no host-visible coherent memory type");
}

void Context::begin_() {
    check(vkResetCommandBuffer(cmd_, 0), "vkResetCommandBuffer");
    VkCommandBufferBeginInfo bi{VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO};
    bi.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
    check(vkBeginCommandBuffer(cmd_, &bi), "vkBeginCommandBuffer");
    vkCmdResetQueryPool(cmd_, timestamps_, 0, 2);
    vkCmdWriteTimestamp(cmd_, VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT, timestamps_, 0);
}

double Context::end_and_wait_() {
    vkCmdWriteTimestamp(cmd_, VK_PIPELINE_STAGE_BOTTOM_OF_PIPE_BIT, timestamps_, 1);
    check(vkEndCommandBuffer(cmd_), "vkEndCommandBuffer");
    VkSubmitInfo si{VK_STRUCTURE_TYPE_SUBMIT_INFO};
    si.commandBufferCount = 1;
    si.pCommandBuffers = &cmd_;
    check(vkResetFences(device_, 1, &fence_), "vkResetFences");
    check(vkQueueSubmit(queue_, 1, &si, fence_), "vkQueueSubmit");
    check(vkWaitForFences(device_, 1, &fence_, VK_TRUE, UINT64_MAX), "vkWaitForFences");
    std::uint64_t ticks[2] = {};
    check(vkGetQueryPoolResults(device_, timestamps_, 0, 2, sizeof ticks, ticks,
                                sizeof(std::uint64_t),
                                VK_QUERY_RESULT_64_BIT | VK_QUERY_RESULT_WAIT_BIT),
          "vkGetQueryPoolResults");
    return static_cast<double>(ticks[1] - ticks[0]) * ns_per_tick_ * 1e-6;
}

Buffer::Buffer(Context& ctx, std::size_t bytes) : device_(ctx.device_), bytes_(bytes) {
    VkBufferCreateInfo bi{VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO};
    bi.size = bytes;
    // Transfer-source too, so a pixel buffer can be copied into a swapchain
    // image, and transfer-destination, so a buffer can be cleared by the
    // device; neither costs anything on a buffer that is never used so.
    bi.usage = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT |
               VK_BUFFER_USAGE_TRANSFER_DST_BIT;
    bi.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
    check(vkCreateBuffer(device_, &bi, nullptr, &buffer_), "vkCreateBuffer");
    VkMemoryRequirements req;
    vkGetBufferMemoryRequirements(device_, buffer_, &req);
    VkMemoryAllocateInfo ai{VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
    ai.allocationSize = req.size;
    ai.memoryTypeIndex = ctx.memory_type_(req.memoryTypeBits);
    check(vkAllocateMemory(device_, &ai, nullptr, &memory_), "vkAllocateMemory");
    check(vkBindBufferMemory(device_, buffer_, memory_, 0), "vkBindBufferMemory");
    check(vkMapMemory(device_, memory_, 0, VK_WHOLE_SIZE, 0, &mapped_), "vkMapMemory");
}

Buffer::Buffer(Buffer&& o) noexcept
    : device_(o.device_), buffer_(o.buffer_), memory_(o.memory_), mapped_(o.mapped_),
      bytes_(o.bytes_) {
    o.buffer_ = VK_NULL_HANDLE;
    o.memory_ = VK_NULL_HANDLE;
    o.mapped_ = nullptr;
}

Buffer::~Buffer() {
    if (!device_) return;
    if (memory_) {
        vkUnmapMemory(device_, memory_);
        vkFreeMemory(device_, memory_, nullptr);
    }
    if (buffer_) vkDestroyBuffer(device_, buffer_, nullptr);
}

Kernel::Kernel(Context& ctx, const char* glsl, const char* tag, std::uint32_t buffers,
               std::uint32_t push_bytes)
    : device_(ctx.device_), buffers_(buffers), push_bytes_(push_bytes) {
    std::vector<VkDescriptorSetLayoutBinding> binds(buffers);
    for (std::uint32_t i = 0; i < buffers; ++i) {
        binds[i].binding = i;
        binds[i].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        binds[i].descriptorCount = 1;
        binds[i].stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
    }
    VkDescriptorSetLayoutCreateInfo li{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
    li.bindingCount = buffers;
    li.pBindings = binds.data();
    check(vkCreateDescriptorSetLayout(device_, &li, nullptr, &set_layout_), "set layout");

    VkPushConstantRange pr{VK_SHADER_STAGE_COMPUTE_BIT, 0, push_bytes};
    VkPipelineLayoutCreateInfo pl{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
    pl.setLayoutCount = 1;
    pl.pSetLayouts = &set_layout_;
    pl.pushConstantRangeCount = push_bytes ? 1 : 0;
    pl.pPushConstantRanges = &pr;
    check(vkCreatePipelineLayout(device_, &pl, nullptr, &layout_), "pipeline layout");

    const auto spv = compile(glsl, tag);
    VkShaderModuleCreateInfo mi{VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO};
    mi.codeSize = spv.size() * sizeof(std::uint32_t);
    mi.pCode = spv.data();
    VkShaderModule mod{};
    check(vkCreateShaderModule(device_, &mi, nullptr, &mod), "shader module");
    VkComputePipelineCreateInfo cp{VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO};
    cp.stage = {VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO};
    cp.stage.stage = VK_SHADER_STAGE_COMPUTE_BIT;
    cp.stage.module = mod;
    cp.stage.pName = "main";
    cp.layout = layout_;
    const VkResult r = vkCreateComputePipelines(device_, VK_NULL_HANDLE, 1, &cp, nullptr, &pipeline_);
    vkDestroyShaderModule(device_, mod, nullptr);
    check(r, "compute pipeline");

    VkDescriptorPoolSize ps{VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, buffers};
    VkDescriptorPoolCreateInfo dp{VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO};
    dp.maxSets = 1;
    dp.poolSizeCount = 1;
    dp.pPoolSizes = &ps;
    check(vkCreateDescriptorPool(device_, &dp, nullptr, &pool_), "descriptor pool");
    VkDescriptorSetAllocateInfo da{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
    da.descriptorPool = pool_;
    da.descriptorSetCount = 1;
    da.pSetLayouts = &set_layout_;
    check(vkAllocateDescriptorSets(device_, &da, &set_), "descriptor set");
}

Kernel::~Kernel() {
    vkDestroyDescriptorPool(device_, pool_, nullptr);
    vkDestroyPipeline(device_, pipeline_, nullptr);
    vkDestroyPipelineLayout(device_, layout_, nullptr);
    vkDestroyDescriptorSetLayout(device_, set_layout_, nullptr);
}

void Kernel::bind(std::span<Buffer* const> buffers) {
    if (buffers.size() != buffers_) throw std::runtime_error("Kernel::bind: wrong buffer count");
    std::vector<VkDescriptorBufferInfo> infos(buffers_);
    std::vector<VkWriteDescriptorSet> writes(buffers_);
    for (std::uint32_t i = 0; i < buffers_; ++i) {
        infos[i] = {buffers[i]->handle(), 0, VK_WHOLE_SIZE};
        writes[i] = {VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET};
        writes[i].dstSet = set_;
        writes[i].dstBinding = i;
        writes[i].descriptorCount = 1;
        writes[i].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        writes[i].pBufferInfo = &infos[i];
    }
    vkUpdateDescriptorSets(device_, buffers_, writes.data(), 0, nullptr);
}

void Kernel::dispatch(VkCommandBuffer cmd, const void* push, std::uint32_t gx,
                      std::uint32_t gy, std::uint32_t gz) const {
    vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline_);
    vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, layout_, 0, 1, &set_, 0, nullptr);
    if (push_bytes_) vkCmdPushConstants(cmd, layout_, VK_SHADER_STAGE_COMPUTE_BIT, 0, push_bytes_, push);
    vkCmdDispatch(cmd, gx, gy, gz);
}

// ── Presenter ────────────────────────────────────────────────────

Presenter::Presenter(Context& ctx, GLFWwindow* window) : ctx_(ctx), window_(window) {
    if (!ctx_.surface_) throw std::runtime_error("Presenter needs a Context made with a window");
    VkDevice dev = ctx_.device_;

    VkCommandBufferAllocateInfo ca{VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO};
    ca.commandPool = ctx_.pool_;
    ca.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
    ca.commandBufferCount = 1;
    check(vkAllocateCommandBuffers(dev, &ca, &cmd_), "vkAllocateCommandBuffers");
    VkSemaphoreCreateInfo si{VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO};
    check(vkCreateSemaphore(dev, &si, nullptr, &acquired_), "vkCreateSemaphore");
    VkFenceCreateInfo fi{VK_STRUCTURE_TYPE_FENCE_CREATE_INFO};
    fi.flags = VK_FENCE_CREATE_SIGNALED_BIT;
    check(vkCreateFence(dev, &fi, nullptr, &fence_), "vkCreateFence");
    VkQueryPoolCreateInfo qp{VK_STRUCTURE_TYPE_QUERY_POOL_CREATE_INFO};
    qp.queryType = VK_QUERY_TYPE_TIMESTAMP;
    qp.queryCount = 2;
    check(vkCreateQueryPool(dev, &qp, nullptr, &timestamps_), "vkCreateQueryPool");

    create_swapchain_();

    // One render pass for ImGui, loading what the copy put there. Its
    // format is the swapchain's and does not change on a resize, so it is
    // made once.
    VkAttachmentDescription att{};
    att.format = format_;
    att.samples = VK_SAMPLE_COUNT_1_BIT;
    att.loadOp = VK_ATTACHMENT_LOAD_OP_LOAD;
    att.storeOp = VK_ATTACHMENT_STORE_OP_STORE;
    att.stencilLoadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE;
    att.stencilStoreOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
    att.initialLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;
    att.finalLayout = VK_IMAGE_LAYOUT_PRESENT_SRC_KHR;
    VkAttachmentReference ref{0, VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL};
    VkSubpassDescription sub{};
    sub.pipelineBindPoint = VK_PIPELINE_BIND_POINT_GRAPHICS;
    sub.colorAttachmentCount = 1;
    sub.pColorAttachments = &ref;
    VkRenderPassCreateInfo rp{VK_STRUCTURE_TYPE_RENDER_PASS_CREATE_INFO};
    rp.attachmentCount = 1;
    rp.pAttachments = &att;
    rp.subpassCount = 1;
    rp.pSubpasses = &sub;
    check(vkCreateRenderPass(dev, &rp, nullptr, &render_pass_), "vkCreateRenderPass");
    create_framebuffers_();

#if SPATIUM_HAS_IMGUI
    VkDescriptorPoolSize ps{VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 16};
    VkDescriptorPoolCreateInfo dp{VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO};
    dp.flags = VK_DESCRIPTOR_POOL_CREATE_FREE_DESCRIPTOR_SET_BIT;
    dp.maxSets = 16;
    dp.poolSizeCount = 1;
    dp.pPoolSizes = &ps;
    check(vkCreateDescriptorPool(dev, &dp, nullptr, &imgui_pool_), "imgui descriptor pool");
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForVulkan(window_, true);
    ImGui_ImplVulkan_InitInfo info{};
    info.Instance = ctx_.instance_;
    info.PhysicalDevice = ctx_.physical_;
    info.Device = dev;
    info.QueueFamily = ctx_.family_;
    info.Queue = ctx_.queue_;
    info.DescriptorPool = imgui_pool_;
    info.RenderPass = render_pass_;
    info.MinImageCount = static_cast<std::uint32_t>(images_.size());
    info.ImageCount = static_cast<std::uint32_t>(images_.size());
    info.MSAASamples = VK_SAMPLE_COUNT_1_BIT;
    ImGui_ImplVulkan_Init(&info);
    ImGui_ImplVulkan_CreateFontsTexture();
#endif
}

Presenter::~Presenter() {
    VkDevice dev = ctx_.device_;
    vkDeviceWaitIdle(dev);
#if SPATIUM_HAS_IMGUI
    ImGui_ImplVulkan_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    vkDestroyDescriptorPool(dev, imgui_pool_, nullptr);
#endif
    destroy_swapchain_();
    vkDestroyRenderPass(dev, render_pass_, nullptr);
    vkDestroyQueryPool(dev, timestamps_, nullptr);
    vkDestroyFence(dev, fence_, nullptr);
    vkDestroySemaphore(dev, acquired_, nullptr);
    vkFreeCommandBuffers(dev, ctx_.pool_, 1, &cmd_);
}

void Presenter::create_swapchain_() {
    VkPhysicalDevice pd = ctx_.physical_;
    VkSurfaceKHR surf = ctx_.surface_;
    VkSurfaceCapabilitiesKHR caps;
    vkGetPhysicalDeviceSurfaceCapabilitiesKHR(pd, surf, &caps);
    if (!(caps.supportedUsageFlags & VK_IMAGE_USAGE_TRANSFER_DST_BIT))
        throw std::runtime_error("swapchain images cannot be copied into on this surface");

    std::uint32_t nf = 0;
    vkGetPhysicalDeviceSurfaceFormatsKHR(pd, surf, &nf, nullptr);
    std::vector<VkSurfaceFormatKHR> fmts(nf);
    vkGetPhysicalDeviceSurfaceFormatsKHR(pd, surf, &nf, fmts.data());
    // UNORM, not SRGB: the kernel writes the same bytes a PNG would hold,
    // and an sRGB image would encode them a second time.
    VkSurfaceFormatKHR chosen = fmts.front();
    for (auto f : fmts)
        if (f.format == VK_FORMAT_R8G8B8A8_UNORM || f.format == VK_FORMAT_B8G8R8A8_UNORM) {
            chosen = f;
            if (f.format == VK_FORMAT_R8G8B8A8_UNORM) break;
        }
    format_ = chosen.format;
    bgra_ = format_ == VK_FORMAT_B8G8R8A8_UNORM;

    int fw = 0, fh = 0;
    glfwGetFramebufferSize(window_, &fw, &fh);
    if (caps.currentExtent.width != UINT32_MAX) {
        extent_ = caps.currentExtent;
    } else {
        extent_.width = std::clamp(static_cast<std::uint32_t>(fw), caps.minImageExtent.width,
                                   caps.maxImageExtent.width);
        extent_.height = std::clamp(static_cast<std::uint32_t>(fh), caps.minImageExtent.height,
                                    caps.maxImageExtent.height);
    }

    std::uint32_t n_images = caps.minImageCount + 1;
    if (caps.maxImageCount && n_images > caps.maxImageCount) n_images = caps.maxImageCount;
    VkSwapchainCreateInfoKHR sc{VK_STRUCTURE_TYPE_SWAPCHAIN_CREATE_INFO_KHR};
    sc.surface = surf;
    sc.minImageCount = n_images;
    sc.imageFormat = format_;
    sc.imageColorSpace = chosen.colorSpace;
    sc.imageExtent = extent_;
    sc.imageArrayLayers = 1;
    sc.imageUsage = VK_IMAGE_USAGE_TRANSFER_DST_BIT | VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT;
    sc.imageSharingMode = VK_SHARING_MODE_EXCLUSIVE;
    sc.preTransform = caps.currentTransform;
    sc.compositeAlpha = VK_COMPOSITE_ALPHA_OPAQUE_BIT_KHR;
    // FIFO is the one mode every implementation must offer, and it paces
    // the loop to the display instead of spinning the GPU on frames no one
    // sees.
    sc.presentMode = VK_PRESENT_MODE_FIFO_KHR;
    sc.clipped = VK_TRUE;
    check(vkCreateSwapchainKHR(ctx_.device_, &sc, nullptr, &swapchain_), "vkCreateSwapchainKHR");

    std::uint32_t n = 0;
    vkGetSwapchainImagesKHR(ctx_.device_, swapchain_, &n, nullptr);
    images_.resize(n);
    vkGetSwapchainImagesKHR(ctx_.device_, swapchain_, &n, images_.data());
    views_.resize(n);
    framebuffers_.assign(n, VK_NULL_HANDLE);
    finished_.resize(n);
    for (std::uint32_t i = 0; i < n; ++i) {
        VkImageViewCreateInfo vi{VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO};
        vi.image = images_[i];
        vi.viewType = VK_IMAGE_VIEW_TYPE_2D;
        vi.format = format_;
        vi.subresourceRange = {VK_IMAGE_ASPECT_COLOR_BIT, 0, 1, 0, 1};
        check(vkCreateImageView(ctx_.device_, &vi, nullptr, &views_[i]), "vkCreateImageView");
        VkSemaphoreCreateInfo si{VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO};
        check(vkCreateSemaphore(ctx_.device_, &si, nullptr, &finished_[i]), "vkCreateSemaphore");
    }
    // The render pass needs the swapchain's format, so on the first call
    // it does not exist yet and the constructor makes the framebuffers
    // once it does; on a resize it exists and they are made here.
    if (render_pass_) create_framebuffers_();
}

void Presenter::create_framebuffers_() {
    for (std::size_t i = 0; i < views_.size(); ++i) {
        VkFramebufferCreateInfo fb{VK_STRUCTURE_TYPE_FRAMEBUFFER_CREATE_INFO};
        fb.renderPass = render_pass_;
        fb.attachmentCount = 1;
        fb.pAttachments = &views_[i];
        fb.width = extent_.width;
        fb.height = extent_.height;
        fb.layers = 1;
        check(vkCreateFramebuffer(ctx_.device_, &fb, nullptr, &framebuffers_[i]),
              "vkCreateFramebuffer");
    }
}

void Presenter::destroy_swapchain_() {
    VkDevice dev = ctx_.device_;
    for (auto f : framebuffers_) if (f) vkDestroyFramebuffer(dev, f, nullptr);
    for (auto v : views_) vkDestroyImageView(dev, v, nullptr);
    for (auto s : finished_) vkDestroySemaphore(dev, s, nullptr);
    framebuffers_.clear();
    views_.clear();
    finished_.clear();
    images_.clear();
    if (swapchain_) vkDestroySwapchainKHR(dev, swapchain_, nullptr);
    swapchain_ = VK_NULL_HANDLE;
}

bool Presenter::begin_frame_() {
    VkDevice dev = ctx_.device_;
    check(vkWaitForFences(dev, 1, &fence_, VK_TRUE, UINT64_MAX), "vkWaitForFences");

    int fw = 0, fh = 0;
    glfwGetFramebufferSize(window_, &fw, &fh);
    if (fw == 0 || fh == 0) return false;   // minimised: nothing to draw into

    const VkResult r = vkAcquireNextImageKHR(dev, swapchain_, UINT64_MAX, acquired_,
                                             VK_NULL_HANDLE, &image_);
    if (r == VK_ERROR_OUT_OF_DATE_KHR) {
        vkDeviceWaitIdle(dev);
        destroy_swapchain_();
        create_swapchain_();
        return false;
    }
    if (r != VK_SUBOPTIMAL_KHR) check(r, "vkAcquireNextImageKHR");

    check(vkResetFences(dev, 1, &fence_), "vkResetFences");
    check(vkResetCommandBuffer(cmd_, 0), "vkResetCommandBuffer");
    VkCommandBufferBeginInfo bi{VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO};
    bi.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
    check(vkBeginCommandBuffer(cmd_, &bi), "vkBeginCommandBuffer");
    vkCmdResetQueryPool(cmd_, timestamps_, 0, 2);
    vkCmdWriteTimestamp(cmd_, VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT, timestamps_, 0);
    return true;
}

bool Presenter::end_frame_(Buffer& pixels, double& gpu_ms) {
    VkDevice dev = ctx_.device_;
    vkCmdWriteTimestamp(cmd_, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, timestamps_, 1);

    // The kernel's writes, visible to the copy.
    VkBufferMemoryBarrier bb{VK_STRUCTURE_TYPE_BUFFER_MEMORY_BARRIER};
    bb.srcAccessMask = VK_ACCESS_SHADER_WRITE_BIT;
    bb.dstAccessMask = VK_ACCESS_TRANSFER_READ_BIT;
    bb.srcQueueFamilyIndex = bb.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    bb.buffer = pixels.handle();
    bb.size = VK_WHOLE_SIZE;
    VkImageMemoryBarrier to_dst{VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER};
    to_dst.srcAccessMask = 0;
    to_dst.dstAccessMask = VK_ACCESS_TRANSFER_WRITE_BIT;
    to_dst.oldLayout = VK_IMAGE_LAYOUT_UNDEFINED;
    to_dst.newLayout = VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL;
    to_dst.srcQueueFamilyIndex = to_dst.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    to_dst.image = images_[image_];
    to_dst.subresourceRange = {VK_IMAGE_ASPECT_COLOR_BIT, 0, 1, 0, 1};
    vkCmdPipelineBarrier(cmd_, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT | VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT,
                         VK_PIPELINE_STAGE_TRANSFER_BIT, 0, 0, nullptr, 1, &bb, 1, &to_dst);

    VkBufferImageCopy region{};
    region.imageSubresource = {VK_IMAGE_ASPECT_COLOR_BIT, 0, 0, 1};
    region.imageExtent = {extent_.width, extent_.height, 1};
    vkCmdCopyBufferToImage(cmd_, pixels.handle(), images_[image_],
                           VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL, 1, &region);

    VkImageMemoryBarrier to_color = to_dst;
    to_color.srcAccessMask = VK_ACCESS_TRANSFER_WRITE_BIT;
    to_color.dstAccessMask = VK_ACCESS_COLOR_ATTACHMENT_READ_BIT | VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT;
    to_color.oldLayout = VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL;
    to_color.newLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;
    vkCmdPipelineBarrier(cmd_, VK_PIPELINE_STAGE_TRANSFER_BIT,
                         VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT, 0, 0, nullptr, 0, nullptr,
                         1, &to_color);

    VkRenderPassBeginInfo rb{VK_STRUCTURE_TYPE_RENDER_PASS_BEGIN_INFO};
    rb.renderPass = render_pass_;
    rb.framebuffer = framebuffers_[image_];
    rb.renderArea.extent = extent_;
    vkCmdBeginRenderPass(cmd_, &rb, VK_SUBPASS_CONTENTS_INLINE);
#if SPATIUM_HAS_IMGUI
    if (ImGui::GetDrawData()) ImGui_ImplVulkan_RenderDrawData(ImGui::GetDrawData(), cmd_);
#endif
    vkCmdEndRenderPass(cmd_);
    check(vkEndCommandBuffer(cmd_), "vkEndCommandBuffer");

    const VkPipelineStageFlags wait = VK_PIPELINE_STAGE_TRANSFER_BIT;
    VkSubmitInfo si{VK_STRUCTURE_TYPE_SUBMIT_INFO};
    si.waitSemaphoreCount = 1;
    si.pWaitSemaphores = &acquired_;
    si.pWaitDstStageMask = &wait;
    si.commandBufferCount = 1;
    si.pCommandBuffers = &cmd_;
    si.signalSemaphoreCount = 1;
    si.pSignalSemaphores = &finished_[image_];
    check(vkQueueSubmit(ctx_.queue_, 1, &si, fence_), "vkQueueSubmit");

    VkPresentInfoKHR pi{VK_STRUCTURE_TYPE_PRESENT_INFO_KHR};
    pi.waitSemaphoreCount = 1;
    pi.pWaitSemaphores = &finished_[image_];
    pi.swapchainCount = 1;
    pi.pSwapchains = &swapchain_;
    pi.pImageIndices = &image_;
    const VkResult r = vkQueuePresentKHR(ctx_.queue_, &pi);

    // Waiting here rather than at the next frame's start costs nothing --
    // one frame is in flight either way -- and it lets the trace's time be
    // read back for the frame that produced it.
    check(vkWaitForFences(dev, 1, &fence_, VK_TRUE, UINT64_MAX), "vkWaitForFences");
    std::uint64_t ticks[2] = {};
    vkGetQueryPoolResults(dev, timestamps_, 0, 2, sizeof ticks, ticks, sizeof(std::uint64_t),
                          VK_QUERY_RESULT_64_BIT | VK_QUERY_RESULT_WAIT_BIT);
    gpu_ms = static_cast<double>(ticks[1] - ticks[0]) * ctx_.ns_per_tick_ * 1e-6;

    if (r == VK_ERROR_OUT_OF_DATE_KHR || r == VK_SUBOPTIMAL_KHR) {
        vkDeviceWaitIdle(dev);
        destroy_swapchain_();
        create_swapchain_();
        return false;
    }
    check(r, "vkQueuePresentKHR");
    return true;
}

}  // namespace spatium::viewer::compute
