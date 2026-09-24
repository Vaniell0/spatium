#include <spatium/viewer/compute.hpp>

#include <shaderc/shaderc.hpp>

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

Context::Context(const char* app_name) {
    VkApplicationInfo ai{VK_STRUCTURE_TYPE_APPLICATION_INFO};
    ai.pApplicationName = app_name;
    ai.apiVersion = VK_API_VERSION_1_3;
    VkInstanceCreateInfo ci{VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO};
    ci.pApplicationInfo = &ai;
    check(vkCreateInstance(&ci, nullptr, &instance_), "vkCreateInstance");

    std::uint32_t count = 0;
    vkEnumeratePhysicalDevices(instance_, &count, nullptr);
    std::vector<VkPhysicalDevice> devs(count);
    vkEnumeratePhysicalDevices(instance_, &count, devs.data());

    const bool allow_cpu = std::getenv("SPATIUM_VK_ALLOW_CPU") != nullptr;
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
            if (!(qs[i].queueFlags & VK_QUEUE_COMPUTE_BIT) || qs[i].timestampValidBits == 0)
                continue;
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
        throw std::runtime_error("no GPU with a compute queue (a CPU implementation is refused "
                                 "unless SPATIUM_VK_ALLOW_CPU is set)");

    float prio = 1.0f;
    VkDeviceQueueCreateInfo qi{VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO};
    qi.queueFamilyIndex = family_;
    qi.queueCount = 1;
    qi.pQueuePriorities = &prio;
    VkDeviceCreateInfo di{VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO};
    di.queueCreateInfoCount = 1;
    di.pQueueCreateInfos = &qi;
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
    bi.usage = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT;
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

}  // namespace spatium::viewer::compute
