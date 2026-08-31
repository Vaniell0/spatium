// Burning Ship fractal — interactive Vulkan window.
//
// Iterates z_{n+1} = (|Re z_n| + i |Im z_n|)^2 + c on the complex plane.
// In real coordinates:
//
//     X_{n+1} = X_n^2 - Y_n^2 + A
//     Y_{n+1} = 2 |X_n Y_n| + B
//
// The escape map is rendered straight on the GPU.  Per-pixel distance
// estimation (mathr.co.uk, 2018) tracks the 2x2 Jacobian J = d(X,Y)/d(A,B)
// alongside z, accounting for the |x|' = sgn(x) discontinuities that make
// the Burning Ship non-analytic.  The estimator
//
//     d = |z|^2 * log|z| / |z * J|
//
// gives the crisp, banding-free outline that defines the iconic look of
// the fractal — escape-time alone produces a grainy image.
//
// Controls:
//   left-drag         pan
//   right-drag        not used (reserved for parameter tweaks)
//   scroll wheel      zoom around the cursor
//   double-click      zoom in 2x around the cursor
//   I / K             increase / decrease max iterations
//   R                 reset to the canonical full-frame view
//   T                 toggle DE / smooth-iteration coloring
//   F11               toggle fullscreen
//   Esc / Q           quit
//
// CLI:
//   --width N         initial window width  (default 1280)
//   --height N        initial window height (default 800)
//   --iter N          starting max iterations (default 256)
//   -h / --help       print this help
//
// The demo links GLFW + Vulkan + shaderc directly, the same set
// SPATIUM_BUILD_VIEWER pulls in.  It does not depend on spatium_viewer.

#include <spatium/algebra/complex.hpp>
#include <spatium/algebra/vector.hpp>
#include <spatium/core/error.hpp>
#include <spatium/io/table.hpp>

#include <vulkan/vulkan.h>
#include <GLFW/glfw3.h>
#include <shaderc/shaderc.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <print>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

using spatium::Complex;
using spatium::Vec;
using spatium::Result;
using spatium::Error;
using spatium::ErrorCode;

namespace {

// ── CLI ────────────────────────────────────────────────────────

struct Args {
    int width = 1280;
    int height = 800;
    int max_iter = 512;
};

void print_help() {
    std::print(
        "burning_ship_demo - interactive Vulkan renderer for the Burning Ship.\n"
        "\n"
        "Usage: burning_ship_demo [options]\n"
        "\n"
        "  --width N        initial window width  (default 1280)\n"
        "  --height N       initial window height (default 800)\n"
        "  --iter N         starting max iterations (default 512)\n"
        "  -h, --help       this message\n"
        "\n"
        "Controls (in window):\n"
        "  drag             pan\n"
        "  scroll           zoom around cursor\n"
        "  double-click     zoom in 2x around cursor\n"
        "  I / K            iter++/iter--  (step of 32)\n"
        "  R                reset to canonical view\n"
        "  T                toggle banded teal (default) <-> distance estimation\n"
        "  F11              fullscreen toggle\n"
        "  Esc, Q           quit\n");
}

bool parse_args(int argc, char** argv, Args& a) {
    auto need = [&](int& i) -> char* {
        if (i + 1 >= argc) {
            std::print(stderr, "{}: missing argument\n", argv[i]);
            return nullptr;
        }
        return argv[++i];
    };
    for (int i = 1; i < argc; ++i) {
        std::string_view s = argv[i];
        if (s == "-h" || s == "--help") { print_help(); std::exit(0); }
        else if (s == "--width")  { auto* v = need(i); if (!v) return false; a.width  = std::atoi(v); }
        else if (s == "--height") { auto* v = need(i); if (!v) return false; a.height = std::atoi(v); }
        else if (s == "--iter")   { auto* v = need(i); if (!v) return false; a.max_iter = std::atoi(v); }
        else { std::print(stderr, "unknown option: {}\n", s); return false; }
    }
    if (a.width < 64 || a.height < 64 || a.max_iter < 16) {
        std::print(stderr, "values out of range\n");
        return false;
    }
    return true;
}

// ── Shaders ────────────────────────────────────────────────────

constexpr const char* VERT_SRC = R"glsl(
#version 450
// Fullscreen triangle (Sascha Willems trick): no vertex buffer needed,
// gl_VertexIndex 0..2 maps to (-1,-1), (3,-1), (-1,3) in clip space, so
// the resulting triangle covers the [-1,1] viewport with overdraw outside.
layout(location = 0) out vec2 v_uv;
void main() {
    v_uv = vec2((gl_VertexIndex << 1) & 2, gl_VertexIndex & 2);
    gl_Position = vec4(v_uv * 2.0 - 1.0, 0.0, 1.0);
}
)glsl";

constexpr const char* FRAG_SRC = R"glsl(
#version 450
layout(location = 0) in vec2 v_uv;          // 0..1 across viewport
layout(location = 0) out vec4 outColor;

// Single-precision push constants — mainstream desktop GPUs do not always
// expose shaderFloat64 in fragment shaders (Iris Xe being a notable example),
// and float is enough for ~10^-5 zoom on the Burning Ship.  Deep zoom would
// need perturbation theory anyway, not just doubles.
layout(push_constant) uniform PC {
    vec2  center;     // 8 bytes
    float zoom;       // 4
    float aspect;     // 4
    int   max_iter;   // 4
    int   color_mode; // 4    0 = distance estimation, 1 = smooth iter
    int   win_w;      // 4
    int   win_h;      // 4
} pc;

vec3 teal_ramp(float t) {
    // Solid black -> deep teal -> cyan -> pale cyan, with a touch of yellow
    // at the very thin escape band so the antenna filaments pop out without
    // taking over the colour story.
    float k = pow(clamp(t, 0.0, 1.0), 0.55);
    vec3 dark   = vec3(0.00, 0.02, 0.03);
    vec3 deep   = vec3(0.00, 0.32, 0.36);
    vec3 cyan   = vec3(0.10, 0.85, 0.95);
    vec3 pale   = vec3(0.78, 1.00, 1.00);
    vec3 accent = vec3(1.00, 0.95, 0.30);

    vec3 c0 = mix(dark, deep, smoothstep(0.00, 0.20, k));
    vec3 c1 = mix(c0,   cyan, smoothstep(0.20, 0.70, k));
    vec3 c2 = mix(c1,   pale, smoothstep(0.70, 0.95, k));
    vec3 c3 = mix(c2, accent, smoothstep(0.95, 1.00, k));
    return c3;
}

float band_intensity(float mu) {
    // Subtle iso-iteration banding: one half-period per integer iteration,
    // gives the "topographic map" feel of the YouTube thumbnail without
    // turning the image into a barber pole.
    return 0.78 + 0.22 * cos(mu * 2.0);
}

void main() {
    // Map fragment to c-plane.  Imaginary axis is flipped so the ship sits
    // upright (canonical orientation in the literature).
    float aspect = pc.aspect;
    float cx = pc.center.x + (v_uv.x - 0.5) * 2.0 * pc.zoom;
    float cy = pc.center.y - (v_uv.y - 0.5) * 2.0 * pc.zoom * aspect;

    // Iteration with full Jacobian for distance estimation.
    // J = | dXdA dXdB |
    //     | dYdA dYdB |
    float X = 0.0, Y = 0.0;
    float dXdA = 0.0, dXdB = 0.0;
    float dYdA = 0.0, dYdB = 0.0;

    int n;
    float r2 = 0.0;
    const float escape = 1.0e6;
    for (n = 0; n < pc.max_iter; ++n) {
        r2 = X * X + Y * Y;
        if (r2 > escape) break;

        float sx = (X >= 0.0) ? 1.0 : -1.0;
        float sy = (Y >= 0.0) ? 1.0 : -1.0;

        // dz/dc update — chain rule with the |x|' = sgn(x) factor.
        // X' = X^2 - Y^2 + A   ->  dX'/dA = 2(X dXdA - Y dYdA) + 1
        // Y' = 2|X Y| + B      ->  dY'/dA = 2 sx sy (X dYdA + dXdA Y)
        float new_dXdA = 2.0 * (X * dXdA - Y * dYdA) + 1.0;
        float new_dXdB = 2.0 * (X * dXdB - Y * dYdB);
        float new_dYdA = 2.0 * sx * sy * (X * dYdA + dXdA * Y);
        float new_dYdB = 2.0 * sx * sy * (X * dYdB + dXdB * Y) + 1.0;

        // Apply abs to z, then square + offset.
        float Xa = abs(X);
        float Ya = abs(Y);
        float newX = Xa * Xa - Ya * Ya + cx;
        float newY = 2.0 * Xa * Ya + cy;

        X = newX; Y = newY;
        dXdA = new_dXdA; dXdB = new_dXdB;
        dYdA = new_dYdA; dYdB = new_dYdB;
    }

    if (n >= pc.max_iter) {
        // Solid black interior — matches the literature's "burning hull".
        outColor = vec4(0.0, 0.0, 0.0, 1.0);
        return;
    }

    // Smooth (continuous) iteration count.
    float mag2 = X * X + Y * Y;
    float mu   = float(n) + 1.0 - log2(0.5 * log(mag2));

    if (pc.color_mode == 0) {
        // Banded teal — main view.  Smooth iter through teal_ramp + iso
        // banding via cos(mu).  Filaments and the antenna spires stay
        // visible at any zoom level because the bands track absolute
        // iteration count, not normalised t.
        float t = mu / float(pc.max_iter);
        vec3 base = teal_ramp(t);
        outColor = vec4(base * band_intensity(mu), 1.0);
    } else {
        // Distance estimation — bright rim, dark interior.  Mathr 2018
        // formula: d = |z|^2 * log|z| / |z * J|.  Better than escape-time
        // for very deep zoom or print stills where banding annoys.
        float mag    = sqrt(mag2);
        float logmag = log(mag);
        float jx = X * dXdA + Y * dYdA;
        float jy = X * dXdB + Y * dYdB;
        float jnorm = sqrt(jx * jx + jy * jy);
        float d = mag2 * logmag / max(jnorm, 1e-30);
        float pix = 2.0 * pc.zoom / float(pc.win_w);
        float s = clamp(d / pix * 1.5, 0.0, 1.0);

        float t = sqrt(clamp(mu / float(pc.max_iter), 0.0, 1.0));
        vec3 base = teal_ramp(t);
        vec3 edge = vec3(0.85, 1.00, 1.00);
        outColor = vec4(mix(edge, base * 0.55, s), 1.0);
    }
}
)glsl";

// ── Burning Ship orbit (CPU-side, via spatium::Complex) ────────
//
// Mirrors the GLSL iteration so we can sanity-check the shader from C++.
// Uses spatium::Complex for the squaring step; the abs-step happens before
// the multiplication, since Complex itself is analytic and Burning Ship
// isn't.  Returns the smooth-iteration count or -1 if the orbit stays
// bounded for max_iter.
float burning_ship_orbit(Vec<float, 2> c, int max_iter) {
    Complex<float> z{0.0f, 0.0f};
    for (int n = 0; n < max_iter; ++n) {
        Complex<float> za{std::abs(z.re), std::abs(z.im)};
        z = za * za + Complex<float>{c[0], c[1]};
        const float r2 = z.magnitude_sq();
        if (r2 > 1.0e6f) {
            return float(n) + 1.0f - std::log2(0.5f * std::log(r2));
        }
    }
    return -1.0f;
}

// ── Vulkan boilerplate ─────────────────────────────────────────

[[noreturn]] void throw_vk(VkResult r, const char* what) {
    Error e{ErrorCode::NotImplemented,        // closest fit — no Vulkan-specific code
            std::string("Vulkan: ") + what
                + " (VkResult " + std::to_string(int(r)) + ")"};
    throw std::runtime_error(e.message);
}

void check(VkResult r, const char* what) {
    if (r != VK_SUCCESS) throw_vk(r, what);
}

std::vector<uint32_t> compile_glsl(const char* src, shaderc_shader_kind kind, const char* tag) {
    shaderc::Compiler compiler;
    shaderc::CompileOptions opts;
    opts.SetTargetEnvironment(shaderc_target_env_vulkan, shaderc_env_version_vulkan_1_3);
    opts.SetOptimizationLevel(shaderc_optimization_level_performance);
    auto res = compiler.CompileGlslToSpv(src, kind, tag, opts);
    if (res.GetCompilationStatus() != shaderc_compilation_status_success) {
        throw std::runtime_error(std::string("Shader compile (") + tag + "): "
            + res.GetErrorMessage());
    }
    return {res.cbegin(), res.cend()};
}

VkShaderModule make_module(VkDevice dev, const std::vector<uint32_t>& spv) {
    VkShaderModuleCreateInfo ci{VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO};
    ci.codeSize = spv.size() * sizeof(uint32_t);
    ci.pCode = spv.data();
    VkShaderModule m{};
    check(vkCreateShaderModule(dev, &ci, nullptr, &m), "vkCreateShaderModule");
    return m;
}

// ── Push constant struct (matches GLSL layout) ─────────────────
//
// `center` is laid out as `vec2` in the shader; on the C++ side we keep it
// as `Vec<float, 2>` so view-state arithmetic uses the spatium overloads
// (norm, +, -, dot — all expression-template paths stay in scope).
struct alignas(8) PushConsts {
    Vec<float, 2> center;
    float zoom;
    float aspect;
    int32_t max_iter;
    int32_t color_mode;
    int32_t win_w;
    int32_t win_h;
};
static_assert(sizeof(PushConsts) == 32, "push constant size mismatch");

// ── Application state ──────────────────────────────────────────

struct App {
    GLFWwindow* window = nullptr;
    int width = 0, height = 0;
    bool resized = false;
    bool fullscreen = false;
    int windowed_x = 0, windowed_y = 0, windowed_w = 0, windowed_h = 0;
    bool should_quit = false;

    VkInstance       instance = VK_NULL_HANDLE;
    VkSurfaceKHR     surface  = VK_NULL_HANDLE;
    VkPhysicalDevice physical = VK_NULL_HANDLE;
    VkDevice         device   = VK_NULL_HANDLE;
    uint32_t         qf_index = 0;
    VkQueue          queue    = VK_NULL_HANDLE;

    VkSwapchainKHR swapchain = VK_NULL_HANDLE;
    VkFormat       sc_format = VK_FORMAT_B8G8R8A8_UNORM;
    VkExtent2D     sc_extent{};
    std::vector<VkImage>       sc_images;
    std::vector<VkImageView>   sc_views;
    std::vector<VkFramebuffer> sc_fbs;

    VkRenderPass     render_pass     = VK_NULL_HANDLE;
    VkPipelineLayout pipeline_layout = VK_NULL_HANDLE;
    VkPipeline       pipeline        = VK_NULL_HANDLE;
    VkShaderModule   vert_mod        = VK_NULL_HANDLE;
    VkShaderModule   frag_mod        = VK_NULL_HANDLE;

    VkCommandPool                cmd_pool = VK_NULL_HANDLE;
    std::vector<VkCommandBuffer> cmd_buffers;
    static constexpr int FRAMES_IN_FLIGHT = 2;
    std::array<VkSemaphore, FRAMES_IN_FLIGHT> sem_acquire{};
    std::array<VkSemaphore, FRAMES_IN_FLIGHT> sem_present{};
    std::array<VkFence,     FRAMES_IN_FLIGHT> fence_inflight{};
    int frame_index = 0;

    // View state
    PushConsts pc{
        .center = Vec<float, 2>{-0.4f, -0.5f},
        .zoom = 1.75f,
        .aspect = 0.625f,
        .max_iter = 256,
        .color_mode = 0,
        .win_w = 0,
        .win_h = 0,
    };

    // Mouse state
    bool dragging = false;
    double drag_x = 0, drag_y = 0;
    double last_click_t = 0;
};

// Forward decl of callback wrappers that read App from window user pointer.
void cb_resize(GLFWwindow* w, int width, int height);
void cb_key(GLFWwindow* w, int key, int sc, int action, int mods);
void cb_mouse(GLFWwindow* w, int button, int action, int mods);
void cb_cursor(GLFWwindow* w, double x, double y);
void cb_scroll(GLFWwindow* w, double xoff, double yoff);

// ── Init: GLFW + window ────────────────────────────────────────

void init_window(App& a, int w, int h) {
    if (!glfwInit()) throw std::runtime_error("glfwInit failed");
    glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API);
    glfwWindowHint(GLFW_RESIZABLE, GLFW_TRUE);
    a.window = glfwCreateWindow(w, h, "burning_ship_demo  -  drag/scroll to explore", nullptr, nullptr);
    if (!a.window) { glfwTerminate(); throw std::runtime_error("glfwCreateWindow failed"); }
    glfwSetWindowUserPointer(a.window, &a);
    glfwSetFramebufferSizeCallback(a.window, cb_resize);
    glfwSetKeyCallback(a.window, cb_key);
    glfwSetMouseButtonCallback(a.window, cb_mouse);
    glfwSetCursorPosCallback(a.window, cb_cursor);
    glfwSetScrollCallback(a.window, cb_scroll);
    glfwGetFramebufferSize(a.window, &a.width, &a.height);
}

// ── Init: instance / surface / device ──────────────────────────

void init_instance(App& a) {
    VkApplicationInfo ai{VK_STRUCTURE_TYPE_APPLICATION_INFO};
    ai.pApplicationName = "burning_ship_demo";
    ai.apiVersion = VK_API_VERSION_1_3;

    uint32_t glfw_count = 0;
    const char** glfw_exts = glfwGetRequiredInstanceExtensions(&glfw_count);
    std::vector<const char*> exts(glfw_exts, glfw_exts + glfw_count);

    VkInstanceCreateInfo ci{VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO};
    ci.pApplicationInfo = &ai;
    ci.enabledExtensionCount = uint32_t(exts.size());
    ci.ppEnabledExtensionNames = exts.data();
    check(vkCreateInstance(&ci, nullptr, &a.instance), "vkCreateInstance");
}

void init_surface(App& a) {
    check(glfwCreateWindowSurface(a.instance, a.window, nullptr, &a.surface),
          "glfwCreateWindowSurface");
}

void init_device(App& a) {
    uint32_t count = 0;
    vkEnumeratePhysicalDevices(a.instance, &count, nullptr);
    if (!count) throw std::runtime_error("no Vulkan physical devices");
    std::vector<VkPhysicalDevice> pdevs(count);
    vkEnumeratePhysicalDevices(a.instance, &count, pdevs.data());

    // Pick first device that supports the surface and a graphics queue.
    for (auto pd : pdevs) {
        uint32_t qcount = 0;
        vkGetPhysicalDeviceQueueFamilyProperties(pd, &qcount, nullptr);
        std::vector<VkQueueFamilyProperties> qprops(qcount);
        vkGetPhysicalDeviceQueueFamilyProperties(pd, &qcount, qprops.data());
        for (uint32_t i = 0; i < qcount; ++i) {
            VkBool32 present = VK_FALSE;
            vkGetPhysicalDeviceSurfaceSupportKHR(pd, i, a.surface, &present);
            if ((qprops[i].queueFlags & VK_QUEUE_GRAPHICS_BIT) && present) {
                a.physical = pd;
                a.qf_index = i;
                break;
            }
        }
        if (a.physical) break;
    }
    if (!a.physical) throw std::runtime_error("no suitable GPU/queue");

    float prio = 1.0f;
    VkDeviceQueueCreateInfo qi{VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO};
    qi.queueFamilyIndex = a.qf_index;
    qi.queueCount = 1;
    qi.pQueuePriorities = &prio;

    const char* devExts[] = {VK_KHR_SWAPCHAIN_EXTENSION_NAME};
    VkPhysicalDeviceFeatures features{};        // no special features needed

    VkDeviceCreateInfo di{VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO};
    di.queueCreateInfoCount = 1;
    di.pQueueCreateInfos = &qi;
    di.enabledExtensionCount = 1;
    di.ppEnabledExtensionNames = devExts;
    di.pEnabledFeatures = &features;
    check(vkCreateDevice(a.physical, &di, nullptr, &a.device), "vkCreateDevice");
    vkGetDeviceQueue(a.device, a.qf_index, 0, &a.queue);
}

// ── Swapchain ──────────────────────────────────────────────────

void create_swapchain(App& a) {
    VkSurfaceCapabilitiesKHR caps{};
    vkGetPhysicalDeviceSurfaceCapabilitiesKHR(a.physical, a.surface, &caps);

    uint32_t fcount = 0;
    vkGetPhysicalDeviceSurfaceFormatsKHR(a.physical, a.surface, &fcount, nullptr);
    std::vector<VkSurfaceFormatKHR> formats(fcount);
    vkGetPhysicalDeviceSurfaceFormatsKHR(a.physical, a.surface, &fcount, formats.data());
    a.sc_format = formats[0].format;
    VkColorSpaceKHR colspace = formats[0].colorSpace;
    for (auto& f : formats) {
        if (f.format == VK_FORMAT_B8G8R8A8_UNORM && f.colorSpace == VK_COLOR_SPACE_SRGB_NONLINEAR_KHR) {
            a.sc_format = f.format;
            colspace = f.colorSpace;
            break;
        }
    }

    a.sc_extent.width  = std::clamp<uint32_t>(uint32_t(a.width),  caps.minImageExtent.width,  caps.maxImageExtent.width);
    a.sc_extent.height = std::clamp<uint32_t>(uint32_t(a.height), caps.minImageExtent.height, caps.maxImageExtent.height);

    uint32_t img_count = caps.minImageCount + 1;
    if (caps.maxImageCount && img_count > caps.maxImageCount) img_count = caps.maxImageCount;

    VkSwapchainCreateInfoKHR ci{VK_STRUCTURE_TYPE_SWAPCHAIN_CREATE_INFO_KHR};
    ci.surface = a.surface;
    ci.minImageCount = img_count;
    ci.imageFormat = a.sc_format;
    ci.imageColorSpace = colspace;
    ci.imageExtent = a.sc_extent;
    ci.imageArrayLayers = 1;
    ci.imageUsage = VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT;
    ci.imageSharingMode = VK_SHARING_MODE_EXCLUSIVE;
    ci.preTransform = caps.currentTransform;
    ci.compositeAlpha = VK_COMPOSITE_ALPHA_OPAQUE_BIT_KHR;
    ci.presentMode = VK_PRESENT_MODE_FIFO_KHR;   // vsync, always supported
    ci.clipped = VK_TRUE;
    check(vkCreateSwapchainKHR(a.device, &ci, nullptr, &a.swapchain), "vkCreateSwapchainKHR");

    uint32_t real_count = 0;
    vkGetSwapchainImagesKHR(a.device, a.swapchain, &real_count, nullptr);
    a.sc_images.resize(real_count);
    vkGetSwapchainImagesKHR(a.device, a.swapchain, &real_count, a.sc_images.data());

    a.sc_views.resize(real_count);
    for (uint32_t i = 0; i < real_count; ++i) {
        VkImageViewCreateInfo vi{VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO};
        vi.image = a.sc_images[i];
        vi.viewType = VK_IMAGE_VIEW_TYPE_2D;
        vi.format = a.sc_format;
        vi.components = {VK_COMPONENT_SWIZZLE_IDENTITY, VK_COMPONENT_SWIZZLE_IDENTITY,
                         VK_COMPONENT_SWIZZLE_IDENTITY, VK_COMPONENT_SWIZZLE_IDENTITY};
        vi.subresourceRange = {VK_IMAGE_ASPECT_COLOR_BIT, 0, 1, 0, 1};
        check(vkCreateImageView(a.device, &vi, nullptr, &a.sc_views[i]), "vkCreateImageView");
    }
}

// ── Render pass ────────────────────────────────────────────────

void create_render_pass(App& a) {
    VkAttachmentDescription att{};
    att.format = a.sc_format;
    att.samples = VK_SAMPLE_COUNT_1_BIT;
    att.loadOp  = VK_ATTACHMENT_LOAD_OP_DONT_CARE;
    att.storeOp = VK_ATTACHMENT_STORE_OP_STORE;
    att.stencilLoadOp  = VK_ATTACHMENT_LOAD_OP_DONT_CARE;
    att.stencilStoreOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
    att.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
    att.finalLayout   = VK_IMAGE_LAYOUT_PRESENT_SRC_KHR;

    VkAttachmentReference ref{};
    ref.attachment = 0;
    ref.layout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;

    VkSubpassDescription sub{};
    sub.pipelineBindPoint = VK_PIPELINE_BIND_POINT_GRAPHICS;
    sub.colorAttachmentCount = 1;
    sub.pColorAttachments = &ref;

    VkSubpassDependency dep{};
    dep.srcSubpass = VK_SUBPASS_EXTERNAL;
    dep.dstSubpass = 0;
    dep.srcStageMask = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT;
    dep.dstStageMask = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT;
    dep.srcAccessMask = 0;
    dep.dstAccessMask = VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT;

    VkRenderPassCreateInfo rpi{VK_STRUCTURE_TYPE_RENDER_PASS_CREATE_INFO};
    rpi.attachmentCount = 1; rpi.pAttachments = &att;
    rpi.subpassCount = 1;    rpi.pSubpasses   = &sub;
    rpi.dependencyCount = 1; rpi.pDependencies = &dep;
    check(vkCreateRenderPass(a.device, &rpi, nullptr, &a.render_pass), "vkCreateRenderPass");
}

void create_framebuffers(App& a) {
    a.sc_fbs.resize(a.sc_views.size());
    for (size_t i = 0; i < a.sc_views.size(); ++i) {
        VkFramebufferCreateInfo fi{VK_STRUCTURE_TYPE_FRAMEBUFFER_CREATE_INFO};
        fi.renderPass = a.render_pass;
        fi.attachmentCount = 1;
        fi.pAttachments = &a.sc_views[i];
        fi.width  = a.sc_extent.width;
        fi.height = a.sc_extent.height;
        fi.layers = 1;
        check(vkCreateFramebuffer(a.device, &fi, nullptr, &a.sc_fbs[i]), "vkCreateFramebuffer");
    }
}

// ── Pipeline ───────────────────────────────────────────────────

void create_pipeline(App& a) {
    auto vspv = compile_glsl(VERT_SRC, shaderc_vertex_shader,   "burning_ship.vert");
    auto fspv = compile_glsl(FRAG_SRC, shaderc_fragment_shader, "burning_ship.frag");
    a.vert_mod = make_module(a.device, vspv);
    a.frag_mod = make_module(a.device, fspv);

    VkPipelineShaderStageCreateInfo stages[2]{};
    stages[0].sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
    stages[0].stage = VK_SHADER_STAGE_VERTEX_BIT;
    stages[0].module = a.vert_mod;
    stages[0].pName = "main";
    stages[1].sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
    stages[1].stage = VK_SHADER_STAGE_FRAGMENT_BIT;
    stages[1].module = a.frag_mod;
    stages[1].pName = "main";

    VkPipelineVertexInputStateCreateInfo vi{VK_STRUCTURE_TYPE_PIPELINE_VERTEX_INPUT_STATE_CREATE_INFO};
    VkPipelineInputAssemblyStateCreateInfo ia{VK_STRUCTURE_TYPE_PIPELINE_INPUT_ASSEMBLY_STATE_CREATE_INFO};
    ia.topology = VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST;

    VkPipelineViewportStateCreateInfo vp{VK_STRUCTURE_TYPE_PIPELINE_VIEWPORT_STATE_CREATE_INFO};
    vp.viewportCount = 1; vp.scissorCount = 1;

    VkPipelineRasterizationStateCreateInfo rs{VK_STRUCTURE_TYPE_PIPELINE_RASTERIZATION_STATE_CREATE_INFO};
    rs.polygonMode = VK_POLYGON_MODE_FILL;
    rs.cullMode = VK_CULL_MODE_NONE;
    rs.frontFace = VK_FRONT_FACE_COUNTER_CLOCKWISE;
    rs.lineWidth = 1.0f;

    VkPipelineMultisampleStateCreateInfo ms{VK_STRUCTURE_TYPE_PIPELINE_MULTISAMPLE_STATE_CREATE_INFO};
    ms.rasterizationSamples = VK_SAMPLE_COUNT_1_BIT;

    VkPipelineColorBlendAttachmentState cba{};
    cba.colorWriteMask = VK_COLOR_COMPONENT_R_BIT | VK_COLOR_COMPONENT_G_BIT
                       | VK_COLOR_COMPONENT_B_BIT | VK_COLOR_COMPONENT_A_BIT;
    cba.blendEnable = VK_FALSE;

    VkPipelineColorBlendStateCreateInfo cb{VK_STRUCTURE_TYPE_PIPELINE_COLOR_BLEND_STATE_CREATE_INFO};
    cb.attachmentCount = 1;
    cb.pAttachments = &cba;

    VkDynamicState dyn_states[]{VK_DYNAMIC_STATE_VIEWPORT, VK_DYNAMIC_STATE_SCISSOR};
    VkPipelineDynamicStateCreateInfo dyn{VK_STRUCTURE_TYPE_PIPELINE_DYNAMIC_STATE_CREATE_INFO};
    dyn.dynamicStateCount = 2;
    dyn.pDynamicStates = dyn_states;

    VkPushConstantRange pcr{};
    pcr.stageFlags = VK_SHADER_STAGE_FRAGMENT_BIT;
    pcr.offset = 0;
    pcr.size = sizeof(PushConsts);

    VkPipelineLayoutCreateInfo plc{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
    plc.pushConstantRangeCount = 1;
    plc.pPushConstantRanges = &pcr;
    check(vkCreatePipelineLayout(a.device, &plc, nullptr, &a.pipeline_layout), "pipeline layout");

    VkGraphicsPipelineCreateInfo gp{VK_STRUCTURE_TYPE_GRAPHICS_PIPELINE_CREATE_INFO};
    gp.stageCount = 2; gp.pStages = stages;
    gp.pVertexInputState = &vi;
    gp.pInputAssemblyState = &ia;
    gp.pViewportState = &vp;
    gp.pRasterizationState = &rs;
    gp.pMultisampleState = &ms;
    gp.pColorBlendState = &cb;
    gp.pDynamicState = &dyn;
    gp.layout = a.pipeline_layout;
    gp.renderPass = a.render_pass;
    gp.subpass = 0;
    check(vkCreateGraphicsPipelines(a.device, VK_NULL_HANDLE, 1, &gp, nullptr, &a.pipeline), "graphics pipeline");
}

// ── Command buffers + sync ─────────────────────────────────────

void create_commands(App& a) {
    VkCommandPoolCreateInfo pi{VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO};
    pi.flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;
    pi.queueFamilyIndex = a.qf_index;
    check(vkCreateCommandPool(a.device, &pi, nullptr, &a.cmd_pool), "vkCreateCommandPool");

    a.cmd_buffers.resize(App::FRAMES_IN_FLIGHT);
    VkCommandBufferAllocateInfo ai{VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO};
    ai.commandPool = a.cmd_pool;
    ai.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
    ai.commandBufferCount = App::FRAMES_IN_FLIGHT;
    check(vkAllocateCommandBuffers(a.device, &ai, a.cmd_buffers.data()), "vkAllocateCommandBuffers");

    VkSemaphoreCreateInfo si{VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO};
    VkFenceCreateInfo fi{VK_STRUCTURE_TYPE_FENCE_CREATE_INFO};
    fi.flags = VK_FENCE_CREATE_SIGNALED_BIT;
    for (int i = 0; i < App::FRAMES_IN_FLIGHT; ++i) {
        check(vkCreateSemaphore(a.device, &si, nullptr, &a.sem_acquire[i]), "sem_acquire");
        check(vkCreateSemaphore(a.device, &si, nullptr, &a.sem_present[i]), "sem_present");
        check(vkCreateFence(a.device, &fi, nullptr, &a.fence_inflight[i]), "fence_inflight");
    }
}

// ── Recreate swapchain on resize ───────────────────────────────

void destroy_swapchain(App& a) {
    for (auto fb : a.sc_fbs)   vkDestroyFramebuffer(a.device, fb, nullptr);
    for (auto v  : a.sc_views) vkDestroyImageView(a.device, v, nullptr);
    a.sc_fbs.clear(); a.sc_views.clear();
    if (a.swapchain) vkDestroySwapchainKHR(a.device, a.swapchain, nullptr);
    a.swapchain = VK_NULL_HANDLE;
}

void recreate_swapchain(App& a) {
    int w = 0, h = 0;
    while (w == 0 || h == 0) {
        glfwGetFramebufferSize(a.window, &w, &h);
        if (w == 0 || h == 0) glfwWaitEvents();
    }
    a.width = w; a.height = h;
    vkDeviceWaitIdle(a.device);
    destroy_swapchain(a);
    create_swapchain(a);
    create_framebuffers(a);
    a.resized = false;
}

// ── Per-frame ──────────────────────────────────────────────────

void update_pc(App& a) {
    a.pc.win_w = a.width;
    a.pc.win_h = a.height;
    a.pc.aspect = float(a.height) / float(a.width);
    if (a.pc.max_iter < 16) a.pc.max_iter = 16;
}

void draw_frame(App& a) {
    int idx = a.frame_index;
    vkWaitForFences(a.device, 1, &a.fence_inflight[idx], VK_TRUE, UINT64_MAX);

    uint32_t image_index = 0;
    VkResult r = vkAcquireNextImageKHR(a.device, a.swapchain, UINT64_MAX,
                                       a.sem_acquire[idx], VK_NULL_HANDLE, &image_index);
    if (r == VK_ERROR_OUT_OF_DATE_KHR) { recreate_swapchain(a); return; }
    if (r != VK_SUCCESS && r != VK_SUBOPTIMAL_KHR) check(r, "vkAcquireNextImageKHR");

    vkResetFences(a.device, 1, &a.fence_inflight[idx]);

    VkCommandBuffer cb = a.cmd_buffers[idx];
    vkResetCommandBuffer(cb, 0);

    VkCommandBufferBeginInfo bi{VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO};
    bi.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
    check(vkBeginCommandBuffer(cb, &bi), "vkBeginCommandBuffer");

    VkClearValue clear{};
    clear.color = {{0.0f, 0.0f, 0.0f, 1.0f}};
    VkRenderPassBeginInfo rpb{VK_STRUCTURE_TYPE_RENDER_PASS_BEGIN_INFO};
    rpb.renderPass = a.render_pass;
    rpb.framebuffer = a.sc_fbs[image_index];
    rpb.renderArea.offset = {0, 0};
    rpb.renderArea.extent = a.sc_extent;
    rpb.clearValueCount = 1;
    rpb.pClearValues = &clear;
    vkCmdBeginRenderPass(cb, &rpb, VK_SUBPASS_CONTENTS_INLINE);

    vkCmdBindPipeline(cb, VK_PIPELINE_BIND_POINT_GRAPHICS, a.pipeline);
    VkViewport vp{0, 0, float(a.sc_extent.width), float(a.sc_extent.height), 0, 1};
    vkCmdSetViewport(cb, 0, 1, &vp);
    VkRect2D sc{{0, 0}, a.sc_extent};
    vkCmdSetScissor(cb, 0, 1, &sc);
    update_pc(a);
    vkCmdPushConstants(cb, a.pipeline_layout, VK_SHADER_STAGE_FRAGMENT_BIT,
                       0, sizeof(PushConsts), &a.pc);
    vkCmdDraw(cb, 3, 1, 0, 0);

    vkCmdEndRenderPass(cb);
    check(vkEndCommandBuffer(cb), "vkEndCommandBuffer");

    VkPipelineStageFlags wait_stage = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT;
    VkSubmitInfo si{VK_STRUCTURE_TYPE_SUBMIT_INFO};
    si.waitSemaphoreCount = 1;
    si.pWaitSemaphores    = &a.sem_acquire[idx];
    si.pWaitDstStageMask  = &wait_stage;
    si.commandBufferCount = 1;
    si.pCommandBuffers    = &cb;
    si.signalSemaphoreCount = 1;
    si.pSignalSemaphores  = &a.sem_present[idx];
    check(vkQueueSubmit(a.queue, 1, &si, a.fence_inflight[idx]), "vkQueueSubmit");

    VkPresentInfoKHR pi{VK_STRUCTURE_TYPE_PRESENT_INFO_KHR};
    pi.waitSemaphoreCount = 1;
    pi.pWaitSemaphores    = &a.sem_present[idx];
    pi.swapchainCount     = 1;
    pi.pSwapchains        = &a.swapchain;
    pi.pImageIndices      = &image_index;
    r = vkQueuePresentKHR(a.queue, &pi);
    if (r == VK_ERROR_OUT_OF_DATE_KHR || r == VK_SUBOPTIMAL_KHR || a.resized) {
        recreate_swapchain(a);
    } else if (r != VK_SUCCESS) {
        check(r, "vkQueuePresentKHR");
    }

    a.frame_index = (a.frame_index + 1) % App::FRAMES_IN_FLIGHT;
}

// ── Cleanup ────────────────────────────────────────────────────

void cleanup(App& a) {
    if (a.device) vkDeviceWaitIdle(a.device);
    for (int i = 0; i < App::FRAMES_IN_FLIGHT; ++i) {
        if (a.sem_acquire[i])    vkDestroySemaphore(a.device, a.sem_acquire[i], nullptr);
        if (a.sem_present[i])    vkDestroySemaphore(a.device, a.sem_present[i], nullptr);
        if (a.fence_inflight[i]) vkDestroyFence(a.device, a.fence_inflight[i], nullptr);
    }
    if (a.cmd_pool)        vkDestroyCommandPool(a.device, a.cmd_pool, nullptr);
    if (a.pipeline)        vkDestroyPipeline(a.device, a.pipeline, nullptr);
    if (a.pipeline_layout) vkDestroyPipelineLayout(a.device, a.pipeline_layout, nullptr);
    if (a.vert_mod)        vkDestroyShaderModule(a.device, a.vert_mod, nullptr);
    if (a.frag_mod)        vkDestroyShaderModule(a.device, a.frag_mod, nullptr);
    destroy_swapchain(a);
    if (a.render_pass) vkDestroyRenderPass(a.device, a.render_pass, nullptr);
    if (a.device)      vkDestroyDevice(a.device, nullptr);
    if (a.surface)     vkDestroySurfaceKHR(a.instance, a.surface, nullptr);
    if (a.instance)    vkDestroyInstance(a.instance, nullptr);
    if (a.window)      glfwDestroyWindow(a.window);
    glfwTerminate();
}

// ── Camera helpers (c-plane <-> screen) ────────────────────────

Vec<float, 2> cursor_to_c(const App& a, double px, double py) {
    double u = px / double(a.width);
    double v = py / double(a.height);
    double aspect = double(a.height) / double(a.width);
    return Vec<float, 2>{
        float(double(a.pc.center[0]) + (u - 0.5) * 2.0 * double(a.pc.zoom)),
        float(double(a.pc.center[1]) - (v - 0.5) * 2.0 * double(a.pc.zoom) * aspect),
    };
}

void zoom_around_cursor(App& a, double px, double py, double factor) {
    Vec<float, 2> before = cursor_to_c(a, px, py);
    double new_zoom = double(a.pc.zoom) * factor;
    if (new_zoom < 1e-6) new_zoom = 1e-6;   // float precision floor
    a.pc.zoom = float(new_zoom);
    Vec<float, 2> after = cursor_to_c(a, px, py);
    a.pc.center = a.pc.center + (before - after);
}

// ── Callbacks ──────────────────────────────────────────────────

void cb_resize(GLFWwindow* w, int width, int height) {
    auto* a = static_cast<App*>(glfwGetWindowUserPointer(w));
    a->resized = true;
    a->width = width; a->height = height;
}

void cb_key(GLFWwindow* w, int key, int /*sc*/, int action, int mods) {
    auto* a = static_cast<App*>(glfwGetWindowUserPointer(w));
    if (action != GLFW_PRESS && action != GLFW_REPEAT) return;
    switch (key) {
    case GLFW_KEY_ESCAPE:
    case GLFW_KEY_Q:
        a->should_quit = true;
        break;
    case GLFW_KEY_R:
        a->pc.center = Vec<float, 2>{-0.4f, -0.5f};
        a->pc.zoom = 1.75f;
        std::print("reset to canonical view\n");
        break;
    case GLFW_KEY_T:
        a->pc.color_mode = 1 - a->pc.color_mode;
        std::print("color mode: {}\n",
            a->pc.color_mode == 0 ? "banded teal" : "distance estimation");
        break;
    case GLFW_KEY_I:
        a->pc.max_iter = std::min(a->pc.max_iter + 32, 8192);
        std::print("max_iter = {}\n", a->pc.max_iter);
        break;
    case GLFW_KEY_K:
        a->pc.max_iter = std::max(a->pc.max_iter - 32, 16);
        std::print("max_iter = {}\n", a->pc.max_iter);
        break;
    case GLFW_KEY_F11: {
        if (!a->fullscreen) {
            glfwGetWindowPos(w, &a->windowed_x, &a->windowed_y);
            glfwGetWindowSize(w, &a->windowed_w, &a->windowed_h);
            GLFWmonitor* mon = glfwGetPrimaryMonitor();
            const GLFWvidmode* mode = glfwGetVideoMode(mon);
            glfwSetWindowMonitor(w, mon, 0, 0, mode->width, mode->height, mode->refreshRate);
            a->fullscreen = true;
        } else {
            glfwSetWindowMonitor(w, nullptr, a->windowed_x, a->windowed_y,
                                 a->windowed_w, a->windowed_h, 0);
            a->fullscreen = false;
        }
        break;
    }
    default: break;
    }
    (void)mods;
}

void cb_mouse(GLFWwindow* w, int button, int action, int /*mods*/) {
    auto* a = static_cast<App*>(glfwGetWindowUserPointer(w));
    if (button != GLFW_MOUSE_BUTTON_LEFT) return;
    if (action == GLFW_PRESS) {
        a->dragging = true;
        glfwGetCursorPos(w, &a->drag_x, &a->drag_y);
        double now = glfwGetTime();
        if (now - a->last_click_t < 0.30) {
            // Double click — zoom 2x toward cursor.
            zoom_around_cursor(*a, a->drag_x, a->drag_y, 0.5);
        }
        a->last_click_t = now;
    } else if (action == GLFW_RELEASE) {
        a->dragging = false;
    }
}

void cb_cursor(GLFWwindow* w, double x, double y) {
    auto* a = static_cast<App*>(glfwGetWindowUserPointer(w));
    if (!a->dragging) return;
    double dx = x - a->drag_x;
    double dy = y - a->drag_y;
    a->drag_x = x; a->drag_y = y;
    // Convert pixel delta to c-plane delta and pan in the opposite direction.
    double aspect = double(a->height) / double(a->width);
    double cdx = (dx / double(a->width))  * 2.0 * double(a->pc.zoom);
    double cdy = (dy / double(a->height)) * 2.0 * double(a->pc.zoom) * aspect;
    a->pc.center = a->pc.center
        + Vec<float, 2>{float(-cdx), float(cdy)};   // image y flipped
}

void cb_scroll(GLFWwindow* w, double /*xoff*/, double yoff) {
    auto* a = static_cast<App*>(glfwGetWindowUserPointer(w));
    double cx, cy;
    glfwGetCursorPos(w, &cx, &cy);
    double factor = (yoff > 0) ? std::pow(0.85, yoff)
                                : std::pow(1.0 / 0.85, -yoff);
    zoom_around_cursor(*a, cx, cy, factor);
}

}  // namespace

int main(int argc, char** argv) {
    Args args;
    if (!parse_args(argc, argv, args)) return 1;

    // CPU sanity check via spatium::Complex — make sure the C++ orbit and
    // the GLSL orbit will agree.  Picks one obviously-interior point (origin),
    // one obviously-exterior point (far outside), and one near the bow of
    // the ship that should escape after a moderate number of iterations.
    {
        spatium::io::section("Burning Ship CPU sanity check (spatium::Complex)");
        spatium::io::Table t{"point", "smooth iter (-1 = bounded)"};
        struct Sample { const char* name; Vec<float, 2> c; } samples[] = {
            {"origin (0, 0)",            {0.0f, 0.0f}},
            {"far escape (3, 3)",        {3.0f, 3.0f}},
            {"near hull (-1.76, -0.03)", {-1.76f, -0.03f}},
        };
        for (auto& s : samples) {
            float n = burning_ship_orbit(s.c, 256);
            t.row(s.name, n < 0 ? std::string{"bounded"}
                                : std::format("{:.2f}", n));
        }
        t.print();
    }

    App a;
    a.pc.max_iter = args.max_iter;

    try {
        init_window(a, args.width, args.height);
        init_instance(a);
        init_surface(a);
        init_device(a);
        create_swapchain(a);
        create_render_pass(a);
        create_framebuffers(a);
        create_pipeline(a);
        create_commands(a);

        std::print("burning_ship_demo: {}x{}, drag to pan, scroll to zoom, "
                   "T to toggle DE/iter, R to reset, Esc to quit\n",
                   a.width, a.height);

        while (!a.should_quit && !glfwWindowShouldClose(a.window)) {
            glfwPollEvents();
            draw_frame(a);
        }

        cleanup(a);
    } catch (const std::exception& e) {
        std::print(stderr, "fatal: {}\n", e.what());
        cleanup(a);
        return 1;
    }
    return 0;
}
