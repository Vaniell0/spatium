// A kernel on the device held against the same computation on the host: the
// compute library alone, no window. Where no GPU is found the cases are
// skipped, not passed; SPATIUM_VK_ALLOW_CPU=1 runs them on a software device.
#include <spatium/viewer/compute.hpp>

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace vc = spatium::viewer::compute;

namespace {

// C = A B, A m x k and B k x n row-major, one invocation per entry of C.
// `b_index` is how the kernel reads B[i][j]; the right one and a wrong one
// are both compiled, so the check below has a case it must fail.
std::string matmul_glsl(const char* b_index) {
    return std::string(R"(#version 450
layout(local_size_x = 16, local_size_y = 16) in;
layout(std430, binding = 0) readonly buffer A_ { float a[]; };
layout(std430, binding = 1) readonly buffer B_ { float b[]; };
layout(std430, binding = 2) writeonly buffer C_ { float c[]; };
layout(push_constant) uniform P { uint m, n, k; } p;
void main() {
    uint row = gl_GlobalInvocationID.y, col = gl_GlobalInvocationID.x;
    if (row >= p.m || col >= p.n) return;
    float s = 0.0;
    for (uint t = 0; t < p.k; ++t) s += a[row * p.k + t] * b[)") +
           b_index + R"(];
    c[row * p.n + col] = s;
}
)";
}

struct Push {
    std::uint32_t m, n, k;
};

std::optional<vc::Context> device() {
    try {
        return std::optional<vc::Context>(std::in_place, "test_compute");
    } catch (const std::runtime_error&) {
        return std::nullopt;
    }
}

// The largest |C_device - C_host| over the bound a sum of k float products may
// differ by in any order of rounding: gamma_k * sum |a||b|, gamma_k = k u / (1 - k u)
// (Higham, Accuracy and Stability, 3.1), here as a ratio -- at most 1 when the
// two agree as float arithmetic allows.
double worst_ratio(vc::Context& ctx, const char* b_index, std::uint32_t m, std::uint32_t n,
                   std::uint32_t k) {
    std::mt19937 rng(7);
    std::uniform_real_distribution<float> u(-1.0f, 1.0f);
    std::vector<float> a(std::size_t(m) * k), b(std::size_t(k) * n);
    for (auto& x : a) x = u(rng);
    for (auto& x : b) x = u(rng);

    auto A = vc::Buffer::from<float>(ctx, a);
    auto B = vc::Buffer::from<float>(ctx, b);
    vc::Buffer C(ctx, std::size_t(m) * n * sizeof(float));
    const std::string src = matmul_glsl(b_index);
    vc::Kernel kernel(ctx, src.c_str(), "matmul", 3, sizeof(Push));
    vc::Buffer* bound[] = {&A, &B, &C};
    kernel.bind(bound);
    const Push push{m, n, k};
    ctx.run([&](VkCommandBuffer cmd) { kernel.dispatch(cmd, &push, (n + 15) / 16, (m + 15) / 16); });

    const auto* c = static_cast<const float*>(C.data());
    const double unit = std::numeric_limits<float>::epsilon() / 2;
    const double gamma = k * unit / (1 - k * unit);
    double worst = 0;
    for (std::uint32_t i = 0; i < m; ++i)
        for (std::uint32_t j = 0; j < n; ++j) {
            float s = 0;
            double mag = 0;
            for (std::uint32_t t = 0; t < k; ++t) {
                s += a[i * k + t] * b[t * n + j];
                mag += std::abs(double(a[i * k + t]) * b[t * n + j]);
            }
            // Both sides round, each within gamma * mag of the exact sum.
            worst = std::max(worst, std::abs(double(c[i * n + j]) - s) / (2 * gamma * mag));
        }
    return worst;
}

}  // namespace

TEST_CASE("A matrix product on the device agrees with the host within float rounding",
          "[compute]") {
    auto ctx = device();
    if (!ctx) SKIP("no Vulkan device (set SPATIUM_VK_ALLOW_CPU=1 for a software one)");
    INFO("device: " << ctx->device_name());
    // Sizes that are not multiples of the 16 x 16 workgroup, so the edge
    // invocations that must write nothing are exercised.
    CHECK(worst_ratio(*ctx, "t * p.n + col", 33, 17, 65) <= 1.0);
    CHECK(worst_ratio(*ctx, "t * p.n + col", 1, 1, 1000) <= 1.0);
}

TEST_CASE("The agreement check fails on a kernel that reads B transposed", "[compute]") {
    auto ctx = device();
    if (!ctx) SKIP("no Vulkan device (set SPATIUM_VK_ALLOW_CPU=1 for a software one)");
    // Square, so the wrong index stays in bounds and only the values are wrong.
    CHECK(worst_ratio(*ctx, "col * p.k + t", 24, 24, 24) > 1.0);
}
