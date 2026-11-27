#pragma once
// A Vulkan call that did not succeed, as an exception naming the call.

#include <vulkan/vulkan.h>

#include <stdexcept>
#include <string>

namespace spatium::viewer::compute {

inline void check(VkResult r, const char* what) {
    if (r != VK_SUCCESS)
        throw std::runtime_error(std::string("Vulkan: ") + what + " failed (" +
                                 std::to_string(static_cast<int>(r)) + ")");
}

}  // namespace spatium::viewer::compute
