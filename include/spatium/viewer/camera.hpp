#pragma once

#include <spatium/algebra/vector.hpp>
#include <array>
#include <cmath>
#include <numbers>

namespace spatium::viewer {

// 4x4 column-major matrix (OpenGL/Vulkan convention)
using Mat4f = std::array<float, 16>;

Mat4f mat4_identity();
Mat4f mat4_multiply(const Mat4f& a, const Mat4f& b);
Mat4f mat4_perspective(float fov_rad, float aspect, float near, float far);
Mat4f mat4_look_at(Vec3f eye, Vec3f center, Vec3f up);

struct Camera {
    float distance = 3.0f;
    float yaw = 0.0f;       // radians around Y
    float pitch = 0.3f;     // radians from XZ plane
    Vec3f target{0.0f, 0.0f, 0.0f};
    float fov = 60.0f;      // degrees
    float near_plane = 0.01f;
    float far_plane = 100.0f;

    Vec3f eye() const;
    Mat4f view() const;
    Mat4f projection(float aspect) const;
    Mat4f view_projection(float aspect) const;

    void orbit(float dx, float dy);
    void zoom(float delta);
    void pan(float dx, float dy);
    void fit(float radius) {
        float half_fov = fov * 0.5f * std::numbers::pi_v<float> / 180.0f;
        distance = radius / std::sin(half_fov) * 1.2f;
        far_plane = distance * 3.0f;
        near_plane = distance * 0.001f;
    }
};

} // namespace spatium::viewer
