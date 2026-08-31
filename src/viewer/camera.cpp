#include <spatium/viewer/camera.hpp>
#include <cmath>
#include <algorithm>

namespace spatium::viewer {

Mat4f mat4_identity() {
    Mat4f m{};
    m[0] = m[5] = m[10] = m[15] = 1.0f;
    return m;
}

Mat4f mat4_multiply(const Mat4f& a, const Mat4f& b) {
    Mat4f r{};
    for (int col = 0; col < 4; ++col)
        for (int row = 0; row < 4; ++row)
            for (int k = 0; k < 4; ++k)
                r[col * 4 + row] += a[k * 4 + row] * b[col * 4 + k];
    return r;
}

Mat4f mat4_perspective(float fov_rad, float aspect, float near, float far) {
    float f = 1.0f / std::tan(fov_rad * 0.5f);
    Mat4f m{};
    m[0] = f / aspect;
    m[5] = -f;  // Vulkan Y-flip
    m[10] = far / (near - far);
    m[11] = -1.0f;
    m[14] = (near * far) / (near - far);
    return m;
}

Mat4f mat4_look_at(Vec3f eye, Vec3f center, Vec3f up) {
    auto f = (center - eye).normalized();
    auto s = f.cross(up).normalized();
    auto u = s.cross(f);

    Mat4f m{};
    m[0] = s[0]; m[4] = s[1]; m[8]  = s[2];
    m[1] = u[0]; m[5] = u[1]; m[9]  = u[2];
    m[2] = -f[0]; m[6] = -f[1]; m[10] = -f[2];
    m[12] = -s.dot(eye);
    m[13] = -u.dot(eye);
    m[14] = f.dot(eye);
    m[15] = 1.0f;
    return m;
}

Vec3f Camera::eye() const {
    float cp = std::cos(pitch);
    return target + Vec3f{
        cp * std::sin(yaw) * distance,
        std::sin(pitch) * distance,
        cp * std::cos(yaw) * distance
    };
}

Mat4f Camera::view() const {
    return mat4_look_at(eye(), target, Vec3f{0.0f, 1.0f, 0.0f});
}

Mat4f Camera::projection(float aspect) const {
    float fov_rad = fov * std::numbers::pi_v<float> / 180.0f;
    return mat4_perspective(fov_rad, aspect, near_plane, far_plane);
}

Mat4f Camera::view_projection(float aspect) const {
    return mat4_multiply(projection(aspect), view());
}

void Camera::orbit(float dx, float dy) {
    yaw -= dx * 0.01f;
    pitch = std::clamp(pitch + dy * 0.01f, -1.5f, 1.5f);
}

void Camera::zoom(float delta) {
    distance = std::max(0.1f, distance - delta * 0.3f);
}

void Camera::pan(float dx, float dy) {
    Vec3f right{std::cos(yaw), 0.0f, -std::sin(yaw)};
    Vec3f up{0.0f, 1.0f, 0.0f};
    target = target + right * (-dx * 0.005f * distance) + up * (dy * 0.005f * distance);
}

} // namespace spatium::viewer
