#pragma once

#include <spatium/viewer/camera.hpp>
#include <spatium/core/concepts.hpp>
#include <spatium/algebra/vector.hpp>
#include <spatium/mesh/mesh.hpp>
#include <spatium/mesh/operations.hpp>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

struct GLFWwindow;
typedef struct VkInstance_T* VkInstance;
typedef struct VkPhysicalDevice_T* VkPhysicalDevice;
typedef struct VkDevice_T* VkDevice;
typedef struct VkSurfaceKHR_T* VkSurfaceKHR;
typedef struct VkSwapchainKHR_T* VkSwapchainKHR;
typedef struct VkRenderPass_T* VkRenderPass;
typedef struct VkPipelineLayout_T* VkPipelineLayout;
typedef struct VkPipeline_T* VkPipeline;
typedef struct VkCommandPool_T* VkCommandPool;
typedef struct VkBuffer_T* VkBuffer;
typedef struct VkDeviceMemory_T* VkDeviceMemory;

namespace spatium::viewer {

struct MeshData {
    std::vector<float> vertices;   // interleaved: x, y, z, nx, ny, nz
    std::vector<uint32_t> indices;
    uint32_t vertex_count = 0;
    uint32_t index_count = 0;
};

// Convert any Mesh<S> to MeshData for rendering.
// Computes proper vertex normals from face cross products (area-weighted).
template<Surface S>
MeshData mesh_to_render_data(const spatium::mesh::Mesh<S>& m) {
    MeshData data;
    data.vertex_count = static_cast<uint32_t>(m.vertices.size());
    data.index_count = static_cast<uint32_t>(m.faces.size() * 3);
    data.vertices.reserve(m.vertices.size() * 6);

    auto normals = spatium::mesh::compute_vertex_normals(m);

    for (std::size_t i = 0; i < m.vertices.size(); ++i) {
        auto& v = m.vertices[i];
        data.vertices.push_back(static_cast<float>(v[0]));
        data.vertices.push_back(static_cast<float>(v[1]));
        data.vertices.push_back((v.data.size() >= 3) ? static_cast<float>(v[2]) : 0.0f);
        data.vertices.push_back(static_cast<float>(normals[i][0]));
        data.vertices.push_back(static_cast<float>(normals[i][1]));
        data.vertices.push_back((normals[i].data.size() >= 3) ? static_cast<float>(normals[i][2]) : 0.0f);
    }

    for (const auto& [a, b, c] : m.faces) {
        data.indices.push_back(a);
        data.indices.push_back(b);
        data.indices.push_back(c);
    }

    return data;
}

// With custom projection morphism — normals computed from projected positions
template<Surface S, typename Proj>
MeshData mesh_to_render_data(const spatium::mesh::Mesh<S>& m, Proj&& proj) {
    // Project vertices first, then compute normals from projected mesh
    using P = Vec<double, 3>;
    std::vector<P> projected;
    projected.reserve(m.vertices.size());
    for (const auto& v : m.vertices)
        projected.push_back(proj(v));

    // Compute normals from projected positions + original faces
    std::vector<P> normals(projected.size());
    for (const auto& [a, b, c] : m.faces) {
        auto e1 = projected[b] - projected[a];
        auto e2 = projected[c] - projected[a];
        P fn{e1.cross(e2)};
        normals[a] = P{normals[a] + fn};
        normals[b] = P{normals[b] + fn};
        normals[c] = P{normals[c] + fn};
    }
    for (auto& n : normals) n = n.normalized();

    MeshData data;
    data.vertex_count = static_cast<uint32_t>(m.vertices.size());
    data.index_count = static_cast<uint32_t>(m.faces.size() * 3);
    data.vertices.reserve(m.vertices.size() * 6);

    for (std::size_t i = 0; i < projected.size(); ++i) {
        data.vertices.push_back(static_cast<float>(projected[i][0]));
        data.vertices.push_back(static_cast<float>(projected[i][1]));
        data.vertices.push_back(static_cast<float>(projected[i][2]));
        data.vertices.push_back(static_cast<float>(normals[i][0]));
        data.vertices.push_back(static_cast<float>(normals[i][1]));
        data.vertices.push_back(static_cast<float>(normals[i][2]));
    }

    for (const auto& [a, b, c] : m.faces) {
        data.indices.push_back(a);
        data.indices.push_back(b);
        data.indices.push_back(c);
    }

    return data;
}

struct PointCloudData {
    std::vector<float> positions;  // x, y, z per point (3 floats/point)
    uint32_t point_count = 0;
};

struct MeshEntry {
    MeshData data;
    Vec4f color{1.0f, 1.0f, 1.0f, 1.0f};
    bool visible{true};
};

struct PointEntry {
    PointCloudData data;
    Vec4f color{1.0f, 1.0f, 1.0f, 1.0f};
    bool visible{true};
};

class App {
public:
    App(const std::string& title, int width = 1280, int height = 720);
    ~App();

    App(const App&) = delete;
    App& operator=(const App&) = delete;

    void add_mesh(MeshData data, Vec4f color = {1.0f, 1.0f, 1.0f, 1.0f});
    void clear_meshes();

    template<Surface S>
    void add_mesh(const spatium::mesh::Mesh<S>& m, Vec4f color = {1.0f, 1.0f, 1.0f, 1.0f}) {
        add_mesh(mesh_to_render_data(m), color);
    }

    template<Surface S, typename Proj>
    void add_mesh(const spatium::mesh::Mesh<S>& m, Proj&& proj, Vec4f color = {1.0f, 1.0f, 1.0f, 1.0f}) {
        add_mesh(mesh_to_render_data(m, std::forward<Proj>(proj)), color);
    }

    void set_mesh_visible(std::size_t index, bool visible);
    std::size_t mesh_count() const { return meshes_.size(); }

    void add_point_cloud(PointCloudData data, Vec4f color = {1.0f, 1.0f, 1.0f, 1.0f});
    void clear_point_clouds();
    void set_point_visible(std::size_t index, bool visible);
    void set_point_color(std::size_t index, Vec4f color);
    std::size_t point_cloud_count() const { return points_.size(); }

    void update_mesh(std::size_t index, MeshData data);
    void update_mesh_vertices(std::size_t index, const MeshData& data);

    using KeyCallback = std::function<void(int key, int action, int mods)>;
    void set_key_callback(KeyCallback cb) { user_key_cb_ = std::move(cb); }

    using FrameCallback = std::function<void()>;
    void set_frame_callback(FrameCallback cb) { frame_cb_ = std::move(cb); }

    using GuiCallback = std::function<void()>;
    void enable_imgui();
    void set_gui_callback(GuiCallback cb) { gui_cb_ = std::move(cb); }

    void fit_camera(float radius);

    // Fit camera to mesh bounding sphere
    template<Surface S>
    void fit_to(const spatium::mesh::Mesh<S>& m, float padding = 1.2f) {
        float max_r2 = 0;
        for (const auto& v : m.vertices) {
            float r2 = 0;
            for (std::size_t d = 0; d < v.data.size() && d < 3; ++d)
                r2 += static_cast<float>(v[d]) * static_cast<float>(v[d]);
            if (r2 > max_r2) max_r2 = r2;
        }
        fit_camera(std::sqrt(max_r2) * padding);
    }
    void save_screenshot(const std::string& path);
    void run();

    Camera camera;
    float point_size_ = 2.0f;
    bool cutaway_enabled_ = false;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;

    int width_, height_;
    bool wireframe_ = true;
    bool imgui_enabled_ = false;
    bool buffers_dirty_ = false;
    std::vector<MeshEntry> meshes_;
    std::vector<PointEntry> points_;
    KeyCallback user_key_cb_;
    FrameCallback frame_cb_;
    GuiCallback gui_cb_;

    void init_vulkan();
    void create_pipeline();
    void create_buffers();
    void draw_frame();
    void cleanup();
    void init_imgui();
    void shutdown_imgui();

    static void key_callback(GLFWwindow* w, int key, int scancode, int action, int mods);
    static void scroll_callback(GLFWwindow* w, double xoff, double yoff);
    static void cursor_callback(GLFWwindow* w, double x, double y);
    static void mouse_button_callback(GLFWwindow* w, int button, int action, int mods);
    static void framebuffer_resize_callback(GLFWwindow* w, int width, int height);

    void recreate_swapchain();
};

} // namespace spatium::viewer
