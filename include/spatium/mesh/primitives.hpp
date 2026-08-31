#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/mesh/mesh.hpp>
#  include <spatium/spaces/sphere.hpp>
#  include <spatium/spaces/euclidean.hpp>
#  include <spatium/spaces/parametric.hpp>
#  include <cmath>
#  include <numbers>
#endif

SPATIUM_EXPORT namespace spatium::mesh {

// Icosahedron inscribed in Sphere<2> — the classic starting mesh
// for spherical subdivision. 12 vertices, 20 faces.

template<Scalar T = double>
Mesh<Sphere<2, T>> icosahedron(const Sphere<2, T>& sphere = {}) {
    using Vec = Vec<T, 3>;
    Mesh<Sphere<2, T>> m;

    using std::sqrt;
    const T phi = (T{1} + sqrt(T{5})) / T{2}; // golden ratio

    auto n = [&](T x, T y, T z) -> Vec {
        return sphere.project(Vec{x, y, z});
    };

    // 12 vertices of icosahedron
    m.vertices = {
        n(-1,  phi, 0), n( 1,  phi, 0), n(-1, -phi, 0), n( 1, -phi, 0),
        n(0, -1,  phi), n(0,  1,  phi), n(0, -1, -phi), n(0,  1, -phi),
        n( phi, 0, -1), n( phi, 0,  1), n(-phi, 0, -1), n(-phi, 0,  1),
    };

    // 20 triangular faces
    m.faces = {
        {0, 11,  5}, {0,  5,  1}, {0,  1,  7}, {0,  7, 10}, {0, 10, 11},
        {1,  5,  9}, {5, 11,  4}, {11, 10,  2}, {10,  7,  6}, {7,  1,  8},
        {3,  9,  4}, {3,  4,  2}, {3,  2,  6}, {3,  6,  8}, {3,  8,  9},
        {4,  9,  5}, {2,  4, 11}, {6,  2, 10}, {8,  6,  7}, {9,  8,  1},
    };

    return m;
}

// Tetrahedron inscribed in Sphere<2> — simpler starting mesh.
// 4 vertices, 4 faces.

template<Scalar T = double>
Mesh<Sphere<2, T>> tetrahedron(const Sphere<2, T>& sphere = {}) {
    Mesh<Sphere<2, T>> m;

    using std::sqrt;
    const T r = sphere.radius;
    const T b = r / T{3};

    m.vertices = {
        sphere.project(Vec<T, 3>{0, 0, r}),
        sphere.project(Vec<T, 3>{sqrt(T{8} / T{9}) * r, 0, -b}),
        sphere.project(Vec<T, 3>{-sqrt(T{2} / T{9}) * r,  sqrt(T{2} / T{3}) * r, -b}),
        sphere.project(Vec<T, 3>{-sqrt(T{2} / T{9}) * r, -sqrt(T{2} / T{3}) * r, -b}),
    };

    m.faces = {{0, 1, 2}, {0, 2, 3}, {0, 3, 1}, {1, 3, 2}};
    return m;
}

// ── Euclidean grid (flat XY plane) ────────────────────────────

template<Scalar T = double>
Mesh<Euclidean<3, T>> grid_mesh(std::size_t nx, std::size_t ny,
                                 T width = T{1}, T height = T{1}) {
    Mesh<Euclidean<3, T>> m;
    m.vertices.reserve((nx + 1) * (ny + 1));

    for (std::size_t j = 0; j <= ny; ++j) {
        T y = static_cast<T>(j) / static_cast<T>(ny) * height;
        for (std::size_t i = 0; i <= nx; ++i) {
            T x = static_cast<T>(i) / static_cast<T>(nx) * width;
            m.vertices.push_back(Vec<T, 3>{x, y, T{0}});
        }
    }

    m.faces.reserve(nx * ny * 2);
    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            auto idx = [&](std::size_t ii, std::size_t jj) -> uint32_t {
                return static_cast<uint32_t>(jj * (nx + 1) + ii);
            };
            uint32_t a = idx(i, j), b = idx(i + 1, j);
            uint32_t c = idx(i + 1, j + 1), d = idx(i, j + 1);
            m.faces.push_back({a, b, c});
            m.faces.push_back({a, c, d});
        }
    }
    return m;
}

// ── UV Sphere (latitude/longitude) ───────────────────────────

template<Scalar T = double>
Mesh<Euclidean<3, T>> uv_sphere_mesh(std::size_t slices, std::size_t stacks,
                                      T radius = T{1}) {
    Mesh<Euclidean<3, T>> m;
    constexpr T pi = std::numbers::pi_v<T>;

    // Vertices: poles + grid
    m.vertices.push_back(Vec<T, 3>{T{0}, T{0}, radius});  // north pole

    for (std::size_t j = 1; j < stacks; ++j) {
        T theta = pi * static_cast<T>(j) / static_cast<T>(stacks);
        for (std::size_t i = 0; i < slices; ++i) {
            T phi = T{2} * pi * static_cast<T>(i) / static_cast<T>(slices);
            m.vertices.push_back(Vec<T, 3>{
                radius * std::sin(theta) * std::cos(phi),
                radius * std::sin(theta) * std::sin(phi),
                radius * std::cos(theta)
            });
        }
    }

    m.vertices.push_back(Vec<T, 3>{T{0}, T{0}, -radius}); // south pole

    // North pole cap
    for (uint32_t i = 0; i < slices; ++i)
        m.faces.push_back({0, 1 + i, 1 + (i + 1) % static_cast<uint32_t>(slices)});

    // Body
    for (std::size_t j = 0; j + 2 < stacks; ++j) {
        uint32_t base = 1 + static_cast<uint32_t>(j * slices);
        uint32_t next = base + static_cast<uint32_t>(slices);
        for (uint32_t i = 0; i < slices; ++i) {
            uint32_t ni = (i + 1) % static_cast<uint32_t>(slices);
            m.faces.push_back({base + i, next + i, next + ni});
            m.faces.push_back({base + i, next + ni, base + ni});
        }
    }

    // South pole cap
    uint32_t south = static_cast<uint32_t>(m.vertices.size() - 1);
    uint32_t last_ring = south - static_cast<uint32_t>(slices);
    for (uint32_t i = 0; i < slices; ++i)
        m.faces.push_back({south, last_ring + (i + 1) % static_cast<uint32_t>(slices), last_ring + i});

    return m;
}

// ── Box mesh (6 faces, 12 triangles) ─────────────────────────

template<Scalar T = double>
Mesh<Euclidean<3, T>> box_mesh(Vec<T, 3> half_extents = {T{0.5}, T{0.5}, T{0.5}}) {
    Mesh<Euclidean<3, T>> m;
    auto h = half_extents;

    m.vertices = {
        { -h[0], -h[1], -h[2] }, {  h[0], -h[1], -h[2] },
        {  h[0],  h[1], -h[2] }, { -h[0],  h[1], -h[2] },
        { -h[0], -h[1],  h[2] }, {  h[0], -h[1],  h[2] },
        {  h[0],  h[1],  h[2] }, { -h[0],  h[1],  h[2] },
    };

    m.faces = {
        {0,2,1}, {0,3,2},  // -Z
        {4,5,6}, {4,6,7},  // +Z
        {0,1,5}, {0,5,4},  // -Y
        {2,3,7}, {2,7,6},  // +Y
        {0,4,7}, {0,7,3},  // -X
        {1,2,6}, {1,6,5},  // +X
    };
    return m;
}

// ── Parametric surface mesh (UV grid) ────────────────────────

template<Scalar T = double>
Mesh<ParametricSurface<T>> parametric_mesh(
    const ParametricSurface<T>& surface,
    std::size_t u_steps, std::size_t v_steps)
{
    Mesh<ParametricSurface<T>> m;
    auto [u_min, u_max, v_min, v_max] = surface.domain();
    bool pu = surface.periodic_u();
    bool pv = surface.periodic_v();

    std::size_t u_verts = pu ? u_steps : u_steps + 1;
    std::size_t v_verts = pv ? v_steps : v_steps + 1;

    T du = (u_max - u_min) / static_cast<T>(u_steps);
    T dv = (v_max - v_min) / static_cast<T>(v_steps);

    // Vertices
    m.vertices.reserve(u_verts * v_verts);
    for (std::size_t j = 0; j < v_verts; ++j) {
        T v = v_min + static_cast<T>(j) * dv;
        for (std::size_t i = 0; i < u_verts; ++i) {
            T u = u_min + static_cast<T>(i) * du;
            m.vertices.push_back(surface(u, v));
        }
    }

    // Faces (two triangles per quad)
    auto idx = [&](std::size_t i, std::size_t j) -> uint32_t {
        return static_cast<uint32_t>((j % v_verts) * u_verts + (i % u_verts));
    };

    for (std::size_t j = 0; j < v_steps; ++j) {
        for (std::size_t i = 0; i < u_steps; ++i) {
            uint32_t a = idx(i, j);
            uint32_t b = idx(i + 1, j);
            uint32_t c = idx(i + 1, j + 1);
            uint32_t d = idx(i, j + 1);
            m.faces.push_back({a, b, c});
            m.faces.push_back({a, c, d});
        }
    }

    return m;
}

} // namespace spatium::mesh

