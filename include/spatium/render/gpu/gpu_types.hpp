#pragma once
// The fp32 std430 layouts a compute shader reads a scene through, on their
// own so that both the packing of a scene and the building of a tree over
// moving instances can include them without including each other.
//
// Each struct is written so that the GLSL struct of the same name has the
// same offsets with no explicit padding: a `vec3` followed by a 4-byte
// scalar packs into one 16-byte slot, and a `vec4` is used everywhere else.
// Sizes are asserted rather than trusted.
//
// Not part of the `spatium.render` module, for the reason
// `render/cooked_scene.hpp` is not.
#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <cstdint>
#  include <cstring>
#endif

namespace spatium::render::gpu {

// GLSL: struct Node { vec3 lo; uint first; vec3 hi; uint count; };
// Same meaning as `BVH::Node`: a leaf has count > 0 and `first` indexes the
// primitive array; an internal node has count == 0, its left child is the
// next node and `first` is the right child.
struct Node {
    float lo[3];
    std::uint32_t first;
    float hi[3];
    std::uint32_t count;
};
static_assert(sizeof(Node) == 32);

// GLSL: struct Triangle { vec4 v0, v1, v2, n0, n1, n2, color_rough, emissive; };
// Vertices and normals in world space; `.w` unused except where named.
struct Triangle {
    float v0[4], v1[4], v2[4];
    float n0[4], n1[4], n2[4];
    float color_rough[4];   // rgb, roughness
    float emissive[4];      // rgb, unused
};
static_assert(sizeof(Triangle) == 128);

// GLSL: struct Quadric { mat4 q; vec4 lo; vec4 hi; };  (q column-major)
// One per shape, shared by every instance of it -- the whole point of an
// instance.
struct Quadric {
    float q[16];      // column-major, so a GLSL mat4 reads it unchanged
    float lo[4];      // clip box; lo.w = 1 when the box closes it into a solid
    float hi[4];
};
static_assert(sizeof(Quadric) == 96);

// GLSL: struct Instance { vec4 r0, r1, r2; vec4 scale_quadric; vec4 color_rough; vec4 emissive_opacity; };
// Rotation rows with the translation in `.w`, so a row is one dot product
// away from a world coordinate. The quadric index travels as the bits of a
// float (`floatBitsToUint` on the device) to keep the struct one layout.
struct Instance {
    float r0[4], r1[4], r2[4];   // rows of R; .w = translation
    float scale_quadric[4];      // scale, quadric index (as bits), unused, unused
    float color_rough[4];
    float emissive_opacity[4];
};
static_assert(sizeof(Instance) == 96);

// GLSL: struct LNode { vec3 lo; uint left; vec3 hi; uint right; };
// Internal: `left` and `right` are node indices. Leaf: `left` is the
// instance index with the top bit set, `right` is unused.
struct LNode {
    float lo[3];
    std::uint32_t left;
    float hi[3];
    std::uint32_t right;
};
static_assert(sizeof(LNode) == 32);

inline constexpr std::uint32_t kLeafBit = 0x80000000u;

namespace detail {

inline std::uint32_t bits_of(float f) {
    std::uint32_t u;
    std::memcpy(&u, &f, sizeof u);
    return u;
}
inline float float_of(std::uint32_t u) {
    float f;
    std::memcpy(&f, &u, sizeof f);
    return f;
}

}  // namespace detail

}  // namespace spatium::render::gpu
