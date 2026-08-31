// Single partition for spatium.render: 6 leaf headers, no mutual includes
// beyond sky.hpp -> color.hpp (handled by include order below, same as
// physics.cppm's mesh/primitives.hpp -> atom_model.hpp ordering).
//
// write_image.hpp deliberately stays HEADER-ONLY, not folded in here. Its
// whole contract is "exactly one translation unit defines
// STB_IMAGE_WRITE_IMPLEMENTATION before including it" (see the header's own
// doc comment) -- every example relies on that. Compiling it into this
// module would bake stb_image_write's externally-linked function bodies
// into spatium_render_module itself; any TU that also defines the macro
// (every example does, and so does tests/test_render.cpp) and links against
// this module in the same binary would then hit duplicate-symbol errors at
// link time. Not a module-vs-header gap to close -- module consumers keep
// `#include`ing write_image.hpp directly, same as eigen_interop.hpp and
// mesh/dec.hpp/primitives.hpp stay header-only for their own documented
// reasons.
module;
// std::vector<int> (Sky::star_buckets is vector<vector<int>>) needs a real
// header here, not just `import std.compat` -- GCC's modules-ts pipeline
// leaves implicit special-member instantiations of std::vector<int>
// unreachable-at-link-time when the only use in the whole program is inside
// a module TU relying purely on the imported std declarations (same class
// of GCC/libstdc++-modules limitation as vector.cppm's <cassert> note).
#include <vector>
export module spatium.render;
import std.compat;
import spatium.core;
import spatium.algebra;
#define SPATIUM_BUILDING_MODULE 1
#include <spatium/render/camera.hpp>
#include <spatium/render/color.hpp>
#include <spatium/render/spectral.hpp>
#include <spatium/render/parallel_for_rows.hpp>
#include <spatium/render/supersample.hpp>
#include <spatium/render/sky.hpp>
