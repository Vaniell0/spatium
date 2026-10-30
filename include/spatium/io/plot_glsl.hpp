#pragma once
// The graph of a function as a WebGL2 fragment shader.
//
// y = f(x), with x the field's U, is drawn where F = f(x) - y is near zero: a pixel is covered by
// |F| / fwidth(F), its distance from the curve in pixels, so the line keeps one width however
// steep the curve is, with no grid of points and no derivative of f -- fwidth is the screen's own.
// Every symbol f reads is a slot (`field_symbols[k]`, named by `symbols`), so moving a parameter
// sets a uniform and the shader is compiled once per form of f.
//
// `plot_values` set to 1 writes F, x and y themselves instead of the colour, into a float
// target; that is what the shader is checked by, against `plot_value` -- the same program run in
// fp32 on the host -- at every pixel's centre.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/error.hpp>
#  include <spatium/io/field.hpp>
#  include <spatium/io/field_glsl.hpp>
#  include <spatium/io/field_pod.hpp>
#  include <format>
#  include <string>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::io::build {

struct PlotShader {
    std::string vertex, fragment;
    std::vector<std::string> symbols;   // the slots, in order: `field_symbols[k]` is `symbols[k]`
};

// A triangle that covers the viewport, made from the vertex index alone.
inline constexpr const char* plot_vertex = R"GLSL(#version 300 es
void main() {
    vec2 p = vec2(float((gl_VertexID & 1) << 2) - 1.0, float((gl_VertexID & 2) << 1) - 1.0);
    gl_Position = vec4(p, 0.0, 1.0);
}
)GLSL";

inline Result<PlotShader> plot_shader(const Field<double>& f) {
    GlslModule m = glsl_es_module();
    if (auto ok = emit_scalar_es(m, f, "plot_f"); !ok) return std::unexpected(ok.error());
    PlotShader out{plot_vertex, {}, m.symbols};
    out.fragment = "#version 300 es\nprecision highp float;\nprecision highp int;\n" + glsl_es_prelude(m) + m.code + R"GLSL(
uniform vec4 plot_view;    // x at the left edge, y at the bottom, x at the right, y at the top
uniform vec2 plot_size;    // the viewport, in pixels
uniform float plot_width;  // the line's half width, in pixels
uniform vec4 plot_ink;
uniform int plot_values;   // 1: F, x and y themselves, for the check
out vec4 plot_out;

void main() {
    vec2 t = gl_FragCoord.xy / plot_size;
    float x = plot_view.x + (plot_view.z - plot_view.x) * t.x;
    float y = plot_view.y + (plot_view.w - plot_view.y) * t.y;
    float F = plot_f(FieldIn(x, 0.0)) - y;
    if (plot_values == 1) { plot_out = vec4(F, x, y, 1.0); return; }
    float d = abs(F) / max(fwidth(F), 1e-30);
    plot_out = vec4(plot_ink.rgb, plot_ink.a * clamp(plot_width + 0.5 - d, 0.0, 1.0));
}
)GLSL";
    return out;
}

// F at one pixel's centre, as the shader computes it, in fp32 on the host: the reference the
// shader's `plot_values` output is held against. `view` and `size` are the uniforms, `slots` the
// values of the shader's symbols.
struct PlotValue {
    float F, x, y;
};

inline PlotValue plot_value(const PodField<float>& f, const float view[4], const float size[2], int px, int py,
                            const std::vector<float>& slots) {
    const float tx = (static_cast<float>(px) + 0.5f) / size[0], ty = (static_cast<float>(py) + 0.5f) / size[1];
    const float x = view[0] + (view[2] - view[0]) * tx, y = view[1] + (view[3] - view[1]) * ty;
    return {interpret(f, x, 0.0f, slots) - y, x, y};
}

} // namespace spatium::io::build
