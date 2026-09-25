#pragma once
// A MetricField written out as GLSL with its derivatives, and the geodesic
// step built on it, for a compute shader to trace rays through the same
// spacetime the host does.
//
// The metric's ten entries share one pool (metric_field.hpp); each op of
// the pool is emitted once as a value and its four partials -- forward-
// mode differentiation, one float and one vec4 per op -- so the device
// gets d g_{mu nu} / d x^kappa exactly as Dual<T> gives the host, with no
// finite differences and no symbolic step in between. Then the Christoffel
// symbols (GLSL's own inverse(mat4)), the geodesic equation and one RK4
// step, transcribed from geodesic.hpp.
//
// The specification on the host is geodesic.hpp itself run in float --
// MetricField evaluates on Dual<float> as it does on Dual<double> -- so a
// device result is checked against the same arithmetic at the same
// precision, and both against double.
//
// Not part of any module: a header of strings, like render/gpu_trace_glsl.hpp.

#include <spatium/io/field.hpp>
#include <spatium/io/field_glsl.hpp>
#include <spatium/physics/relativity/metric_field.hpp>
#include <format>
#include <string>

namespace spatium::physics::relativity {

// `void <name>(vec4 x, out float g[10], out vec4 dg[10])`: the ten entries
// in MetricField::kEntries order, and dg[e][k] = d g_e / d x^k.
inline std::string emit_metric_glsl(const MetricField<double>& metric, const std::string& name) {
    using io::build::Op;
    const auto& pool = metric.pool();
    std::string body;
    const auto f = [](double v) { return io::build::detail::glsl_float(v); };
    for (std::size_t i = 0; i < pool.size(); ++i) {
        const auto& n = pool[i];
        const auto a = std::format("v{}", n.a), da = std::format("d{}", n.a);
        const auto b = std::format("v{}", n.b), db = std::format("d{}", n.b);
        std::string v, d;
        switch (static_cast<Op>(n.code)) {
            case Op::Const: v = f(n.value); d = "vec4(0.0)"; break;
            case Op::Coord:
                v = std::format("x[{}]", n.k);
                d = std::format("vec4({}, {}, {}, {})", n.k == 0 ? "1.0" : "0.0", n.k == 1 ? "1.0" : "0.0",
                                n.k == 2 ? "1.0" : "0.0", n.k == 3 ? "1.0" : "0.0");
                break;
            case Op::Add: v = a + " + " + b; d = da + " + " + db; break;
            case Op::Sub: v = a + " - " + b; d = da + " - " + db; break;
            case Op::Mul: v = a + " * " + b; d = std::format("{} * {} + {} * {}", da, b, a, db); break;
            case Op::Div:
                v = a + " / " + b;
                d = std::format("({} * {} - {} * {}) / ({} * {})", da, b, a, db, b, b);
                break;
            case Op::Min:
                v = std::format("({0} < {1} ? {0} : {1})", a, b);
                d = std::format("({0} < {1} ? {2} : {3})", a, b, da, db);
                break;
            case Op::Max:
                v = std::format("({0} < {1} ? {1} : {0})", a, b);
                d = std::format("({0} < {1} ? {3} : {2})", a, b, da, db);
                break;
            case Op::Sin: v = "sin(" + a + ")"; d = std::format("cos({}) * {}", a, da); break;
            case Op::Cos: v = "cos(" + a + ")"; d = std::format("-sin({}) * {}", a, da); break;
            case Op::Sqrt: v = "sqrt(" + a + ")"; d = std::format("{} / (2.0 * sqrt({}))", da, a); break;
            case Op::Less: v = std::format("({} < {} ? 1.0 : 0.0)", a, b); d = "vec4(0.0)"; break;
            default: v = "0.0"; d = "vec4(0.0)"; break;   // refused by MetricField::make
        }
        body += std::format("    float v{0} = {1};\n    vec4 d{0} = {2};\n", i, v, d);
    }
    // The roots, recovered from the same place operator() reads them.
    std::string out;
    for (std::size_t e = 0; e < 10; ++e)
        out += std::format("    g[{0}] = v{1};\n    dg[{0}] = d{1};\n", e, metric.root(e));
    return std::format("void {}(vec4 x, out float g[10], out vec4 dg[10]) {{\n{}{}}}\n", name, body, out);
}

// The geodesic equation and one RK4 step, over a metric function emitted by
// emit_metric_glsl() under `metric_name`. Defines
//   void geodesic_rhs(vec4 x, vec4 u, out vec4 dx, out vec4 du)
//   void geodesic_step(inout vec4 x, inout vec4 u, float dl)
inline std::string geodesic_glsl(const std::string& metric_name) {
    return std::format(R"GLSL(
// Entry index of (i, j) in MetricField::kEntries.
int metric_entry(int i, int j) {{
    int a = min(i, j), b = max(i, j);
    return a == 0 ? b : (a == 1 ? 3 + b : (a == 2 ? 5 + b : 9));
}}

void geodesic_rhs(vec4 x, vec4 u, out vec4 dx, out vec4 du) {{
    float g[10]; vec4 dg[10];
    {0}(x, g, dg);
    mat4 gm;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) gm[i][j] = g[metric_entry(i, j)];
    mat4 gi = inverse(gm);
    // Gamma^l_{{mn}} u^m u^n = 1/2 g^{{ls}} (2 d_m g_{{sn}} - d_s g_{{mn}}) u^m u^n
    vec4 w;   // w_s = (2 d_m g_{{sn}} - d_s g_{{mn}}) u^m u^n
    for (int s = 0; s < 4; ++s) {{
        float acc = 0.0;
        for (int m = 0; m < 4; ++m)
            for (int n = 0; n < 4; ++n)
                acc += (2.0 * dg[metric_entry(s, n)][m] - dg[metric_entry(m, n)][s]) * u[m] * u[n];
        w[s] = acc;
    }}
    dx = u;
    du = -0.5 * (gi * w);
}}

void geodesic_step(inout vec4 x, inout vec4 u, float dl) {{
    vec4 k1x, k1u, k2x, k2u, k3x, k3u, k4x, k4u;
    geodesic_rhs(x, u, k1x, k1u);
    geodesic_rhs(x + 0.5 * dl * k1x, u + 0.5 * dl * k1u, k2x, k2u);
    geodesic_rhs(x + 0.5 * dl * k2x, u + 0.5 * dl * k2u, k3x, k3u);
    geodesic_rhs(x + dl * k3x, u + dl * k3u, k4x, k4u);
    x += dl / 6.0 * (k1x + 2.0 * k2x + 2.0 * k3x + k4x);
    u += dl / 6.0 * (k1u + 2.0 * k2u + 2.0 * k3u + k4u);
}}
)GLSL", metric_name);
}

} // namespace spatium::physics::relativity
