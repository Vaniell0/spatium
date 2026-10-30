#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/error.hpp>
#  include <spatium/io/build.hpp>
#  include <spatium/io/document.hpp>
#  include <cmath>
#  include <cstdint>
#  include <format>
#  include <stdexcept>
#  include <string>
#  include <string_view>
#  include <utility>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::io::build {

// ── A document's shapes, built as a scene ────────────────────────
//
// The same `Trace` a program builds by calls, built from a document's lines in their order: a
// value is a number the scene names (`let`), a parameter the number chosen for it (`choose`), a
// root the root of its equation (`let` with a bracket), and a shape the call its operation names,
// its numbers evaluated with every name above bound and its fields handed over as written, so the
// scene says them as the document did. An integral is not a number a scene takes, and is refused.
// A builder's refusal comes back naming the line's definition.

struct DocumentScene {
    Trace<double> trace;
    std::vector<std::pair<std::string, std::size_t>> shapes;   // each shape's name and its node

    // the last shape, which is what a document's scene is
    std::size_t root() const { return shapes.back().second; }
};

namespace detail {

inline double scene_number(const Trace<double>& trace, const Written& w, const char* where) {
    const double x = trace.bound(w.field, where)(0.0, 0.0);
    if (!std::isfinite(x)) throw std::invalid_argument(std::format("{}: '{}' is not a finite number", where, w.text));
    return x;
}

// a count or a seed: a whole number from 0 to `most`
inline std::uint64_t scene_whole(const Trace<double>& trace, const Written& w, std::uint64_t most, const char* where) {
    const double x = scene_number(trace, w, where);
    if (x < 0 || x > static_cast<double>(most) || x != std::floor(x))
        throw std::invalid_argument(std::format("{}: '{}' is {}, not a whole number from 0 to {}", where, w.text, x, most));
    return static_cast<std::uint64_t>(x);
}

} // namespace detail

// `chosen` gives each parameter its value; a parameter left out is refused, since a scene is
// drawn at one value of it.
inline Result<DocumentScene> to_trace(const Document& doc,
                                      const std::vector<std::pair<std::string, double>>& chosen = {}) {
    if (auto ok = check(doc); !ok) return std::unexpected(ok.error());
    DocumentScene out;
    Trace<double>& t = out.trace;
    const auto refuse = [](const Statement& d, std::string_view why) {
        return std::unexpected(Error(ErrorCode::InvalidArgument, std::format("{}: {}", d.name, why)));
    };
    const auto node_of = [&](const Written& w) {
        for (const auto& [name, index] : out.shapes)
            if (name == w.text) return Handle<double>{&t, index};
        return Handle<double>{&t, 0};   // not reached: check() has held every shape argument a shape above
    };
    for (const Statement& d : doc.statements) {
        try {
            switch (d.kind) {
                case StatementKind::Value: t.let(d.name, detail::scene_number(t, d.value, "value")); break;
                case StatementKind::Parameter: {
                    const auto it = std::ranges::find(chosen, d.name, &std::pair<std::string, double>::first);
                    if (it == chosen.end())
                        return refuse(d, std::format("a parameter on [{}, {}], drawn at a value it is not given", d.lo.text, d.hi.text));
                    t.choose(d.name, it->second, detail::scene_number(t, d.lo, "lo"), detail::scene_number(t, d.hi, "hi"));
                    break;
                }
                case StatementKind::Root: {
                    if (!has_bracket(d)) return refuse(d, "a root with no bracket, which a scene cannot find; give it one, on [a, b]");
                    // an equation written `= 0` is handed over as its left side, so the scene says it as written
                    const bool zero = d.rhs.field.free_symbols().empty() && d.rhs.field(0.0, 0.0) == 0.0;
                    t.let(d.name, zero ? d.value.field : d.value.field - d.rhs.field, detail::scene_number(t, d.lo, "lo"),
                          detail::scene_number(t, d.hi, "hi"));
                    break;
                }
                case StatementKind::Integral: return refuse(d, "an integral, which is not a number a scene takes");
                case StatementKind::Shape: {
                    const auto& a = d.args;
                    const auto n = [&](std::size_t k) { return detail::scene_number(t, a[k], d.operation.c_str()); };
                    Handle<double> h{&t, 0};
                    if (d.operation == "torus") h = t.torus(n(0), n(1));
                    else if (d.operation == "cylinder") h = t.cylinder(n(0), n(1));
                    else if (d.operation == "sphere") h = t.sphere(n(0));
                    else if (d.operation == "flake") h = a.size() > 3 ? t.flake({n(0), n(1), n(2)}, n(3)) : t.flake({n(0), n(1), n(2)});
                    else if (d.operation == "offset") h = t.offset(node_of(a[0]), a[1].field);
                    else if (d.operation == "scatter") {
                        const auto count = detail::scene_whole(t, a[2], 1'000'000'000, "scatter");
                        h = a.size() > 3 ? t.scatter(node_of(a[0]), node_of(a[1]), count,
                                                     static_cast<std::uint32_t>(detail::scene_whole(t, a[3], 0xffff'ffff, "scatter")))
                                         : t.scatter(node_of(a[0]), node_of(a[1]), count);
                    } else if (d.operation == "union") {
                        std::vector<Handle<double>> children;
                        for (const auto& w : a) children.push_back(node_of(w));
                        h = t.compose(children);
                    } else {
                        return refuse(d, std::format("'{}' is not an operation a shape is made by", d.operation));
                    }
                    out.shapes.emplace_back(d.name, h.index);
                    break;
                }
            }
        } catch (const std::invalid_argument& e) {
            return refuse(d, e.what());
        }
    }
    if (out.shapes.empty()) return std::unexpected(Error(ErrorCode::InvalidArgument, "the document has no shape"));
    return out;
}

} // namespace spatium::io::build
