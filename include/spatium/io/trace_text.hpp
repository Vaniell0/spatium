#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/io/build.hpp>
#  include <spatium/io/field_text.hpp>
#  include <bit>
#  include <cstdint>
#  include <format>
#  include <string>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::io::build {

// ── A scene as LaTeX ──────────────────────────────────────────────
//
// What a page shows of a scene: its named numbers, then the nodes the root reaches, one aligned
// line each. A number given as it is reads r = 1; a chosen one reads p = 1 with its interval; one
// given by an equation reads as the equation was written, with its bracket and the root found on
// it, so the page says both what was asked and what it came to. A node is its operation over the nodes before it: a factory with its
// numbers, an offset with its thickness as written (names and all), a scatter, a union. A
// factory's numbers are printed as numbers, since a factory keeps its parameters as values; a
// thickness that cannot be written (an opaque leaf) is refused, like `to_latex` of a field.

namespace detail {

inline const char* factory_name(std::uint64_t id) {
    switch (id) {
        case 1: return "torus";
        case 2: return "cylinder";
        case 3: return "sphere";
        case 4: return "flake";
        default: return "shape";
    }
}

inline std::string node_name(std::size_t i) { return std::format("S_{{{}}}", i); }

} // namespace detail

inline Result<std::string> to_latex(const Trace<double>& trace, std::size_t root) {
    using detail::TextStyle;
    std::vector<std::string> lines;
    for (const Definition<double>& d : trace.definitions()) {
        const std::string name = detail::name_latex(d.name);
        const std::string value = detail::number_text(d.value, TextStyle::Latex);
        if (d.chosen) {
            lines.push_back(std::format("{} &= {},\\quad {} \\in \\left[{}, {}\\right]", name, value, name,
                                        detail::number_text(d.lo, TextStyle::Latex),
                                        detail::number_text(d.hi, TextStyle::Latex)));
            continue;
        }
        if (!d.equation) {
            lines.push_back(std::format("{} &= {}", name, value));
            continue;
        }
        const auto eq = to_latex(*d.equation);
        if (!eq) return std::unexpected(eq.error());
        lines.push_back(std::format("& {} = 0,\\quad {} \\in \\left[{}, {}\\right] \\;\\Rightarrow\\; {} = {}", *eq,
                                    name, detail::number_text(d.lo, TextStyle::Latex),
                                    detail::number_text(d.hi, TextStyle::Latex), name, value));
    }

    // the nodes the root reaches, in the order they were made
    std::vector<char> reached(trace.size(), 0);
    std::vector<std::size_t> stack{root};
    while (!stack.empty()) {
        const std::size_t i = stack.back();
        stack.pop_back();
        if (i >= trace.size() || reached[i]) continue;
        reached[i] = 1;
        const auto& n = trace.node(i);
        if (n.kind == Kind::Offset) stack.push_back(n.base);
        if (n.kind == Kind::Scatter) stack.push_back(n.item), stack.push_back(n.target);
        if (n.kind == Kind::Compose) stack.insert(stack.end(), n.children.begin(), n.children.end());
    }
    for (std::size_t i = 0; i < trace.size(); ++i) {
        if (!reached[i]) continue;
        const auto& n = trace.node(i);
        std::string rhs;
        switch (n.kind) {
            case Kind::Space: {
                if (n.recipe.empty()) { rhs = "\\operatorname{chart}"; break; }
                std::string args;
                for (std::size_t k = 1; k < n.recipe.size(); ++k)
                    args += (k > 1 ? ", " : "") +
                            detail::number_text(std::bit_cast<double>(n.recipe[k]), TextStyle::Latex);
                rhs = std::format("\\operatorname{{{}}}\\left({}\\right)", detail::factory_name(n.recipe[0]), args);
                break;
            }
            case Kind::Literal: rhs = "\\operatorname{mesh}"; break;
            case Kind::Offset: {
                const auto h = to_latex(n.thickness_written);
                if (!h) return std::unexpected(h.error());
                rhs = std::format("\\operatorname{{offset}}\\left({}, {}\\right)", detail::node_name(n.base), *h);
                break;
            }
            case Kind::Scatter:
                rhs = std::format("\\operatorname{{scatter}}\\left({}, {}, {}\\right)", detail::node_name(n.item),
                                  detail::node_name(n.target), n.count);
                break;
            case Kind::Compose:
                for (std::size_t k = 0; k < n.children.size(); ++k)
                    rhs += (k ? " \\cup " : "") + detail::node_name(n.children[k]);
                if (rhs.empty()) rhs = "\\varnothing";
                break;
        }
        if (n.transform) rhs += "\\;\\text{(moving)}";
        lines.push_back(std::format("{} &= {}", detail::node_name(i), rhs));
    }

    std::string out = "\\begin{aligned}\n";
    for (std::size_t k = 0; k < lines.size(); ++k) out += lines[k] + (k + 1 < lines.size() ? " \\\\\n" : "\n");
    return out + "\\end{aligned}";
}

} // namespace spatium::io::build
