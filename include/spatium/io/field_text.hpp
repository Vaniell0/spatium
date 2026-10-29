#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/io/field.hpp>
#  include <spatium/io/field_pod.hpp>
#  include <format>
#  include <string>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::io::build {

// ── A field as text ───────────────────────────────────────────────
//
// What a person reads of a field: its lowered ops written out in the usual notation, the
// inputs under names the caller gives (u and v by default), a constant in the shortest digits
// that read back to the same double. A synthesised answer is a field, and this is how it is
// read. A field that does not lower (an opaque leaf) has no text and says why.

inline Result<std::string> to_text(const Field<double>& f, std::string u = "u", std::string v = "v") {
    const auto pod = lower(f);
    if (!pod)
        return std::unexpected(pod.error());
    std::vector<std::string> text(pod->ops.size());
    for (std::size_t i = 0; i < pod->ops.size(); ++i) {
        const auto& n = pod->ops[i];
        const auto& a = text[n.a];
        const auto& b = text[n.b];
        switch (static_cast<Op>(n.code)) {
            case Op::Const: text[i] = std::format("{}", n.value); break;
            case Op::U:     text[i] = u; break;
            case Op::V:     text[i] = v; break;
            case Op::Add:   text[i] = std::format("({} + {})", a, b); break;
            case Op::Sub:   text[i] = std::format("({} - {})", a, b); break;
            case Op::Mul:   text[i] = std::format("{} * {}", a, b); break;
            case Op::Div:   text[i] = std::format("{} / {}", a, b); break;
            case Op::Min:   text[i] = std::format("min({}, {})", a, b); break;
            case Op::Max:   text[i] = std::format("max({}, {})", a, b); break;
            case Op::Less:  text[i] = std::format("less({}, {})", a, b); break;
            case Op::Sin: case Op::Cos: case Op::Sqrt: case Op::Log: case Op::Exp: {
                std::string name = op_name(static_cast<Op>(n.code));
                name[0] = static_cast<char>(name[0] - 'A' + 'a');
                text[i] = std::format("{}({})", name, a);
                break;
            }
            default:
                text[i] = std::format("{}[{}]", op_name(static_cast<Op>(n.code)), n.k);
                break;
        }
    }
    return text.empty() ? std::string("0") : text.back();
}

} // namespace spatium::io::build
