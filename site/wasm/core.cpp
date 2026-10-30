// The core the page calls, built to WebAssembly: what a person types goes in as text, and the
// library's answer comes out as the JSON the page reads -- the same door (`rsc::respond`) the
// gallery is written through, so a live answer and a card cannot differ. Nothing is composed on
// the page.
//
//   respond(text)      the answer: evaluated, solved in closed form with its steps, or synthesised
//   evaluate(request)  a definition of a document given its value: the document as text or JSON,
//                      the name it asks for and the parameters it chooses, one JSON object
//   preview(text)      the question as LaTeX while it is typed, or the parse error and its column
//
// Built by build.sh with Emscripten; check.mjs holds it against the host build of the same calls.

#include <rsc/spatium/respond.hpp>
#include <spatium/io/field_parse.hpp>
#include <spatium/io/field_text.hpp>
#include <spatium/io/json.hpp>

#include <string>
#include <string_view>

namespace {

// the returned text lives until the next call, which is what the page's binding reads before
std::string last;

const char* hand_back(std::string s) {
    last = std::move(s);
    return last.c_str();
}

// The question as it is typed, as LaTeX: the expression, or each side of an equation, or the
// part before " for " of a question with a parameter. A parse error says its column, which the
// page underlines.
spatium::io::JsonValue preview(std::string_view text) {
    namespace bd = spatium::io::build;
    using Object = spatium::io::JsonValue::ObjectType;
    // the leading word is blanked rather than cut, so the column an error names is the column
    // in the text as typed
    std::string body(text);
    for (std::string_view word : {"solve ", "integrate "})
        if (body.starts_with(word)) body.replace(0, word.size(), word.size(), ' ');
    if (const auto at = body.find(" for "); at != std::string::npos) body.resize(at);
    bd::FieldReader r(body, "", "");
    auto lhs = r.expression();
    if (!lhs) return Object{{"error", lhs.error().message}};
    std::string latex = *bd::to_latex(*lhs);
    if (r.take('=')) {
        if (r.done()) return Object{{"latex", latex + " ="}};
        auto rhs = r.expression();
        if (!rhs) return Object{{"error", rhs.error().message}};
        latex += " = " + *bd::to_latex(*rhs);
    }
    if (!r.done()) return Object{{"error", r.fail("an operator or the end").error().message}};
    return Object{{"latex", latex}};
}

} // namespace

extern "C" {

const char* spatium_respond(const char* text, double tolerance) {
    return hand_back(rsc::respond(text, tolerance > 0 ? tolerance : 1e-8).dump());
}

const char* spatium_evaluate(const char* request) { return hand_back(rsc::evaluate_request(request).dump()); }

const char* spatium_preview(const char* text) { return hand_back(preview(text).dump()); }

}

// The host build of the same calls, for the check: one question a line on stdin, one answer a
// line on stdout, in the order they came.
#ifndef __EMSCRIPTEN__
#include <iostream>
int main() {
    for (std::string line; std::getline(std::cin, line);) {
        if (line.starts_with("preview ")) std::cout << spatium_preview(line.c_str() + 8) << "\n";
        else if (line.starts_with("evaluate ")) std::cout << spatium_evaluate(line.c_str() + 9) << "\n";
        else std::cout << spatium_respond(line.c_str(), 0) << "\n";
    }
}
#endif
