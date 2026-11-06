// The core the page calls, built to WebAssembly: what a person types goes in as text, and the
// library's answer comes out as the JSON the page reads -- the same door (`spatium::answer::respond`) the
// gallery is written through, so a live answer and a card cannot differ. Nothing is composed on
// the page.
//
//   respond(text)      the answer: evaluated, solved, or synthesised; several lines are a
//                      document and its question in the last
//   rows(request)      the same answer for the rows of a LaTeX editor, each a formula as MathJSON
//   evaluate(request)  a definition of a document given its value: the document as text or JSON,
//                      the name it asks for and the parameters it chooses, one JSON object
//   plot(request)      a value drawn as the graph of a function of one parameter: a WebGL2
//                      shader whose slots are the other parameters, and with a grid, F at each
//                      pixel as the shader's check mode computes it, in fp32
//   zeta(request)      the zeta page's data: Z(t), the zeros proved, Gram points, zeta on a
//                      grid, and the notes, from the proved zeta (rsc/spatium/zeta_request.hpp)
//
// Built by build.sh with Emscripten; check.mjs holds it against the host build of the same calls.

#include <spatium/answer/respond.hpp>
#include <rsc/spatium/zeta_request.hpp>
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

} // namespace

extern "C" {

const char* spatium_respond(const char* text, double tolerance) {
    return hand_back(spatium::answer::respond(text, {.tolerance = tolerance > 0 ? tolerance : 1e-8}).dump());
}

const char* spatium_rows(const char* request) { return hand_back(spatium::answer::respond_rows(request).dump()); }

const char* spatium_evaluate(const char* request) { return hand_back(spatium::answer::evaluate_request(request).dump()); }

const char* spatium_plot(const char* request) { return hand_back(spatium::answer::plot_request(request).dump()); }

const char* spatium_zeta(const char* request) { return hand_back(rsc::zeta_request(request).dump()); }


}

// The host build of the same calls, for the check: one question a line on stdin, one answer a
// line on stdout, in the order they came.
#ifndef __EMSCRIPTEN__
#include <iostream>
int main() {
    for (std::string line; std::getline(std::cin, line);) {
        if (line.starts_with("evaluate ")) std::cout << spatium_evaluate(line.c_str() + 9) << "\n";
        else if (line.starts_with("rows ")) std::cout << spatium_rows(line.c_str() + 5) << "\n";
        else if (line.starts_with("plot ")) std::cout << spatium_plot(line.c_str() + 5) << "\n";
        else if (line.starts_with("zeta ")) std::cout << spatium_zeta(line.c_str() + 5) << "\n";
        else std::cout << spatium_respond(line.c_str(), 0) << "\n";
    }
}
#endif
