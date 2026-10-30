// The gallery of solved questions, as the page reads it: `site/gallery/gallery.json`.
//
// Each card is a question typed as text and what the library made of it, through `respond`, the
// door a live page uses too -- an expression evaluated, an equation solved in closed form with
// the steps it was solved by, or a question with a parameter whose answer was synthesised and
// certified, or refused with the reason -- and a scene written as LaTeX. Nothing
// on a card is typed by hand: the file is this program's output, and `--check` (a ctest) fails
// when the file and the library disagree. An error is written as the power of ten above it, so
// the file does not move with the last digits of a libm.
//
//   gallery --write site/gallery/gallery.json
//   gallery --check site/gallery/gallery.json

#include <rsc/spatium/respond.hpp>
#include <rsc/spatium/synthesis.hpp>
#include <spatium/algebra/quadrature.hpp>
#include <spatium/geometry/make.hpp>
#include <spatium/geometry/transform.hpp>
#include <spatium/io/document_trace.hpp>
#include <spatium/io/field_text.hpp>
#include <spatium/io/json.hpp>
#include <spatium/io/trace_text.hpp>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

namespace bd = spatium::io::build;
using spatium::io::JsonValue;
using Object = JsonValue::ObjectType;
using Array = JsonValue::ArrayType;

// a card is the library's own answer to its question, with an id the page links to; at the
// tolerance a live answer has, since the question's text is all the page sends
JsonValue card(std::string id, std::string question) {
    auto out = rsc::respond(question);
    Object card{{"id", std::move(id)}};
    for (auto& [key, value] : out.as_object()) {
        // a route's cost counts the evaluations it took, and a root found to the last bit takes as
        // many as the libm's last digits lead it to; the card keeps the routes, the live answer the counts
        if (key == "candidates") {
            Array routes;
            for (const auto& c : value.as_array()) {
                Object kept;
                for (const auto& [k, v] : c.as_object())
                    if (k != "cost") kept.emplace_back(k, v);
                routes.push_back(kept);
            }
            card.emplace_back(key, routes);
            continue;
        }
        if (key != "answer") {
            card.emplace_back(key, value);
            continue;
        }
        // a Chebyshev piece's coefficients carry the last digits of the libm that computed the
        // reference; the card keeps its form and degree, the live answer keeps the numbers
        Array pieces;
        for (const auto& piece : value.as_array()) {
            Object kept;
            for (const auto& [k, v] : piece.as_object())
                if (k != "coefficients") kept.emplace_back(k, v);
            pieces.push_back(kept);
        }
        card.emplace_back(key, pieces);
    }
    return card;
}

std::string gallery() {
    std::vector<JsonValue> cards;
    cards.push_back(card("arithmetic", "2(3 + 4)^2 - 1/4"));
    cards.push_back(card("quadratic", "2x^2 - 3x - 5 = 0"));
    cards.push_back(card("complex", "x^2 + 2x + 5 = 0"));
    cards.push_back(card("sqrt", "solve x^2 - p = 0 for x in [0, 10], p in [0.5, 4]"));
    cards.push_back(card("power", "integrate exp(p * log(x)) for x in [0, 1], p in [0, 3]"));
    cards.push_back(card("gauss", "integrate exp(-p * x^2) for x in [0, 1], p in [0.1, 4]"));
    // functions known by a name the library does not have: the sine integral, the error function
    // and the Fresnel integral, each an integral up to the parameter, answered as pieces
    cards.push_back(card("si", "integrate sin(t)/t for t in [0, x], x in [0.5, 6]"));
    cards.push_back(card("erf", "integrate 2 / sqrt(pi) * exp(-t^2) for t in [0, x], x in [0.25, 3]"));
    cards.push_back(card("fresnel", "integrate cos(pi * t^2 / 2) for t in [0, x], x in [0.25, 2]"));
    cards.push_back(card("refused", "solve x^2 + p = 0 for x in [0, 10], p in [0.5, 4]"));
    {
        // the volume under a triangle down to z = 0 as a function of one vertex's height, which
        // the library computes as an integral and synthesis writes as 1/6 p + 1
        using namespace spatium;
        using namespace spatium::geometry;
        const auto volume = [](double q) {
            const Triangle3 top = tri(Vec3{1, 0, 2}, Vec3{0, 1, 4}, Vec3{0, 0, q});
            const auto drop_z = scale(1, 1, 0);
            const auto& [a, b, c] = top.vertices;
            const double area = tri(drop_z(a), drop_z(b), drop_z(c)).area();
            const auto height = [&](double u, double v) { return a[2] + u * (b[2] - a[2]) + v * (c[2] - a[2]); };
            const auto across = [&](double u) {
                return integrate([&](double v) { return height(u, v); }, Finite<double>{0, 1 - u}).value;
            };
            return 2 * area * integrate(across, Finite<double>{0, 1}).value;
        };
        const auto s = rsc::synthesize(volume, {.lo = 0, .hi = 10, .tolerance = 1e-12, .relative = false});
        rsc::Answer a;
        a.pieces = {s};
        a.answered = s.accepted;
        Object prism{{"id", "prism"},
                     {"question", "the volume under the triangle (1, 0, 2), (0, 1, 4), (0, 0, p) down to z = 0, p in [0, 10]"},
                     {"route", "synthesis"},
                     // asked in code, not typed: the page opens it as its card, not in the question line
                     {"typed", false}};
        for (auto& kv : rsc::answer_json(a, "p")) prism.push_back(kv);
        cards.push_back(prism);
    }
    {
        // a scene written as a document, which names its numbers: the shell's thickness is the root
        // of an equation in them
        const auto doc = bd::read_document("R = 2\n"
                                           "r = 0.5\n"
                                           "h: (r + h)^2 - 2 * r^2 = 0 on [0, 1]\n"
                                           "S0 = torus(R, r)\n"
                                           "S1 = offset(S0, h)\n");
        if (!doc) throw std::runtime_error(doc.error().message);
        const auto built = bd::to_trace(*doc);
        if (!built) throw std::runtime_error(built.error().message);
        const bd::DocumentScene& scene = *built;
        cards.push_back(Object{{"id", "scene"}, {"route", "scene"}, {"latex", *bd::to_latex(scene.trace, scene.root())}});
    }
    // one card per line, so a change to one card is a change to one line
    std::string out = "{\"cards\": [\n";
    for (std::size_t i = 0; i < cards.size(); ++i) out += "  " + cards[i].dump() + (i + 1 < cards.size() ? ",\n" : "\n");
    return out + "]}\n";
}

} // namespace

int main(int argc, char** argv) {
    const std::string_view mode = argc > 2 ? argv[1] : "";
    if (mode != "--write" && mode != "--check") {
        std::fprintf(stderr, "usage: gallery --write|--check <path to gallery.json>\n");
        return 2;
    }
    const std::string now = gallery();
    if (mode == "--write") {
        std::ofstream(argv[2]) << now;
        return 0;
    }
    std::ifstream in(argv[2]);
    std::stringstream was;
    was << in.rdbuf();
    if (was.str() == now) return 0;
    // the first line that differs, so the failure says which card moved
    std::istringstream a(was.str()), b(now);
    std::string la, lb;
    for (int line = 1;; ++line) {
        const bool more_a = static_cast<bool>(std::getline(a, la)), more_b = static_cast<bool>(std::getline(b, lb));
        if (!more_a && !more_b) break;
        if (la != lb) {
            std::fprintf(stderr, "gallery.json is out of date at line %d\n  file:    %s\n  library: %s\n"
                                 "regenerate: gallery --write %s\n", line, la.c_str(), lb.c_str(), argv[2]);
            return 1;
        }
        la.clear(), lb.clear();
    }
    return 1;
}
