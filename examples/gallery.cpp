// The gallery of solved questions, as the page reads it: `site/gallery/gallery.json`.
//
// Each card is a question typed as text and what the library made of it -- an equation solved in
// closed form with the steps it was solved by, or a question with a parameter whose answer was
// synthesised and certified, or refused with the reason -- and a scene written as LaTeX. Nothing
// on a card is typed by hand: the file is this program's output, and `--check` (a ctest) fails
// when the file and the library disagree. An error is written as the power of ten above it, so
// the file does not move with the last digits of a libm.
//
//   gallery --write site/gallery/gallery.json
//   gallery --check site/gallery/gallery.json

#include <rsc/spatium/closed_form.hpp>
#include <rsc/spatium/query.hpp>
#include <rsc/spatium/synthesis.hpp>
#include <spatium/algebra/quadrature.hpp>
#include <spatium/geometry/make.hpp>
#include <spatium/geometry/transform.hpp>
#include <spatium/io/field_text.hpp>
#include <spatium/io/json.hpp>
#include <spatium/io/trace_text.hpp>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace {

namespace bd = spatium::io::build;
using spatium::io::JsonValue;
using Object = JsonValue::ObjectType;
using Array = JsonValue::ArrayType;

double decade_above(double e) { return e > 0 ? std::pow(10.0, std::ceil(std::log10(e))) : 0.0; }

JsonValue steps_json(const std::vector<bd::Step>& steps) {
    Array out;
    for (const auto& s : steps) out.push_back(Object{{"label", s.label}, {"latex", s.latex}});
    return out;
}

JsonValue closed_form_card(std::string id, std::string question) {
    const auto solved = rsc::solve_closed_form(question);
    Object card{{"id", id}, {"kind", "closed form"}, {"question", question}};
    if (!solved) {
        card.emplace_back("refused", solved.error().message);
        return card;
    }
    card.emplace_back("unknown", solved->unknown);
    card.emplace_back("steps", steps_json(solved->steps));
    return card;
}

// a question with a parameter: the answer's pieces, written out when there are few of them
JsonValue synthesis_card(std::string id, std::string question, const rsc::Answer& a, const std::string& parameter) {
    Object card{{"id", id}, {"kind", "synthesis"}, {"question", question}, {"answered", a.answered}};
    if (!a.refused.empty()) card.emplace_back("refused", a.refused);
    double worst = 0;
    for (const auto& s : a.pieces) worst = std::max(worst, s.error);
    card.emplace_back("pieces", a.pieces.size());
    card.emplace_back("error_below", decade_above(worst));
    if (!a.pieces.empty() && a.pieces.size() <= 2) {
        Array pieces;
        for (const auto& s : a.pieces)
            pieces.push_back(Object{{"lo", s.lo}, {"hi", s.hi}, {"latex", *bd::to_latex(s.expr, parameter)},
                                    {"text", *bd::to_text(s.expr, parameter)}, {"accepted", s.accepted}});
        card.emplace_back("answer", pieces);
    }
    return card;
}

JsonValue asked_card(std::string id, std::string question, double tolerance) {
    auto q = rsc::ask(question);
    if (!q) return Object{{"id", id}, {"kind", "synthesis"}, {"question", question}, {"refused", q.error().message}};
    q->tolerance = tolerance;
    return synthesis_card(std::move(id), std::move(question), rsc::answer(*q), q->parameter);
}

std::string gallery() {
    std::vector<JsonValue> cards;
    cards.push_back(closed_form_card("quadratic", "2x^2 - 3x - 5 = 0"));
    cards.push_back(closed_form_card("complex", "x^2 + 2x + 5 = 0"));
    cards.push_back(asked_card("sqrt", "solve x^2 - p = 0 for x in [0, 10], p in [0.5, 4]", 1e-10));
    cards.push_back(asked_card("power", "integrate exp(p * log(x)) for x in [0, 1], p in [0, 3]", 1e-9));
    cards.push_back(asked_card("gauss", "integrate exp(-p * x^2) for x in [0, 1], p in [0.1, 4]", 1e-6));
    cards.push_back(asked_card("refused", "solve x^2 + p = 0 for x in [0, 10], p in [0.5, 4]", 1e-8));
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
        cards.push_back(synthesis_card("prism",
                                       "the volume under the triangle (1, 0, 2), (0, 1, 4), (0, 0, p) down to z = 0, "
                                       "p in [0, 10]",
                                       a, "p"));
    }
    {
        // a scene that names its numbers: the shell's thickness is the root of an equation in them
        bd::Trace<double> t;
        const auto R = t.let("R", 2.0), r = t.let("r", 0.5);
        const auto rr = r.field(), h = bd::Field<double>::symbol("h");
        const auto hh = t.let("h", (rr + h) * (rr + h) - bd::Field<double>(2.0) * rr * rr, 0.0, 1.0);
        const auto shell = t.offset(t.torus(R, r), hh.field());
        cards.push_back(Object{{"id", "scene"}, {"kind", "scene"}, {"latex", *bd::to_latex(t, shell.index)}});
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
