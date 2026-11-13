#pragma once

// What the zeta page asks the core, as one JSON object, and the library's answer as JSON. Every
// number comes from zeta in `algebra/zeta.hpp`: proved, in intervals, where it is a zero, a Gram
// point or Z(t), whose drawn value is the middle of an interval narrower than a pixel; the colour
// map is the same sum in doubles, a picture.
//
//   {"what": "hardy", "t0", "t1", "n"}            Z(t) at n points: {"t": [...], "z": [...]}
//   {"what": "zeros", "t0", "t1"}                 the zeros proved: {"zeros": [{"lo", "hi", "text",
//                                                 "digits"}], "counted", "expected", "complete"}
//   {"what": "gram", "t0", "t1"}                  {"gram": [{"n", "t", "good"}]}, good null if undecided
//   {"what": "grid", "s0": [re, im], "s1": [re, im], "w", "h"}
//                                                 zeta on a w x h grid, rows from s0's im up:
//                                                 {"re": [...], "im": [...]}, null where unbounded
//   {"what": "notes", "t0", "t1"}                 what the zeros and Gram points there show, in words

#include <spatium/algebra/zeta.hpp>
#include <spatium/core/interval.hpp>
#include <spatium/io/json.hpp>

#include <algorithm>
#include <cmath>
#include <format>
#include <numbers>
#include <string>
#include <string_view>
#include <vector>

namespace spatium::answer {

namespace detail {

using spatium::io::JsonValue;
using Object = JsonValue::ObjectType;
using Array = JsonValue::ArrayType;
using I = spatium::Interval<double>;

inline JsonValue finite_or_null(double x) { return std::isfinite(x) ? JsonValue(x) : JsonValue(nullptr); }

inline JsonValue zeros_json(const spatium::CriticalZeros<double>& z) {
    Array zeros;
    for (const I& r : z.zeros) {
        const auto p = spatium::proved(r);
        zeros.push_back(Object{{"lo", r.lo}, {"hi", r.hi}, {"text", p.text}, {"digits", p.digits}});
    }
    return Object{{"zeros", std::move(zeros)},
                  {"counted", z.counted.is_empty() ? JsonValue(nullptr) : JsonValue(Array{z.counted.lo, z.counted.hi})},
                  {"expected", z.expected < 0 ? JsonValue(nullptr) : JsonValue(static_cast<double>(z.expected))},
                  {"complete", z.complete},
                  {"evaluations", z.evaluations}};
}

inline JsonValue gram_json(const std::vector<spatium::GramPoint<double>>& g) {
    Array out;
    for (const auto& p : g)
        out.push_back(Object{{"n", static_cast<double>(p.n)}, {"t", p.t.mid()},
                             {"good", p.good == spatium::Truth::Undecided ? JsonValue(nullptr) : JsonValue(p.good == spatium::Truth::Yes)}});
    return Object{{"gram", std::move(out)}};
}

inline std::string number(double x, int digits) { return std::format("{:.{}f}", x, digits); }

// Below t = 0 everything is the mirror of what is above: zeta(conj s) = conj zeta(s), so Z(-t) = Z(t),
// a zero at 1/2 - it for each at 1/2 + it, and theta odd, so g_{-n} = -g_n.
inline I flip(const I& r) { return r.is_empty() ? r : I(-r.hi, -r.lo); }

inline spatium::CriticalZeros<double> zeros_on(double t0, double t1) {
    if (t0 >= 0) return spatium::critical_zeros(t0, t1);
    const auto mirrored = [](spatium::CriticalZeros<double> z) {
        std::ranges::reverse(z.zeros);
        for (I& r : z.zeros) r = flip(r);
        z.counted = flip(z.counted);
        return z;
    };
    if (t1 <= 0) return mirrored(spatium::critical_zeros(-t1, -t0));
    auto below = mirrored(spatium::critical_zeros(0.0, -t0));
    const auto above = spatium::critical_zeros(0.0, t1);
    below.zeros.insert(below.zeros.end(), above.zeros.begin(), above.zeros.end());
    const bool counted = below.expected >= 0 && above.expected >= 0;
    below.counted = counted ? I(below.counted.lo, above.counted.hi) : I::empty();
    below.expected = counted ? below.expected + above.expected : -1;
    below.complete = below.complete && above.complete;
    below.evaluations += above.evaluations;
    return below;
}

inline std::vector<spatium::GramPoint<double>> gram_on(double t0, double t1) {
    std::vector<spatium::GramPoint<double>> out;
    if (t0 < 0)
        for (auto p : spatium::gram_points(std::max(0.0, -t1), -t0)) out.insert(out.begin(), {-p.n, flip(p.t), p.good});
    if (t1 > 0)
        for (const auto& p : spatium::gram_points(std::max(0.0, t0), t1)) out.push_back(p);
    return out;
}

// The notes are what this window's zeros and Gram points show, each with where it is.
inline JsonValue notes(double t0, double t1) {
    Array out;
    const auto note = [&](std::string title, std::string text, std::string strength, JsonValue anchor) {
        out.push_back(Object{{"title", std::move(title)}, {"text", std::move(text)}, {"strength", std::move(strength)}, {"anchor", std::move(anchor)}});
    };
    const auto z = zeros_on(t0, t1);
    std::size_t inside = 0;
    for (const I& r : z.zeros) inside += t0 <= r.mid() && r.mid() <= t1;
    if (z.complete)
        note("Every zero here is on the line",
             std::format("The {} zeros on [{}, {}] are all there are. N(T) counts {} zeros of \\(\\zeta\\) with ordinate between {} and {}, "
                         "and Z changes sign {} times there, so each is simple and has real part \\(\\tfrac12\\).",
                         inside, number(t0, 2), number(t1, 2), z.expected, number(z.counted.lo, 2), number(z.counted.hi, 2), z.expected),
             "proved here", Object{{"t0", z.counted.lo}, {"t1", z.counted.hi}});
    else
        note("Zeros found, their count not proved",
             std::format("Z changes sign {} times here. N(T) was not decided at both ends, so these are not proved to be all.", z.zeros.size()),
             "computed here", Object{{"t0", t0}, {"t1", t1}});
    const auto g = gram_on(t0, t1);
    for (const auto& p : g)
        if (p.good == spatium::Truth::No)
            note(std::format("Gram's law fails at Gram point {}", p.n),
                 std::format("At \\(g_{{{}}} \\approx {}\\), \\((-1)^{{{}}} Z(g) < 0\\): the zero is not in its Gram interval.", p.n, number(p.t.mid(), 4), p.n),
                 "proved here", Object{{"t", p.t.mid()}});
    // two zeros much closer than the mean spacing, as Lehmer's pair near 7005
    for (std::size_t i = 1; i < z.zeros.size(); ++i) {
        const double a = z.zeros[i - 1].mid(), b = z.zeros[i].mid(), mean = spatium::zeta_detail::spacing((a + b) / 2);
        if (b - a < mean / 4)
            note("Two zeros close together",
                 std::format("Zeros at {} and {} are {} apart, a {} of the mean spacing there; Z barely crosses between them.",
                             number(a, 6), number(b, 6), number(b - a, 4), number((b - a) / mean, 3)),
                 "proved here", Object{{"t0", a}, {"t1", b}});
    }
    // the zeros of zeta' here, right of the line as Speiser's theorem ties to the Riemann hypothesis
    const auto critical = spatium::zeta_prime_zeros({-1, std::max(t0, 1.0)}, {3, t1});
    if (!critical.empty()) {
        double leftmost = critical.front().re;
        for (const auto& z : critical) leftmost = std::min(leftmost, z.re);
        note("Zeros of zeta' found here",
             std::format("Newton's method on \\(\\zeta'\\) from a grid of starts finds {} zeros of \\(\\zeta'\\) with \\(-1 \\le \\sigma \\le 3\\) here; "
                         "the leftmost has \\(\\sigma \\approx {}\\), right of the critical line.",
                         critical.size(), number(leftmost, 4)),
             "computed here", Object{{"s", Array{critical.front().re, critical.front().im}}});
        note("Where the zeros of zeta' lie",
             "The Riemann hypothesis holds exactly when \\(\\zeta'\\) has no zero with \\(0 < \\sigma < \\tfrac12\\).", "theorem",
             Object{{"source", "Speiser 1934; Levinson and Montgomery 1974"}});
    }
    // the moments over [0, t1], beside their asymptotics, where that integral is short
    if (0 < t1 && t1 <= 1000) {
        const double T = t1, pi = std::numbers::pi, gamma = std::numbers::egamma;
        const auto m1 = spatium::zeta_moment(1, T), m2 = spatium::zeta_moment(2, T);
        const double a1 = T * std::log(T / (2 * pi)) + (2 * gamma - 1) * T, a2 = T * std::pow(std::log(T), 4) / (2 * pi * pi);
        note("Moments of zeta on the line",
             std::format("\\(M_1({0}) = {1}\\), \\(M_1 / (T\\log\\tfrac{{T}}{{2\\pi}} + (2\\gamma - 1)T) = {2}\\); "
                         "\\(M_2({0}) = {3}\\), \\(M_2 / \\tfrac{{T\\log^4 T}}{{2\\pi^2}} = {4}\\), still far from its limit 1 at this height.",
                         number(T, 0), number(m1.value, 4), number(m1.value / a1, 4), number(m2.value, 2), number(m2.value / a2, 4)),
             "computed here", Object{{"t0", 0.0}, {"t1", T}});
        note("The moments for k at least 3",
             "\\(M_k(T) \\sim c_k\\, T \\log^{k^2} T\\) with \\(c_1, c_2\\) known (Hardy and Littlewood, Ingham); for \\(k \\ge 3\\) the "
             "constant is conjectured (Conrey and Ghosh, Keating and Snaith), not proved.",
             "conjecture", Object{{"source", "Keating and Snaith 2000"}});
    }
    if (t0 < 0)
        note("Below t = 0, the mirror",
             "\\(\\zeta(\\bar s) = \\overline{\\zeta(s)}\\), so \\(Z(-t) = Z(t)\\): each zero \\(\\tfrac12 + it\\) has its twin \\(\\tfrac12 - it\\), "
             "and the zeros and Gram points below 0 are those above, reflected.",
             "theorem", Object{});
    note("Z is real on the critical line",
         "\\(Z(t) = e^{i\\theta(t)}\\zeta(\\tfrac12 + it)\\) and \\(|Z(t)| = |\\zeta(\\tfrac12 + it)|\\): a zero of \\(\\zeta\\) on the line is a change of Z's sign.",
         "theorem", Object{});
    return Object{{"notes", std::move(out)}};
}

} // namespace detail

inline spatium::io::JsonValue zeta_request(std::string_view request) {
    using namespace detail;
    const auto refuse = [](std::string why) { return JsonValue(Object{{"refused", std::move(why)}}); };
    const auto json = spatium::io::parse(request);
    if (!json || !json->is_object()) return refuse("the request is not a JSON object");
    const std::string what = json->string_or("what", "");
    const double t0 = json->number_or("t0", 0), t1 = json->number_or("t1", 50);
    if (what == "hardy" || what == "zeros" || what == "gram" || what == "notes")
        if (!(t0 < t1) || t1 > 1e5 || t0 < -1e5) return refuse("the window is not -100000 <= t0 < t1 <= 100000");
    if (what == "hardy") {
        const int n = static_cast<int>(json->number_or("n", 512));
        if (n < 2 || n > 8192) return refuse("n is not between 2 and 8192");
        Array ts, zs;
        for (int i = 0; i < n; ++i) {
            const double t = t0 + (t1 - t0) * i / (n - 1);
            const I z = spatium::hardy_z(I(t));
            ts.push_back(t);
            zs.push_back(finite_or_null(z.mid()));
        }
        return Object{{"t", std::move(ts)}, {"z", std::move(zs)}};
    }
    if (what == "zeros") return zeros_json(zeros_on(t0, t1));
    if (what == "gram") return gram_json(gram_on(t0, t1));
    if (what == "notes") return notes(t0, t1);
    if (what == "grid") {
        const auto corner = [&](const char* key) -> std::pair<double, double> {
            const JsonValue* c = json->find(key);
            if (!c || !c->is_array() || c->as_array().size() != 2) return {NAN, NAN};
            return {c->as_array()[0].as_number(NAN), c->as_array()[1].as_number(NAN)};
        };
        const auto [x0, y0] = corner("s0");
        const auto [x1, y1] = corner("s1");
        const int w = static_cast<int>(json->number_or("w", 64)), h = static_cast<int>(json->number_or("h", 64));
        if (!(x0 < x1 && y0 < y1) || w < 1 || h < 1 || w * h > 1 << 18) return refuse("the grid is not s0 < s1 with w x h at most 262144");
        Array re, im;
        for (int j = 0; j < h; ++j)
            for (int i = 0; i < w; ++i) {
                const double x = x0 + (x1 - x0) * (i + 0.5) / w, y = y0 + (y1 - y0) * (j + 0.5) / h;
                const auto z = spatium::zeta(spatium::Complex<double>{x, y});
                re.push_back(finite_or_null(z.re));
                im.push_back(finite_or_null(z.im));
            }
        return Object{{"re", std::move(re)}, {"im", std::move(im)}};
    }
    return refuse("what is none of hardy, zeros, gram, grid, notes");
}

} // namespace spatium::answer
