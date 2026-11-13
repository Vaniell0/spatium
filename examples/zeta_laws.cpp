// A law of the zeta function found from its numbers alone, by the library's synthesis of a law from
// noisy data (spatium/answer/law.hpp): the shortest description of the data, its coefficients
// read as named constants only where they could not have landed there by chance.
//
//   zeta_laws moment <blocks> <noise power>   the moment integral over [0, T] at each block's end,
//                                             from what zeta_moments printed (C D conjecture
//                                             reality ...); the noise of each point taken as
//                                             T^power, 0.25 for the second moment
//   zeta_laws zeros <table>                   the count of zeros up to T, from a table of the
//                                             zeros' heights one a line (Odlyzko's zeros1)
//
// For the second moment to T = 5e6 it finds T log T + (2 gamma - 1 - log 2pi) T, Hardy and
// Littlewood's law in Ingham's form; for the count of zeros to 7e4 the leading terms of Riemann
// and von Mangoldt's (T / 2pi) log(T / 2pi) - T / 2pi, its 7/8 left as digits, which the jumps
// of the count leave too uncertain to name.

#include <spatium/answer/law.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace {

void print(const spatium::answer::Law& l, std::size_t points, double top) {
    std::printf("%zu points, to T = %.4g\nthe law: %s\n", points, top, l.text.c_str());
    for (std::size_t j = 0; j < l.terms.size(); ++j) {
        std::printf("  %-14s %.10g +- %.2g", spatium::io::build::to_text(l.terms[j], "T").value_or("?").c_str(),
                    l.coefficients[j], l.errors[j]);
        if (l.named[j]) {
            const spatium::answer::NamedConstant& k = *l.named[j];
            std::printf("   = %s, %.10g (%.1f bits; by chance: %.2g)", k.text.c_str(), k.value, k.cost, k.by_chance);
        }
        std::printf("\n");
    }
    std::printf("residual rms, in units of the noise: %.3g\n", l.rms);
}

} // namespace

int main(int argc, char** argv) {
    const std::string_view mode = argc > 2 ? argv[1] : "";
    std::vector<double> x, y, sigma;
    if (mode == "moment" && argc > 3) {
        std::ifstream in(argv[2]);
        const double power = std::stod(argv[3]);
        std::map<double, double> block;   // its end, its integral; blocks come in any order
        for (std::string line; std::getline(in, line);) {
            std::istringstream s(line);
            double c, d, conjecture, reality;
            if (s >> c >> d >> conjecture >> reality) block[d] = reality;
        }
        double sum = 0;
        for (const auto& [end, value] : block) sum += value, x.push_back(end), y.push_back(sum), sigma.push_back(std::pow(end, power));
    } else if (mode == "zeros") {
        std::ifstream in(argv[2]);
        std::vector<double> zeros;
        for (double g; in >> g;) zeros.push_back(g);
        std::ranges::sort(zeros);
        for (double t = 100; !zeros.empty() && t < zeros.back(); t *= 1.05) {
            x.push_back(t);
            y.push_back(static_cast<double>(std::ranges::upper_bound(zeros, t) - zeros.begin()));
            sigma.push_back(1.0);
        }
    } else {
        std::fprintf(stderr, "usage: zeta_laws moment <blocks> <noise power> | zeta_laws zeros <table>\n");
        return 2;
    }
    if (x.size() < 8) {
        std::fprintf(stderr, "too few points: %zu\n", x.size());
        return 1;
    }
    // terms of up to five nodes over T with its logarithm and root, sums of up to three of them
    const spatium::answer::LawSpec spec{.max_size = 5, .max_terms = 3};
    print(spatium::answer::synthesize_law(x, y, sigma, spec), x.size(), x.back());
}
