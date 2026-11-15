# What a scientific calculator does, against what the question language does: each feature of the
# Casio fx-991EX (as its manual and a walk through it list them) asked of the library's own door,
# the host build of the page's core, as a person would type it -- a line of the language, or the
# MathJSON row the page's formula editor writes for its notation. A feature is
#
#   language   the door answers it;
#   library    the door does not, and the library has it: a header and the name it declares there,
#              checked to be in that header;
#   absent     neither.
#
#   c++ -std=c++23 -O2 -I include -I . site/wasm/core.cpp -o core_host
#   python3 scripts/calculator_inventory.py ./core_host
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# (group, feature, what is typed: a line, "rows <json>" or "evaluate <json>", or None;
#  (header under include/spatium, name it declares) where the library has it, or None)
FEATURES = [
    ("numbers", "fractions, exact and as a decimal", "1/3 + 1/6", None),
    ("numbers", "decimals exactly (ten times 0.1 is 1)", "0.1 * 10", None),
    ("numbers", "mixed numbers", 'rows {"rows": [["InvisibleOperator", 1, ["Divide", 13, 20]]]}', None),
    ("numbers", "percent", 'rows {"rows": [["Multiply", {"num": "50e-2"}, 8]]}', None),
    ("numbers", "roots in simplest form (sqrt 24 + sqrt 150 = 7 sqrt 6)", "sqrt(24) + sqrt(150)", None),
    ("numbers", "pi in exact results (2 pi)", "pi * 2", None),
    ("numbers", "powers and n-th roots", "8^(1/3)", None),
    ("numbers", "engineering prefixes, Fix/Sci/Norm, rounding", None, None),
    ("numbers", "base n and logic (hex, bin, and, or)", None, None),
    ("numbers", "degrees, minutes, seconds", None, None),
    ("numbers", "variables and the last answer", 'rows {"rows": [{"text": "a = 2"}, {"text": "a^2 + 1"}]}', None),
    ("numbers", "substitute a value (CALC)",
     'evaluate {"document": "p in [0, 4]\\nf = p^2 + 1", "name": "f", "chosen": {"p": 2}}', None),
    ("functions", "absolute value", "abs(-3)", None),
    ("functions", "natural logarithm", "log(10)", None),
    ("functions", "logarithm in base 10 and base b", 'rows {"rows": [["Log", 100, 10]]}', None),
    ("functions", "sine and cosine, exact at pi/6, pi/4, pi/3", "sin(pi / 6)", None),
    ("functions", "angles in degrees", 'rows {"rows": [["Sin", ["Degrees", 30]]]}', ("algebra/literals.hpp", "_deg")),
    ("functions", "inverse trigonometric functions", 'rows {"rows": [["Arcsin", {"num": "0.5"}]]}', ("algebra/dual.hpp", "asin")),
    ("functions", "hyperbolic functions", 'rows {"rows": [["Sinh", 1]]}', ("algebra/dual.hpp", "sinh")),
    ("functions", "factorial", 'rows {"rows": [["Factorial", 5]]}', ("discrete/combinatorics.hpp", "factorial")),
    ("functions", "combinations and permutations", 'rows {"rows": [["Binomial", 5, 2]]}',
     ("discrete/combinatorics.hpp", "binomial_coefficient")),
    ("functions", "random numbers", None, None),
    ("functions", "prime factorisation, gcd, lcm", None, None),
    ("functions", "polar and rectangular coordinates", None, ("algebra/complex.hpp", "from_polar")),
    ("calculus", "derivative at a point", "d/dx x^3 at 2", None),
    ("calculus", "definite integral", "I = integral t^2 dt on [0, 1]", None),
    ("calculus", "sum of a series", 'rows {"rows": [["Sum", ["Power", "k", 2], ["Tuple", "k", 1, 10]]]}', None),
    ("equations", "a root of any equation, numerically (SOLVE)", "x: cos(x) = x on [0, 1]", None),
    ("equations", "quadratic, in closed form with its complex roots", "x^2 + 2x + 5 = 0", None),
    ("equations", "cubic and quartic", "x^4 - 5x^2 + 4 = 0", None),
    ("equations", "linear systems of two to four unknowns", 'rows {"rows": [{"text": "x + 2y = 5"}, {"text": "3x - y = 1"}]}', None),
    ("equations", "systems that are not linear", "x, y: x^2 + y^2 = 4, x - y = 1", None),
    ("equations", "inequalities, as the set where they hold", "x: x^2 - 4 > 0", None),
    ("equations", "proportions", "x: x / 3 = 4 / 6", None),
    ("equations", "the vertex of a parabola", None, None),
    ("complex", "complex arithmetic: abs, arg, conjugate",
     'rows {"rows": [["Multiply", ["Add", 1, ["Multiply", 2, "ImaginaryUnit"]], ["Subtract", 3, "ImaginaryUnit"]]]}',
     ("algebra/complex.hpp", "conjugate")),
    ("matrices", "determinant and inverse, to 4 x 4", None, ("algebra/matrix.hpp", "determinant")),
    ("vectors", "dot and cross product, angle, unit vector", None, ("algebra/functions.hpp", "cross")),
    ("statistics", "measures of a sample (mean, deviation, quartiles)", None, None),
    ("statistics", "regressions", None, ("answer/synthesis.hpp", "least_squares")),
    ("statistics", "distributions (normal, binomial, Poisson)", None, None),
    ("statistics", "a spreadsheet of cells", None, None),
    ("statistics", "a table of a function's values", None, None),
    ("units", "physical constants", None, ("physics/mechanics/units.hpp", "h_planck")),
    ("units", "unit conversions", None, ("physics/mechanics/units.hpp", "Quantity")),
]


def declares(header, name):
    text = (ROOT / "include" / "spatium" / header).read_text()
    return re.search(r"(?<![A-Za-z0-9_])" + re.escape(name) + r"(?![A-Za-z0-9_])", text) is not None


def main(core):
    asked = [f[2] for f in FEATURES if f[2]]
    run = subprocess.run([core], input="\n".join(asked) + "\n", capture_output=True, text=True, check=True)
    answers = dict(zip(asked, run.stdout.splitlines()))
    rows, counts = [], {"language": 0, "library": 0, "absent": 0}
    for group, feature, typed, where in FEATURES:
        answered = False
        why = ""
        if typed:
            a = json.loads(answers[typed])
            answered = "refused" not in a
            why = a.get("refused", "")
        if answered:
            kind = "language"
        elif where and declares(*where):
            kind = "library"
        elif where:
            raise SystemExit(f"{feature}: {where[0]} does not declare {where[1]}")
        else:
            kind = "absent"
        counts[kind] += 1
        rows.append((group, feature, kind, f"{where[0]}: {where[1]}" if where and kind == "library" else "", why[:90]))
    print("| group | feature | | where in the library | what the door says |\n|---|---|---|---|---|")
    for r in rows:
        print("| " + " | ".join(c.replace("|", "\\|") for c in r) + " |")
    print(f"\n{len(rows)} features: {counts['language']} in the language, {counts['library']} in the library only, "
          f"{counts['absent']} absent")


if __name__ == "__main__":
    main(sys.argv[1])
