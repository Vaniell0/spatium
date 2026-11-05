# The language

Written by `examples/site_reference.cpp` from the library's own tables, as the site's reference is; do not edit. Every example here is answered by the library, and its formula is the one it writes.

## Writing a question

In the formula editor, what is typed becomes the formula shown, which reads as the line of the language beside it.

| type | shows | reads as | |
|---|---|---|---|
| `sqrt` | $\sqrt{2}$ | `sqrt(2)` | a square root |
| `/` | $\frac{1}{3}$ | `1 / 3` | a fraction |
| `^` | $2^{10}$ | `2^10` | a power; the arrow keys leave it |
| `pi` | $\pi$ | `pi` | the number pi |
| `<=` | $x^{2} \le 4$ | `x: x^2 <= 4` | an inequality; >= gives the other |
| `in` | $x^{2} = 2, x \in \left[0, 2\right]$ | `x: x^2 = 2 on [0, 2]` | where to look for the roots |
| `Shift+Enter` | $x + y = 3$<br>$x - y = 1$ | `x, y: x + y = 3, x - y = 1` | a new row; rows of equations are a system |

## Questions

| form | what it asks | example | answer |
|---|---|---|---|
| $f$ | an expression, evaluated, its digits proved | $2\left(3 + 4\right)^{2} - \frac{1}{4}$ | Exact |
| $f = g$ | an equation in one unknown: in closed form with its steps where it is a quadratic, and every real root proved | $2x^{2} - 3x - 5 = 0$ | Exact, closed form, checked by Vieta's formulas |
| $f = g,\ x \in [a, b]$ | every root on the interval, each proved alone in its own | $\sin\left(x\right) = 0.5, x \in \left[0, 10\right]$ | Proved, 4 roots, no other on the box |
| $f > g$ | an inequality: where it holds, by the interval method | $x^{2} - 4 > 0$ | Proved, by the interval method |
| $f = g \\ h = k$ | a system: one row per equation, an equation per unknown | $x + y = 3$<br>$x - y = 1$ | Exact, by elimination in fractions |
| $p \in [c, d]$ | a parameter: a row below that reads it is answered as a function of p, synthesised and certified | $p \in \left[0.5, 4\right]$<br>$x^{2} - p = 0, x \in \left[0, 10\right]$ | Relative error below 10⁻¹⁵, certified on 2001 points |
| $I = \int_a^b f\,dt$ | an integral; a limit may be a parameter above | $x \in \left[0.5, 6\right]$<br>$\mathrm{Si} = \int_{0}^{x} \frac{\sin\left(t\right)}{t}\,dt$ | Relative error below 10⁻¹², certified on 2001 points |

## Functions

The editor's notation, the language's text, where each has a value, and an example.

| written | in the language | | defined for | example |
|---|---|---|---|---|
| $\sin x$ | `sin` | sine, of an angle in radians | every x | `sin(pi / 6)` |
| $\cos x$ | `cos` | cosine, of an angle in radians | every x | `cos(pi / 3)` |
| $\tan x$ | `sin(pi / 4) / cos(pi / 4)` | tangent, sin x / cos x | every x where cos x is not 0 | `sin(pi / 4) / cos(pi / 4)` |
| $\exp x$ | `exp` | exponential, e raised to x | every x | `exp(1)` |
| $\ln x$ | `log` | natural logarithm, the inverse of exp | x > 0 | `log(2)` |
| $\lg x$ | `log(100) / log(10)` | logarithm in base 10; \log x alone is the same | x > 0 | `log(100) / log(10)` |
| $\log_{b} x$ | `log(81) / log(3)` | logarithm in base b | x > 0, b > 0, b not 1 | `log(81) / log(3)` |
| $\log_{2} x$ | `log(8) / log(2)` | logarithm in base 2 | x > 0 | `log(8) / log(2)` |
| $\sqrt{x}$ | `sqrt` | square root | x >= 0 | `sqrt(2)` |
| $\sqrt[n]{x}$ | `8^(1/3)` | the n-th root | x > 0 | `8^(1/3)` |
| $\left\|x\right\|$ | `abs` | absolute value, the size of x | every x | `abs(-3)` |
| $\min\left(x, y\right)$ | `min` | the smaller of two numbers | every x, y | `min(2, 3)` |
| $\max\left(x, y\right)$ | `max` | the larger of two numbers | every x, y | `max(2, 3)` |
|  | `less` | 1 where the first is below the second, else 0: a step that composes by multiplication | every x, y | `less(1, 2)` |

## Operators

Tightest last.

| written | | example |
|---|---|---|
| `a + b, a - b` | sum and difference, read left to right: a - b - c is (a - b) - c | `7 - 2 - 1` |
| `a * b, a / b` | product and quotient, tighter than a sum, read left to right | `12 / 3 / 2` |
| `2x, 2(a + b), (a)(b)` | a product may leave out its * after a number or a closing parenthesis, not after a name | `2(3 + 4)` |
| `-a` | minus, tighter than a product and looser than a power: -x^2 is -(x^2) | `-3^2` |
| `a^b` | a power: a whole n from 1 to 64 is a product of n factors; any other b, a^b for a > 0 (a^0 is 1) | `2^0.5` |
| `(a)` | grouping | `(1 + 2) * 3` |
| `a = b` | an equation in one unknown, solved for it | `x^2 - 4 = 0` |

## Constants

| name | |
|---|---|
| `pi` | the double nearest pi |
| `inf` | infinity |
| `nan` | not a number: what a function gives outside its domain |

## The language of the canvas

A document is lines, a definition each; a line reads the names the lines above it define.

| line | | example | as a formula |
|---|---|---|---|
| `name = f` | a value, named; a line below reads it by its name | `r = 0.5<br>A = pi * r^2` | $r = 0.5$<br>$A = \pi \cdot r^{2}$ |
| `x: f = g` | the roots of an equation in x | `x: 2x^2 - 3x - 5 = 0` | $2x^{2} - 3x - 5 = 0$ |
| `x: f = g on [a, b]` | every root on the interval | `x: sin(x) = 0.5 on [0, 10]` | $\sin\left(x\right) = 0.5, x \in \left[0, 10\right]$ |
| `x: f > g` | where an inequality holds; also <, <= and >= | `x: x^2 - 4 > 0` | $x^{2} - 4 > 0$ |
| `x, y: f = g, h = k` | a system: an equation for each unknown | `x, y: x + y = 3, x - y = 1` | $x + y = 3$<br>$x - y = 1$ |
| `p in [a, b]` | a parameter: a line that reads it is answered as a function of it | `p in [0.5, 4]<br>x: x^2 - p = 0 on [0, 10]` | $p \in \left[0.5, 4\right]$<br>$x^{2} - p = 0, x \in \left[0, 10\right]$ |
| `I = integral f dt on [a, b]` | an integral in the variable after d | `I = integral t^2 dt on [0, 1]` | $I = \int_{0}^{1} t^{2}\,dt$ |
| `S = torus(R, r)` | a shape, from the operations below; another reads it by name | `S0 = torus(2, 0.5)<br>S1 = offset(S0, 0.1)` |  |
| `# words` | a comment, to the end of the line | `r = 2  # the radius` | $r = 2$ |

Shapes:

| operation | |
|---|---|
| `torus(nn)` | a torus: the radius to the tube's centre, and the tube's |
| `cylinder(nn)` | a closed cylinder on the z axis: its radius and its height |
| `sphere(n)` | a sphere at the origin: its radius |
| `flake(nnnn)` | a shallow cap in a slab: the slab's half extents, and how much wider its sphere is |
| `offset(sf)` | a closed surface pushed along its normal by a thickness |
| `scatter(ssnn)` | copies of a shape on a surface: the shape, the surface, how many, and a seed |
| `union(s)` | the shapes together |

`n` is a number, `s` a shape.

What the library says of a document it does not take:

| | document | told |
|---|---|---|
| a name read before a line defines it | `A = pi * r^2` | line 1: A reads 'r', which no line above defines |
| a name defined twice | `r = 1<br>r = 2` | line 2: 'r' is defined twice |
| a cycle: a line reads only the lines above it | `a = b + 1<br>b = a - 1` | line 1: a reads 'b', which no line above defines |

## Answers

How strongly an answer holds, strongest first.

| | promises |
|---|---|
| Exact | The value itself, in fractions or in closed form: nothing in it is rounded. |
| Every root proved, and no other | Each root lies in an interval that holds exactly one, and the rest of the interval or box is proved to hold none. |
| Proved to their digits | Every digit shown is proved: the true value lies between two bounds that share them. |
| Certified on a grid | A function of a parameter, its error below the bound at every point of a grid; between the points it is not proved. |
| Error estimated | The error is estimated, not bounded: likely that close, not proved. |
| Computed | Computed in doubles; nothing is proved about its digits. |
| Not answered, and why | No answer, and the reason the library gives for it. |
