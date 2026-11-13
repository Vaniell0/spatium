# The sixth moment of zeta to T = 2·10⁸

The integral of |ζ(½ + it)|⁶ over [0, T], computed to T = 2·10⁸ — 85 times past the last table of it we
know, Table 5.1.2 of Conrey, Farmer, Keating, Rubinstein and Snaith, *Integral moments of L-functions*
(Proc. London Math. Soc. 2005), which stops at 2 350 000 — and set against their conjecture

  ∫₀ᵀ |ζ(½ + it)|⁶ dt ≈ ∫₀ᵀ P₃(log(t / 2π)) dt,

P₃ the degree-9 polynomial of their (5.1.6). Hiary and Odlyzko (Math. Comp. 2012) tested the same
conjecture much higher, near 10²², but on windows: the integral from 0 they call out of reach at such
heights. This is the integral from 0.

The numbers: [`sixth-moment-blocks.tsv`](sixth-moment-blocks.tsv), 4000 blocks [C, D] of length 50 000,
each the integral of |ζ|⁶ there ("reality") beside the integral of P₃ ("conjecture"). Everything below is
recomputed from that file by [`scripts/external/zeta/analyse.py`](../../scripts/external/zeta/analyse.py).

## What it shows

| T | ∫₀ᵀ \|ζ\|⁶ / ∫₀ᵀ P₃ | E(T) = difference | E / √T |
|---|---|---|---|
| 2.35·10⁶ | 1.0000184 | 6.1·10⁷ | 4.0·10⁴ |
| 10⁷ | 0.9997625 | −6.9·10⁹ | −2.2·10⁶ |
| 10⁸ | 1.0000049 | 3.9·10⁹ | 3.9·10⁵ |
| 2·10⁸ | 0.9999789 | −4.5·10¹⁰ | −3.2·10⁶ |

- **The conjecture holds to 2·10⁸**, to 2·10⁻⁵ of the integral and closer at most heights.
- **The remainder grows faster than √T.** The rms of E/√T over each decade of T is 4.2·10⁵, 1.1·10⁶
  and 2.9·10⁶ from 2·10⁵ to 2·10⁸: log|E| against log T has slope 0.90. CFKRS expect a remainder
  O(T^{1/2+ε}); √T·log⁶T has local slope 0.5 + 6/log T ≈ 0.88 over this range, so the data agree with
  that, and three decades cannot tell it from a power 0.9.
- **The random-matrix leading term alone is far off here.** Keating and Snaith's leading coefficient
  42·a₃/9! ≈ 5.7·10⁻⁶ is the coefficient of log⁹; at these heights that term is 2–5 % of P₃, the rest
  being the lower coefficients, and it would be half of it only near t ≈ 7.5·10⁴¹.
- **The form cannot be read off the numbers.** Fitting the cumulative integral by ∫ Q(log(t/2π)) for a
  polynomial Q of each degree, the noise taken as s·√T: a degree of at least 5 is needed, and BIC picks
  5, not 9 — the higher coefficients are not resolved at these heights. P₃ itself, with no free
  parameter, beats every fit by 18.4 in BIC: the conjecture is the shorter description of the data.
  (Block by block the same test misleads: neighbouring blocks' differences are E(D) − E(C), strongly
  anticorrelated, not independent noise.)

## Table 5.1.2 of CFKRS

Their conjecture column agrees with ours to 10⁻¹² in all 29 blocks. Their "reality" column, computed
with Mathematica, is lower than ours in every one of the 29, by 1.9·10⁻⁷ to 4.2·10⁻⁶: one sign, so not
rounding. Three routes here agree on it:

| route to ∫₀^{50 000} \|ζ\|⁶ | value |
|---|---|
| CFKRS Table 5.1.2 | 7 231 005 642.3 |
| adaptive Gauss–Kronrod (`zeta_moment`: Euler–Maclaurin below t = 1000, Riemann–Siegel above) | 7 231 006 986.4067 |
| Riemann–Siegel Z, Gauss–Legendre with turned phases (`zeta_moment_fixed`) | agrees to 1.2·10⁻¹⁰ on every block to 2.35·10⁶ |
| Euler–Maclaurin ζ in doubles, Gauss–Legendre 40 points a unit, written apart | 7 231 006 986.4031 |

So their value of that block is low by 1344, and the whole [0, 2 350 000] by 2.9·10⁹ (ours
3 317 498 943 727, theirs 3 317 496 016 045); the ratio to the conjecture, 1.0000184 against their
1.0000176, is unchanged in what it says.

## What each number is worth

Computed, with an error budget; not proved in interval arithmetic.

| source of error | how it is known | size, relative to a block |
|---|---|---|
| where the Riemann–Siegel formula is cut | proved: Gabcke's \|R₂(t)\| ≤ 0.011 t^{−7/4} for t ≥ 200 (as Gourdon 2004 quotes it), plus the C₃ and C₄ terms the sum adds, `riemann_siegel_truncation`; carried into the moment as `error` | 3·10⁻⁸ at t ~ 10³, 6·10⁻¹⁸ at 2·10⁸ |
| rounding | Z against mpmath at 30 digits on random points: 6·10⁻¹¹ at 10³, 10⁻¹³ at 10⁵, 10⁻¹¹ at 10⁷, 9·10⁻¹⁰ at 2.6·10⁸, phases reduced in long double; 24 blocks spread over the range recomputed so against the file | ≤ 1.2·10⁻⁸ |
| the quadrature | the same rule at twice the points on one unit of every block | ≤ 4·10⁻⁸ on a unit, less on a block |
| below t = 1000 | ζ by Euler–Maclaurin with its remainder bounded, adaptive quadrature at 10⁻¹⁰ | ≤ 10⁻¹⁰ |

The file was computed before the phases moved to long double; with phases in doubles Z was off by up
to 6·10⁻⁶ at t = 2.6·10⁸, but the errors cancel over a block, and the 24 blocks recomputed differ from
the file by at most 1.2·10⁻⁸. The effects above are 10⁻⁵ in the cumulative ratio and 10⁻² in a block's.

Z itself is held apart from the moments: by the proved Z of Euler–Maclaurin in intervals (to 10⁻¹¹ on
[10³, 3·10⁵]), by mpmath, and by Odlyzko's tables of zeros
([`scripts/external/zeta/odlyzko_check.cpp`](../../scripts/external/zeta/odlyzko_check.cpp)): the first
100 000 zeros all found and none extra, at most 2.99·10⁻⁹ from the table (its stated accuracy 3·10⁻⁹);
the 9998 from zero number 10¹² + 2 on, at t ≈ 2.7·10¹¹, likewise, at most 2.53·10⁻⁹ (stated 10⁻⁸).

## Reproducing it

Spatium is header-only here: the headers and a C++23 compiler (GCC 15 tested) are all it takes, no
CMake, no dependency.

```bash
git clone https://github.com/Vaniell0/spatium && cd spatium
g++ -std=c++23 -O3 -march=native -pthread -Iinclude examples/zeta_moments.cpp -o zeta_moments
./zeta_moments 3 0 2350000 50000 12 > table.txt        # k, from, to, block, threads: CFKRS's whole table
./zeta_moments 3 0 200000000 50000 12 > blocks.txt     # the file here
python3 scripts/external/zeta/analyse.py docs/zeta/sixth-moment-blocks.tsv
```

In the CMake build it is the target `zeta_moments`.

`zeta_moments` prints a line per block as it finishes, the highest first: C, D, conjecture, reality,
their ratio, evaluations, seconds, the rule against itself at twice the points on the block's middle
unit, and the bound on the cut of the formula (for a block that starts below t = 1000, where the
adaptive rule runs, its own error estimate). CFKRS's whole table, [0, 2 350 000], takes 42 s. On an Intel i5-1235U, 12 threads, [0, 2·10⁸] takes about
6 hours; a block at 10⁸ about 18 s on one core alone.

How it is fast: on each unit of t, Gauss–Legendre with k·log(t/2π) + 16 points, which is past the band
|Z|^{2k} oscillates in; every node steps by exactly 1 from one unit to the next, so each phase t·log n
of the Riemann–Siegel sum goes on by a product with n^{−i} instead of a cosine, taken afresh every 64
units; sixteen units in one pass over the terms, so the phases stay in cache; C₀…C₄ from a table in p.
