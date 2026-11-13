# What docs/zeta/sixth-moment.md says about the sixth moment, recomputed from the blocks it gives:
#
#   python3 scripts/external/zeta/analyse.py docs/zeta/sixth-moment-blocks.tsv
#
# 1. Table 5.1.2 of Conrey, Farmer, Keating, Rubinstein and Snaith ("Integral moments of L-functions",
#    Proc. London Math. Soc. 2005) against the same blocks here: the conjecture column and the
#    "reality" column, each as a relative difference.
# 2. The ratio of the moment to the conjecture over [0, T], and the remainder E(T) = moment -
#    conjecture against sqrt(T), decade by decade.
# 3. The form found blind: the cumulative moment at T a factor 1.15 apart, fitted by the integral of
#    a polynomial Q in log(t / 2 pi) of each degree, the noise taken as s sqrt(T); BIC picks the
#    degree, and P_3 with no free parameter is set against the best fit.
import math
import sys

P3 = [0.9165155076378930590178543, 6.59391302064975810465713392, 10.274330830703446134183009522,
      8.60731914578120675614834763629, 3.984385094823534724747964073429, 1.0459251779054883439385323798059,
      0.14840073080150272680851401518774, 0.011072455215246998350410400826667,
      0.00040502133088411440331215332025984, 0.000005708527034652788398376841445252313]  # lowest power first

# Table 5.1.2 of CFKRS: [C, D], conjecture, reality
CFKRS = [(0, 50000, 7236872972.7, 7231005642.3), (50000, 100000, 15696470555.3, 15723919113.6),
         (100000, 150000, 21568672884.1, 21536840937.9), (150000, 200000, 26381397608.2, 26246250354.1),
         (200000, 250000, 30556177136.5, 30692229217.8), (250000, 300000, 34290291841.0, 34414329738.9),
         (300000, 350000, 37695829854.3, 37683495193.0), (350000, 400000, 40843941365.7, 40566252008.5),
         (400000, 450000, 43783216365.2, 43907511751.1), (450000, 500000, 46548617846.7, 46531247056.9),
         (500000, 550000, 49166313161.9, 49136264678.2), (550000, 600000, 51656498739.2, 51744796875.0),
         (600000, 650000, 54035153255.1, 53962410634.2), (650000, 700000, 56315178564.8, 56541799179.3),
         (700000, 750000, 58507171421.6, 58365383245.2), (750000, 800000, 60619962488.2, 60870809317.1),
         (800000, 850000, 62661003164.6, 62765220708.6), (850000, 900000, 64636649728.0, 64227164326.1),
         (900000, 950000, 66552376294.2, 65994874052.2), (950000, 1000000, 68412937271.4, 68961125079.8),
         (1000000, 1050000, 70222493232.7, 70233393177.0), (1050000, 1100000, 71984709805.4, 72919426905.7),
         (1100000, 1150000, 73702836332.4, 72567024812.4), (1150000, 1200000, 75379769148.4, 76267763314.7),
         (1200000, 1250000, 77018102997.5, 76750297112.6), (1250000, 1300000, 78620173202.6, 78315210623.9),
         (1300000, 1350000, 80188090542.5, 80320710380.9), (1350000, 1400000, 81723770322.2, 80767881132.6),
         (1400000, 1450000, 83228956776.3, 83782957374.3)]
CFKRS_WHOLE = (3317437762612.4, 3317496016044.9)   # [0, 2350000]


def int_power(T, j):
    """The integral of log(t / 2 pi)^j over [0, T]: with t = 2 pi e^x, 2 pi e^x sum_i (-1)^(j-i) j!/i! x^i."""
    x = math.log(T / (2 * math.pi))
    return 2 * math.pi * math.exp(x) * sum((-1) ** (j - i) * math.factorial(j) / math.factorial(i) * x ** i for i in range(j + 1))


def least_squares(cols, rhs):
    """Modified Gram-Schmidt: no normal equations, whose condition is the square of the basis's."""
    q = [c[:] for c in cols]
    m = len(q)
    r = [[0.0] * m for _ in range(m)]
    for j in range(m):
        for i in range(j):
            r[i][j] = sum(a * b for a, b in zip(q[i], q[j]))
            q[j] = [a - r[i][j] * b for a, b in zip(q[j], q[i])]
        r[j][j] = math.sqrt(sum(a * a for a in q[j]))
        q[j] = [a / r[j][j] for a in q[j]]
    z = [sum(a * b for a, b in zip(q[j], rhs)) for j in range(m)]
    c = [0.0] * m
    for j in reversed(range(m)):
        c[j] = (z[j] - sum(r[j][i] * c[i] for i in range(j + 1, m))) / r[j][j]
    return c


def main(path):
    blocks = {}
    for line in open(path):
        if line.startswith("#") or line.startswith("C\t"):
            continue
        c, d, conj, real = line.split()
        blocks[(int(c), int(d))] = (float(conj), float(real))
    keys = sorted(blocks)

    print("1. CFKRS Table 5.1.2 against these blocks (relative differences, theirs to ours)")
    for c, d, conj, real in CFKRS:
        oc, orl = blocks[(c, d)]
        print(f"   [{c}, {d}]  conjecture {(conj - oc) / oc: .1e}  reality {(real - orl) / orl: .2e}")
    conj = sum(blocks[k][0] for k in keys if k[1] <= 2350000)
    real = sum(blocks[k][1] for k in keys if k[1] <= 2350000)
    print(f"   [0, 2350000]  ours {real:.1f} / {conj:.1f} = {real / conj:.7f}; "
          f"theirs {CFKRS_WHOLE[1]} / {CFKRS_WHOLE[0]} = {CFKRS_WHOLE[1] / CFKRS_WHOLE[0]:.7f}")

    print("2. The moment over [0, T] against the conjecture")
    cum, conj, ends = 0.0, 0.0, []
    for k in keys:
        conj += blocks[k][0]
        cum += blocks[k][1]
        ends.append((k[1], cum, conj))
    at = {T: (m, c) for T, m, c in ends}
    for T in (2350000, 10 ** 7, 10 ** 8, 2 * 10 ** 8):
        m, c = at[T]
        print(f"   T = {T:.3g}: ratio {m / c:.7f}, E = {m - c:.4g}, E / sqrt T = {(m - c) / math.sqrt(T):.3g}")
    for lo, hi in ((2e5, 2e6), (2e6, 2e7), (2e7, 2e8 + 1)):
        v = [(m - c) / math.sqrt(T) for T, m, c in ends if lo <= T < hi]
        print(f"   T in [{lo:.0e}, {hi:.0e}): rms of E / sqrt T {math.sqrt(sum(x * x for x in v) / len(v)):.3g} over {len(v)} ends")

    print("3. The form, blind: the cumulative moment fitted by int_0^T Q(log(t / 2 pi)) dt, noise s sqrt(T)")
    pick, last = [], 0
    for T, m, c in ends:
        if T >= 2e5 and T >= last * 1.15:
            pick.append((T, m, c))
            last = T
    if pick[-1] != ends[-1]:
        pick.append(ends[-1])
    n = len(pick)
    w = [1 / math.sqrt(T) for T, _, _ in pick]
    y = [m * wi for (T, m, _), wi in zip(pick, w)]
    rss0 = sum(((m - c) * wi) ** 2 for (T, m, c), wi in zip(pick, w))
    best = None
    for deg in range(12):
        cols = [[int_power(T, j) * wi for (T, _, _), wi in zip(pick, w)] for j in range(deg + 1)]
        # the powers of log are nearly collinear over the range: scale each column before the fit
        scale = [math.sqrt(sum(v * v for v in col)) for col in cols]
        coef = least_squares([[v / s for v in col] for col, s in zip(cols, scale)], y)
        rss = sum((sum(coef[j] * cols[j][i] / scale[j] for j in range(deg + 1)) - y[i]) ** 2 for i in range(n))
        bic = n * math.log(rss / n) + (deg + 1) * math.log(n)
        print(f"   degree {deg:2d}: rms of (moment - fit) / sqrt T {math.sqrt(rss / n):.3g}, BIC {bic:.1f}")
        if best is None or bic < best[1]:
            best = (deg, bic)
    bic0 = n * math.log(rss0 / n)
    print(f"   BIC picks degree {best[0]}; P_3 with no free parameter: BIC {bic0:.1f}, "
          f"{'better' if bic0 < best[1] else 'worse'} by {abs(bic0 - best[1]):.1f}")


if __name__ == "__main__":
    main(sys.argv[1])
