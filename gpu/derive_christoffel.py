#!/usr/bin/env python3
"""Symbolic derivation of the closed-form Christoffel symbols hand-
transcribed into christoffel_closed_form.hpp.

Not part of the build -- a provenance/reproducibility script. Run it
whenever christoffel_closed_form.hpp's formulas are in question:
    nix-shell -p "python3.withPackages(ps: [ps.sympy])" --run \
        'python3 gpu/derive_christoffel.py'

It (1) derives Gamma^lambda_{mu nu} = 1/2 g^{lambda sigma}(d_mu g_{sigma
nu} + d_nu g_{sigma mu} - d_sigma g_{mu nu}) directly from the exact
metric forms in schwarzschild.hpp/kerr.hpp, (2) checks that Kerr at a=0
reduces term-by-term to the Schwarzschild expressions (independent of
christoffel_closed_form.hpp -- this is the derivation's own internal
consistency check), and (3) prints each nonzero component with Sigma/
Delta substituted in, in the same form transcribed into the header.
The header's own correctness against the CPU engine is a separate,
numeric check: tests/test_relativity.cpp's "closed-form Christoffel"
cases compare christoffel_closed_form.hpp against geodesic.hpp's
Dual<T>-exact christoffel() directly.
"""
import sympy as sp

t, r, th, ph, M, a = sp.symbols('t r theta phi M a', real=True)
coords = [t, r, th, ph]
names = ['t', 'r', 'th', 'ph']


def christoffel(g):
    ginv = g.inv()
    n = 4
    dg = [[[sp.diff(g[mu, nu], coords[kappa]) for nu in range(n)] for mu in range(n)] for kappa in range(n)]
    Gamma = {}
    for lam in range(n):
        for mu in range(n):
            for nu in range(mu, n):
                s = sp.Integer(0)
                for sigma in range(n):
                    s += ginv[lam, sigma] * (dg[mu][sigma][nu] + dg[nu][sigma][mu] - dg[sigma][mu][nu])
                val = sp.simplify(sp.Rational(1, 2) * s)
                if val != 0:
                    Gamma[(lam, mu, nu)] = val
    return Gamma


f = 1 - 2 * M / r
g_schw = sp.diag(-f, 1 / f, r**2, r**2 * sp.sin(th)**2)
Gamma_schw = christoffel(g_schw)

sigma_s, delta_s = sp.symbols('sigma delta')
sigma_expr = r**2 + a**2 * sp.cos(th)**2
delta_expr = r**2 - 2 * M * r + a**2

g_kerr = sp.zeros(4, 4)
g_kerr[0, 0] = -(1 - 2 * M * r / sigma_expr)
g_kerr[0, 3] = g_kerr[3, 0] = -2 * M * r * a * sp.sin(th)**2 / sigma_expr
g_kerr[1, 1] = sigma_expr / delta_expr
g_kerr[2, 2] = sigma_expr
g_kerr[3, 3] = (r**2 + a**2 + 2 * M * r * a**2 * sp.sin(th)**2 / sigma_expr) * sp.sin(th)**2
Gamma_kerr = christoffel(g_kerr)

print(f"Schwarzschild: {len(Gamma_schw)} nonzero (mu<=nu) components")
print(f"Kerr:          {len(Gamma_kerr)} nonzero (mu<=nu) components")

print("\n=== sanity check: Kerr|a=0 == Schwarzschild, term by term ===")
ok = True
for key, val in Gamma_kerr.items():
    val0 = sp.simplify(val.subs(a, 0))
    ref = Gamma_schw.get(key, sp.Integer(0))
    if sp.simplify(val0 - ref) != 0:
        lam, mu, nu = key
        print(f"MISMATCH at Gamma^{names[lam]}_{names[mu]}{names[nu]}: {val0} != {ref}")
        ok = False
print("PASSED" if ok else "FAILED")

print("\n=== Kerr components (Sigma/Delta substituted) ===")
for key in sorted(Gamma_kerr.keys()):
    lam, mu, nu = key
    e = sp.together(Gamma_kerr[key]).subs(sigma_expr, sigma_s).subs(delta_expr, delta_s)
    print(f"Gamma^{names[lam]}_{names[mu]}{names[nu]} = {sp.simplify(e)}")
