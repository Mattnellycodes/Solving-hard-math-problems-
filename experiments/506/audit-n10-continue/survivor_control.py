"""Positive control for survivor_exact_real.py: same equations and labelling, but with the
saturation rho(rho^2-1) != 0 dropped and rho = 1, x^2 + y^2 = 1 imposed (all ten points and O on the
unit circle): a genuine (degenerate) solution exists, so the Gröbner basis must NOT be {1}."""
import json, sympy as sp
from lib import mask, bits, check_structure
S, R = list(range(5)), list(range(5, 10))
rec = json.load(open("runs/survivors.json"))[0]
L = [mask(l) for l in rec["L"]]
x, y, rho, c, s = sp.symbols("x y rho c s")
rel = [4 * c**2 - 2 * c - 1, s**2 + c**2 - 1]
def cs(k):
    k %= 10
    return sp.expand(sp.chebyshevt(k, c)), (sp.Integer(0) if k == 0 else sp.expand(sp.chebyshevu(k - 1, c) * s))
lab = [0, 2, 8, 4, 6, 0, 2, 8, 4, 6]
P = {}
for p in range(10):
    C, Sn = cs(lab[p]); P[p] = (C, Sn) if p < 5 else (rho * C, rho * Sn)
eqs = []
for l in L:
    a, b, d = bits(l)
    eqs.append(sp.expand(sp.Matrix([[P[p][0]**2 + P[p][1]**2, P[p][0], P[p][1], 1] for p in (a, b, d)] + [[x**2 + y**2, x, y, 1]]).det()))
G = sp.groebner(eqs + rel + [rho - 1, x**2 + y**2 - 1], x, y, rho, c, s, order="grevlex", domain="QQ")
print("control (rho = 1, |O| = 1, no saturation): GB == {1}?", list(G.exprs) == [1], "(must be False)")
assert list(G.exprs) != [1]
print("control passed")
