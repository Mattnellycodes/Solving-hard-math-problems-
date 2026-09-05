"""Positive controls for two5_exact.exact_kill (the machinery must NOT return {1} when a real
solution exists).
 (a) concentric regular pentagons S at angles 2pi*k/5, R at angles 2pi*(2k+1)/10, with the 5 real
     lines x = const through s_k and r_{k}, r_{k-1}... at rho = 1/cos(pi/5): O = inf case must
     be satisfiable  -> Gröbner basis != {1}.
 (b) the survivor's line set with the saturation rho(rho^2-1) != 0 removed: rho = 1 (all points
     on one circle) with O on that circle satisfies everything -> basis != {1}."""
import json
import sympy as sp
from fractions import Fraction
from itertools import combinations
from lib import mask, bits
import two5_exact as T

# (a): labelling angles*10: S = 0,2,4,6,8 ; R = 1,3,5,7,9  (R at odd tenths)
x = [Fraction(k, 5) for k in range(5)] + [Fraction(2 * k + 1, 10) for k in range(5)]
# 4-blocks: s_i + s_j = r_k + r_l  (mod 1)
F4 = []
for i, j in combinations(range(5), 2):
    for k, l in combinations(range(5, 10), 2):
        if (x[i] + x[j] - x[k] - x[l]) % 1 == 0:
            F4.append(mask((i, j, k, l)))
print("(a) 4-blocks from the angle rule:", len(F4))
F = [mask(range(5)), mask(range(5, 10))] + F4
# vertical lines: s_k at angle k/5 (x-coordinate cos(2pi k/5)), r at angles (2k+1)/10 and (2k-1)/10
# have x-coordinate rho*cos(pi/5)... the line through s_0=(1,0) and r_{1/10}, r_{9/10} (x = rho cos(pi/5) = 1).
# by rotation symmetry: s_k with r at angles k/5 +- 1/10.
L = []
for k in range(5):
    a = Fraction(k, 5)
    rs = [p for p in range(5, 10) if (x[p] - a) % 1 in (Fraction(1, 10), Fraction(9, 10))]
    L.append(mask([k] + rs))
print("(a) lines:", [bits(l) for l in L])
fin_ok, inf_ok = T.exact_kill(F, L, x, verbose=True)
print("(a) O=inf case returned 'no solution'?", inf_ok, " (must be False)")
assert not inf_ok

# (b) survivor without saturation
surv = json.load(open("runs/survivors.json"))[0]
Fs = [mask(b) for b in surv["F"]]; Ls = [mask(l) for l in surv["L"]]
info = T.angle_labellings(Fs[2:])
xs = info["labellings"][0]
# rebuild the finite-O ideal without saturation (copy of exact_kill's construction)
zeta, rho, z0, w0, t = sp.symbols("zeta rho z0 w0 t")
e = 10
Phi = sp.cyclotomic_poly(e, zeta)
n_of = [int(xi * e) % e for xi in xs]
z, zb = {}, {}
for p in range(10):
    pw, pwb = zeta ** n_of[p], zeta ** ((e - n_of[p]) % e)
    z[p], zb[p] = (pw, pwb) if p < 5 else (rho * pw, rho * pwb)
eqs = []
for l in Ls:
    a, b, c = bits(l)
    rows = [[z[p] * zb[p], z[p], zb[p], 1] for p in (a, b, c)] + [[z0 * w0, z0, w0, 1]]
    eqs.append(sp.rem(sp.expand(sp.Matrix(rows).det()), Phi, zeta))
G = sp.groebner(eqs + [Phi, rho - 1, z0 * w0 - 1], z0, w0, rho, zeta, order="grevlex", domain="QQ")
print("(b) survivor, rho = 1 and |O| = 1 imposed, no saturation: basis == {1}?", list(G.exprs) == [1], "(must be False)")
assert list(G.exprs) != [1]
print("controls passed")
