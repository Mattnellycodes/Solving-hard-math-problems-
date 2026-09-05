"""Sanity test of two5_exact on a REALISABLE configuration: two concentric regular pentagons rotated by
36 degrees with radius ratio rho = phi^2 = (3+sqrt5)/2 (then each chord s s' of the unit pentagon
passes through two vertices of the outer pentagon: five 4-point lines).  The O = oo test must
return an ideal containing rho^2 - 3 rho + 1 (i.e. NOT (1)); the finite-O equations must vanish
numerically at a constructed solution."""
import math, cmath, itertools
import sympy as sp
from fractions import Fraction
import two5_exact as T

rho = (3 + math.sqrt(5)) / 2
pts = [cmath.exp(2j * math.pi * k / 5) for k in range(5)] + [rho * cmath.exp(2j * math.pi * (k + 0.5) / 5) for k in range(5)]
# numeric structure extraction (tolerance)
def concyc(a, b, c, d):
    cr = (a - c) * (b - d) / ((a - d) * (b - c))
    return abs(cr.imag) < 1e-9
def collin(a, b, c):
    return abs(((b - a) * (c - a).conjugate()).imag) < 1e-9
F = []
for Q in itertools.combinations(range(10), 4):
    a, b, c, d = [pts[i] for i in Q]
    if concyc(a, b, c, d):
        F.append(list(Q))
# merge to maximal blocks
blocks = set()
for T3 in itertools.combinations(range(10), 3):
    a, b, c = [pts[i] for i in T3]
    blk = set(T3)
    for q in range(10):
        if q not in blk and (concyc(a, b, c, pts[q]) if not collin(a, b, c) else collin(a, b, pts[q])):
            blk.add(q)
    if len(blk) >= 4:
        blocks.add(frozenset(blk))
F = [sorted(b) for b in blocks]
L = [sorted(b) for b in blocks if collin(*[pts[i] for i in list(b)[:3]])]
# 3-point lines
for T3 in itertools.combinations(range(10), 3):
    if collin(*[pts[i] for i in T3]) and not any(set(T3) <= b for b in blocks):
        L.append(list(T3))
print('F sizes', sorted(len(b) for b in F), 'L =', L)
phi = tuple(Fraction(k, 5) for k in range(5)) + tuple((Fraction(2 * k + 1, 10)) for k in range(5))
T.numeric_check(F, phi, rho=rho)
Gf, Gi = T.exact_test(F, L, phi, 10)
print('O = oo basis:', list(Gi.exprs))
print('finite O basis:', list(Gf.exprs))
# finite-O numeric check on one circle through O: take O on the circle through pts 0,1,5
import numpy as np
a, b, c = pts[0], pts[1], pts[5]
# circumcentre
ax, ay = a.real, a.imag; bx, by = b.real, b.imag; cx, cy = c.real, c.imag
d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
ux = ((ax**2 + ay**2) * (by - cy) + (bx**2 + by**2) * (cy - ay) + (cx**2 + cy**2) * (ay - by)) / d
uy = ((ax**2 + ay**2) * (cx - bx) + (bx**2 + by**2) * (ax - cx) + (cx**2 + cy**2) * (bx - ax)) / d
rad = abs(a - complex(ux, uy))
O = complex(ux, uy) + rad * cmath.exp(0.7j)
z, rr, u, v, t = sp.symbols('z rho u v t')
Gf2, _ = T.exact_test(F, [[0, 1, 5]], phi, 10)
# evaluate the single finite-O equation numerically
Nn = 20
zeta = cmath.exp(2j * math.pi / Nn)
eqs = [e for e in Gf2.exprs]
# rebuild the raw equation instead (before groebner) for the numeric check
import sympy
Phi = sp.cyclotomic_poly(Nn, z)
ii = z ** 5
coords = []
for p in range(10):
    k = phi[p]; a_ = int(k * 10) * 2 % Nn
    x = (z ** a_ + z ** ((Nn - a_) % Nn)) / 2; y = -ii * (z ** a_ - z ** ((Nn - a_) % Nn)) / 2
    if p >= 5: x, y = rr * x, rr * y
    coords.append((x, y))
rows = [[1 if p < 5 else rr ** 2, coords[p][0], coords[p][1], 1] for p in (0, 1, 5)] + [[u ** 2 + v ** 2, u, v, 1]]
val = sp.Matrix(rows).det().subs({z: zeta, rr: rho, u: O.real, v: O.imag})
print('finite-O determinant at a constructed solution:', complex(val))
val2 = sp.Matrix(rows).det().subs({z: zeta, rr: rho, u: O.real + 0.3, v: O.imag})
print('... at a non-solution:', complex(val2))
