"""Independent numerical sanity check of the exact verdict for the theorem-only survivor (concentric regular
pentagons + 10 three-point lines): search for rho > 0 and a point O (or O = infinity) lying on the circles of all
10 lines by random-start least squares; report the best residual (should stay far from 0)."""
import json, itertools, cmath, math, sys
import numpy as np
from scipy.optimize import least_squares
sys.path.insert(0, '/home/user/Solving-hard-math-problems-/experiments/506/n10-continue')
post = json.load(open('runs/post_7.json')); rec = post[4]
F = [frozenset(B) for B in rec['blocks']]; L = [frozenset(S) for S in rec['lines']]
fives = sorted([B for B in F if len(B) == 5], key=sorted); S, R = sorted(fives[0]), sorted(fives[1])
idx = {p: i for i, p in enumerate(S + R)}
labellings = [(0, 2, 6, 4, 8, 1, 7, 3, 9, 5), (0, 2, 6, 4, 8, 6, 2, 8, 4, 0), (0, 4, 2, 8, 6, 2, 4, 6, 8, 0), (0, 4, 2, 8, 6, 7, 9, 1, 3, 5),
              (0, 6, 8, 2, 4, 3, 1, 9, 7, 5), (0, 6, 8, 2, 4, 8, 6, 4, 2, 0), (0, 8, 4, 6, 2, 4, 8, 2, 6, 0), (0, 8, 4, 6, 2, 9, 3, 7, 1, 5)]
def pts(lab, rho):
    z = {}
    for p in S: z[p] = cmath.exp(2j * math.pi * lab[idx[p]] / 10)
    for p in R: z[p] = rho * cmath.exp(2j * math.pi * lab[idx[p]] / 10)
    return z
def concyc(a, b, c, d):
    # normalised imaginary part of the cross-ratio (a,b;c,d)
    cr = ((a - c) * (b - d)) / ((a - d) * (b - c))
    return cr.imag / (abs(cr) + 1)
def resid(v, lab, infinite):
    rho = v[0]
    z = pts(lab, rho)
    out = []
    for Sl in L:
        p, q, r = sorted(Sl)[:3]
        if infinite:
            d = (z[q] - z[p]) * (z[r] - z[p]).conjugate(); out.append(d.imag / (abs(d) + 1))
        else:
            O = complex(v[1], v[2]); out.append(concyc(O, z[p], z[q], z[r]))
    return np.array(out)
rng = np.random.default_rng(1)
for k, lab in enumerate(labellings):
    z = pts(lab, 2.0)
    assert all(abs(concyc(*[z[p] for p in sorted(B)])) < 1e-12 for B in F if len(B) == 4)   # 4-blocks concyclic
    best_f = min(least_squares(resid, [rng.uniform(0.05, 5), *rng.normal(size=2) * 3], args=(lab, False)).cost for _ in range(300))
    best_i = min(least_squares(resid, [rng.uniform(0.05, 5)], args=(lab, True)).cost for _ in range(100))
    print(f"labelling {k}: best residual (finite O) {math.sqrt(2*best_f):.3e}, (O = infinity) {math.sqrt(2*best_i):.3e}")
