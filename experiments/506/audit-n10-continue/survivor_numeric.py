"""Independent numerical sanity check of the survivor in the concentric normal form: unknowns
(rho, O=(x,y)); residual = sum over the 10 three-point lines of the normalised concyclicity
determinant of (a, b, c, O).  Random restarts; report the best residual with rho away from 0, +-1,
and separately the O = infinity residual (collinearity) as a function of rho on a fine grid."""
import json, random, math
import numpy as np
from scipy.optimize import least_squares
from lib import mask, bits
S, R = list(range(5)), list(range(5, 10))
rec = json.load(open("runs/survivors.json"))[0]
F = [mask(b) for b in rec["F"]]; L = [list(l) for l in rec["L"]]
labs = [([0, 2, 8, 4, 6], [0, 2, 8, 4, 6]), ([0, 2, 8, 4, 6], [5, 7, 3, 9, 1]), ([0, 4, 6, 8, 2], [0, 4, 6, 8, 2]), ([0, 4, 6, 8, 2], [5, 9, 1, 3, 7]),
        ([0, 6, 4, 2, 8], [0, 6, 4, 2, 8]), ([0, 6, 4, 2, 8], [5, 1, 9, 7, 3]), ([0, 8, 2, 6, 4], [0, 8, 2, 6, 4]), ([0, 8, 2, 6, 4], [5, 3, 7, 1, 9])]


def pts(lab, rho):
    a = lab[0] + lab[1]
    return [np.array([math.cos(a[p] * math.pi / 5), math.sin(a[p] * math.pi / 5)]) * (1 if p < 5 else rho) for p in range(10)]


def cocirc(A, B, C, O):
    # scale-invariant residual: imaginary part of the cross ratio (a,b;c,O); it vanishes iff the four
    # points are concyclic or collinear and does NOT tend to 0 as O -> infinity
    a, b, c, o = (complex(P[0], P[1]) for P in (A, B, C, O))
    cr = (a - c) * (b - o) / ((a - o) * (b - c))
    return cr.imag


def res(v, lab):
    rho, x, y = v
    P = pts(lab, rho); O = np.array([x, y])
    return np.array([cocirc(P[a], P[b], P[c], O) for a, b, c in L])


random.seed(0); np.random.seed(0)
for lab in labs:
    best = (1e9, None)
    for trial in range(300):
        v0 = [random.choice([random.uniform(0.05, 0.9), random.uniform(1.1, 6)]), random.uniform(-4, 4), random.uniform(-4, 4)]
        r = least_squares(res, v0, args=(lab,), max_nfev=400)
        rho = r.x[0]
        if abs(rho) < 0.03 or abs(abs(rho) - 1) < 0.03:
            continue
        val = float(np.sum(r.fun ** 2))
        if val < best[0]:
            best = (val, r.x)
    # O = infinity: collinearity residual on a grid of rho
    def coll(rho):
        P = pts(lab, rho)
        return sum(abs(np.linalg.det(np.array([[P[p][0], P[p][1], 1.0] for p in tr]))) for tr in L)
    grid = [rho for rho in np.concatenate([np.linspace(0.05, 0.95, 400), np.linspace(1.05, 30, 3000)])]
    cm = min((coll(rho), rho) for rho in grid)
    print(f"labelling S={lab[0]} R={lab[1]}: best finite-O residual (rho away from 0,+-1) = {best[0]:.3e} at (rho,x,y)={best[1]}; "
          f"min O=inf collinearity sum on grid = {cm[0]:.3e} at rho={cm[1]:.3f}")
