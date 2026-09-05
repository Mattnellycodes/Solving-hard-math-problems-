"""Numerical realisation attempts (Levenberg-Marquardt, random restarts) for abstract structures.
(a) A 'Mobius structure' on n points in R^2: each block must be concyclic-or-collinear, encoded by
    the 4x4 determinants det[x^2+y^2, x, y, 1] of all quadruples of the block (all must vanish).
    Three points are fixed by a Mobius transformation ((0,0),(1,0),(0,1)); the rest are unknowns.
(b) A 'line structure' on m points in R^2: each line = collinear triple determinants vanish;
    four points fixed as an affine frame (only when no three of them are on a structure line).
Distinctness is enforced by a barrier term  eps / |p-q|^2  (scaled), reported separately.
A residual ~1e-12 with all points distinct = numerical realisation (would refute the report)."""
import numpy as np, itertools, sys
from scipy.optimize import least_squares
rng = np.random.default_rng(12345)

def mobius_residuals(X, blocks, fixed, n):
    P = np.zeros((n, 2)); P[list(fixed.keys())] = list(fixed.values())
    free = [i for i in range(n) if i not in fixed]
    P[free] = X.reshape(-1, 2)
    res = []
    for b in blocks:
        for q in itertools.combinations(sorted(b), 4):
            M = np.array([[P[i, 0] ** 2 + P[i, 1] ** 2, P[i, 0], P[i, 1], 1.0] for i in q])
            res.append(np.linalg.det(M))
    # separation barrier (soft)
    for i, j in itertools.combinations(range(n), 2):
        d2 = np.sum((P[i] - P[j]) ** 2)
        res.append(1e-3 / (d2 + 1e-9))
    return np.array(res)

def line_residuals(X, lines, fixed, m):
    P = np.zeros((m, 2)); P[list(fixed.keys())] = list(fixed.values())
    free = [i for i in range(m) if i not in fixed]
    P[free] = X.reshape(-1, 2)
    res = []
    for l in lines:
        for t in itertools.combinations(sorted(l), 3):
            M = np.array([[P[i, 0], P[i, 1], 1.0] for i in t])
            res.append(np.linalg.det(M))
    for i, j in itertools.combinations(range(m), 2):
        d2 = np.sum((P[i] - P[j]) ** 2)
        res.append(1e-3 / (d2 + 1e-9))
    return np.array(res)

def attempt(resfun, nfree, restarts=200, scale=3.0):
    best = None
    for r in range(restarts):
        x0 = rng.normal(0, scale, 2 * nfree)
        sol = least_squares(resfun, x0, method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
        f = sol.fun
        inc = np.sum(f[:-(len(f) - (len(f) - 0))] ** 2) if False else None
        if best is None or sol.cost < best.cost: best = sol
    return best

def report(name, sol, nres_incidence):
    f = sol.fun
    inc = np.sqrt(np.sum(f[:nres_incidence] ** 2)); bar = np.sqrt(np.sum(f[nres_incidence:] ** 2))
    print(f"{name}: best total cost={sol.cost:.3e}  incidence-residual={inc:.3e}  barrier={bar:.3e}")
    return inc

if __name__ == "__main__":
    # ---- positive control 1: the n=8 cube structure (SQS(8) minus a parallel class) is realisable
    SQS8 = [frozenset(b) for b in [(1,2,5,6),(3,4,5,6),(1,3,5,7),(2,4,5,7),(2,3,6,7),(1,4,6,7),(2,3,5,8),(1,4,5,8),
             (1,3,6,8),(2,4,6,8),(1,2,7,8),(3,4,7,8),(1,2,3,4),(5,6,7,8)]]
    cube = [frozenset(i - 1 for i in b) for b in SQS8[:12]]
    fixed = {0: (0.0, 0.0), 1: (1.0, 0.0), 2: (0.0, 1.0)}
    ninc = sum(1 for b in cube for _ in itertools.combinations(b, 4))
    sol = attempt(lambda X: mobius_residuals(X, cube, fixed, 8), 5, restarts=60)
    report("control: cube structure n=8 (realisable)", sol, ninc)
    # ---- positive control 2: Pappus 9_3 line structure is realisable
    pappus = [frozenset(l) for l in [(0,1,2),(3,4,5),(6,7,8),(0,4,8),(0,5,7),(1,3,8),(1,5,6),(2,3,7),(2,4,6)]]
    fixedL = {0: (0.0, 0.0), 1: (1.0, 0.0), 3: (0.0, 1.0), 4: (1.0, 1.0)}
    assert not any(len(set(fixedL) & l) >= 3 for l in pappus)
    ninc = sum(1 for l in pappus for _ in itertools.combinations(l, 3))
    sol = attempt(lambda X: line_residuals(X, pappus, fixedL, 9), 5, restarts=60)
    report("control: Pappus 9_3 lines (realisable)", sol, ninc)
    # ---- target 1: Mobius-Kantor (8_3) lines {i,i+1,i+3}
    MK = [frozenset({i, (i + 1) % 8, (i + 3) % 8}) for i in range(8)]
    fixedMK = {0: (0.0, 0.0), 1: (1.0, 0.0), 2: (0.0, 1.0), 5: (1.0, 1.0)}
    assert not any(len(set(fixedMK) & l) >= 3 for l in MK)
    ninc = sum(1 for l in MK for _ in itertools.combinations(l, 3))
    sol = attempt(lambda X: line_residuals(X, MK, fixedMK, 8), 4, restarts=400)
    report("target: Mobius-Kantor (8_3) lines", sol, ninc)
    # ---- target 2: the 9-point candidate structure from the enumeration (two 5-blocks + 12 four-blocks)
    cand = [frozenset(b) for b in [(0,1,2,3,4),(0,5,6,7,8)]] + [frozenset(b) for b in SQS8[:12]]
    ninc = sum(1 for b in cand for _ in itertools.combinations(b, 4))
    sol = attempt(lambda X: mobius_residuals(X, cand, fixed, 9), 6, restarts=400)
    inc = report("target: 9-point candidate (2x5-block + 12 four-blocks)", sol, ninc)
    # ---- target 3: same, but only the blocks through one degree-7 point (the Fano-derived part) -- should fail
    p = 1
    sub = [b for b in cand if p in b]
    ninc = sum(1 for b in sub for _ in itertools.combinations(b, 4))
    sol = attempt(lambda X: mobius_residuals(X, sub, fixed, 9), 6, restarts=200)
    report("target: only the 7 blocks through point 1 (Fano-derived)", sol, ninc)
    # ---- target 4: candidate minus the blocks through point 1 (is the rest realisable?)
    sub2 = [b for b in cand if p not in b]
    ninc = sum(1 for b in sub2 for _ in itertools.combinations(b, 4))
    sol = attempt(lambda X: mobius_residuals(X, sub2, fixed, 9), 6, restarts=100)
    report("info: candidate minus blocks through point 1", sol, ninc)
