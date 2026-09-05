"""Numerical realisation attempt for an abstract structure (F, L) on 10 points: find real points with
p0 = (0,0), p1 = (1,0) (a similarity normalisation; similarities are the Möbius maps fixing infinity)
such that every block of F is concyclic (or collinear) and every line of L is collinear.
Least squares from random starts (scipy).  A solution is accepted only if the residual is ~0, the
points are distinct, and the incidence structure computed from the coordinates contains (F, L).
"""
import itertools, math, sys, json
import numpy as np
from scipy.optimize import least_squares


def equations(F, L):
    eqs = []   # list of ('c', a, b, c, d) or ('l', a, b, c)
    for B in F:
        B = sorted(B); a, b, c = B[:3]
        for d in B[3:]:
            eqs.append(('c', a, b, c, d))
    for S in L:
        S = sorted(S); a, b = S[:2]
        for c in S[2:]:
            eqs.append(('l', a, b, c))
    return eqs


def residuals(v, eqs):
    P = np.zeros((10, 2)); P[1, 0] = 1.0
    P[2:] = v.reshape(8, 2)
    out = []
    for e in eqs:
        if e[0] == 'c':
            a, b, c, d = e[1:]
            M = np.array([[P[i, 0] ** 2 + P[i, 1] ** 2, P[i, 0], P[i, 1], 1.0] for i in (a, b, c, d)])
            out.append(np.linalg.det(M))
        else:
            a, b, c = e[1:]
            out.append((P[b, 0] - P[a, 0]) * (P[c, 1] - P[a, 1]) - (P[b, 1] - P[a, 1]) * (P[c, 0] - P[a, 0]))
    return np.array(out)


def structure_of(P, tol=1e-7):
    """blocks (circles or lines) through >= 3 points, computed with tolerance; returns (F, L) with F all
    blocks of size >= 4 and L all lines of size >= 3, plus a flag if two points nearly coincide."""
    n = len(P)
    coincide = any(np.hypot(*(P[i] - P[j])) < 1e-6 for i, j in itertools.combinations(range(n), 2))
    blocks = {}
    for i, j, k in itertools.combinations(range(n), 3):
        a, b, c = P[i], P[j], P[k]
        d = 2 * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]))
        if abs(d) < tol:
            A = b[1] - a[1]; Bc = a[0] - b[0]; C = A * a[0] + Bc * a[1]
            nr = math.hypot(A, Bc); A, Bc, C = A / nr, Bc / nr, C / nr
            if A < -1e-9 or (abs(A) < 1e-9 and Bc < 0):
                A, Bc, C = -A, -Bc, -C
            key = ('L', round(A, 5), round(Bc, 5), round(C, 5))
        else:
            a2 = a @ a; b2 = b @ b; c2 = c @ c
            ux = (a2 * (b[1] - c[1]) + b2 * (c[1] - a[1]) + c2 * (a[1] - b[1])) / d
            uy = (a2 * (c[0] - b[0]) + b2 * (a[0] - c[0]) + c2 * (b[0] - a[0])) / d
            r = math.hypot(a[0] - ux, a[1] - uy)
            key = ('C', round(ux, 5), round(uy, 5), round(r, 5))
        blocks.setdefault(key, set()).update((i, j, k))
    Fr = [frozenset(s) for k, s in blocks.items() if len(s) >= 4]
    Lr = [frozenset(s) for k, s in blocks.items() if k[0] == 'L']
    return Fr, Lr, coincide


def try_realise(F, L, tries=200, seed=0, verbose=False):
    F = [frozenset(B) for B in F]; L = [frozenset(S) for S in L]
    eqs = equations(F, L)
    rng = np.random.default_rng(seed)
    best = None
    for t in range(tries):
        v0 = rng.normal(size=16) * rng.choice([0.5, 1.0, 2.0])
        try:
            res = least_squares(residuals, v0, args=(eqs,), method=('lm' if len(eqs) >= 16 else 'trf'), xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=4000)
        except Exception:
            continue
        r = np.max(np.abs(res.fun))
        P = np.zeros((10, 2)); P[1, 0] = 1.0; P[2:] = res.x.reshape(8, 2)
        Fr, Lr, coincide = structure_of(P)
        if best is None or r < best[0]:
            best = (r, P, coincide)
        if r < 1e-9 and not coincide and not any(len(B) == 10 for B in Fr):
            contains = all(any(B <= B2 for B2 in Fr) for B in F) and all(any(S <= S2 for S2 in Lr) for S in L)
            if contains:
                extra = (sorted(sorted(B) for B in Fr) != sorted(sorted(B) for B in F)) or (sorted(sorted(S) for S in Lr) != sorted(sorted(S) for S in L))
                if verbose:
                    print(f"  REALISED at try {t}: residual {r:.2e}, extra incidences: {extra}")
                return {'realised': True, 'points': P.tolist(), 'extra': extra, 'Fr': [sorted(B) for B in Fr], 'Lr': [sorted(S) for S in Lr]}
    return {'realised': False, 'best_residual': None if best is None else float(best[0]), 'coincide': None if best is None else bool(best[2])}


if __name__ == '__main__':
    # sanity test: the 8-point record configuration (17 circles) padded is not 10 points; instead test the
    # antipodal 10-point structure: 9 points on a circle in antipodal pairs (+ one) and the centre.
    F = [list(range(1, 10))]
    L = [[0, 1, 2], [0, 3, 4], [0, 5, 6], [0, 7, 8]]
    print(try_realise(F, L, tries=50, verbose=True)['realised'])
