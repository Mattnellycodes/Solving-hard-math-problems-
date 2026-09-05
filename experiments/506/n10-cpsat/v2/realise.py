"""Numerical realisation of an abstract structure (F, L) on 10 points: find real points p0 = (0,0), p1 = (1,0)
(similarity normalisation), p2..p9 free, such that every block of F is concyclic-or-collinear and every line of L
is collinear.  Random-start least squares (scipy); a hit is accepted only if the residual is ~ 0, the points are
distinct, and the exact incidence structure recomputed from the coordinates CONTAINS (F, L).  Own code."""
import itertools, math
import numpy as np
from scipy.optimize import least_squares


def equations(F, L):
    eqs = []
    for B in F:
        B = sorted(B)
        for d in B[3:]:
            eqs.append(('c', B[0], B[1], B[2], d))
    for S in L:
        S = sorted(S)
        for c in S[2:]:
            eqs.append(('l', S[0], S[1], c))
    return eqs


def residuals(v, eqs):
    P = np.zeros((10, 2)); P[1, 0] = 1.0; P[2:] = v.reshape(8, 2)
    out = np.empty(len(eqs))
    for i, e in enumerate(eqs):
        if e[0] == 'c':
            a, b, c, d = e[1:]
            M = np.array([[P[k, 0] ** 2 + P[k, 1] ** 2, P[k, 0], P[k, 1], 1.0] for k in (a, b, c, d)])
            out[i] = np.linalg.det(M)
        else:
            a, b, c = e[1:]
            out[i] = (P[b, 0] - P[a, 0]) * (P[c, 1] - P[a, 1]) - (P[b, 1] - P[a, 1]) * (P[c, 0] - P[a, 0])
    return out


def structure_of(P, tol=1e-8):
    """rich blocks (size >= 4) and lines (size >= 3) of the point set P (numerical, tolerance tol on
    normalised incidence residuals)."""
    n = len(P)
    keys = {}
    for i, j, k in itertools.combinations(range(n), 3):
        a, b, c = P[i], P[j], P[k]
        d = 2 * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]))
        scale = max(np.linalg.norm(b - a), np.linalg.norm(c - a), np.linalg.norm(c - b))
        if abs(d) < tol * scale ** 2:
            A = b[1] - a[1]; Bc = a[0] - b[0]; nr = math.hypot(A, Bc); A, Bc = A / nr, Bc / nr
            C = A * a[0] + Bc * a[1]
            if A < -1e-9 or (abs(A) < 1e-9 and Bc < 0): A, Bc, C = -A, -Bc, -C
            key = ('L', round(A, 6), round(Bc, 6), round(C, 6))
        else:
            a2 = a @ a; b2 = b @ b; c2 = c @ c
            ux = (a2 * (b[1] - c[1]) + b2 * (c[1] - a[1]) + c2 * (a[1] - b[1])) / d
            uy = (a2 * (c[0] - b[0]) + b2 * (a[0] - c[0]) + c2 * (b[0] - a[0])) / d
            r = math.hypot(a[0] - ux, a[1] - uy)
            key = ('C', round(ux, 6), round(uy, 6), round(r, 6))
        keys.setdefault(key, set()).update((i, j, k))
    Fr = [frozenset(s) for kk, s in keys.items() if len(s) >= 4]
    Lr = [frozenset(s) for kk, s in keys.items() if kk[0] == 'L']
    return Fr, Lr


def try_realise(F, L, tries=300, seed=0, verbose=False, scale_choices=(0.5, 1.0, 2.0)):
    F = [frozenset(B) for B in F]; L = [frozenset(S) for S in L]
    eqs = equations(F, L)
    rng = np.random.default_rng(seed)
    best = None
    for t in range(tries):
        v0 = rng.normal(size=16) * rng.choice(scale_choices)
        try:
            res = least_squares(residuals, v0, args=(eqs,), method='lm' if len(eqs) >= 16 else 'trf',
                                xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=5000)
        except Exception:
            continue
        P = np.zeros((10, 2)); P[1, 0] = 1.0; P[2:] = res.x.reshape(8, 2)
        mind = min(np.linalg.norm(P[i] - P[j]) for i, j in itertools.combinations(range(10), 2))
        maxd = max(np.linalg.norm(P[i] - P[j]) for i, j in itertools.combinations(range(10), 2))
        r = np.max(np.abs(res.fun)) / maxd ** 4
        if best is None or (mind > 1e-4 and r < best[0]):
            best = (r, P, mind)
        if r < 1e-12 and mind > 1e-4 * maxd:
            Fr, Lr = structure_of(P)
            if any(len(B) == 10 for B in Fr):
                continue
            contains = all(any(B <= B2 for B2 in Fr) for B in F) and all(any(S <= S2 for S2 in Lr) for S in L)
            if contains:
                exact = (sorted(sorted(B) for B in Fr) == sorted(sorted(B) for B in F)) and (sorted(sorted(S) for S in Lr) == sorted(sorted(S) for S in L))
                if verbose:
                    print(f"   realised at try {t}: residual {r:.1e}, exact structure: {exact}; extra blocks {[sorted(B) for B in Fr if not any(B <= B0 for B0 in F)]}")
                return {'realised': True, 'exact': exact, 'points': P.tolist(), 'Fr': [sorted(B) for B in Fr], 'Lr': [sorted(S) for S in Lr]}
    return {'realised': False, 'best_residual': None if best is None else float(best[0]), 'min_dist': None if best is None else float(best[2])}


if __name__ == '__main__':
    F = [list(range(1, 10))]; L = [[0, 1, 2], [0, 3, 4], [0, 5, 6], [0, 7, 8]]
    print("antipodal structure:", try_realise(F, L, tries=50, verbose=True)['realised'])
