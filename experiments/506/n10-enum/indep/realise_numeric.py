"""Numerical realisability probe (evidence only, not a proof) for a structure (F, L) on n points:
find planar points p_0..p_{n-1} (p_0 = (0,0), p_1 = (1,0) by similarity) such that every block of F
not in L is concyclic, every line of L is collinear, and points are pairwise distinct.
Least-squares from random starts (scipy).  Reports the best residual found."""
import sys, json, math, random
import numpy as np
from scipy.optimize import least_squares


def residuals_factory(n, blocks, lines, eps=0.05, wpen=3.0):
    circ_eqs = []   # (a,b,c,d) concyclic
    for B in blocks:
        if any(set(B) <= set(L) for L in lines):
            continue
        a, b, c = B[0], B[1], B[2]
        for d in B[3:]:
            circ_eqs.append((a, b, c, d))
    col_eqs = []    # (a,b,c) collinear
    for L in lines:
        a, b = L[0], L[1]
        for c in L[2:]:
            col_eqs.append((a, b, c))
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]

    def unpack(x):
        P = np.zeros((n, 2)); P[1, 0] = 1.0
        P[2:] = x.reshape(n - 2, 2)
        return P

    def res(x):
        P = unpack(x)
        out = []
        for a, b, c, d in circ_eqs:
            M = np.array([[P[i, 0], P[i, 1], P[i, 0] ** 2 + P[i, 1] ** 2, 1.0] for i in (a, b, c, d)])
            out.append(np.linalg.det(M))
        for a, b, c in col_eqs:
            out.append((P[b, 0] - P[a, 0]) * (P[c, 1] - P[a, 1]) - (P[b, 1] - P[a, 1]) * (P[c, 0] - P[a, 0]))
        for i, j in pairs:
            d = np.hypot(P[i, 0] - P[j, 0], P[i, 1] - P[j, 1])
            out.append(wpen * max(0.0, eps - d))
        # exactness: a point outside a block must not lie on the block's circle (Moebius-invariant
        # cross-ratio test), a point outside a line must not lie on it
        Z = P[:, 0] + 1j * P[:, 1]
        for B in blocks:
            a, b, c = B[0], B[1], B[2]
            for p in range(n):
                if p in B:
                    continue
                den = (Z[p] - Z[c]) * (Z[b] - Z[a])
                if abs(den) < 1e-12:
                    out.append(wpen); continue
                cr = ((Z[p] - Z[a]) * (Z[b] - Z[c])) / den
                out.append(wpen * max(0.0, 0.03 - abs(cr.imag) / (1.0 + abs(cr))))
        for L in lines:
            a, b = L[0], L[1]
            for p in range(n):
                if p in L:
                    continue
                col = (P[b, 0] - P[a, 0]) * (P[p, 1] - P[a, 1]) - (P[b, 1] - P[a, 1]) * (P[p, 0] - P[a, 0])
                scale = np.hypot(P[b, 0] - P[a, 0], P[b, 1] - P[a, 1]) * np.hypot(P[p, 0] - P[a, 0], P[p, 1] - P[a, 1]) + 1e-12
                out.append(wpen * max(0.0, 0.03 - abs(col) / scale))
        # non-degeneracy: the first three points of every circle-block must not be collinear
        for a, b, c, d in circ_eqs:
            col = (P[b, 0] - P[a, 0]) * (P[c, 1] - P[a, 1]) - (P[b, 1] - P[a, 1]) * (P[c, 0] - P[a, 0])
            out.append(wpen * max(0.0, 0.02 - abs(col)))
        # keep points bounded
        for i in range(2, n):
            r = np.hypot(P[i, 0], P[i, 1])
            out.append(0.1 * max(0.0, r - 6.0))
        return np.array(out)
    return res, unpack, len(circ_eqs), len(col_eqs)


def probe(n, blocks, lines, starts=200, seed=0, verbose=False):
    rng = random.Random(seed)
    res, unpack, nc, nl = residuals_factory(n, blocks, lines)
    best = (math.inf, None)
    for s in range(starts):
        x0 = np.array([rng.uniform(-2, 2) for _ in range(2 * (n - 2))])
        try:
            sol = least_squares(res, x0, method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
        except Exception as e:
            continue
        r = float(np.sum(sol.fun ** 2))
        if r < best[0]:
            best = (r, unpack(sol.x))
            if verbose:
                print(f"  start {s}: residual {r:.3e}", flush=True)
    return best, nc, nl


if __name__ == "__main__":
    spec = json.load(open(sys.argv[1]))
    starts = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    (r, P), nc, nl = probe(spec['n'], spec['blocks'], spec['lines'], starts=starts, verbose=True)
    print(f"equations: {nc} concyclicity + {nl} collinearity; best residual {r:.3e}")
    if P is not None:
        print("points:", [(round(float(x), 6), round(float(y), 6)) for x, y in P])
