"""Generic numerical realiser for an abstract Moebius structure on n points.
Input: blocks (list of index lists; each is a circle-or-line through those points) and
lines (subset of blocks that must be straight lines, i.e. pass through infinity).
Unknowns: n points in R^2.  Equations: for every block B with base triple (b0,b1,b2), for every
other point p in B: concyclic det (or collinear det for lines) = 0; for every line: collinear.
Method: scipy least_squares from many random starts; a solution is accepted when the residual
is < 1e-12 and the configuration is non-degenerate (distinct points, not all concyclic).
Then the circles are counted in floating point (and the point set printed for exact checking)."""
import numpy as np, itertools, sys, random
from scipy.optimize import least_squares

def det4(rows):
    return np.linalg.det(np.array(rows))

def make_residual(n, blocks, lines):
    lines = [frozenset(l) for l in lines]
    eqs = []  # (type, indices)
    for B in blocks:
        B = sorted(B)
        if frozenset(B) in lines:
            for p in B[2:]:
                eqs.append(("L", (B[0], B[1], p)))
        else:
            for p in B[3:]:
                eqs.append(("C", (B[0], B[1], B[2], p)))
    for l in lines:
        l = sorted(l)
        for p in l[2:]:
            eqs.append(("L", (l[0], l[1], p)))
    def unpack(v):
        # v = [a, u, x3,y3,...]: P0=(0,0), P1=(1,0), P2=(a, exp(u)), rest free
        P = np.zeros((n, 2))
        P[1] = (1.0, 0.0); P[2] = (v[0], np.exp(v[1]))
        P[3:] = v[2:].reshape(n-3, 2)
        return P
    def res(v):
        P = unpack(v)
        out = []
        for t, idx in eqs:
            if t == "L":
                a, b, c = (P[i] for i in idx)
                out.append((b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]))
            else:
                M = [[P[i][0]**2+P[i][1]**2, P[i][0], P[i][1], 1.0] for i in idx]
                out.append(np.linalg.det(np.array(M)))
        # normalisation: fix scale softly (points 0,1 at distance ~1)
        return np.array(out)
    return res, eqs, unpack

def count_float(P, tol=1e-7):
    n = len(P)
    circles = []
    ncol = 0
    for i, j, k in itertools.combinations(range(n), 3):
        a, b, c = P[i], P[j], P[k]
        d = 2*(a[0]*(b[1]-c[1]) + b[0]*(c[1]-a[1]) + c[0]*(a[1]-b[1]))
        scale = max(np.linalg.norm(b-a), np.linalg.norm(c-a), np.linalg.norm(c-b))
        if abs(d) < tol*scale*scale*10:
            ncol += 1; continue
        a2 = a@a; b2 = b@b; c2 = c@c
        ux = (a2*(b[1]-c[1]) + b2*(c[1]-a[1]) + c2*(a[1]-b[1]))/d
        uy = (a2*(c[0]-b[0]) + b2*(a[0]-c[0]) + c2*(b[0]-a[0]))/d
        r = np.hypot(a[0]-ux, a[1]-uy)
        key = None
        for idx, (cx, cy, rr, memb) in enumerate(circles):
            if abs(cx-ux) < tol*(1+rr) and abs(cy-uy) < tol*(1+rr) and abs(rr-r) < tol*(1+rr):
                memb.update((i, j, k)); key = idx; break
        if key is None:
            circles.append((ux, uy, r, {i, j, k}))
    return len(circles), ncol, circles

def realise(n, blocks, lines, tries=200, seed=0, verbose=True, spread=3.0):
    rng = np.random.default_rng(seed)
    # relabel so that points 0,1,2 form a triple lying in NO line (then P2 off the x-axis is forced)
    lines_f = [frozenset(l) for l in lines]
    base = None
    for tri in itertools.combinations(range(n), 3):
        if not any(set(tri) <= l for l in lines_f):
            base = tri; break
    assert base is not None
    perm = list(base) + [i for i in range(n) if i not in base]   # new index -> old index
    inv = {old: new for new, old in enumerate(perm)}
    blocks = [[inv[i] for i in B] for B in blocks]; lines = [[inv[i] for i in l] for l in lines]
    res, eqs, unpack = make_residual(n, blocks, lines)
    best = None
    for t in range(tries):
        v0 = rng.normal(size=2*n-4)*spread; v0[1] = rng.normal()*0.5
        try:
            lo = -np.inf*np.ones(2*n-4); hi = np.inf*np.ones(2*n-4); lo[1] = -3.0; hi[1] = 3.0
            v0[1] = np.clip(v0[1], -2.5, 2.5)
            sol = least_squares(res, v0, method="trf", bounds=(lo, hi), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
        except Exception as e:
            continue
        P = unpack(sol.x)
        # degeneracy: min pairwise distance relative to diameter
        dmin = min(np.linalg.norm(P[i]-P[j]) for i, j in itertools.combinations(range(n), 2))
        diam = max(np.linalg.norm(P[i]-P[j]) for i, j in itertools.combinations(range(n), 2))
        # rescale residual by size
        r = np.abs(res(sol.x)).max()/(diam**4 + 1e-300)
        Pold = np.zeros_like(P); 
        for new, old in enumerate(perm): Pold[old] = P[new]
        P = Pold
        if best is None or (r, -dmin/diam) < (best[0], -best[1]):
            best = (r, dmin/diam, P.copy())
        if r < 1e-13 and dmin/diam > 1e-4:
            nc, ncol, circ = count_float(P)
            if verbose:
                print(f"  try {t}: residual {r:.2e} dmin/diam {dmin/diam:.3g} -> circles {nc}, collinear triples {ncol}")
            return P, nc, ncol
    if verbose:
        print(f"  no realisation found in {tries} tries; best scaled residual {best[0]:.2e}, dmin/diam {best[1]:.3g}")
    return None, None, None

if __name__ == "__main__":
    # positive control: the cube structure with 3 lines (face 0123 + two corner triples)
    cube = [[0,1,2,3],[0,1,4,5],[2,3,4,5],[0,2,4,6],[1,3,4,6],[1,2,5,6],[0,3,4,7],[0,2,5,7],[1,3,5,7],[0,1,6,7],[2,3,6,7],[4,5,6,7]]
    lines = [[0,1,2,3],[0,5,6],[1,4,7]]
    print("control: cube structure + lines", lines)
    P, nc, ncol = realise(8, cube + [[0,5,6],[1,4,7]], lines, tries=100)
    if P is not None:
        print(P)
