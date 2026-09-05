"""Numerical realisation attempts for abstract Mobius structures, written from scratch.
Input structure on n points: blocks (size>=4, concyclic-or-collinear) and lines (size>=3).
We work in the Mobius picture: Q = P + {infinity} (index n).  Every block of Q must be concyclic
(lines are blocks through n).  Normalisation: three points a,b,c of Q that form a 3-block of Q
(no block of size>=4 contains them) are fixed at (0,0),(1,0),(0,1) using the Mobius group.
Non-degeneracy: an extra variable t with t*det(a,b,c,d) - 1 = 0 for one further point d, so the
all-concyclic collapse is impossible; coincident points are rejected afterwards.
Residuals: for every block B, det[x,y,x^2+y^2,1] over (B[0],B[1],B[2],B[k]), k>=3.
Success: residual < 1e-10, min pairwise distance > 1e-4; then the Euclidean circle count of
P is |blocks(Q)| - deg(n) computed numerically and reported.
"""
import numpy as np, sys, json, itertools
from scipy.optimize import least_squares

def cdet(a, b, c, d):
    M = np.array([[p[0], p[1], p[0]*p[0] + p[1]*p[1], 1.0] for p in (a, b, c, d)])
    return np.linalg.det(M)

def setup(n, blocks, lines):
    N = n + 1
    Qblocks = [tuple(B) for B in blocks] + [tuple(L) + (n,) for L in lines]
    # choose a 3-block triple of Q as frame
    frame = None
    for t in itertools.combinations(range(N), 3):
        if not any(set(t) <= set(B) for B in Qblocks):
            frame = t; break
    assert frame is not None
    d = [p for p in range(N) if p not in frame][0]
    free = [p for p in range(N) if p not in frame]
    return N, Qblocks, frame, d, free

def make_residual(N, Qblocks, frame, d, free):
    fixed = {frame[0]: (0.0, 0.0), frame[1]: (1.0, 0.0), frame[2]: (0.0, 1.0)}
    pairs = np.array(list(itertools.combinations(range(N), 2)))
    nfree = 2 * len(free)
    quads = np.array([(B[0], B[1], B[2], k) for B in Qblocks for k in B[3:]] + [(frame[0], frame[1], frame[2], d)])
    base = np.zeros((N, 2))
    for p, xy in fixed.items(): base[p] = xy
    def unpack(v):
        P = base.copy()
        P[free] = v[:nfree].reshape(-1, 2)
        return P, v[nfree], v[nfree + 1:]
    def dets(P):
        X = P[quads]                       # (Q,4,2)
        M = np.concatenate([X, (X**2).sum(axis=2, keepdims=True), np.ones((len(quads), 4, 1))], axis=2)
        return np.linalg.det(M)
    def res(v):
        P, t, tp = unpack(v)
        dd = dets(P)
        dist2 = ((P[pairs[:, 0]] - P[pairs[:, 1]])**2).sum(axis=1)
        return np.concatenate([dd[:-1], [t * dd[-1] - 1.0], tp * dist2 - 1.0])
    def init(rng):
        v = np.zeros(nfree + 1 + len(pairs))
        v[:nfree] = rng.normal(size=nfree) * rng.choice([0.3, 1.0, 3.0])
        P, _, _ = unpack(v)
        v[nfree] = 1.0 / dets(P)[-1]
        v[nfree + 1:] = 1.0 / ((P[pairs[:, 0]] - P[pairs[:, 1]])**2).sum(axis=1)
        return v
    return res, unpack, init

def blocks_float(P, tol=1e-6):
    """all blocks (circles incl. lines, via stereographic lifting: use circle centres, lines flagged) of a float set."""
    N = len(P); cl = []
    for i, j, k in itertools.combinations(range(N), 3):
        a, b, c = P[i], P[j], P[k]
        det = (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])
        if abs(det) < tol * (np.linalg.norm(b-a) * np.linalg.norm(c-a)):
            key = ('L', None)
        else:
            s = [np.dot(p, p) for p in (a, b, c)]; D = 2*det
            ux = ((s[1]-s[0])*(c[1]-a[1]) - (s[2]-s[0])*(b[1]-a[1])) / D
            uy = ((b[0]-a[0])*(s[2]-s[0]) - (c[0]-a[0])*(s[1]-s[0])) / D
            key = ('C', (ux, uy))
        for entry in cl:
            if key[0] == 'L' and entry[0] == 'L' and len(entry[1] & {i, j, k}) >= 2:
                entry[1].update((i, j, k)); break
            if key[0] == 'C' and entry[0] == 'C' and abs(entry[2][0]-ux) < 1e-6*(1+abs(ux)) and abs(entry[2][1]-uy) < 1e-6*(1+abs(uy)):
                entry[1].update((i, j, k)); break
        else:
            cl.append([key[0], {i, j, k}, key[1]])
    return [(e[0], sorted(e[1])) for e in cl]

def attempt(n, blocks, lines, tries=200, seed=0, verbose=True, maxnfev=600):
    rng = np.random.default_rng(seed)
    N, Qblocks, frame, d, free = setup(n, blocks, lines)
    res, unpack, init = make_residual(N, Qblocks, frame, d, free)
    nvar = 2*len(free) + 1 + N*(N-1)//2
    best = (np.inf, 0, None)
    for tr in range(tries):
        v0 = init(rng)
        nres = len(res(v0)); meth = 'lm' if nres >= nvar else 'trf'
        try:
            r = least_squares(res, v0, method=meth, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=maxnfev)
        except Exception as e:
            print('   exception', e); continue
        P, t, _ = unpack(r.x)
        dmin = min(np.linalg.norm(P[i]-P[j]) for i, j in itertools.combinations(range(N), 2))
        val = np.max(np.abs(r.fun))
        if dmin > 1e-4 and val < best[0]:
            best = (val, dmin, P.copy())
        if val < 1e-10 and dmin > 1e-4:
            bl = blocks_float(P)
            deg_inf = sum(1 for kind, s in bl if n in s)
            euclid = len(bl) - deg_inf
            if verbose:
                print(f'  try {tr}: residual {val:.2e} dmin {dmin:.3g} -> |B(Q)|={len(bl)} deg(inf)={deg_inf} EUCLIDEAN CIRCLES = {euclid}')
            return dict(success=True, residual=float(val), dmin=float(dmin), Q=P.tolist(), circles=euclid,
                        blocks_found=[(k, s) for k, s in bl if len(s) >= 4 or k == 'L'])
    if verbose:
        print(f'  no realisation in {tries} tries; best residual {best[0]:.2e} (dmin {best[1]:.3g})')
    return dict(success=False, residual=float(best[0]), dmin=float(best[1]), Q=None if best[2] is None else best[2].tolist())

if __name__ == '__main__':
    recs = json.load(open(sys.argv[1]))
    tries = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    for i, r in enumerate(recs):
        n = r['n'] if 'n' in r else max(max(b) for b in r['blocks'] + r['lines']) + 1
        lines = [tuple(l) for l in r['lines']]
        blocks = [tuple(b) for b in r['blocks'] if set(b) not in [set(l) for l in lines]]
        print(f'structure {i}: n={n} D+l={r.get("D_plus_l")} target circles={r.get("circles")} blocks>=4: {len(r["blocks"])} lines: {len(lines)}')
        r['realise'] = attempt(n, blocks, lines, tries=tries)
    json.dump(recs, open(sys.argv[1].replace('.json', '_realised.json'), 'w'), default=str)
