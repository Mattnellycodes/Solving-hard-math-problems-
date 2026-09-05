"""Try to realise a prescribed Moebius block structure numerically: find n points such that
each listed quadruple is concyclic or collinear (4x4 determinant = 0), via least squares from
random starts; then evaluate the count with inversion-centre optimisation.
Usage: python3 realise_design.py  (runs the built-in designs)"""
import sys, os, itertools, math, time
import numpy as np
from scipy.optimize import least_squares
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import symcore as S


def resid(x, quads, n):
    P = x.reshape(n, 2)
    Q = P[quads]
    X = Q[..., 0]; Y = Q[..., 1]
    M = np.stack([X, Y, X * X + Y * Y, np.ones_like(X)], -1)
    return np.linalg.det(M)


def realise(n, blocks, tries=300, seed=0, log=print, fix=None):
    """blocks: list of point-index lists (size >= 4); each 4-subset of a block must be coplanar-on-sphere."""
    quads = []
    for b in blocks:
        for q in itertools.combinations(sorted(b), 4): quads.append(q)
    quads = np.array(sorted(set(quads)))
    qset = set(map(tuple, quads))
    nonq = np.array([q for q in itertools.combinations(range(n), 4) if q not in qset])
    rng = np.random.default_rng(seed)
    sols = []
    for t in range(tries):
        x0 = rng.normal(size=2 * n)
        # penalty against collapse: add weak repulsion residuals
        def f(x):
            r = resid(x, quads, n)
            P = x.reshape(n, 2)
            d = P[:, None, :] - P[None, :, :]
            dist = np.sqrt((d * d).sum(-1) + 1e-12)
            iu = np.triu_indices(n, 1)
            rep = 0.02 / dist[iu]          # keeps points apart
            spread = 0.05 * (np.abs(P).max() - 2.0) if np.abs(P).max() > 2 else 0.0
            rn = resid(x, nonq, n) if len(nonq) else np.zeros(0)
            anti = 0.01 / (np.abs(rn) + 0.02)    # keeps non-design quadruples non-concyclic
            return np.concatenate([r, rep - rep.mean(), [spread], anti])
        res = least_squares(f, x0, method='lm', max_nfev=4000)
        r = resid(res.x, quads, n)
        P = res.x.reshape(n, 2)
        # scale-normalise residual
        sc = np.abs(P).max()
        if np.abs(r).max() / sc ** 4 < 1e-8:
            # polish without repulsion
            res2 = least_squares(lambda x: resid(x, quads, n), res.x, method='trf', max_nfev=2000, xtol=1e-15, ftol=1e-15, gtol=1e-15)
            P = res2.x.reshape(n, 2)
            d = P[:, None, :] - P[None, :, :]
            dmin = np.sqrt((d * d).sum(-1) + np.eye(n) * 1e9).min()
            if dmin < 1e-3 * np.abs(P).max(): continue
            ev = S.evaluate_fast(P)
            if ev is None: continue
            sols.append((ev['best'], ev['nblocks'], ev['deg'], P))
            log(f"  try {t}: realised; best count {ev['best']} nblocks {ev['nblocks']} deg {ev['deg']} sizes {ev['sizes']}")
            if len(sols) >= 5: break
    return sols


if __name__ == '__main__':
    # n=7: 2-(7,4,2) design = complements of Fano lines
    fano7 = [[0,1,2,3],[0,1,4,5],[2,3,4,5],[0,2,4,6],[1,3,4,6],[1,2,5,6],[0,3,5,6]]
    print("n=7, 7 four-blocks (complements of Fano lines):")
    s = realise(7, fano7, tries=400)
    print("  solutions found:", len(s))
    # n=7: 6 of the 7 blocks (should be realisable: triangle+midpoints+centroid gives 11)
    print("n=7, 6 four-blocks:")
    s = realise(7, fano7[:6], tries=60)
    print("  solutions found:", len(s), "best", min([x[0] for x in s], default=None))
    # n=8: 12 blocks (AG(3,2) minus a parallel class)
    ag = [[0,1,2,3],[0,1,4,5],[2,3,4,5],[0,2,4,6],[1,3,4,6],[1,2,5,6],[0,3,4,7],[0,2,5,7],[1,3,5,7],[0,1,6,7],[2,3,6,7],[4,5,6,7]]
    print("n=8, 12 four-blocks:")
    s = realise(8, ag, tries=150)
    print("  solutions found:", len(s), "best", min([x[0] for x in s], default=None))
    # n=8: 13 blocks (12 + one tetrahedron) -- expected unrealisable
    print("n=8, 13 four-blocks:")
    s = realise(8, ag + [[0,3,5,6]], tries=150)
    print("  solutions found:", len(s))
