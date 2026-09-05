"""Robust float evaluation: normalise a configuration (centre, unit RMS radius), compute the Moebius
block structure at two tolerances; accept only if they agree.  Returns (best count, euclid count, |B|, centre)."""
import numpy as np, itertools, math
import polyhedra as ph

def normalise(P):
    P = np.array(P, float); P = P - P.mean(axis=0); s = np.sqrt((P**2).sum(axis=1).mean()); return P/s

def evaluate(P, tols=(1e-7, 1e-10)):
    P = normalise(P)
    if min(np.linalg.norm(P[i]-P[j]) for i, j in itertools.combinations(range(len(P)), 2)) < 1e-6: return None
    out = []
    for tol in tols:
        ph_tol = tol
        B = ph.blocks_planar(P, tol=tol)
        if any(len(b[1]) == len(P) for b in B): return None  # all concyclic/collinear
        nb = len(B); nl = sum(1 for b in B if b[0] == "L")
        best = (nb - nl, None, 0)
        for p, q in itertools.combinations(B, 2):
            for pt in ph.intersections(p, q):
                if np.linalg.norm(pt) > 1e4: continue
                if any(np.linalg.norm(pt - P[t]) < 1e-5 for t in range(len(P))): continue
                deg = sum(1 for b in B if ph.on_block(b, pt[0], pt[1], tol=tol*10))
                if nb - deg < best[0]: best = (nb - deg, pt, deg)
        out.append((best[0], nb - nl, nb, best[1]))
    if out[0][:3] != out[1][:3]:
        return ("UNSTABLE", out)
    return out[0]
