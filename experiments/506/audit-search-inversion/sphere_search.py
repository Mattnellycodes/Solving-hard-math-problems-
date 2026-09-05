"""Independent exhaustive search on an integer lattice sphere x^2+y^2+z^2=N (written for the audit).

For every n-subset S of the integer points U and every projection centre Q on the sphere
(Q not in S), circles(proj_Q S) = #{planes with >= 3 points of S} - #{such planes containing Q}.
The optimum over Q is attained either at a point on <= 1 block (lines <= 1) or at an
intersection point of two block-circles; the latter are enumerated exactly (rational points
and conjugate pairs in Q(sqrt(D))).  Output: minimum circle count over all (S, Q) and all
subsets attaining <= target.
"""
import sys, time
from fractions import Fraction as Fr
from itertools import combinations
from math import comb, isqrt
import numpy as np
from sphere_audit import lattice_points, rich_planes

def is_square_fr(x):
    if x < 0: return False
    p, q = x.numerator, x.denominator
    return isqrt(p) ** 2 == p and isqrt(q) ** 2 == q

def line_sphere(k1, k2, N):
    """Intersection of planes k1,k2 with |x|^2=N. Returns list of ('rat', Q) or ('irr', p0, u, a,b,c)."""
    a1, b1, c1, d1 = k1; a2, b2, c2, d2 = k2
    u = (b1 * c2 - c1 * b2, c1 * a2 - a1 * c2, a1 * b2 - b1 * a2)
    if u == (0, 0, 0):
        return []  # parallel planes
    # point p0 on the line: solve two equations with one coordinate set to 0
    p0 = None
    for drop in range(3):
        idx = [i for i in range(3) if i != drop]
        m11, m12 = k1[idx[0]], k1[idx[1]]; m21, m22 = k2[idx[0]], k2[idx[1]]
        det = m11 * m22 - m12 * m21
        if det != 0:
            s1 = Fr(d1 * m22 - m12 * d2, det); s2 = Fr(m11 * d2 - d1 * m21, det)
            p0 = [Fr(0)] * 3; p0[idx[0]] = s1; p0[idx[1]] = s2
            break
    a = Fr(sum(x * x for x in u)); b = 2 * sum(p * x for p, x in zip(p0, u)); c = sum(p * p for p in p0) - N
    disc = b * b - 4 * a * c
    if disc < 0:
        return []
    if is_square_fr(disc):
        s = Fr(isqrt(disc.numerator), isqrt(disc.denominator))
        out = []
        for sg in (1, -1):
            t = (-b + sg * s) / (2 * a)
            out.append(('rat', tuple(p + t * x for p, x in zip(p0, u))))
        return out
    return [('irr', tuple(p0), u)]

def plane_contains(pl, cand):
    a, b, c, d = pl
    if cand[0] == 'rat':
        Q = cand[1]
        return a * Q[0] + b * Q[1] + c * Q[2] == d
    p0, u = cand[1], cand[2]
    # irrational point on line p0 + t u: contained iff the whole line is in the plane
    return (a * p0[0] + b * p0[1] + c * p0[2] == d) and (a * u[0] + b * u[1] + c * u[2] == 0)

def build_universe(N):
    U = lattice_points(N)
    idx = {p: i for i, p in enumerate(U)}
    planes = rich_planes(U)
    plist = list(planes.keys())
    pmask = np.array([sum(1 << i for i in planes[pl]) for pl in plist], dtype=np.uint32)
    return U, idx, plist, pmask

def build_centres(U, idx, plist, N):
    """All candidate centres lying on >= 2 universe planes, as (plane-bitmask, universe point index or -1).
    Plane membership is tested exactly with integer numpy arithmetic (coordinates cleared of denominators)."""
    from math import lcm
    P = len(plist)
    A = np.array([pl[:3] for pl in plist], dtype=np.int64)
    D = np.array([pl[3] for pl in plist], dtype=np.int64)
    def mask_of(sel):
        return int.from_bytes(np.packbits(sel, bitorder='little').tobytes(), 'little')
    cands = {}
    for i, j in combinations(range(P), 2):
        for cand in line_sphere(plist[i], plist[j], N):
            if cand[0] == 'rat':
                key = cand
                if key in cands: continue
                Q = cand[1]
                den = lcm(*(q.denominator for q in Q))
                qn = np.array([int(q * den) for q in Q], dtype=np.int64)
                sel = (A @ qn) == D * den
                up = idx.get(tuple(int(q) for q in Q), -1) if den == 1 else -1
            else:
                key = ('irr', cand[1], cand[2])
                if key in cands: continue
                p0, u = cand[1], cand[2]
                den = lcm(*(q.denominator for q in p0))
                pn = np.array([int(q * den) for q in p0], dtype=np.int64)
                un = np.array(u, dtype=np.int64)
                sel = ((A @ pn) == D * den) & ((A @ un) == 0)
                up = -1
            cands[key] = (mask_of(sel), up)
    uniq = {}
    for m, up in cands.values():
        uniq[(m, up)] = 1
    return list(uniq.keys())

_CACHE = {}
def search(N, n, target, chunk=20000, verbose=True):
    t0 = time.time()
    if N not in _CACHE:
        U, idx, plist, pmask = build_universe(N)
        cands = build_centres(U, idx, plist, N)
        _CACHE[N] = (U, idx, plist, pmask, cands)
    U, idx, plist, pmask, cands = _CACHE[N]
    P = len(plist)
    if verbose:
        print(f"N={N}: |U|={len(U)}, planes={P}, centre candidates={len(cands)} ({time.time()-t0:.1f}s)", flush=True)
    # index candidates by plane pairs for fast lookup
    from collections import defaultdict
    bypair = defaultdict(list)
    for ci, (m, up) in enumerate(cands):
        bits = np.nonzero(np.unpackbits(np.frombuffer(m.to_bytes((P + 7) // 8, 'little'), dtype=np.uint8), bitorder='little')[:P])[0].tolist()
        for a, b in combinations(bits, 2):
            bypair[(a, b)].append(ci)
    # enumerate subsets
    m = len(U)
    masks = np.fromiter((sum(1 << i for i in c) for c in combinations(range(m), n)), dtype=np.uint32, count=comb(m, n))
    if verbose:
        print(f"  subsets: {len(masks)} ({time.time()-t0:.1f}s)", flush=True)
    maxlines = comb(n, 2) // 3  # safe: lines pairwise share <= 1 point, each has >= 3 points
    surv = []
    Rmin = 10 ** 9
    for s in range(0, len(masks), chunk):
        ch = masks[s:s + chunk]
        cnt = np.bitwise_count(ch[:, None] & pmask[None, :])
        R = (cnt >= 3).sum(axis=1)
        Rmin = min(Rmin, int(R.min()))
        sel = np.nonzero(R <= target + maxlines)[0]
        for k in sel:
            surv.append((int(ch[k]), int(R[k])))
    if verbose:
        print(f"  min #blocks over subsets: {Rmin}; survivors with R <= {target+maxlines}: {len(surv)} ({time.time()-t0:.1f}s)", flush=True)
    best = 10 ** 9; hits = []
    for S, R in surv:
        blocks = [k for k in range(P) if bin(int(pmask[k]) & S).count('1') >= 3]
        assert len(blocks) == R
        bm = sum(1 << k for k in blocks)
        L = 1 if R >= 1 else 0
        bestQ = None
        seen = set()
        for a, b in combinations(blocks, 2):
            for ci in bypair.get((a, b), ()):
                if ci in seen: continue
                seen.add(ci)
                cm, up = cands[ci]
                if up >= 0 and (S >> up) & 1:
                    continue  # centre would be a point of S
                l = bin(cm & bm).count('1')
                if l > L:
                    L = l; bestQ = ci
        degenerate = any(bin(int(pmask[k]) & S).count('1') == n for k in blocks)
        c = R - L
        if not degenerate:
            best = min(best, c)
            if c <= target:
                hits.append((S, R, L, c, bestQ))
    if verbose:
        print(f"  RESULT N={N} n={n}: min circles = {best}; #subsets with <= {target}: {len(hits)} ({time.time()-t0:.1f}s)", flush=True)
    return best, hits, U, plist, cands

if __name__ == "__main__":
    N = int(sys.argv[1])
    for spec in sys.argv[2:]:
        n, target = map(int, spec.split(':'))
        best, hits, U, plist, cands = search(N, n, target)
        from collections import Counter
        print("hit structures (R, L, circles):", Counter((h[1], h[2], h[3]) for h in hits), flush=True)
        for S, R, L, c, q in hits[:3]:
            pts = [U[i] for i in range(len(U)) if S >> i & 1]
            print("  example:", pts, "R=", R, "L=", L, "circles=", c, flush=True)
