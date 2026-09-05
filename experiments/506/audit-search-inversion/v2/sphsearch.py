"""Audit v2: exhaustive subset + projection-centre search on the lattice sphere x^2+y^2+z^2 = N.

For S (n lattice points on the sphere) and a centre Q on the sphere, Q not in S, the planar
image has  circles = R(S) - L(S,Q)  where R = #planes containing >= 3 points of S and
L = #such planes containing Q.  L >= 2 forces Q onto the intersection of two block planes
with the sphere, so all such Q are enumerated exactly: rational points, or an irrational
conjugate pair (a rational plane contains one iff it contains the whole line).
Every subset (mask) is scanned with numpy; the exact best L is computed for all subsets that
could reach the target.  Rational-centre hits are re-verified by exact inversion + coplanar
count (s2lat_check.count_coplanar).
usage: python3 sphsearch.py N n1 n2 ...   (target = f(n) reported together with min)
"""
import sys, time
from fractions import Fraction as Fr
from itertools import combinations
from math import comb, gcd, isqrt, lcm
from collections import Counter, defaultdict
import numpy as np
from s2lat_check import invert, count_coplanar, cross, sub, dot, prim

def f(n): return comb(n - 1, 2) + 1 - (n - 1) // 2

def lattice_points(N):
    r = isqrt(N)
    return sorted((x, y, z) for x in range(-r, r + 1) for y in range(-r, r + 1) for z in range(-r, r + 1)
                  if x * x + y * y + z * z == N)

def plane_key(X, Y, Z):
    nrm = cross(sub(Y, X), sub(Z, X))
    if not any(nrm): return None
    a, b, c = prim(nrm)
    return (a, b, c, a * X[0] + b * X[1] + c * X[2])

def rich_planes(U):
    pl = defaultdict(set)
    for i, j, k in combinations(range(len(U)), 3):
        key = plane_key(U[i], U[j], U[k])
        if key is not None: pl[key].update((i, j, k))
    return pl

def line_of(k1, k2):
    """intersection line of two planes: (point p0 (Fractions), direction u (ints)) or None if parallel"""
    u = cross(k1[:3], k2[:3])
    if not any(u): return None
    for drop in range(3):
        idx = [i for i in range(3) if i != drop]
        m11, m12, m21, m22 = k1[idx[0]], k1[idx[1]], k2[idx[0]], k2[idx[1]]
        det = m11 * m22 - m12 * m21
        if det:
            p0 = [Fr(0)] * 3
            p0[idx[0]] = Fr(k1[3] * m22 - m12 * k2[3], det)
            p0[idx[1]] = Fr(m11 * k2[3] - k1[3] * m21, det)
            return tuple(p0), u
    return None

def sqrt_fr(x):
    p, q = x.numerator, x.denominator
    sp, sq = isqrt(p), isqrt(q)
    return Fr(sp, sq) if sp * sp == p and sq * sq == q else None

def centre_candidates(plist, N):
    """returns dict key -> ('rat', Q) or ('irr', p0, u); keys are canonical"""
    out = {}
    for k1, k2 in combinations(plist, 2):
        ln = line_of(k1, k2)
        if ln is None: continue
        p0, u = ln
        A = Fr(dot(u, u)); B = 2 * dot(p0, u); C = dot(p0, p0) - N
        disc = B * B - 4 * A * C
        if disc < 0: continue
        s = sqrt_fr(disc)
        if s is None:
            # canonical line key: primitive direction, foot from origin
            d = prim(u)
            foot = tuple(x - Fr(dot(p0, d), dot(d, d)) * di for x, di in zip(p0, d))
            out[('irr', d, foot)] = ('irr', foot, d)
        else:
            for sg in (1, -1):
                t = (-B + sg * s) / (2 * A)
                Q = tuple(p + t * x for p, x in zip(p0, u))
                out[('rat', Q)] = ('rat', Q)
    return list(out.values())

def contains(pl, cand):
    a, b, c, d = pl
    if cand[0] == 'rat':
        Q = cand[1]; return a * Q[0] + b * Q[1] + c * Q[2] == d
    p0, u = cand[1], cand[2]
    return a * p0[0] + b * p0[1] + c * p0[2] == d and a * u[0] + b * u[1] + c * u[2] == 0

def pack(bits, P):
    W = (P + 63) // 64
    w = [0] * W
    for b in bits: w[b // 64] |= 1 << (b % 64)
    return w  # unsigned 64-bit words

_UNIV = {}
def universe(N):
    if N in _UNIV: return _UNIV[N]
    U = lattice_points(N); m = len(U); uidx = {p: i for i, p in enumerate(U)}
    planes = rich_planes(U); plist = list(planes); P = len(plist)
    pmask = np.array([sum(1 << i for i in planes[k]) for k in plist], dtype=np.int64)
    cands = centre_candidates(plist, N)
    W = (P + 63) // 64
    A = np.array([k[:3] for k in plist], dtype=np.int64); Dv = np.array([k[3] for k in plist], dtype=np.int64)
    cmat = np.zeros((len(cands), W), dtype=np.uint64)
    cpoint = np.full(len(cands), -1, dtype=np.int64)   # universe index if the candidate is a lattice point
    for ci, cand in enumerate(cands):
        if cand[0] == 'rat':
            Q = cand[1]; den = lcm(*(q.denominator for q in Q))
            qn = np.array([int(q * den) for q in Q], dtype=np.int64)
            sel = (A @ qn) == Dv * den
            if den == 1: cpoint[ci] = uidx.get(tuple(int(q) for q in Q), -1)
        else:
            p0, u = cand[1], cand[2]; den = lcm(*(q.denominator for q in p0))
            pn = np.array([int(q * den) for q in p0], dtype=np.int64)
            sel = ((A @ pn) == Dv * den) & ((A @ np.array(u, dtype=np.int64)) == 0)
        cmat[ci] = pack(np.nonzero(sel)[0].tolist(), P)
    # spot-check the vectorised membership against the exact predicate on 200 random candidates
    import random
    for ci in random.sample(range(len(cands)), min(200, len(cands))):
        bits = [k for k in range(P) if contains(plist[k], cands[ci])]
        assert list(cmat[ci]) == pack(bits, P)
    _UNIV[N] = (U, m, uidx, plist, P, pmask, cands, cmat, cpoint, W)
    return _UNIV[N]

def search(N, n, target=None):
    t0 = time.time()
    if target is None: target = f(n)
    U, m, uidx, plist, P, pmask, cands, cmat, cpoint, W = universe(N)
    Lcap = int(np.bitwise_count(cmat).sum(axis=1).max())
    print(f"N={N} n={n}: |U|={m} rich planes={P} centre candidates={len(cands)} (rational {sum(c[0]=='rat' for c in cands)}) "
          f"max planes through a candidate={Lcap} [{time.time()-t0:.1f}s]", flush=True)
    # scan all subsets
    Rmin = 10 ** 9; surv = []
    chunk = 1 << 15
    combs = combinations(range(m), n)
    total = comb(m, n)
    done = 0
    while done < total:
        take = min(chunk, total - done)
        masks = np.fromiter((sum(1 << i for i in c) for _, c in zip(range(take), combs)), dtype=np.int64, count=take)
        done += take
        cnt = np.bitwise_count(masks[:, None] & pmask[None, :])
        R = (cnt >= 3).sum(axis=1)
        Rmin = min(Rmin, int(R.min()))
        surv.extend((int(masks[k]), int(R[k])) for k in np.nonzero(R <= max(target, Rmin - 1) + min(Lcap, comb(n, 2) // 3))[0])
    # two blocks through Q meet the sphere only in Q and one further point, so they share at most one
    # point of S: the blocks through Q are 'lines' pairwise sharing <= 1 point, each with >= 3 points of S,
    # hence L(S,Q) <= floor(C(n,2)/3).  (Lcap, the max number of rich planes through a candidate, is only informative.)
    Lb_max = min(Lcap, comb(n, 2) // 3)
    thr = max(target, Rmin - 1) + Lb_max
    surv = [(S, R) for S, R in surv if R <= thr]
    print(f"  scanned {done} subsets; min R={Rmin}; survivors (R <= max(target,Rmin-1)+min(Lcap,C(n,2)//3)={thr}): {len(surv)} [{time.time()-t0:.1f}s]", flush=True)
    best = 10 ** 9; hits = []; bestex = None
    for S, R in surv:
        bl = [k for k in range(P) if bin(int(pmask[k]) & S).count('1') >= 3]
        assert len(bl) == R
        if any(bin(int(pmask[k]) & S).count('1') == n for k in bl):
            continue  # degenerate: all of S on one circle
        bvec = np.array(pack(bl, P), dtype=np.uint64)
        L = np.bitwise_count(cmat & bvec[None, :]).sum(axis=1)
        # exclude candidates that are points of S
        bad = (cpoint >= 0) & (((np.int64(S) >> np.maximum(cpoint, 0)) & 1) == 1)
        L = np.where(bad, 0, L)
        assert int(L.max()) <= min(R, comb(n, 2) // 3)
        ci = int(L.argmax()); Lb = max(1, int(L[ci]))
        c = R - Lb
        if c < best: best, bestex = c, (S, R, Lb, ci)
        if c <= target: hits.append((S, R, Lb, ci))
    print(f"  RESULT N={N} n={n}: min circles over universe = {best} (f(n)={f(n)}, gap {best - f(n):+d}); hits <= {target}: {len(hits)} [{time.time()-t0:.1f}s]", flush=True)
    print("  hit structures Counter[(R, L, circles)]:", Counter((R, L, R - L) for _, R, L, _ in hits), flush=True)
    # verify: recount up to 5 hits (or the best example) by exact inversion if the centre is rational
    todo = hits[:5] if hits else ([bestex] if bestex else [])
    for S, R, L, ci in todo:
        pts = [U[i] for i in range(m) if S >> i & 1]
        cand = cands[ci]
        if cand[0] == 'rat':
            img = [invert(p, cand[1]) for p in pts]
            nc, nl, csz, lsz = count_coplanar(img)
            print(f"   verify pts={pts} Q={tuple(str(q) for q in cand[1])}: planar count circles={nc} lines={nl} (plane count R-L={R-L}) {'OK' if nc == R-L else 'MISMATCH'} max circle {csz[0]}, lines {lsz}", flush=True)
        else:
            print(f"   pts={pts} irrational centre on line foot={cand[1]} dir={cand[2]}: R-L={R-L} (not re-projected)", flush=True)
    return best

if __name__ == "__main__":
    N = int(sys.argv[1])
    for a in sys.argv[2:]:
        search(N, int(a))
