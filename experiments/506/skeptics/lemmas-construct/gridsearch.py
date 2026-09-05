"""Construction attack: exhaustive search over k-subsets Q of a rational 'universe' U (grid and/or
rational circle points).  For each Q we compute (exactly, via precomputed triple->block ids)
   |B(Q)| = number of circles+lines with >=3 points of Q,
   best(Q) = |B(Q)| - max_{q in Q} deg(q)  (Euclidean count after inverting about q, n = k-1),
   euclid(Q) = |B(Q)| - #lines                (Euclidean count of Q itself, n = k).
Vectorised with numpy.  Reports minima found and the witnesses.
Usage: python3 gridsearch.py universe k   (universe: gridA, gridAxB, gridA+circle)
"""
import sys, itertools, numpy as np, json
from fractions import Fraction as Fr
from math import gcd

def universe(name):
    pts = []
    if name.startswith('grid'):
        spec = name[4:].split('+')[0]
        if 'x' in spec:
            a, b = map(int, spec.split('x'))
        else:
            a = b = int(spec)
        pts = [(Fr(x), Fr(y)) for x in range(a) for y in range(b)]
        if '+circle' in name:
            # rational points on circle centred at the grid centre with radius chosen to hit many rational points
            cx, cy = Fr(a - 1, 2), Fr(b - 1, 2)
            R2 = Fr(25, 4)   # radius 5/2: points (cx+-5/2, cy), (cx, cy+-5/2), (cx+-3/2, cy+-2), (cx+-2, cy+-3/2)
            for dx, dy in [(Fr(5,2),0),(0,Fr(5,2)),(Fr(3,2),2),(2,Fr(3,2))]:
                for sx in (1,-1):
                    for sy in (1,-1):
                        p = (cx + sx*dx, cy + sy*dy)
                        if p not in pts: pts.append(p)
    elif name == 'antipodal':
        # unit circle rational points in antipodal pairs + centre + some grid
        pts = [(Fr(0),Fr(0))]
        for t in range(1, 7):
            x = Fr(1 - t*t, 1 + t*t); y = Fr(2*t, 1 + t*t)
            pts += [(x, y), (-x, -y)]
        pts += [(Fr(1),Fr(0)),(Fr(-1),Fr(0)),(Fr(0),Fr(1)),(Fr(0),Fr(-1))]
        pts += [(Fr(2),Fr(0)),(Fr(-2),Fr(0)),(Fr(0),Fr(2)),(Fr(0),Fr(-2)),(Fr(1),Fr(1)),(Fr(-1),Fr(-1)),(Fr(1),Fr(-1)),(Fr(-1),Fr(1))]
    out = []
    for p in pts:
        if p not in out: out.append(p)
    return out

def block_key(a, b, c):
    (x1, y1), (x2, y2), (x3, y3) = a, b, c
    det = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
    if det == 0:
        A, B = y2 - y1, x1 - x2; C = -(A * x1 + B * y1)
        den = 1
        for q in (A, B, C): den = den * q.denominator // gcd(den, q.denominator)
        A, B, C = int(A*den), int(B*den), int(C*den)
        g = gcd(gcd(abs(A), abs(B)), abs(C)); A, B, C = A//g, B//g, C//g
        if A < 0 or (A == 0 and B < 0): A, B, C = -A, -B, -C
        return ('L', A, B, C)
    s1 = x1*x1 + y1*y1; s2 = x2*x2 + y2*y2; s3 = x3*x3 + y3*y3; D = 2*det
    ux = ((s2 - s1)*(y3 - y1) - (s3 - s1)*(y2 - y1)) / D
    uy = ((x2 - x1)*(s3 - s1) - (x3 - x1)*(s2 - s1)) / D
    return ('C', ux, uy, (x1 - ux)**2 + (y1 - uy)**2)

def triple_table(pts):
    m = len(pts); ids = {}; T = np.zeros((m, m, m), dtype=np.int32); isline = []
    for i, j, k in itertools.combinations(range(m), 3):
        key = block_key(pts[i], pts[j], pts[k])
        if key not in ids:
            ids[key] = len(ids); isline.append(key[0] == 'L')
        t = ids[key]
        for p in itertools.permutations((i, j, k)): T[p] = t
    return T, np.array(isline)

def search(pts, k, chunk=100000, want_best=True):
    m = len(pts); T, isline = triple_table(pts)
    trip_idx = np.array(list(itertools.combinations(range(k), 3)))
    best_min, euclid_min = 10**9, 10**9; wit_best, wit_euclid = [], []
    combos = itertools.combinations(range(m), k)
    total = 0
    while True:
        batch = list(itertools.islice(combos, chunk))
        if not batch: break
        S = np.array(batch, dtype=np.int32); total += len(S)
        ids = T[S[:, trip_idx[:, 0]], S[:, trip_idx[:, 1]], S[:, trip_idx[:, 2]]]   # (B, C(k,3))
        srt = np.sort(ids, axis=1)
        nb = 1 + np.sum(srt[:, 1:] != srt[:, :-1], axis=1)
        # euclid count of Q itself: blocks minus lines (distinct line ids)
        lid = np.where(isline[ids], ids, -1)
        ls = np.sort(lid, axis=1)
        nl = np.sum((ls[:, 1:] != ls[:, :-1]) & (ls[:, 1:] >= 0), axis=1) + (ls[:, 0] >= 0)
        eu = nb - nl
        # exclude degenerate: all on one block  (nb == 1)
        ok = nb > 1
        eu = np.where(ok, eu, 10**9)
        mn = eu.min()
        if mn < euclid_min: euclid_min, wit_euclid = int(mn), [batch[i] for i in np.where(eu == mn)[0][:5]]
        elif mn == euclid_min: wit_euclid += [batch[i] for i in np.where(eu == mn)[0][:2]]
        if want_best:
            bs = np.full(len(S), 10**9, dtype=np.int64)
            for pos in range(k):
                mask = (trip_idx == pos).any(axis=1)
                sub = np.sort(ids[:, mask], axis=1)
                deg = 1 + np.sum(sub[:, 1:] != sub[:, :-1], axis=1)
                rest = ids[:, ~mask]
                degenerate_q = (rest.max(axis=1) == rest.min(axis=1))   # Q minus q on one block
                val = np.where(ok & ~degenerate_q, nb - deg, 10**9)
                bs = np.minimum(bs, val)
            mn = bs.min()
            if mn < best_min: best_min, wit_best = int(mn), [batch[i] for i in np.where(bs == mn)[0][:5]]
            elif mn == best_min: wit_best += [batch[i] for i in np.where(bs == mn)[0][:2]]
    return total, euclid_min, wit_euclid, best_min, wit_best

if __name__ == '__main__':
    name, k = sys.argv[1], int(sys.argv[2])
    pts = universe(name)
    print(f'universe {name}: {len(pts)} points; k={k}')
    total, em, we, bm, wb = search(pts, k)
    print(f'subsets scanned: {total}')
    print(f'min Euclidean circles of a {k}-subset itself: {em}   (formula f({k}) = {(k-1)*(k-2)//2 + 1 - (k-1)//2})')
    print('   witnesses:', [[(str(pts[i][0]), str(pts[i][1])) for i in w] for w in we[:3]])
    print(f'min over ({k})-subsets Q and q in Q of |B(Q)|-deg(q), i.e. best n={k-1} count reachable by inversion: {bm}   (formula f({k-1}) = {(k-2)*(k-3)//2 + 1 - (k-2)//2})')
    print('   witnesses:', [[(str(pts[i][0]), str(pts[i][1])) for i in w] for w in wb[:3]])
