"""Stage 2 (exact) of the audit of the 5x5-grid exhaustive claim.

For every candidate subset S (|B(S)| <= 43, from stage 1) we compute exactly
    m(S) = min over inversion centres O not in S of  #circles(inv_O(S)) = |B(S)| - max_O deg_S(O),
where deg_S(O) = number of blocks of S through O.  Every O with deg >= 2 is an intersection point
of two blocks of S, so it suffices to intersect blocks pairwise.  Intersection points are computed
in exact arithmetic in Q(sqrt(Delta)):  a point is (xp + xq*sqrt(Delta), yp + yq*sqrt(Delta)) / r
with integers, Delta a positive non-square integer (or Delta = 0 for a rational point, xq = yq = 0).
Membership of such a point on a block A(x^2+y^2)+Bx+Cy+D=0 is "rational part = 0 and surd part = 0",
which is exact because sqrt(Delta) is irrational.
Degree <= 1 centres give |B(S)| or |B(S)|-1; in pure mode the centre at infinity has degree L(S).
"""
import sys, time, json
from itertools import combinations
from math import gcd, isqrt, comb
import numpy as np
sys.path.insert(0, "/home/user/Solving-hard-math-problems-/experiments/506/audit-search-local")
from count_exact import block_key
from grid5_exhaustive import line_key, GRID, N


def blocks_of(pts, with_inf):
    blocks = {}
    for a, b, c in combinations(range(len(pts)), 3):
        k = block_key(pts[a], pts[b], pts[c])
        blocks.setdefault(k, set()).update((a, b, c))
    if with_inf:
        for a, b in combinations(range(len(pts)), 2):
            k = line_key(pts[a], pts[b])
            blocks.setdefault(k, set()).update((a, b))
    return blocks


def solve_line_circle(Bl, Cl, Dl, A, B, C, D):
    """Points of line Bl x + Cl y + Dl = 0 on circle A(x^2+y^2)+Bx+Cy+D=0. Returns list of
    (xp, xq, yp, yq, r, Delta)."""
    if Cl != 0:
        a = A * (Cl * Cl + Bl * Bl)
        b = 2 * A * Bl * Dl + B * Cl * Cl - C * Bl * Cl
        c = A * Dl * Dl - C * Dl * Cl + D * Cl * Cl
        # x = (-b +- s)/(2a), y = -(Bl x + Dl)/Cl ; common denominator r = 2 a Cl
        def pt(sign, s, Delta):
            return (-b * Cl + sign * s * Cl, sign * Cl if Delta else 0,
                    Bl * b - 2 * a * Dl - sign * s * Bl, -sign * Bl if Delta else 0,
                    2 * a * Cl, Delta)
    else:
        a = A * Bl * Bl
        b = C * Bl * Bl
        c = A * Dl * Dl - B * Dl * Bl + D * Bl * Bl
        # x = -Dl/Bl, y = (-b +- s)/(2a) ; common denominator r = 2 a Bl
        def pt(sign, s, Delta):
            return (-2 * a * Dl, 0,
                    -b * Bl + sign * s * Bl, sign * Bl if Delta else 0,
                    2 * a * Bl, Delta)
    assert a != 0
    Delta = b * b - 4 * a * c
    if Delta < 0:
        return []
    if Delta == 0:
        return [pt(1, 0, 0)]
    s = isqrt(Delta)
    if s * s == Delta:
        return [pt(1, s, 0), pt(-1, s, 0)]
    return [pt(1, 0, Delta), pt(-1, 0, Delta)]


def intersect(k1, k2):
    A1, B1, C1, D1 = k1
    A2, B2, C2, D2 = k2
    if A1 == 0 and A2 == 0:
        det = B1 * C2 - B2 * C1
        if det == 0:
            return []
        return [(C1 * D2 - C2 * D1, 0, B2 * D1 - B1 * D2, 0, det, 0)]
    if A1 == 0:
        Bl, Cl, Dl = B1, C1, D1
        circ = k2
    elif A2 == 0:
        Bl, Cl, Dl = B2, C2, D2
        circ = k1
    else:
        Bl, Cl, Dl = A2 * B1 - A1 * B2, A2 * C1 - A1 * C2, A2 * D1 - A1 * D2
        if Bl == 0 and Cl == 0:
            return []          # concentric distinct circles
        circ = k1
    return solve_line_circle(Bl, Cl, Dl, *circ)


def on_block_surd(key, P):
    A, B, C, D = key
    xp, xq, yp, yq, r, Delta = P
    rat = A * (xp * xp + yp * yp + Delta * (xq * xq + yq * yq)) + r * (B * xp + C * yp) + D * r * r
    if rat != 0:
        return False
    sur = 2 * A * (xp * xq + yp * yq) + r * (B * xq + C * yq)
    return sur == 0


def in_S(P, pts):
    xp, xq, yp, yq, r, Delta = P
    if xq != 0 or yq != 0:
        return False
    return any(xp == gx * r and yp == gy * r for gx, gy in pts)


def analyse_subset(pts, with_inf):
    blocks = blocks_of(pts, with_inf)
    keys = list(blocks)
    NB = len(keys)
    L = sum(1 for k in keys if k[0] == 0)
    maxdeg, argbest = 0, None
    for i, j in combinations(range(NB), 2):
        for P in intersect(keys[i], keys[j]):
            # sanity: P lies on both blocks
            assert on_block_surd(keys[i], P) and on_block_surd(keys[j], P), (keys[i], keys[j], P)
            if in_S(P, pts):
                continue
            deg = sum(1 for k in keys if on_block_surd(k, P))
            assert deg >= 2
            if deg > maxdeg:
                maxdeg, argbest = deg, P
    return NB, L, maxdeg, argbest, blocks


def main(mode, nb_max):
    with_inf = (mode != "pure")
    cand = np.load(f"cand_{mode}.npy")
    t0 = time.time()
    rows = []
    for idx, (mask, NB1, L1, E1) in enumerate(cand.tolist()):
        if NB1 > nb_max:
            continue
        pts = [GRID[t] for t in range(N) if (mask >> t) & 1]
        NB, L, maxdeg, P, blocks = analyse_subset(pts, with_inf)
        assert NB == NB1 and L == L1, (pts, NB, NB1, L, L1)
        best = max(maxdeg, L) if not with_inf else maxdeg
        val = NB - best
        sizes = sorted((len(v) for v in blocks.values()), reverse=True)
        rows.append({"pts": pts, "NB": NB, "L": L, "maxdeg_finite": maxdeg, "best_deg": best,
                     "min_circles": val, "block_sizes": sizes, "argbest": P})
        if idx % 100 == 0:
            print(f"  {idx}/{len(cand)} t={time.time()-t0:.0f}s", flush=True)
    from collections import Counter
    hist = Counter(r["min_circles"] for r in rows)
    print(f"[{mode}] subsets treated exactly: {len(rows)} (|B(S)| <= {nb_max}); t={time.time()-t0:.0f}s")
    print(f"[{mode}] histogram of min over all centres: {sorted(hist.items())}")
    print(f"[{mode}] histogram of (NB, best_deg): {sorted(Counter((r['NB'], r['best_deg']) for r in rows).items())}")
    json.dump(rows, open(f"stage2_{mode}.json", "w"), default=str)
    return rows


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 43)
