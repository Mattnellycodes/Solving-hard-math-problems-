"""Independent exact check of a list of sympy coordinate strings: pure sympy radicals, no number
field machinery.  For every triple: exact circumcircle (centre, r^2) simplified with radsimp;
grouping by exact equality tested with sympy simplify(...)==0 on differences; distinctness of
groups certified numerically at 60 digits (a difference > 1e-40 cannot be an exact zero for
the small-height algebraic numbers involved).
Usage: python3 independent_check.py records.json [n]
"""
import sys, json, itertools
import sympy as sp
from collections import defaultdict


def circ(a, b, c):
    (ax, ay), (bx, by), (cx, cy) = a, b, c
    d = sp.simplify(2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by)))
    if d == 0: return None
    a2 = ax**2 + ay**2; b2 = bx**2 + by**2; c2 = cx**2 + cy**2
    ux = sp.radsimp(sp.simplify((a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d))
    uy = sp.radsimp(sp.simplify((a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d))
    r2 = sp.radsimp(sp.simplify((ax - ux)**2 + (ay - uy)**2))
    return (ux, uy, r2)


def check(coords):
    pts = [(sp.sympify(x), sp.sympify(y)) for x, y in coords]
    n = len(pts)
    for i, j in itertools.combinations(range(n), 2):
        assert sp.simplify(pts[i][0] - pts[j][0]) != 0 or sp.simplify(pts[i][1] - pts[j][1]) != 0, "duplicate"
    groups = []          # list of (circle, set)
    numkeys = {}
    coll = 0
    for i, j, k in itertools.combinations(range(n), 3):
        c = circ(pts[i], pts[j], pts[k])
        if c is None:
            coll += 1; continue
        key = tuple(round(float(v.evalf(60)), 30) for v in c)
        # numeric bucket, then exact confirmation
        found = None
        for g, (cc, s, nk) in enumerate(groups):
            if all(abs(float((cc[t] - c[t]).evalf(60))) < 1e-40 for t in range(3)):
                # confirm exactly
                assert all(sp.simplify(cc[t] - c[t]) == 0 for t in range(3)), "exact mismatch in bucket"
                found = g; break
        if found is None:
            groups.append((c, {i, j, k}, key))
        else:
            groups[found][1].update((i, j, k))
    # certify distinctness pairwise numerically
    for g1, g2 in itertools.combinations(range(len(groups)), 2):
        c1, c2 = groups[g1][0], groups[g2][0]
        assert max(abs(float((c1[t] - c2[t]).evalf(60))) for t in range(3)) > 1e-40
    degenerate = any(len(s) == n for _, s, _ in groups)
    return len(groups), coll, sorted((len(s) for _, s, _ in groups), reverse=True), degenerate


if __name__ == '__main__':
    d = json.load(open(sys.argv[1]))
    only = int(sys.argv[2]) if len(sys.argv) > 2 else None
    for rec in d:
        if only and rec['n'] != only: continue
        nc, coll, sizes, degen = check(rec['coordinates_exact'])
        ok = (nc == rec['circles']) and not degen
        print(f"n={rec['n']} claimed={rec['circles']} independent={nc} collinear={coll} sizes={sizes} degenerate={degen} -> {'OK' if ok else 'MISMATCH'}")
