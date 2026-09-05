#!/usr/bin/env python3
"""Positive control for hcheck: the block structure (F, L) of real point sets (computed here from
coordinates with exact rational arithmetic) must pass every tier at every point, including the
10/11-point derived structures with the point at infinity."""
from fractions import Fraction as Fr
from itertools import combinations
from hcheck import check_structure, cumulative

def collinear(a, b, c):
    return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]) == 0

def concyclic(a, b, c, d):
    # determinant test (points in R^2)
    rows = [[p[0], p[1], p[0]*p[0]+p[1]*p[1], 1] for p in (a, b, c, d)]
    def det(m):
        if len(m) == 1: return m[0][0]
        return sum((-1)**j * m[0][j] * det([r[:j]+r[j+1:] for r in m[1:]]) for j in range(len(m)))
    return det(rows) == 0

def structure(P):
    n = len(P)
    blocks = set(); lines = set()
    for t in combinations(range(n), 3):
        a, b, c = (P[i] for i in t)
        if collinear(a, b, c):
            L = frozenset(i for i in range(n) if collinear(a, b, P[i]))
            lines.add(L)
        else:
            B = frozenset(i for i in range(n) if i in t or concyclic(a, b, c, P[i]))
            blocks.add(B)
    rich = [sorted(B) for B in blocks | lines if len(B) >= 4]
    return rich, [sorted(L) for L in lines], len([B for B in blocks])  # circles = non-line blocks

tests = {
 "n=8 record (17 circles)": [(0,0),(0,5),(0,10),(0,15),(3,6),(3,9),(5,5),(5,10)],
 "n=8 two squares (18)": [(1,0),(0,1),(-1,0),(0,-1),(2,0),(0,2),(-2,0),(0,-2)],
 "n=9 antipodal (25)": [(1,0),(-1,0),(0,1),(0,-1),(Fr(3,5),Fr(4,5)),(Fr(-3,5),Fr(-4,5)),(Fr(4,5),Fr(3,5)),(Fr(-4,5),Fr(-3,5)),(0,0)],
 "n=10 antipodal (33)": [(1,0),(-1,0),(0,1),(0,-1),(Fr(3,5),Fr(4,5)),(Fr(-3,5),Fr(-4,5)),(Fr(4,5),Fr(3,5)),(Fr(-4,5),Fr(-3,5)),(Fr(5,13),Fr(12,13)),(0,0)],
 "n=10 record-8 + 2 pts": [(0,0),(0,5),(0,10),(0,15),(3,6),(3,9),(5,5),(5,10),(7,3),(-4,2)],
}
for name, P in tests.items():
    P = [(Fr(x), Fr(y)) for x, y in P]
    rich, lines, ncirc = structure(P)
    n = len(P)
    rep = check_structure(n, rich, [lines])
    cum = cumulative(rep)
    ok = bool(cum['C'])
    print(f"{name}: n={n} rich blocks={len(rich)} lines={lines} circles={ncirc} -> "
          f"passes all tiers: {ok}" + ("" if ok else f"  VIOLATIONS: {rep['F_only']} {rep['line_sets']}"))
