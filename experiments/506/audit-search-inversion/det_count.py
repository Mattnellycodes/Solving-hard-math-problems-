"""Second independent recount: blocks via the 4x4 concyclicity determinant (points lifted to the paraboloid).
Four points are concyclic or collinear iff det[[x,y,x^2+y^2,1]] = 0; a triple is collinear iff the 3x3
orientation determinant vanishes. Integer arithmetic only."""
import json
from itertools import combinations
from math import comb
def det3(m):
    return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1]) - m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
            + m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))
def det4(m):
    s = 0
    for c in range(4):
        minor = [[m[r][cc] for cc in range(4) if cc != c] for r in range(1, 4)]
        s += (-1) ** c * m[0][c] * det3(minor)
    return s
def lift(p): return [p[0], p[1], p[0]*p[0]+p[1]*p[1], 1]
def count(pts):
    n = len(pts); circles = set(); lines = set()
    for t in combinations(range(n), 3):
        col = det3([[pts[i][0], pts[i][1], 1] for i in t]) == 0
        blk = set(t)
        for m in range(n):
            if m in blk: continue
            if col:
                if det3([[pts[i][0], pts[i][1], 1] for i in (t[0], t[1], m)]) == 0: blk.add(m)
            else:
                if det4([lift(pts[i]) for i in (*t, m)]) == 0: blk.add(m)
        (lines if col else circles).add(frozenset(blk))
    assert sum(comb(len(b), 3) for b in circles | lines) == comb(n, 3)
    return len(circles), len(lines), max(map(len, circles | lines))
recs = json.load(open('/home/user/Solving-hard-math-problems-/experiments/506/search-inversion/runs/antipodal_certified.json'))
for k, v in sorted(recs.items(), key=lambda kv: int(kv[0])):
    pts = [tuple(p) for p in v['points']]; n = len(pts)
    c, l, mx = count(pts)
    print(f"n={n:2d} circles={c:3d} lines={l} largest block={mx} f(n)={comb(n-1,2)+1-(n-1)//2} match={c == v['circles'] == comb(n-1,2)+1-(n-1)//2 and mx < n}")
