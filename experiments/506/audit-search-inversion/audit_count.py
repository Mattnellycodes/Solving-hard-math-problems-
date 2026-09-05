"""Independent exact counter for Erdos #506 (written from scratch for the audit).

circles(P) = number of distinct circles passing through >= 3 points of P.
Exact arithmetic over Fractions. A circle is identified by (cx, cy, r^2) of the circumcircle
of a non-collinear triple; a line by its normalised integer coefficients.
Also reports: distinct points, all-collinear / all-concyclic degeneracy, block sizes.
"""
from fractions import Fraction as Fr
from itertools import combinations
from collections import defaultdict
from math import gcd, comb

def norm_line(p, q):
    (x1, y1), (x2, y2) = p, q
    a, b, c = y2 - y1, x1 - x2, -(x1 * (y2 - y1) + y1 * (x1 - x2))
    # a x + b y + c = 0, with Fractions -> scale to integers
    den = 1
    for v in (a, b, c):
        den = den * v.denominator // gcd(den, v.denominator)
    a, b, c = int(a * den), int(b * den), int(c * den)
    g = gcd(gcd(abs(a), abs(b)), abs(c))
    a, b, c = a // g, b // g, c // g
    if a < 0 or (a == 0 and b < 0):
        a, b, c = -a, -b, -c
    return (a, b, c)

def circumcircle(p, q, r):
    (ax, ay), (bx, by), (cx, cy) = p, q, r
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    if d == 0:
        return None
    a2, b2, c2 = ax * ax + ay * ay, bx * bx + by * by, cx * cx + cy * cy
    ux = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d
    uy = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d
    r2 = (ax - ux) ** 2 + (ay - uy) ** 2
    return (ux, uy, r2)

def analyse(points):
    pts = [(Fr(x), Fr(y)) for x, y in points]
    n = len(pts)
    assert len(set(pts)) == n, "repeated point"
    circles = defaultdict(set)
    lines = defaultdict(set)
    for i, j, k in combinations(range(n), 3):
        cc = circumcircle(pts[i], pts[j], pts[k])
        if cc is None:
            lines[norm_line(pts[i], pts[j])].update((i, j, k))
        else:
            circles[cc].update((i, j, k))
    # consistency: every triple in exactly one block
    assert sum(comb(len(s), 3) for s in circles.values()) + sum(comb(len(s), 3) for s in lines.values()) == comb(n, 3)
    csizes = sorted((len(s) for s in circles.values()), reverse=True)
    lsizes = sorted((len(s) for s in lines.values()), reverse=True)
    degenerate = (csizes and csizes[0] == n) or (lsizes and lsizes[0] == n)
    return dict(n=n, circles=len(circles), lines=len(lines), circle_sizes=csizes, line_sizes=lsizes,
                degenerate=bool(degenerate))

def f(n):
    return comb(n - 1, 2) + 1 - (n - 1) // 2

if __name__ == "__main__":
    import json, sys
    recs = json.load(open('/home/user/Solving-hard-math-problems-/experiments/506/search-inversion/runs/antipodal_certified.json'))
    ok = True
    for k, v in sorted(recs.items(), key=lambda kv: int(kv[0])):
        pts = [tuple(p) for p in v['points']]
        r = analyse(pts)
        n = r['n']
        on_circle = [p for p in pts if p[0] ** 2 + p[1] ** 2 == 625]
        print(f"n={n:2d} circles={r['circles']:3d} f(n)={f(n):3d} lines={r['lines']} "
              f"circle_sizes(top3)={r['circle_sizes'][:3]} #3-circles={r['circle_sizes'].count(3)} line_sizes={r['line_sizes']} "
              f"degenerate={r['degenerate']} on_x2+y2=625: {len(on_circle)} centre_in={(0,0) in pts} claimed={v['circles']}")
        ok &= (r['circles'] == f(n) == v['circles']) and not r['degenerate']
    print("ALL RECORDS MATCH f(n):", ok)
