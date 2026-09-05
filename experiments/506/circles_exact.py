"""Exact counter for the number of circles determined by a finite planar point set.

A circle is "determined" by P if it passes through at least three points of P.
Collinear triples determine no circle. All arithmetic is exact (fractions.Fraction),
so the count is a certificate for the given rational coordinates.

Usage: python3 scripts/circles.py            (runs built-in examples)
       from scripts.circles import count_circles
"""
from fractions import Fraction as Fr
from itertools import combinations
from collections import defaultdict


def circumcircle(a, b, c):
    """Return (cx, cy, r2) of the circle through a,b,c, or None if collinear. Exact."""
    ax, ay = a; bx, by = b; cx, cy = c
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    if d == 0:
        return None
    a2 = ax * ax + ay * ay; b2 = bx * bx + by * by; c2 = cx * cx + cy * cy
    ux = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d
    uy = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d
    r2 = (ax - ux) ** 2 + (ay - uy) ** 2
    return (ux, uy, r2)


def analyze(points):
    """Return dict with circles (map circle->set of point indices), collinear triple count."""
    pts = [(Fr(x), Fr(y)) for x, y in points]
    assert len(set(pts)) == len(pts), "duplicate points"
    circles = defaultdict(set)
    collinear = 0
    for i, j, k in combinations(range(len(pts)), 3):
        cc = circumcircle(pts[i], pts[j], pts[k])
        if cc is None:
            collinear += 1
        else:
            circles[cc].update((i, j, k))
    return {"circles": dict(circles), "num_circles": len(circles), "collinear_triples": collinear}


def count_circles(points):
    return analyze(points)["num_circles"]


def formula(n):
    """Elliott / Purdy-Smith value, proven minimal for n > 393."""
    return (n - 1) * (n - 2) // 2 + 1 - (n - 1) // 2


def all_on_one_circle_or_line(points):
    pts = [(Fr(x), Fr(y)) for x, y in points]
    n = len(pts)
    if n < 3:
        return True
    # find first non-collinear triple
    for i, j, k in combinations(range(n), 3):
        cc = circumcircle(pts[i], pts[j], pts[k])
        if cc is not None:
            ux, uy, r2 = cc
            return all((x - ux) ** 2 + (y - uy) ** 2 == r2 for x, y in pts)
    return True  # all collinear


if __name__ == "__main__":
    # 1. Formula configuration: n-1 points on a circle in antipodal pairs + centre.
    def formula_config(n):
        # rational points on unit circle via Pythagorean parametrisation, in antipodal pairs
        pts = [(Fr(0), Fr(0))]
        m = n - 1
        t = 1
        while len(pts) < n:
            x = Fr(1 - t * t, 1 + t * t); y = Fr(2 * t, 1 + t * t)
            for p in [(x, y), (-x, -y)]:
                if len(pts) < n and p not in pts:
                    pts.append(p)
            t += 1
        return pts
    for n in range(4, 13):
        c = count_circles(formula_config(n))
        print(f"n={n:2d} formula={formula(n):3d} antipodal-config={c:3d}")
    # 2. Two concentric axis-aligned squares (n=8).
    sq = [(1, 0), (0, 1), (-1, 0), (0, -1), (2, 0), (0, 2), (-2, 0), (0, -2)]
    res = analyze(sq)
    print("two concentric squares n=8:", res["num_circles"], "circles,", res["collinear_triples"], "collinear triples; formula", formula(8))
    print("degenerate?", all_on_one_circle_or_line(sq))
    sizes = sorted(len(s) for s in res["circles"].values())
    print("points per circle:", sizes)
