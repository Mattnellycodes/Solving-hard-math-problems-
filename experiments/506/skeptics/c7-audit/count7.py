"""Independent exact block/circle counter (skeptic c7-audit).

Different from circles_exact.py / certify.py: a circle is identified by the exact rational
coefficients (a, b, c) of x^2 + y^2 + a x + b y + c = 0 (Fraction-normalised, hence canonical);
a line by its normalised equation.  We then recompute the Moebius identity
    circles = C(n,3) - D - ell,   D = sum_{blocks of size>=4} (C(k,3)-1),  ell = #lines(>=3 pts)
from the raw block list, as an internal consistency check of the reformulation.
"""
from fractions import Fraction as Fr
from itertools import combinations
from math import comb
import sys


def circle_coeffs(p, q, r):
    (x1, y1), (x2, y2), (x3, y3) = p, q, r
    # x^2+y^2+ax+by+c=0 through the three points; subtract the first equation from the others
    A1, B1, R1 = x2 - x1, y2 - y1, -(x2 * x2 + y2 * y2 - x1 * x1 - y1 * y1)
    A2, B2, R2 = x3 - x1, y3 - y1, -(x3 * x3 + y3 * y3 - x1 * x1 - y1 * y1)
    det = A1 * B2 - A2 * B1
    if det == 0:
        return None
    a = Fr(R1 * B2 - R2 * B1) / det
    b = Fr(A1 * R2 - A2 * R1) / det
    c = -(x1 * x1 + y1 * y1) - a * x1 - b * y1
    return ("C", a, b, c)


def line_coeffs(p, q):
    (x1, y1), (x2, y2) = p, q
    A, B = y2 - y1, x1 - x2
    C = A * x1 + B * y1
    if A != 0:
        return ("L", Fr(1), Fr(B) / A, Fr(C) / A)
    return ("L", Fr(0), Fr(1), Fr(C) / B)


def blocks_of(points):
    pts = [(Fr(x), Fr(y)) for x, y in points]
    assert len(set(pts)) == len(pts), "repeated point"
    blocks = {}
    for i, j, k in combinations(range(len(pts)), 3):
        key = circle_coeffs(pts[i], pts[j], pts[k])
        if key is None:
            key = line_coeffs(pts[i], pts[j])
        blocks.setdefault(key, set()).update((i, j, k))
    return blocks


def report(points, name):
    n = len(points)
    blocks = blocks_of(points)
    circles = {k: v for k, v in blocks.items() if k[0] == "C"}
    lines = {k: v for k, v in blocks.items() if k[0] == "L"}
    degenerate = any(len(v) == n for v in blocks.values())
    D = sum(comb(len(v), 3) - 1 for v in blocks.values() if len(v) >= 4)
    ell = len(lines)
    # sanity: every triple in exactly one block
    assert sum(comb(len(v), 3) for v in blocks.values()) == comb(n, 3)
    # sanity: blocks pairwise share <= 2 points, lines pairwise <= 1
    bl = list(blocks.values())
    for a, b in combinations(bl, 2):
        assert len(a & b) <= 2
    for a, b in combinations(list(lines.values()), 2):
        assert len(a & b) <= 1
    print(f"== {name}: n={n}")
    print(f"   circles(>=3 pts) = {len(circles)}   lines(>=3 pts) = {ell}   degenerate = {degenerate}")
    print(f"   circle sizes = {sorted(len(v) for v in circles.values())}")
    print(f"   line sizes   = {sorted(len(v) for v in lines.values())}")
    print(f"   D = {D}, ell = {ell}, C(n,3)-D-ell = {comb(n,3) - D - ell}  (must equal #circles)")
    assert comb(n, 3) - D - ell == len(circles)
    print("   4-blocks (circles):", sorted(sorted(v) for v in circles.values() if len(v) >= 4))
    print("   lines:", sorted(sorted(v) for v in lines.values()))
    return len(circles), blocks


if __name__ == "__main__":
    P6 = [(0, 0), (20, 0), (5, 15), (5, 0), (2, 6), (10, 10)]
    P7 = P6 + [(5, 5)]
    report(P6, "theory n=6 construction")
    c7, b7 = report(P7, "theory n=7 construction")
    # Check the geometric description: (5,5) should be the orthocentre of triangle (0,0),(20,0),(5,15)
    A, B, C, H = (Fr(0), Fr(0)), (Fr(20), Fr(0)), (Fr(5), Fr(15)), (Fr(5), Fr(5))
    def dot(u, v): return u[0] * v[0] + u[1] * v[1]
    def sub(u, v): return (u[0] - v[0], u[1] - v[1])
    print("   H orthocentre check: (H-A).(B-C) =", dot(sub(H, A), sub(B, C)), " (H-B).(A-C) =", dot(sub(H, B), sub(A, C)))
    print("RESULT n=7 circles =", c7)
