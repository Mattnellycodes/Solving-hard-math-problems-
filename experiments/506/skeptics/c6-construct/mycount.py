"""Independent exact circle counter (skeptic c6-construct).

A circle/line is identified by the primitive integer vector (a, d, e, f) of its equation
a (x^2 + c y^2) + d x + e y + f = 0 (a = 0 for a line), obtained by exact Fraction arithmetic.
The quadratic form x^2 + c y^2 lets the same code handle the square lattice (c = 1) and the
triangular lattice (c = 3, with (X, Y) = (2a + b, b) for Eisenstein point a + b*omega).
"""
from fractions import Fraction as Fr
from itertools import combinations
from math import gcd, comb


def _primitive(vec):
    den = 1
    for q in vec:
        den = den * q.denominator // gcd(den, q.denominator)
    ints = [int(q * den) for q in vec]
    g = 0
    for v in ints:
        g = gcd(g, abs(v))
    ints = [v // g for v in ints]
    for v in ints:
        if v != 0:
            if v < 0:
                ints = [-w for w in ints]
            break
    return tuple(ints)


def block_key(p, q, r, c=1):
    """Key of the unique circle or line through three distinct points (Fractions)."""
    (x1, y1), (x2, y2), (x3, y3) = p, q, r
    m = (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
    if m == 0:  # collinear: line d x + e y + f = 0
        d, e = y2 - y1, x1 - x2
        f = -(d * x1 + e * y1)
        return ("L",) + _primitive((d, e, f))
    # circle: N + d x + e y + f = 0 with N = x^2 + c y^2; solve 3x3 by Cramer
    n1, n2, n3 = x1 * x1 + c * y1 * y1, x2 * x2 + c * y2 * y2, x3 * x3 + c * y3 * y3
    # rows: [x_i, y_i, 1] * (d, e, f) = -n_i
    det = m
    dd = (-n1 * (y2 - y3) - n2 * (y3 - y1) - n3 * (y1 - y2)) / det
    ee = (x1 * (-n2 + n3) + x2 * (-n3 + n1) + x3 * (-n1 + n2)) / det
    ff = (x1 * (y2 * (-n3) - y3 * (-n2)) - y1 * (x2 * (-n3) - x3 * (-n2)) + (x2 * y3 - x3 * y2) * (-n1)) / det
    return ("C",) + _primitive((Fr(1), dd, ee, ff))


def analyse(points, c=1):
    pts = [(Fr(x), Fr(y)) for x, y in points]
    n = len(pts)
    assert len(set(pts)) == n, "duplicate points"
    blocks = {}
    for i, j, k in combinations(range(n), 3):
        blocks.setdefault(block_key(pts[i], pts[j], pts[k], c), set()).update((i, j, k))
    circles = {k: v for k, v in blocks.items() if k[0] == "C"}
    lines = {k: v for k, v in blocks.items() if k[0] == "L"}
    # sanity: every block key must be consistent (each point on it satisfies the equation)
    for key, members in blocks.items():
        for i in members:
            x, y = pts[i]
            if key[0] == "C":
                _, a, d, e, f = key
                assert a * (x * x + c * y * y) + d * x + e * y + f == 0
            else:
                _, d, e, f = key
                assert d * x + e * y + f == 0
    degenerate = any(len(v) == n for v in blocks.values())
    return {
        "n": n,
        "circles": len(circles),
        "lines": len(lines),
        "circle_sizes": sorted(len(v) for v in circles.values()),
        "line_sizes": sorted(len(v) for v in lines.values()),
        "degenerate": degenerate,
        "blocks": blocks,
    }


def f(n):
    return comb(n - 1, 2) + 1 - (n - 1) // 2


if __name__ == "__main__":
    P6 = [(0, 0), (20, 0), (5, 15), (5, 0), (2, 6), (10, 10)]
    r = analyse(P6)
    print("n=6 claimed construction:", {k: v for k, v in r.items() if k != "blocks"})
    for key, members in r["blocks"].items():
        print("   ", key, sorted(members))
    # cross-check formula (1): circles = C(n,3) - D - l
    D = sum(comb(len(v), 3) - 1 for v in r["blocks"].values() if len(v) >= 4)
    print("   C(6,3) - D - l =", comb(6, 3) - D - r["lines"])
