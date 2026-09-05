"""Independent exact block/circle counter (hostile-referee re-implementation).

Representation deliberately different from circles_exact.py / certify.py:
  * a circle is identified by the exact rational coefficients (D, E, F) of x^2 + y^2 + D x + E y + F = 0
    obtained by Gaussian elimination over Fractions;
  * a line is identified by its normalised equation (first non-zero coefficient of (a, b) scaled to 1).
Every triple of distinct points is assigned to exactly one block (line or circle).  We then report
   circles = #blocks that are circles,  and check the Mobius identity  circles = C(n,3) - D - ell.
"""
from fractions import Fraction as Fr
from itertools import combinations
from math import comb
import sys


def solve3(M, v):
    """Solve 3x3 rational system M z = v by Gaussian elimination with partial pivoting; None if singular."""
    A = [row[:] + [v[i]] for i, row in enumerate(M)]
    n = 3
    for c in range(n):
        piv = next((r for r in range(c, n) if A[r][c] != 0), None)
        if piv is None:
            return None
        A[c], A[piv] = A[piv], A[c]
        p = A[c][c]
        A[c] = [x / p for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return [A[i][n] for i in range(n)]


def block_key(p, q, r):
    (x1, y1), (x2, y2), (x3, y3) = p, q, r
    if (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1) == 0:      # collinear -> line
        a, b = y2 - y1, x1 - x2
        c = a * x1 + b * y1
        if a != 0:
            return ("line", Fr(1), b / a, c / a)
        return ("line", Fr(0), Fr(1), c / b)
    sol = solve3([[x1, y1, Fr(1)], [x2, y2, Fr(1)], [x3, y3, Fr(1)]],
                 [-(x1 * x1 + y1 * y1), -(x2 * x2 + y2 * y2), -(x3 * x3 + y3 * y3)])
    assert sol is not None
    return ("circle",) + tuple(sol)


def analyse(points):
    pts = [(Fr(x), Fr(y)) for x, y in points]
    n = len(pts)
    assert len(set(pts)) == n, "points not distinct"
    blocks = {}
    for i, j, k in combinations(range(n), 3):
        blocks.setdefault(block_key(pts[i], pts[j], pts[k]), set()).update((i, j, k))
    # sanity: every block's points really satisfy its equation, and every triple is in exactly one block
    for key, S in blocks.items():
        for i in S:
            x, y = pts[i]
            if key[0] == "line":
                assert key[1] * x + key[2] * y == key[3]
            else:
                assert x * x + y * y + key[1] * x + key[2] * y + key[3] == 0
    assert sum(comb(len(S), 3) for S in blocks.values()) == comb(n, 3)
    circles = {k: S for k, S in blocks.items() if k[0] == "circle"}
    lines = {k: S for k, S in blocks.items() if k[0] == "line"}
    degenerate = any(len(S) == n for S in blocks.values())
    D = sum(comb(len(S), 3) - 1 for S in blocks.values() if len(S) >= 4)
    ell = len(lines)
    return dict(n=n, circles=len(circles), lines=ell, D=D,
                mobius_identity=(comb(n, 3) - D - ell == len(circles)),
                circle_sizes=sorted(len(S) for S in circles.values()),
                line_sizes=sorted(len(S) for S in lines.values()),
                circle_blocks=sorted(sorted(S) for S in circles.values()),
                line_blocks=sorted(sorted(S) for S in lines.values()),
                degenerate=degenerate)


if __name__ == "__main__":
    P6 = [(0, 0), (20, 0), (5, 15), (5, 0), (2, 6), (10, 10)]
    P7 = P6 + [(5, 5)]
    P8 = [(0, 0), (0, 5), (0, 10), (0, 15), (3, 6), (3, 9), (5, 5), (5, 10)]
    for name, P in [("n=6 claimed 8", P6), ("n=7 claimed 11", P7), ("n=8 claimed 17", P8)]:
        r = analyse(P)
        print(name, ":", {k: v for k, v in r.items() if k not in ("circle_blocks", "line_blocks")})
        print("   circle blocks:", r["circle_blocks"])
        print("   line blocks:  ", r["line_blocks"])
