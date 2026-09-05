"""Independent exact circle counter (written from scratch for the c9 audit).
A circle/line through 3 points is represented by the primitive integer vector (A,B,C,D) of
A(x^2+y^2) + Bx + Cy + D = 0 (A=0 <=> line).  Rational coordinates are cleared to integers first,
so everything is exact integer arithmetic (no Fractions needed after scaling)."""
from fractions import Fraction as Fr
from itertools import combinations
from math import gcd, comb

def to_int_points(points):
    pts = [(Fr(x), Fr(y)) for x, y in points]
    den = 1
    for x, y in pts:
        den = den * x.denominator // gcd(den, x.denominator)
        den = den * y.denominator // gcd(den, y.denominator)
    return [(int(x * den), int(y * den)) for x, y in pts]

def block_key(p, q, r):
    """primitive (A,B,C,D) of the circle/line through p,q,r (None if two coincide)."""
    # Solve for the 3x4 nullspace by cofactors of the matrix rows (x^2+y^2, x, y, 1)
    rows = [(x * x + y * y, x, y, 1) for x, y in (p, q, r)]
    def det3(cols):
        m = [[rows[i][c] for c in cols] for i in range(3)]
        return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))
    A = det3((1, 2, 3)); B = -det3((0, 2, 3)); C = det3((0, 1, 3)); D = -det3((0, 1, 2))
    g = gcd(gcd(abs(A), abs(B)), gcd(abs(C), abs(D)))
    if g == 0:
        return None
    A, B, C, D = A // g, B // g, C // g, D // g
    for v in (A, B, C, D):
        if v != 0:
            if v < 0: A, B, C, D = -A, -B, -C, -D
            break
    return (A, B, C, D)

def blocks(points):
    P = to_int_points(points)
    assert len(set(P)) == len(P), "duplicate points"
    bl = {}
    for i, j, k in combinations(range(len(P)), 3):
        key = block_key(P[i], P[j], P[k])
        bl.setdefault(key, set()).update((i, j, k))
    return bl

def analyse(points):
    bl = blocks(points)
    n = len(points)
    circ = {k: v for k, v in bl.items() if k[0] != 0}
    lines = {k: v for k, v in bl.items() if k[0] == 0}
    degenerate = any(len(v) == n for v in bl.values())
    return {"n": n, "circles": len(circ), "lines": len(lines), "degenerate": degenerate,
            "circle_sizes": sorted((len(v) for v in circ.values()), reverse=True),
            "line_sizes": sorted((len(v) for v in lines.values()), reverse=True),
            "blocks": bl}

def formula(n):
    return comb(n - 1, 2) + 1 - (n - 1) // 2

if __name__ == "__main__":
    # antipodal configuration: 8 points on unit circle in 4 antipodal pairs + centre
    anti = [(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1), (Fr(3, 5), Fr(4, 5)), (Fr(-3, 5), Fr(-4, 5)),
            (Fr(3, 5), Fr(-4, 5)), (Fr(-3, 5), Fr(4, 5))]
    r = analyse(anti)
    print("antipodal n=9:", {k: v for k, v in r.items() if k != "blocks"}, "formula", formula(9))
    # generic antipodal: chords through a non-centre point
    # 8 points on unit circle (rational), 9th point = intersection of 4 chords? use centre variant with
    # 4 chords through (0,0) that are NOT diameters is impossible on a circle centred at 0; instead test
    # 8 concyclic points with a 9th point on 4 chords for a non-centre point via the theory's claim later.
    for name, P in {"n=6 rec": [(0,0),(20,0),(5,15),(5,0),(2,6),(10,10)],
                    "n=7 rec": [(0,0),(20,0),(5,15),(5,0),(2,6),(10,10),(5,5)],
                    "n=8 rec": [(0,0),(0,5),(0,10),(0,15),(3,6),(3,9),(5,5),(5,10)]}.items():
        r = analyse(P); print(name, {k: v for k, v in r.items() if k != "blocks"}, "formula", formula(len(P)))
