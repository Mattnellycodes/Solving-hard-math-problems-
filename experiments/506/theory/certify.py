"""Self-contained exact certificate checker for the record configurations of Erdos #506.
Independent of circles_exact.py: circles are identified by the exact rational triple (cx, cy, r^2)
computed by solving the 2x2 linear system for the centre; lines by an exact normalised equation.
Prints, for each configuration: n, #circles (>=3 points), #lines (>=3 points), block sizes,
and the formula value f(n) = C(n-1,2)+1-floor((n-1)/2)."""
from fractions import Fraction as Fr
from itertools import combinations
from math import comb, gcd

def f(n): return comb(n - 1, 2) + 1 - (n - 1) // 2

def circle_key(a, b, c):
    (x1, y1), (x2, y2), (x3, y3) = a, b, c
    # centre (u,v): 2(x2-x1)u + 2(y2-y1)v = x2^2-x1^2+y2^2-y1^2, same with 3
    A11, A12, B1 = 2 * (x2 - x1), 2 * (y2 - y1), x2 * x2 - x1 * x1 + y2 * y2 - y1 * y1
    A21, A22, B2 = 2 * (x3 - x1), 2 * (y3 - y1), x3 * x3 - x1 * x1 + y3 * y3 - y1 * y1
    det = A11 * A22 - A12 * A21
    if det == 0:
        return None
    u = (B1 * A22 - B2 * A12) / det; v = (A11 * B2 - A21 * B1) / det
    return ("C", u, v, (x1 - u) ** 2 + (y1 - v) ** 2)

def line_key(a, b):
    (x1, y1), (x2, y2) = a, b
    A, B, C = y2 - y1, x1 - x2, (y2 - y1) * x1 + (x1 - x2) * y1   # A x + B y = C
    # normalise to primitive integers with positive leading coefficient
    den = 1
    for q in (A, B, C): den = den * q.denominator // gcd(den, q.denominator)
    A, B, C = int(A * den), int(B * den), int(C * den)
    g = gcd(gcd(abs(A), abs(B)), abs(C)); A, B, C = A // g, B // g, C // g
    if A < 0 or (A == 0 and B < 0): A, B, C = -A, -B, -C
    return ("L", A, B, C)

def analyse(pts):
    pts = [(Fr(x), Fr(y)) for x, y in pts]
    n = len(pts); assert len(set(pts)) == n
    blocks = {}
    for i, j, k in combinations(range(n), 3):
        key = circle_key(pts[i], pts[j], pts[k])
        if key is None:
            key = line_key(pts[i], pts[j])
        blocks.setdefault(key, set()).update((i, j, k))
    circles = [b for k, b in blocks.items() if k[0] == "C"]
    lines = [b for k, b in blocks.items() if k[0] == "L"]
    degenerate = any(len(b) == n for b in blocks.values())
    return len(circles), len(lines), sorted(len(b) for b in circles), sorted(len(b) for b in lines), degenerate

RECORDS = {
 "n=6 orthic (triangle + feet of altitudes)": [(0,0),(4,0),(1,3),(1,0),(Fr(2,5),Fr(6,5)),(2,2)],
 "n=7 orthocentric (triangle + feet + orthocentre)": [(0,0),(4,0),(1,3),(1,0),(Fr(2,5),Fr(6,5)),(2,2),(1,1)],
 "n=8 inverted orthocentric system": [(0,Fr(3,2)),(Fr(1,2),Fr(1,2)),(Fr(3,10),Fr(9,10)),(0,Fr(1,2)),(Fr(3,10),Fr(3,5)),(Fr(1,2),1),(0,0),(0,1)],
 "n=8 two concentric squares (18, old record)": [(1,0),(0,1),(-1,0),(0,-1),(2,0),(0,2),(-2,0),(0,-2)],
 "n=9 antipodal (formula)": [(0,0),(1,0),(-1,0),(0,1),(0,-1),(Fr(3,5),Fr(4,5)),(Fr(-3,5),Fr(-4,5)),(Fr(3,5),Fr(-4,5)),(Fr(-3,5),Fr(4,5))],
}
if __name__ == "__main__":
    for name, P in RECORDS.items():
        c, l, cs, ls, deg = analyse(P)
        print(f"{name}: n={len(P)} circles={c} (formula {f(len(P))}), 3+-lines={l}, circle sizes={cs}, line sizes={ls}, degenerate={deg}")
