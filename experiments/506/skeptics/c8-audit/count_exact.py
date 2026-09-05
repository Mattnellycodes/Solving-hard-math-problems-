"""Independent exact circle counter (integer arithmetic only, no Fractions).
A circle/line through 3 integer points is identified by the primitive integer coefficient
vector of  A(x^2+y^2) + Bx + Cy + D = 0  (A = 0 <=> line), obtained from 3x3 minors of the
classical 4x4 circle determinant.  Every triple is assigned to exactly one key; a key is a
circle iff A != 0.  Cross-check: for each key we recompute the full incidence set by
evaluating the equation at all points."""
from itertools import combinations
from math import gcd
from functools import reduce
import sys

def det3(m):
    return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
          - m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
          + m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))

def norm(v):
    g = reduce(gcd, [abs(a) for a in v if a], 0)
    v = tuple(a // g for a in v)
    for a in v:
        if a:
            return v if a > 0 else tuple(-b for b in v)
    raise ValueError("zero vector")

def key_of(p, q, r):
    pts = [p, q, r]
    A = det3([[x, y, 1] for x, y in pts])
    if A == 0:
        (x1, y1), (x2, y2) = p, q
        return ('L',) + norm((0, y2 - y1, x1 - x2, x2 * y1 - x1 * y2))
    B = -det3([[x * x + y * y, y, 1] for x, y in pts])
    C = det3([[x * x + y * y, x, 1] for x, y in pts])
    D = -det3([[x * x + y * y, x, y] for x, y in pts])
    k = ('C',) + norm((A, B, C, D))
    for x, y in pts:  # sanity: the three points satisfy the equation
        assert k[1] * (x * x + y * y) + k[2] * x + k[3] * y + k[4] == 0
    return k

def on(k, pt):
    x, y = pt
    return k[1] * (x * x + y * y) + k[2] * x + k[3] * y + k[4] == 0

def analyze(P):
    P = [tuple(int(c) for c in p) for p in P]
    assert len(set(P)) == len(P)
    groups = {}
    for i, j, k in combinations(range(len(P)), 3):
        key = key_of(P[i], P[j], P[k])
        groups.setdefault(key, set()).update((i, j, k))
    # cross-check incidence sets by direct evaluation
    for key, S in groups.items():
        full = {i for i, p in enumerate(P) if on(key, p)}
        assert full == S, (key, S, full)
    circles = {k: v for k, v in groups.items() if k[0] == 'C'}
    lines = {k: v for k, v in groups.items() if k[0] == 'L'}
    return circles, lines

if __name__ == "__main__":
    P = [(0,0),(0,5),(0,10),(0,15),(3,6),(3,9),(5,5),(5,10)]
    circles, lines = analyze(P)
    n = len(P)
    print("points:", P)
    print("#circles (>=3 points):", len(circles))
    print("circle sizes:", sorted(len(v) for v in circles.values()))
    print("#lines (>=3 points):", len(lines), "line sizes:", sorted(len(v) for v in lines.values()))
    for k, v in lines.items():
        print("   line", k[1:], "points", sorted(P[i] for i in v))
    big = [v for v in list(circles.values()) + list(lines.values()) if len(v) >= 4]
    D = sum((len(v) * (len(v) - 1) * (len(v) - 2)) // 6 - 1 for v in big)
    print("blocks of size>=4:", len(big), " D =", D, " ell =", len(lines), " C(8,3)-D-ell =", 56 - D - len(lines))
    print("all on one circle/line?", any(len(v) == n for v in list(circles.values()) + list(lines.values())))
    # degrees of points in 4-blocks
    deg = [sum(1 for v in big if i in v) for i in range(n)]
    print("degrees in big blocks:", deg)
    # verify that two blocks share <= 2 points and lines share <= 1
    allb = list(circles.values()) + list(lines.values())
    print("max |B∩B'| =", max(len(a & b) for a, b in combinations(allb, 2)),
          " max |L∩L'| =", max((len(a & b) for a, b in combinations(list(lines.values()), 2)), default=0))
    # the same set scaled / translated should give the same count (sanity)
    Q = [(3*x + 7, 3*y - 4) for x, y in P]
    print("affine copy count:", len(analyze(Q)[0]))
    # sanity on known configurations: two concentric squares -> 18; antipodal 9 -> 25
    sq = [(1,0),(0,1),(-1,0),(0,-1),(2,0),(0,2),(-2,0),(0,-2)]
    print("two concentric squares (expected 18):", len(analyze(sq)[0]))
    ant = [(0,0),(5,0),(-5,0),(0,5),(0,-5),(3,4),(-3,-4),(3,-4),(-3,4)]
    print("antipodal n=9 (expected 25):", len(analyze(ant)[0]))
