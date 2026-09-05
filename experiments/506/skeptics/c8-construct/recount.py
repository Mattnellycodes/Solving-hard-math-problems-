"""Independent exact recount of circles for a rational point set.
Method: for every triple, compute the circle in normalised INTEGER form
(a x^2 + a y^2 + b x + c y + d = 0 with a>0, gcd=1) using integer determinants;
distinct keys = distinct circles.  Also cross-check via quadruple concyclicity
determinants (no circle solving at all)."""
from fractions import Fraction as Fr
from itertools import combinations
from math import gcd
from functools import reduce

def circle_key(p, q, r):
    # general circle a(x^2+y^2)+bx+cy+d=0 through p,q,r: solve via cofactors of
    # | x^2+y^2  x  y  1 | rows.
    rows = [(x*x+y*y, x, y, 1) for (x, y) in (p, q, r)]
    def det3(m):
        return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
              - m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
              + m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))
    cols = [0, 1, 2, 3]
    coef = []
    for j in range(4):
        sub = [[rows[i][k] for k in cols if k != j] for i in range(3)]
        coef.append((-1)**j * det3(sub))
    a, b, c, d = coef
    if a == 0:
        return None  # collinear
    g = reduce(gcd, [abs(a), abs(b), abs(c), abs(d)])
    if a < 0:
        g = -g
    return (a//g, b//g, c//g, d//g)

def count(points):
    # scale to integers
    pts = [(Fr(x), Fr(y)) for x, y in points]
    den = 1
    for x, y in pts:
        den = den * x.denominator // gcd(den, x.denominator)
        den = den * y.denominator // gcd(den, y.denominator)
    ipts = [(int(x*den), int(y*den)) for x, y in pts]
    assert len(set(ipts)) == len(ipts)
    circles = {}
    lines = 0
    for t in combinations(range(len(ipts)), 3):
        k = circle_key(*[ipts[i] for i in t])
        if k is None:
            lines += 1
        else:
            circles.setdefault(k, set()).update(t)
    return circles, lines, ipts

def concyclic(p, q, r, s):
    m = [(x*x+y*y, x, y, 1) for (x, y) in (p, q, r, s)]
    # 4x4 determinant
    def det(m):
        if len(m) == 1: return m[0][0]
        return sum((-1)**j * m[0][j] * det([row[:j]+row[j+1:] for row in m[1:]]) for j in range(len(m)))
    return det(m) == 0

def collinear(p, q, r):
    return (q[0]-p[0])*(r[1]-p[1]) - (q[1]-p[1])*(r[0]-p[0]) == 0

if __name__ == "__main__":
    P = [(0,0),(0,5),(0,10),(0,15),(3,6),(3,9),(5,5),(5,10)]
    circles, lines, ipts = count(P)
    print("points:", ipts)
    print("collinear triples:", lines)
    print("number of circles:", len(circles))
    sizes = sorted(len(s) for s in circles.values())
    print("circle sizes:", sizes)
    for k, s in sorted(circles.items(), key=lambda kv: -len(kv[1])):
        print("  ", k, sorted(s))
    # Cross-check 2: build blocks purely from concyclicity/collinearity tests on quadruples
    n = len(ipts)
    # union-find-free approach: for each triple, the block = triple + all s with concyclic
    blocks = set()
    ntri_col = 0
    for t in combinations(range(n), 3):
        p, q, r = [ipts[i] for i in t]
        if collinear(p, q, r):
            ntri_col += 1
            blk = set(t) | {s for s in range(n) if s not in t and collinear(p, q, ipts[s])}
            blocks.add(("L", frozenset(blk)))
        else:
            blk = set(t) | {s for s in range(n) if s not in t and concyclic(p, q, r, ipts[s])}
            blocks.add(("C", frozenset(blk)))
    ncirc = sum(1 for k, b in blocks if k == "C")
    nlin = sum(1 for k, b in blocks if k == "L")
    print("cross-check: circles =", ncirc, " lines(>=3 pts) =", nlin, " line blocks:", [sorted(b) for k, b in blocks if k == "L"])
    # Also check not all on one circle / line
    print("all concyclic?", any(len(s) == n for s in circles.values()))
