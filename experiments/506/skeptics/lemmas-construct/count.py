"""Independent exact circle counter (written from scratch by the lemmas-construct skeptic).
A block is identified by an exact canonical key: for circles the pair (centre, r^2) in Fractions;
for lines the primitive integer equation.  Also computes the Mobius decomposition
(blocks of size>=4, lines) and the D + ell bookkeeping to test identity (i) of the report."""
from fractions import Fraction as F
from itertools import combinations
from math import comb, gcd

def blocks_of(pts):
    pts = [(F(x), F(y)) for x, y in pts]
    n = len(pts)
    assert len(set(pts)) == n, "duplicate points"
    key2set = {}
    for i, j, k in combinations(range(n), 3):
        (x1, y1), (x2, y2), (x3, y3) = pts[i], pts[j], pts[k]
        # twice signed area
        det = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
        if det == 0:
            a, b = y2 - y1, x1 - x2
            c = -(a * x1 + b * y1)
            den = 1
            for q in (a, b, c):
                den = den * q.denominator // gcd(den, q.denominator)
            a, b, c = int(a * den), int(b * den), int(c * den)
            g = gcd(gcd(abs(a), abs(b)), abs(c))
            a, b, c = a // g, b // g, c // g
            if a < 0 or (a == 0 and b < 0):
                a, b, c = -a, -b, -c
            key = ('L', a, b, c)
        else:
            # centre via perpendicular bisectors
            s1 = x1 * x1 + y1 * y1; s2 = x2 * x2 + y2 * y2; s3 = x3 * x3 + y3 * y3
            D = 2 * det
            ux = ((s2 - s1) * (y3 - y1) - (s3 - s1) * (y2 - y1)) / D
            uy = ((x2 - x1) * (s3 - s1) - (x3 - x1) * (s2 - s1)) / D
            r2 = (x1 - ux) ** 2 + (y1 - uy) ** 2
            key = ('C', ux, uy, r2)
        key2set.setdefault(key, set()).update((i, j, k))
    return key2set

def analyse(pts):
    n = len(pts)
    bl = blocks_of(pts)
    circles = [frozenset(s) for k, s in bl.items() if k[0] == 'C']
    lines = [frozenset(s) for k, s in bl.items() if k[0] == 'L']
    degenerate = any(len(s) == n for s in bl.values())
    D = sum(comb(len(s), 3) - 1 for s in bl.values() if len(s) >= 4)
    ell = len(lines)
    identity_ok = (len(circles) == comb(n, 3) - D - ell)
    return dict(n=n, circles=len(circles), lines=ell, D=D, identity_ok=identity_ok,
                circle_sizes=sorted(len(s) for s in circles), line_sizes=sorted(len(s) for s in lines),
                degenerate=degenerate)

def f(n):
    return comb(n - 1, 2) + 1 - (n - 1) // 2

if __name__ == '__main__':
    R = {
     'n6': [(0,0),(20,0),(5,15),(5,0),(2,6),(10,10)],
     'n7': [(0,0),(20,0),(5,15),(5,0),(2,6),(10,10),(5,5)],
     'n8': [(0,0),(0,5),(0,10),(0,15),(3,6),(3,9),(5,5),(5,10)],
     'n8 squares': [(1,0),(0,1),(-1,0),(0,-1),(2,0),(0,2),(-2,0),(0,-2)],
     'n9 antipodal': [(0,0),(1,0),(-1,0),(0,1),(0,-1),(F(3,5),F(4,5)),(F(-3,5),F(-4,5)),(F(3,5),F(-4,5)),(F(-3,5),F(4,5))],
    }
    for name, P in R.items():
        a = analyse(P)
        print(name, 'f(n)=', f(a['n']), a)
