"""Independent exact circle/line counter (written from scratch for the c9 audit).
Representation: the unique (up to scale) 4-vector (a,b,c,d) with a(x^2+y^2)+bx+cy+d=0 through
three points = cofactor vector of the 3x4 matrix with rows (x^2+y^2, x, y, 1).  a==0 <=> collinear.
Normalised to a primitive integer vector with positive leading entry -> hashable key.
Checks: (i) antipodal 9-point configuration has 25 circles; (ii) Mobius identity
circles = C(n,3) - D - ell on several sets (including random rational sets)."""
from fractions import Fraction as Fr
from itertools import combinations
from math import comb, gcd
import random

def det3(M):
    return (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
            - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
            + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))

def key3(p, q, r):
    rows = [(x*x+y*y, x, y, Fr(1)) for x, y in (p, q, r)]
    v = []
    for j in range(4):
        M = [[row[k] for k in range(4) if k != j] for row in rows]
        v.append((-1)**j * det3(M))
    den = 1
    for t in v:
        den = den*t.denominator//gcd(den, t.denominator)
    w = [int(t*den) for t in v]
    g = 0
    for t in w: g = gcd(g, abs(t))
    assert g != 0, "zero cofactor vector (points not distinct?)"
    w = [t//g for t in w]
    for t in w:
        if t != 0:
            if t < 0: w = [-u for u in w]
            break
    return tuple(w)

def blocks_of(points):
    pts = [(Fr(x), Fr(y)) for x, y in points]
    assert len(set(pts)) == len(pts)
    blk = {}
    for i, j, k in combinations(range(len(pts)), 3):
        blk.setdefault(key3(pts[i], pts[j], pts[k]), set()).update((i, j, k))
    return blk

def summary(points):
    n = len(points)
    blk = blocks_of(points)
    circles = [s for k, s in blk.items() if k[0] != 0]
    lines = [s for k, s in blk.items() if k[0] == 0]
    D = sum(comb(len(s), 3) - 1 for s in blk.values() if len(s) >= 4)
    ell = len(lines)
    degenerate = any(len(s) == n for s in blk.values())
    return dict(n=n, circles=len(circles), lines=ell, D=D, identity=(len(circles) == comb(n, 3) - D - ell),
                circle_sizes=sorted(len(s) for s in circles), line_sizes=sorted(len(s) for s in lines),
                degenerate=degenerate)

if __name__ == "__main__":
    antipodal = [(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1), (Fr(3, 5), Fr(4, 5)), (Fr(-3, 5), Fr(-4, 5)),
                 (Fr(3, 5), Fr(-4, 5)), (Fr(-3, 5), Fr(4, 5))]
    print("antipodal n=9:", summary(antipodal))
    # a second 25-point witness: 8 points on a circle, centre NOT used, 9th point on 4 chords (non-diameters)?
    # simple sanity sets
    sets = {
        "n=8 record": [(0,0),(0,5),(0,10),(0,15),(3,6),(3,9),(5,5),(5,10)],
        "n=7 record": [(0,0),(20,0),(5,15),(5,0),(2,6),(10,10),(5,5)],
        "n=6 record": [(0,0),(20,0),(5,15),(5,0),(2,6),(10,10)],
        "two squares n=8": [(1,0),(0,1),(-1,0),(0,-1),(2,0),(0,2),(-2,0),(0,-2)],
    }
    for name, P in sets.items():
        print(name, summary(P))
    random.seed(7)
    bad = 0
    for trial in range(300):
        n = random.choice([6, 7, 8, 9])
        P = set()
        while len(P) < n:
            P.add((Fr(random.randint(-3, 3)), Fr(random.randint(-3, 3))))
        s = summary(list(P))
        if not s["identity"]:
            bad += 1; print("IDENTITY FAILS", P, s)
    print("random small-grid sets checked: 300, identity failures:", bad)
