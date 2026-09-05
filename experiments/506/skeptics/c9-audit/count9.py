"""Exact circle counter written from scratch for the c(9) audit (upper bound c(9) <= 25).

A circle/line through three points p,q,r is identified by the primitive integer vector
(A,B,C,Dd) of the equation A(x^2+y^2) + Bx + Cy + Dd = 0 (A = 0 <=> line), obtained as the
cofactors of the 3x4 matrix with rows (x^2+y^2, x, y, 1).  Exact rational arithmetic throughout.
Also verifies the Moebius identity  circles = C(n,3) - D - ell  on the tested sets and on random sets.
"""
from fractions import Fraction as Fr
from itertools import combinations
from math import comb, gcd
import random, sys

def det3(a, b, c):
    return (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]))

def block_key(p, q, r):
    rows = [(x * x + y * y, x, y, Fr(1)) for (x, y) in (p, q, r)]
    cof = []
    for j in range(4):
        sub = [tuple(v for k, v in enumerate(row) if k != j) for row in rows]
        cof.append((-1) ** j * det3(*sub))
    den = 1
    for v in cof:
        den = den * v.denominator // gcd(den, v.denominator)
    w = [int(v * den) for v in cof]
    g = 0
    for v in w:
        g = gcd(g, abs(v))
    assert g > 0, "points not distinct"
    w = [v // g for v in w]
    lead = next(v for v in w if v != 0)
    if lead < 0:
        w = [-v for v in w]
    return tuple(w)

def blocks(points):
    pts = [(Fr(x), Fr(y)) for (x, y) in points]
    assert len(set(pts)) == len(pts), "duplicate points"
    B = {}
    for i, j, k in combinations(range(len(pts)), 3):
        B.setdefault(block_key(pts[i], pts[j], pts[k]), set()).update((i, j, k))
    return B

def report(points, name=""):
    n = len(points)
    B = blocks(points)
    circles = {k: s for k, s in B.items() if k[0] != 0}
    lines = {k: s for k, s in B.items() if k[0] == 0}
    D = sum(comb(len(s), 3) - 1 for s in B.values() if len(s) >= 4)
    ell = len(lines)
    ident = (len(circles) == comb(n, 3) - D - ell)
    degenerate = any(len(s) == n for s in B.values())
    print(f"{name}: n={n} circles={len(circles)} lines(>=3 pts)={ell} D={D} "
          f"identity C(n,3)-D-ell holds={ident} degenerate(all on one circle/line)={degenerate}")
    print(f"   circle sizes={sorted(len(s) for s in circles.values())} line sizes={sorted(len(s) for s in lines.values())}")
    return len(circles), degenerate

if __name__ == "__main__":
    # antipodal configuration: 4 antipodal pairs on the unit circle + centre
    pyth = [(Fr(1), Fr(0)), (Fr(0), Fr(1)), (Fr(3, 5), Fr(4, 5)), (Fr(-5, 13), Fr(12, 13))]
    anti = [(Fr(0), Fr(0))] + [p for (x, y) in pyth for p in ((x, y), (-x, -y))]
    c1, d1 = report(anti, "antipodal 9-point set")
    # non-symmetric witness: 8 points on a circle, 9th point on 4 chords which are NOT diameters
    # (chords through an interior point q=(1/3,0): each rational point p on the unit circle gives the
    # second intersection of line pq with the circle)
    def second_intersection(p, q):
        (px, py), (qx, qy) = p, q
        dx, dy = px - qx, py - qy
        # points q + t(d): |q + t d|^2 = 1 -> (d.d) t^2 + 2 (q.d) t + (q.q - 1) = 0, t=1 is a root
        a = dx * dx + dy * dy; b = 2 * (qx * dx + qy * dy)
        t2 = -b / a - 1           # sum of roots = -b/a, one root is 1
        return (qx + t2 * dx, qy + t2 * dy)
    q = (Fr(1, 3), Fr(0))
    base = [(Fr(1), Fr(0)), (Fr(3, 5), Fr(4, 5)), (Fr(-5, 13), Fr(12, 13)), (Fr(8, 17), Fr(-15, 17))]
    chords = [q] + [p for b in base for p in (b, second_intersection(b, q))]
    c2, d2 = report(chords, "8 concyclic points + interior point on 4 non-diameter chords")
    assert all(x * x + y * y == 1 for (x, y) in chords[1:]), "not on the unit circle"
    assert not d1 and not d2
    print("upper bound c(9) <=", min(c1, c2))
    # random sanity check of the Moebius identity
    random.seed(20260905)
    bad = 0
    for trial in range(200):
        n = random.choice([6, 7, 8, 9])
        P = set()
        while len(P) < n:
            P.add((Fr(random.randint(-3, 3)), Fr(random.randint(-3, 3))))
        B = blocks(list(P))
        circles = sum(1 for k in B if k[0] != 0)
        D = sum(comb(len(s), 3) - 1 for s in B.values() if len(s) >= 4)
        ell = sum(1 for k in B if k[0] == 0)
        if circles != comb(n, 3) - D - ell:
            bad += 1
    print("random grid sets: 200 tested, identity failures:", bad)
