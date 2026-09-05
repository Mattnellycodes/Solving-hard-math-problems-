"""Independent exact audit of lattice-sphere configurations (written from scratch).

S = set of integer points on x^2+y^2+z^2 = N, Q = a point of the sphere not in S.
Stereographic projection from Q is a bijection sphere\{Q} -> plane taking circles on the sphere
(= plane sections) to circles (Q not on the plane) or lines (Q on the plane).  Hence
   circles(proj_Q(S)) = #{planes through >= 3 points of S, Q not on plane},
   lines = #{planes through >= 3 points of S, Q on plane}.
Degenerate iff S is coplanar.  All arithmetic is exact (integers / Fractions / quadratic surds).
"""
from itertools import combinations, product
from math import gcd, comb, isqrt
from fractions import Fraction as Fr
from collections import defaultdict

def lattice_points(N):
    r = isqrt(N)
    return [(x, y, z) for x in range(-r, r + 1) for y in range(-r, r + 1) for z in range(-r, r + 1)
            if x * x + y * y + z * z == N]

def plane_of(p, q, r):
    ux, uy, uz = q[0] - p[0], q[1] - p[1], q[2] - p[2]
    vx, vy, vz = r[0] - p[0], r[1] - p[1], r[2] - p[2]
    a, b, c = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
    if a == b == c == 0:
        return None  # collinear (cannot happen for 3 distinct sphere points)
    d = a * p[0] + b * p[1] + c * p[2]
    g = gcd(gcd(abs(a), abs(b)), gcd(abs(c), abs(d)))
    a, b, c, d = a // g, b // g, c // g, d // g
    if (a, b, c) < (0, 0, 0):
        a, b, c, d = -a, -b, -c, -d
    return (a, b, c, d)

def rich_planes(S):
    planes = defaultdict(set)
    for i, j, k in combinations(range(len(S)), 3):
        pl = plane_of(S[i], S[j], S[k])
        assert pl is not None
        planes[pl].update((i, j, k))
    assert sum(comb(len(v), 3) for v in planes.values()) == comb(len(S), 3)
    return planes

def on_plane(pl, Q):
    a, b, c, d = pl
    return a * Q[0] + b * Q[1] + c * Q[2] == d

def count_from_Q(S, Q, planes=None):
    planes = planes or rich_planes(S)
    assert Q not in S
    assert sum(q * q for q in Q) == sum(s * s for s in S[0])
    circ = [v for pl, v in planes.items() if not on_plane(pl, Q)]
    lin = [v for pl, v in planes.items() if on_plane(pl, Q)]
    degenerate = max(len(v) for v in planes.values()) == len(S)
    return dict(circles=len(circ), lines=len(lin), circle_sizes=sorted(map(len, circ), reverse=True),
                line_sizes=sorted(map(len, lin), reverse=True), degenerate=degenerate)

if __name__ == "__main__":
    S5 = [(0,-2,-1),(0,-2,1),(0,-1,-2),(0,-1,2),(0,1,-2),(0,1,2),(0,2,-1),(0,2,1),(2,1,0)]
    S6 = [(-2,1,-1),(-1,-2,-1),(-1,-2,1),(-1,-1,-2),(-1,-1,2),(-1,1,-2),(-1,1,2),(-1,2,-1),(-1,2,1)]
    for name, S, N in (("S2lat-5", S5, 5), ("S2lat-6", S6, 6)):
        assert all(sum(c * c for c in p) == N for p in S) and len(set(S)) == 9
        planes = rich_planes(S)
        print(name, "rich planes:", len(planes), "sizes:", sorted((len(v) for v in planes.values()), reverse=True))
        # every integer point of the sphere as projection centre
        best = None
        for Q in lattice_points(N):
            if Q in S:
                continue
            r = count_from_Q(S, Q, planes)
            if best is None or r['circles'] < best[1]['circles']:
                best = (Q, r)
            if r['circles'] <= 25:
                print("  Q =", Q, r)
        print("  best integer centre:", best)
