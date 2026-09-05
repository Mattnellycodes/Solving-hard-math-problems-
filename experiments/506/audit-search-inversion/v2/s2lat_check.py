"""Independent check of the two lattice-sphere n=9 sets (S2lat-5, S2lat-6).

Method: invert R^3 about the rational sphere point Q (power 1): X' = Q + (X-Q)/|X-Q|^2.
The sphere through Q maps to a plane; the images are rational points of that plane.
Then count circles/lines in that plane using 3D exact arithmetic (circumcentre + radius^2
keys; line keys by primitive direction + foot of perpendicular from the origin).
This never uses plane-counting on the sphere, so it is a different computation from
the 3D 'integer plane counting' certificate of search-inversion.
"""
from fractions import Fraction as Fr
from itertools import combinations
from math import gcd, comb
from collections import defaultdict

def sub(a, b): return tuple(x - y for x, y in zip(a, b))
def dot(a, b): return sum(x * y for x, y in zip(a, b))
def cross(a, b): return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])

def invert(X, Q):
    d = sub(X, Q); s = dot(d, d)
    return tuple(q + Fr(v) / s for q, v in zip(Q, d))

def prim(v):
    den = 1
    from math import lcm
    for x in v: den = lcm(den, Fr(x).denominator)
    w = [int(Fr(x) * den) for x in v]
    g = 0
    for x in w: g = gcd(g, abs(x))
    w = [x // g for x in w]
    for x in w:
        if x != 0:
            if x < 0: w = [-y for y in w]
            break
    return tuple(w)

def solve3(M, rhs):
    # Cramer's rule over Fractions
    def det(m):
        return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1]) - m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
                + m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))
    D = det(M); assert D != 0
    out = []
    for c in range(3):
        Mc = [list(row) for row in M]
        for r in range(3): Mc[r][c] = rhs[r]
        out.append(det(Mc) / D)
    return tuple(out)

def count_coplanar(pts):
    n = len(pts)
    normal = None
    for X, Y, Z in combinations(pts, 3):
        c = cross(sub(Y, X), sub(Z, X))
        if any(c): normal = c; break
    assert normal is not None
    for X in pts: assert dot(sub(X, pts[0]), normal) == 0, "not coplanar"
    circles = defaultdict(set); lines = defaultdict(set)
    for i, j, k in combinations(range(n), 3):
        X, Y, Z = pts[i], pts[j], pts[k]
        c = cross(sub(Y, X), sub(Z, X))
        if not any(c):
            d = prim(sub(Y, X))
            foot = tuple(x - Fr(dot(X, d), dot(d, d)) * di for x, di in zip(X, d))
            lines[(d, foot)].update((i, j, k))
        else:
            M = [[2 * v for v in sub(Y, X)], [2 * v for v in sub(Z, X)], list(normal)]
            rhs = [dot(Y, Y) - dot(X, X), dot(Z, Z) - dot(X, X), dot(normal, X)]
            C = solve3(M, rhs)
            r2 = dot(sub(X, C), sub(X, C))
            assert r2 == dot(sub(Y, C), sub(Y, C)) == dot(sub(Z, C), sub(Z, C))
            circles[(C, r2)].update((i, j, k))
    assert sum(comb(len(s), 3) for s in list(circles.values()) + list(lines.values())) == comb(n, 3)
    csz = sorted((len(s) for s in circles.values()), reverse=True)
    lsz = sorted((len(s) for s in lines.values()), reverse=True)
    return len(circles), len(lines), csz, lsz

def check(name, N, S, Q):
    assert all(dot(p, p) == N for p in S) and dot(Q, Q) == N and Q not in S and len(set(S)) == len(S)
    img = [invert(p, Q) for p in S]
    nc, nl, csz, lsz = count_coplanar(img)
    deg = (csz and csz[0] == len(S)) or (lsz and lsz[0] == len(S))
    print(f"{name}: N={N} n={len(S)} centre={Q}: circles={nc} lines={nl} circle_sizes={csz} line_sizes={lsz} degenerate={bool(deg)}")
    return nc

S5 = [(0,-2,-1),(0,-2,1),(0,-1,-2),(0,-1,2),(0,1,-2),(0,1,2),(0,2,-1),(0,2,1),(2,1,0)]
S6 = [(-2,1,-1),(-1,-2,-1),(-1,-2,1),(-1,-1,-2),(-1,-1,2),(-1,1,-2),(-1,1,2),(-1,2,-1),(-1,2,1)]
if __name__ == "__main__":
    # certified.json exact centres, and the 'description' centres, for both sets
    check("S2lat-5 (certified centre)", 5, S5, (2, 0, -1))
    check("S2lat-5 (antipode of (2,1,0))", 5, S5, (-2, -1, 0))
    check("S2lat-6 (certified centre)", 6, S6, (-2, -1, 1))
    check("S2lat-6 (description centre)", 6, S6, (-2, -1, -1))
    # exhaustive over all lattice centres on the sphere not in S, for reference
    for name, N, S in (("S2lat-5", 5, S5), ("S2lat-6", 6, S6)):
        r = int(N ** 0.5) + 1
        best = None
        for x in range(-r, r + 1):
            for y in range(-r, r + 1):
                for z in range(-r, r + 1):
                    if x*x + y*y + z*z == N and (x, y, z) not in S:
                        img = [invert(p, (x, y, z)) for p in S]
                        nc, nl, csz, lsz = count_coplanar(img)
                        if best is None or nc < best[0]: best = (nc, (x, y, z), nl)
        print(f"{name}: best over lattice centres: circles={best[0]} at Q={best[1]} (lines={best[2]})")
