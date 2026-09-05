"""Integer points on spheres x^2 + y^2 + z^2 = N as Moebius universes, with an EXACT certificate.

A finite set V of points on a sphere S is a Moebius configuration; its blocks are the planes
containing >= 3 points of V (a plane meets S in a circle).  Stereographic projection from a sphere
point Q not in V maps the blocks through Q to straight lines and every other block to a circle, so
the projected planar set determines exactly

        circles = #(planes through >= 3 points of V) - #(such planes through Q).

For integer V and rational Q both numbers are computed with exact integer/rational arithmetic
(coplanarity = vanishing 3x3 determinant), which certifies the circle count of the planar set
without ever writing down its (in general irrational) planar coordinates.  Planar coordinates in
Q(sqrt(N)) can be produced with `planar_coords` and re-certified with exact_verify.certify.
"""
import itertools, math
from fractions import Fraction as Fr
from collections import defaultdict
import numpy as np


def integer_sphere_points(N):
    """All integer (x,y,z) with x^2+y^2+z^2 = N."""
    r = math.isqrt(N)
    out = []
    for x in range(-r, r + 1):
        for y in range(-r, r + 1):
            z2 = N - x * x - y * y
            if z2 < 0:
                continue
            z = math.isqrt(z2)
            if z * z == z2:
                if z == 0:
                    out.append((x, y, 0))
                else:
                    out.append((x, y, z)); out.append((x, y, -z))
    return sorted(set(out))


def _plane_key(p, q, r):
    """Canonical integer key of the plane through three non-collinear rational points."""
    ux, uy, uz = (q[0] - p[0], q[1] - p[1], q[2] - p[2])
    vx, vy, vz = (r[0] - p[0], r[1] - p[1], r[2] - p[2])
    nx = uy * vz - uz * vy; ny = uz * vx - ux * vz; nz = ux * vy - uy * vx   # cross product
    if nx == 0 and ny == 0 and nz == 0:
        return None
    d = nx * p[0] + ny * p[1] + nz * p[2]
    # normalise: rational -> integers, divide by gcd, fix sign
    vals = [Fr(nx), Fr(ny), Fr(nz), Fr(d)]
    den = 1
    for v in vals:
        den = den * v.denominator // math.gcd(den, v.denominator)
    ints = [int(v * den) for v in vals]
    g = 0
    for v in ints:
        g = math.gcd(g, abs(v))
    ints = [v // g for v in ints]
    for v in ints[:3]:
        if v != 0:
            if v < 0:
                ints = [-w for w in ints]
            break
    return tuple(ints)


def planes_exact(points):
    """All planes through >= 3 of the (rational) points, as {plane_key: set(indices)}.
    Assumes no three points collinear (true on a sphere)."""
    pts = [tuple(Fr(c) for c in p) for p in points]
    assert len(set(pts)) == len(pts), "duplicate points"
    planes = defaultdict(set)
    n = len(pts)
    for i, j, k in itertools.combinations(range(n), 3):
        key = _plane_key(pts[i], pts[j], pts[k])
        assert key is not None, "collinear triple on a sphere?!"
        planes[key].update((i, j, k))
    return dict(planes)


def on_plane(key, Q):
    a, b, c, d = key
    return a * Q[0] + b * Q[1] + c * Q[2] == d


def certify_sphere(points, Q):
    """Exact circle count of the stereographic projection of `points` (rational, on a common
    sphere) from the rational sphere point Q (not in the set).  Returns dict."""
    pts = [tuple(Fr(c) for c in p) for p in points]
    Q = tuple(Fr(c) for c in Q)
    assert Q not in pts, "projection centre coincides with a point"
    # sanity: all on one sphere through Q?  (centre-free check: |p|^2 equal after common centre)
    planes = planes_exact(pts)
    n = len(pts)
    degenerate = any(len(s) == n for s in planes.values())
    through = [k for k in planes if on_plane(k, Q)]
    circles = len(planes) - len(through)
    sizes_c = sorted(len(planes[k]) for k in planes if not on_plane(k, Q))
    sizes_l = sorted(len(planes[k]) for k in through)
    return dict(n=n, blocks=len(planes), lines=len(through), circles=circles,
                circle_sizes=sizes_c, line_sizes=sizes_l, degenerate=degenerate)


def sphere_check(points, Q=None):
    """Verify that all points (and Q) lie on one sphere centred at the origin; return radius^2."""
    r2 = set(sum(Fr(c) ** 2 for c in p) for p in points)
    if Q is not None:
        r2.add(sum(Fr(c) ** 2 for c in Q))
    assert len(r2) == 1, f"points not on one origin-centred sphere: {r2}"
    return r2.pop()


def planar_coords(points, Q, N=None):
    """Exact planar coordinates (sympy, in Q(sqrt(N)) at worst) of the stereographic projection of
    the origin-centred sphere points from Q: inversion in the sphere of centre Q and radius^2 2N
    maps the sphere onto the plane tangent at -Q; we then express the image in an orthonormal
    frame of that plane.  Returns (coords, radicand) where coords are sympy pairs."""
    import sympy as sp
    if N is None:
        N = sphere_check(points, Q)
    N = sp.Rational(N)
    Qv = sp.Matrix([sp.Rational(c) for c in Q])
    # orthonormal frame of the plane orthogonal to Q
    q = Qv / sp.sqrt(N)
    # pick a vector not parallel to q
    a = sp.Matrix([1, 0, 0]) if abs(float(q[0])) < 0.9 else sp.Matrix([0, 1, 0])
    e1 = a - (a.dot(q)) * q; e1 = e1 / sp.sqrt(e1.dot(e1))
    e2 = q.cross(e1)
    out = []
    for p in points:
        X = sp.Matrix([sp.Rational(c) for c in p])
        D = X - Qv
        Y = Qv + (2 * N / D.dot(D)) * D          # image on the tangent plane at -Q
        out.append((sp.nsimplify(sp.simplify(e1.dot(Y))), sp.nsimplify(sp.simplify(e2.dot(Y)))))
    return out


# ---------------------------------------------------------------- universe builders
def lattice_sphere_universe(N):
    P = integer_sphere_points(N)
    V = np.array(P, float) / math.sqrt(N)
    labels = [f"({x},{y},{z})" for x, y, z in P]
    return (f"S2lat-{N}", V, labels, {'points3d': P, 'N': N})


def lattice_sphere_universe_union(Ns):
    """Points of several spheres cannot be on one sphere; this is NOT a Moebius universe.
    (Kept only to document the negative: not used.)"""
    raise NotImplementedError


if __name__ == "__main__":
    # certificate sanity: cube (N=3) projected from a degree-2 point (0,0,sqrt3) is irrational; use
    # instead N = 9: 30 points; pick Q = (3,0,0)?  That is a configuration point.  Use N = 9 subset
    # test with Q on the sphere: rational points of x^2+y^2+z^2 = 9: e.g. (1,2,2).
    for N in (1, 2, 3, 5, 6, 9, 11, 14, 17, 18, 21, 26, 27, 29, 30, 33, 35, 38, 41):
        P = integer_sphere_points(N)
        pl = planes_exact(P)
        sizes = defaultdict(int)
        for s in pl.values():
            sizes[len(s)] += 1
        print(f"N={N:3d}: {len(P):3d} points, {len(pl):5d} planes, sizes {dict(sorted(sizes.items()))}")
    # cube from N=3, Q=(1,1,1)?? that's a vertex.  Test: N=9 subset = cube (±1,±2,±2)?  not a cube.
    # Test the certificate against exact_verify on a rational example: N = 9, S = 8 points, Q rational.
    P = integer_sphere_points(9)
    S = [p for p in P if p[0] == 3 or p[2] == -3 or (abs(p[0]) == 1 and p[1] == 2)]
    Q = (Fr(3), Fr(0), Fr(0))
    S = [p for p in S if tuple(Fr(c) for c in p) != Q][:9]
    print("test subset", S, "Q", Q, certify_sphere(S, Q))
    from exact_verify import certify
    coords = planar_coords(S, Q, N=9)
    print("planar coords:", coords)
    print("re-certified:", {k: v for k, v in certify(coords).items() if k != 'coords'})
