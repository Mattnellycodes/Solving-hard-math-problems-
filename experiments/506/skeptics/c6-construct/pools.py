"""Generate structured point pools (doubles) for the floating-point 6-subset search.
Families: triangle centres/orthocentric systems (+ inversions), complete quadrilaterals
(+ Miquel inversion), complete quadrangles (cyclic and generic), regular polygons on concentric
circles, projections of Platonic/Archimedean solids.  Each pool is written as pool_<name>.txt.
"""
import numpy as np, itertools, math, os
rng = np.random.default_rng(506)
OUT = os.path.dirname(os.path.abspath(__file__))


def dedupe(pts, eps=1e-9):
    out = []
    for p in pts:
        p = np.asarray(p, float)
        if not np.all(np.isfinite(p)) or np.max(np.abs(p)) > 1e4:
            continue
        if all(np.linalg.norm(p - q) > eps for q in out):
            out.append(p)
    return out


def write(name, pts, tol=1e-7):
    pts = dedupe(pts)
    with open(os.path.join(OUT, f"pool_{name}.txt"), "w") as f:
        f.write(f"{len(pts)} 7 {tol}\n")
        for p in pts:
            f.write(f"{p[0]!r} {p[1]!r}\n")
    print(f"pool_{name}: {len(pts)} points, C(N,6)={math.comb(len(pts),6)}")


def invert(pts, c, r2=1.0):
    c = np.asarray(c, float); out = []
    for p in pts:
        d = p - c; n = d @ d
        if n > 1e-12:
            out.append(c + r2 * d / n)
    return out


def line_inter(p1, p2, p3, p4):
    d1 = p2 - p1; d2 = p4 - p3
    den = d1[0] * d2[1] - d1[1] * d2[0]
    if abs(den) < 1e-12:
        return None
    t = ((p3 - p1)[0] * d2[1] - (p3 - p1)[1] * d2[0]) / den
    return p1 + t * d1


def circumcentre(a, b, c):
    d = 2 * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]))
    if abs(d) < 1e-12:
        return None
    a2, b2, c2 = a @ a, b @ b, c @ c
    return np.array([(a2 * (b[1] - c[1]) + b2 * (c[1] - a[1]) + c2 * (a[1] - b[1])) / d,
                     (a2 * (c[0] - b[0]) + b2 * (a[0] - c[0]) + c2 * (b[0] - a[0])) / d])


def foot(p, a, b):
    d = b - a; t = ((p - a) @ d) / (d @ d); return a + t * d


def triangle_pool(A, B, C):
    A, B, C = map(lambda v: np.asarray(v, float), (A, B, C))
    O = circumcentre(A, B, C); G = (A + B + C) / 3; H = 3 * G - 2 * O
    D, E, F = foot(A, B, C), foot(B, C, A), foot(C, A, B)
    Ma, Mb, Mc = (B + C) / 2, (C + A) / 2, (A + B) / 2
    Nn = (O + H) / 2
    pts = [A, B, C, O, G, H, D, E, F, Ma, Mb, Mc, Nn, (A + H) / 2, (B + H) / 2, (C + H) / 2,
           2 * Ma - A, 2 * Mb - B, 2 * Mc - C,              # antipodes of A,B,C on circumcircle
           2 * D - H, 2 * E - H, 2 * F - H,                 # reflections of H in sides
           2 * Ma - H, 2 * Mb - H, 2 * Mc - H]              # reflections of H in midpoints
    a, b, c = np.linalg.norm(B - C), np.linalg.norm(C - A), np.linalg.norm(A - B)
    I = (a * A + b * B + c * C) / (a + b + c)
    Ia = (-a * A + b * B + c * C) / (-a + b + c); Ib = (a * A - b * B + c * C) / (a - b + c); Ic = (a * A + b * B - c * C) / (a + b - c)
    pts += [I, Ia, Ib, Ic, foot(I, B, C), foot(I, C, A), foot(I, A, B)]
    K = (a * a * A + b * b * B + c * c * C) / (a * a + b * b + c * c)   # symmedian point
    pts += [K, 2 * O - H, 2 * G - O, (O + G) / 2, 2 * Nn - O, 2 * A - O, 2 * B - O, 2 * C - O]
    pts += [2 * O - A, 2 * O - B, 2 * O - C]
    return pts


TRIANGLES = {
    "tri_random": (rng.normal(size=2), rng.normal(size=2), rng.normal(size=2)),
    "tri_345": ((0, 0), (4, 0), (0, 3)),
    "tri_isos40": ((0, 0), (1, 0), (0.5, 0.5 / math.tan(math.radians(20)))),
    "tri_306090": ((0, 0), (math.sqrt(3), 0), (0, 1)),
    "tri_equi": ((0, 0), (1, 0), (0.5, math.sqrt(3) / 2)),
    "tri_report_b1": ((0, 3), (1, 0), (3, 0)),
    "tri_obtuse": ((0, 0), (3, 0), (-1, 1)),
    "tri_random2": (rng.normal(size=2), rng.normal(size=2), rng.normal(size=2)),
}
for name, (A, B, C) in TRIANGLES.items():
    base = triangle_pool(A, B, C)
    write(name, base)
    A, B, C = map(lambda v: np.asarray(v, float), (A, B, C))
    O = circumcentre(A, B, C); H = 3 * (A + B + C) / 3 - 2 * O
    for cname, cen in (("H", H), ("O", O), ("A", A), ("N", (O + H) / 2), ("D", foot(A, B, C))):
        write(f"{name}_inv{cname}", base + invert(base, cen))


def quadrilateral_pool(lines):
    """lines: list of 4 (point, direction)."""
    P = {}
    for i, j in itertools.combinations(range(4), 2):
        p1, d1 = lines[i]; p2, d2 = lines[j]
        q = line_inter(p1, p1 + d1, p2, p2 + d2)
        if q is not None:
            P[(i, j)] = q
    pts = list(P.values())
    # diagonal midpoints (Gauss line)
    for (i, j), (k, l) in (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))):
        if (i, j) in P and (k, l) in P:
            pts.append((P[(i, j)] + P[(k, l)]) / 2)
    # the 4 triangles: circumcentres, orthocentres, centroids, incentres; Miquel point
    for omit in range(4):
        idx = [i for i in range(4) if i != omit]
        tri = [P.get((min(a, b), max(a, b))) for a, b in itertools.combinations(idx, 2)]
        if any(t is None for t in tri):
            continue
        A, B, C = tri
        O = circumcentre(A, B, C)
        if O is None:
            continue
        G = (A + B + C) / 3; H = 3 * G - 2 * O
        a, b, c = np.linalg.norm(B - C), np.linalg.norm(C - A), np.linalg.norm(A - B)
        pts += [O, G, H, (a * A + b * B + c * C) / (a + b + c), (O + H) / 2]
    return pts, P


def rand_line():
    return (rng.normal(size=2), np.array([math.cos(t := rng.uniform(0, math.pi)), math.sin(t)]))


QUADS = {
    "quad_random": [rand_line() for _ in range(4)],
    "quad_random2": [rand_line() for _ in range(4)],
    "quad_60": [(np.array([0., 0.]), np.array([1., 0.])), (np.array([0., 0.]), np.array([0.5, math.sqrt(3) / 2])),
                (np.array([1., 0.]), np.array([-0.5, math.sqrt(3) / 2])), (np.array([0.3, 0.9]), np.array([1., 0.2]))],
    "quad_sym": [(np.array([0., 0.]), np.array([1., 1.])), (np.array([0., 0.]), np.array([1., -1.])),
                 (np.array([2., 0.]), np.array([1., 2.])), (np.array([2., 0.]), np.array([1., -2.]))],
    "quad_orthic": [(np.array([0., 0.]), np.array([1., 0.])), (np.array([0., 0.]), np.array([1., 3.])),
                    (np.array([4., 0.]), np.array([-3., 3.])), (np.array([1., 0.]), np.array([0., 1.]))],
}
for name, L in QUADS.items():
    pts, P = quadrilateral_pool(L)
    write(name, pts)
    # Miquel point: second intersection of circumcircles of two of the triangles -> approximate by
    # circumcircle intersection; simpler: invert about each vertex and about a random point
    for k, v in list(P.items())[:3]:
        write(f"{name}_inv{k[0]}{k[1]}", pts + invert(pts, v))


def quadrangle_pool(Q):
    Q = [np.asarray(q, float) for q in Q]
    pts = list(Q)
    diag = []
    for (a, b), (c, d) in (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))):
        x = line_inter(Q[a], Q[b], Q[c], Q[d])
        if x is not None:
            diag.append(x)
    pts += diag
    for a, b in itertools.combinations(range(4), 2):
        pts.append((Q[a] + Q[b]) / 2)
    for tri in itertools.combinations(range(4), 3):
        A, B, C = (Q[i] for i in tri)
        O = circumcentre(A, B, C)
        if O is None:
            continue
        G = (A + B + C) / 3; pts += [O, G, 3 * G - 2 * O]
    if len(diag) == 3:
        pts += [(diag[0] + diag[1]) / 2, (diag[1] + diag[2]) / 2, (diag[0] + diag[2]) / 2, sum(diag) / 3]
    return pts


ang = rng.uniform(0, 2 * math.pi, 4)
cyc = [(math.cos(t), math.sin(t)) for t in ang]
# cyclic quadrangle with perpendicular diagonals: angles t1,t2,t3,t4 with t1+t3 = t2+t4 + pi (chords perpendicular)
t1, t2, t3 = ang[:3]; t4 = t1 + t3 - t2 - math.pi
QUADRANGLES = {
    "qr_random": [rng.normal(size=2) for _ in range(4)],
    "qr_cyclic": cyc + [(0, 0)],
    "qr_cyclic_perp": [(math.cos(t), math.sin(t)) for t in (t1, t2, t3, t4)] + [(0, 0)],
    "qr_square": [(0, 0), (1, 0), (1, 1), (0, 1)],
    "qr_orthocentric": [(0, 3), (1, 0), (3, 0), (0, -1)],
}
for name, Q in QUADRANGLES.items():
    extra = Q[4:]; Q4 = Q[:4]
    pts = quadrangle_pool(Q4) + [np.asarray(e, float) for e in extra]
    write(name, pts)
    write(f"{name}_inv0", pts + invert(pts, Q4[0]))
    write(f"{name}_invc", pts + invert(pts, np.mean(np.array(Q4, float), axis=0)))

# regular polygons on concentric circles
def ngon(n, r, phase=0.0):
    return [(r * math.cos(phase + 2 * math.pi * k / n), r * math.sin(phase + 2 * math.pi * k / n)) for k in range(n)]
write("poly12", ngon(24, 1) + ngon(12, 2) + ngon(12, 0.5) + ngon(12, 3, math.pi / 12) + [(0, 0)])
phi = (1 + math.sqrt(5)) / 2
write("poly10", ngon(10, 1) + ngon(10, phi) + ngon(10, 1 / phi) + ngon(10, phi * phi, math.pi / 10) + [(0, 0)])
write("poly8", ngon(16, 1) + ngon(8, math.sqrt(2)) + ngon(8, 2) + ngon(8, 1 / math.sqrt(2), math.pi / 8) + ngon(8, 1 + math.sqrt(2)) + [(0, 0)])
write("poly7", ngon(14, 1) + ngon(7, 2) + ngon(7, 1 / (2 * math.cos(math.pi / 7))) + ngon(7, 2 * math.cos(math.pi / 7), math.pi / 7) + [(0, 0)])

# polyhedra projections
def project(V, dirs):
    pts = []
    for d in dirs:
        d = np.asarray(d, float); d = d / np.linalg.norm(d)
        u = np.cross(d, [1, 0, 0] if abs(d[0]) < 0.9 else [0, 1, 0]); u /= np.linalg.norm(u); v = np.cross(d, u)
        pts += [(np.dot(p, u), np.dot(p, v)) for p in V]
    return pts


cube = [np.array(s) for s in itertools.product((-1, 1), repeat=3)]
octa = [np.array(p) for p in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))]
ico = [np.array(p) for p in itertools.chain.from_iterable(
    [((0, s1, s2 * phi), (s1, s2 * phi, 0), (s2 * phi, 0, s1)) for s1 in (-1, 1) for s2 in (-1, 1)])]
dode = cube + [np.array(p) for p in itertools.chain.from_iterable(
    [((0, s1 / phi, s2 * phi), (s1 / phi, s2 * phi, 0), (s2 * phi, 0, s1 / phi)) for s1 in (-1, 1) for s2 in (-1, 1)])]
cubocta = [np.array(p) for p in set(itertools.permutations((0, 1, 1))) | set(itertools.permutations((0, -1, 1))) | set(itertools.permutations((0, -1, -1))) | set(itertools.permutations((0, 1, -1)))]
dirs = [(0, 0, 1), (1, 1, 0), (1, 1, 1), (1, 2, 3), rng.normal(size=3)]
write("poly_cube_octa", project(cube + octa, dirs))
write("poly_ico", project(ico, dirs))
write("poly_dode", project(dode, dirs[:3]))
write("poly_cubocta", project(cubocta + [np.zeros(3)], dirs))
