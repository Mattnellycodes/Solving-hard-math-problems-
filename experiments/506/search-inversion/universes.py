"""Rich Moebius configurations ("universes") for the subset + inversion search.

Each builder returns (name, planar_points, labels, exact) where planar_points is a list of (x, y)
pairs (Fractions / floats) or 'inf', and exact is None or a dict {'planar': [(sympy, sympy) | 'inf'],
'radicands': (...)} that allows exact certification.
"""
from fractions import Fraction as Fr
import math, itertools
import sympy
from sympy import Rational, sqrt, nsimplify


def _dedupe(pts, labels, tol=1e-9):
    out, lab = [], []
    for p, l in zip(pts, labels):
        dup = False
        for q in out:
            if p == 'inf' or q == 'inf':
                if p == q:
                    dup = True; break
                continue
            if abs(float(p[0]) - float(q[0])) < tol and abs(float(p[1]) - float(q[1])) < tol:
                dup = True; break
        if not dup:
            out.append(p); lab.append(l)
    return out, lab


def _issq(q):
    """exact square root of a Fraction if it is a perfect square, else None."""
    if isinstance(q, Fr):
        if q < 0:
            return None
        n, d = q.numerator, q.denominator
        rn, rd = math.isqrt(n), math.isqrt(d)
        if rn * rn == n and rd * rd == d:
            return Fr(rn, rd)
        return None
    return math.sqrt(q)


# ---------------------------------------------------------------- triangle universe
def triangle_universe(A, B, C, incentre=True, name=None, extra=True):
    """All the classical special points of triangle ABC (exact if A,B,C rational and, for the
    incentre family, the triangle is Heronian)."""
    def V(p): return (p[0], p[1])
    def add(p, q): return (p[0] + q[0], p[1] + q[1])
    def sub(p, q): return (p[0] - q[0], p[1] - q[1])
    def mul(p, s): return (p[0] * s, p[1] * s)
    def dot(p, q): return p[0] * q[0] + p[1] * q[1]
    exact_ok = all(isinstance(c, Fr) for P in (A, B, C) for c in P)
    A, B, C = V(A), V(B), V(C)
    pts, lab = [], []
    def put(p, l):
        pts.append(p); lab.append(l)
    put(A, 'A'); put(B, 'B'); put(C, 'C')
    # circumcentre
    ax, ay = A; bx, by = B; cx, cy = C
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    a2 = ax * ax + ay * ay; b2 = bx * bx + by * by; c2 = cx * cx + cy * cy
    O = ((a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d, (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d)
    H = sub(add(add(A, B), C), mul(O, 2))
    G = mul(add(add(A, B), C), Fr(1, 3) if exact_ok else 1 / 3)
    N = mul(add(O, H), Fr(1, 2) if exact_ok else 0.5)
    half = Fr(1, 2) if exact_ok else 0.5
    put(H, 'H'); put(O, 'O'); put(N, 'N'); put(G, 'G')
    verts = {'A': A, 'B': B, 'C': C}
    names = ['A', 'B', 'C']
    for i, X in enumerate(names):
        Y, Z = names[(i + 1) % 3], names[(i + 2) % 3]
        P, Q, R = verts[X], verts[Y], verts[Z]
        # foot of altitude from X on YZ
        QR = sub(R, Q)
        foot = add(Q, mul(QR, dot(sub(P, Q), QR) / dot(QR, QR)))
        put(foot, f'D{X.lower()}')
        M = mul(add(Q, R), half); put(M, f'M{X.lower()}')
        put(mul(add(P, H), half), f'E{X.lower()}')           # Euler point
        put(sub(mul(foot, 2), H), f'H{X.lower()}')           # reflection of H in YZ (on circumcircle)
        put(sub(mul(O, 2), P), f"{X}'")                       # antipode of X
        put(sub(add(Q, R), O), f'O{X.lower()}')              # reflection of O in YZ
    if incentre:
        a = _issq(dot(sub(B, C), sub(B, C))); b = _issq(dot(sub(A, C), sub(A, C))); c = _issq(dot(sub(A, B), sub(A, B)))
        if a is None or b is None or c is None:
            # not Heronian: use floats for these points
            a = math.dist(B, C); b = math.dist(A, C); c = math.dist(A, B)
            exact_ok = False
            A = (float(A[0]), float(A[1])); B = (float(B[0]), float(B[1])); C = (float(C[0]), float(C[1]))
        s = (a + b + c) / 2
        I = mul(add(add(mul(A, a), mul(B, b)), mul(C, c)), 1 / (a + b + c))
        Ia = mul(add(add(mul(A, -a), mul(B, b)), mul(C, c)), 1 / (-a + b + c))
        Ib = mul(add(add(mul(A, a), mul(B, -b)), mul(C, c)), 1 / (a - b + c))
        Ic = mul(add(add(mul(A, a), mul(B, b)), mul(C, -c)), 1 / (a + b - c))
        put(I, 'I'); put(Ia, 'Ia'); put(Ib, 'Ib'); put(Ic, 'Ic')
        for X, Ix, Iy, Iz in (('a', Ia, Ib, Ic), ('b', Ib, Ic, Ia), ('c', Ic, Ia, Ib)):
            put(mul(add(I, Ix), half), f'W{X}')        # arc midpoint (midpoint of I I_x)
            put(mul(add(Iy, Iz), half), f"W{X}'")      # opposite arc midpoint
        # feet of bisectors and incircle contact points
        for X, P, Q, R, la, lb, lc in (('a', A, B, C, a, b, c), ('b', B, C, A, b, c, a), ('c', C, A, B, c, a, b)):
            # bisector from P meets QR at (lb*R + lc*Q)/(lb+lc)?  side lengths: |PR| = lb? use: L = (|PR| Q + |PQ| R)/(|PR|+|PQ|)
            # with P=A: |AC| = b, |AB| = c: L = (b*B + c*C)/(b+c)
            # generic: for vertex P with opposite side QR, the other sides are PQ and PR.
            PQ = math.sqrt(float(dot(sub(P, Q), sub(P, Q)))); PR = math.sqrt(float(dot(sub(P, R), sub(P, R))))
            PQe = _issq(dot(sub(P, Q), sub(P, Q))) if exact_ok else PQ
            PRe = _issq(dot(sub(P, R), sub(P, R))) if exact_ok else PR
            L = mul(add(mul(Q, PRe), mul(R, PQe)), 1 / (PRe + PQe))
            put(L, f'L{X}')
            # contact point of the incircle on QR: distance from Q equals s - |PR|... (tangent lengths: from Q: s - PR)
            QR = sub(R, Q); lenQR = _issq(dot(QR, QR)) if exact_ok else math.sqrt(dot(QR, QR))
            T = add(Q, mul(QR, (s - PRe) / lenQR))
            put(T, f'T{X}')
    put('inf', 'inf')
    pts, lab = _dedupe(pts, lab)
    exact = None
    if exact_ok:
        exact = {'planar': [p if p == 'inf' else (Rational(p[0].numerator, p[0].denominator), Rational(p[1].numerator, p[1].denominator)) for p in pts],
                 'radicands': ()}
    return (name or f"triangle{A},{B},{C}", pts, lab, exact)


HERONIAN = {
    '13-14-15': ((5, 12), (0, 0), (14, 0)),
    '3-4-5': ((0, 3), (0, 0), (4, 0)),
    '5-5-6': ((3, 4), (0, 0), (6, 0)),
    '5-5-8': ((4, 3), (0, 0), (8, 0)),
    '10-13-13': ((5, 12), (0, 0), (10, 0)),
    '7-15-20': ((Fr(28, 5), Fr(21, 5)), (0, 0), (20, 0)),
    '9-10-17': ((Fr(135, 17), Fr(72, 17)), (0, 0), (17, 0)),
    '4-13-15': ((Fr(12, 5), Fr(16, 5)), (0, 0), (15, 0)),
    '3-25-26': ((Fr(15, 13), Fr(36, 13)), (0, 0), (26, 0)),
    '25-25-14': ((7, 24), (0, 0), (14, 0)),
    '17-17-16': ((8, 15), (0, 0), (16, 0)),
    '6-25-29': ((Fr(126, 29), Fr(120, 29)), (0, 0), (29, 0)),
}


def heronian_universes():
    out = []
    for key, (A, B, C) in HERONIAN.items():
        A = tuple(Fr(c) for c in A); B = tuple(Fr(c) for c in B); C = tuple(Fr(c) for c in C)
        out.append(triangle_universe(A, B, C, incentre=True, name=f"tri-{key}"))
    return out


def special_triangle_universes():
    out = []
    # the b-family of the n=8 record (rational, no incentre points)
    for b in (Fr(1), Fr(1, 2), Fr(3, 2), Fr(2, 3)):
        A = (Fr(0), Fr(3)); B = (b, Fr(0)); C = (3 / b, Fr(0))
        out.append(triangle_universe(A, B, C, incentre=False, name=f"tri-bfam-b={b}"))
    # right isosceles, equilateral, 30-60-90 (floats)
    out.append(triangle_universe((0.0, 1.0), (0.0, 0.0), (1.0, 0.0), incentre=True, name="tri-right-isosceles"))
    out.append(triangle_universe((0.5, math.sqrt(3) / 2), (0.0, 0.0), (1.0, 0.0), incentre=True, name="tri-equilateral"))
    out.append(triangle_universe((0.0, math.sqrt(3)), (0.0, 0.0), (1.0, 0.0), incentre=True, name="tri-30-60-90"))
    out.append(triangle_universe((1.0, 2.0), (0.0, 0.0), (3.0, 0.0), incentre=True, name="tri-generic-float"))
    return out


# ---------------------------------------------------------------- grids
def grid_universe(a, b, half=False, inf=True, name=None):
    pts, lab = [], []
    for i in range(a):
        for j in range(b):
            pts.append((Fr(i), Fr(j))); lab.append(f"({i},{j})")
    if half:
        for i in range(a - 1):
            for j in range(b - 1):
                pts.append((Fr(2 * i + 1, 2), Fr(2 * j + 1, 2))); lab.append(f"({i}.5,{j}.5)")
    if inf:
        pts.append('inf'); lab.append('inf')
    exact = {'planar': [p if p == 'inf' else (Rational(p[0].numerator, p[0].denominator), Rational(p[1].numerator, p[1].denominator)) for p in pts], 'radicands': ()}
    return (name or f"grid{a}x{b}{'+half' if half else ''}{'+inf' if inf else ''}", pts, lab, exact)


# ---------------------------------------------------------------- concentric polygons
def polygon_universe(m, radii, rots, centre=True, inf=True, name=None):
    """Union of regular m-gons of the given radii and rotations (in units of pi/m)."""
    pts, lab = [], []
    for li, (r, rot) in enumerate(zip(radii, rots)):
        for k in range(m):
            th = 2 * math.pi * k / m + rot * math.pi / m
            pts.append((r * math.cos(th), r * math.sin(th))); lab.append(f"P{li}_{k}")
    if centre:
        pts.append((0.0, 0.0)); lab.append('C')
    if inf:
        pts.append('inf'); lab.append('inf')
    pts, lab = _dedupe(pts, lab)
    return (name or f"poly{m}-r{[round(r, 4) for r in radii]}-rot{rots}", pts, lab, None)


# ---------------------------------------------------------------- polyhedra (sphere points)
def polyhedron_universe(name):
    """Returns (name, sphere_points (list of 3-vectors), labels)."""
    phi = (1 + math.sqrt(5)) / 2
    P = []
    if name == 'cube':
        P = list(itertools.product([-1, 1], repeat=3))
    elif name == 'octahedron':
        P = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
    elif name == 'cuboctahedron':
        P = [p for p in itertools.product([-1, 0, 1], repeat=3) if sum(abs(c) for c in p) == 2]
    elif name == 'cube+octahedron':
        P = list(itertools.product([-1, 1], repeat=3)) + [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
    elif name == 'cube+octa+cubocta':   # all 26 directions of {-1,0,1}^3
        P = [p for p in itertools.product([-1, 0, 1], repeat=3) if any(p)]
    elif name == 'icosahedron':
        P = [(0, s1, s2 * phi) for s1 in (-1, 1) for s2 in (-1, 1)] + [(s1, s2 * phi, 0) for s1 in (-1, 1) for s2 in (-1, 1)] + [(s2 * phi, 0, s1) for s1 in (-1, 1) for s2 in (-1, 1)]
    elif name == 'dodecahedron':
        P = list(itertools.product([-1, 1], repeat=3)) + [(0, s1 / phi, s2 * phi) for s1 in (-1, 1) for s2 in (-1, 1)] + [(s1 / phi, s2 * phi, 0) for s1 in (-1, 1) for s2 in (-1, 1)] + [(s2 * phi, 0, s1 / phi) for s1 in (-1, 1) for s2 in (-1, 1)]
    elif name == 'icosidodecahedron':
        P = [(0, 0, s * phi) for s in (-1, 1)] + [(0, s * phi, 0) for s in (-1, 1)] + [(s * phi, 0, 0) for s in (-1, 1)]
        P += [(s1 / 2, s2 * phi / 2, s3 * phi * phi / 2) for s1 in (-1, 1) for s2 in (-1, 1) for s3 in (-1, 1)]
        P += [(s1 * phi / 2, s2 * phi * phi / 2, s3 / 2) for s1 in (-1, 1) for s2 in (-1, 1) for s3 in (-1, 1)]
        P += [(s1 * phi * phi / 2, s2 / 2, s3 * phi / 2) for s1 in (-1, 1) for s2 in (-1, 1) for s3 in (-1, 1)]
    elif name == 'icosa+dodeca':
        _, a, _ = polyhedron_universe('icosahedron'); _, b, _ = polyhedron_universe('dodecahedron')
        P = a + b
    elif name == 'truncated-tetrahedron':
        P = [p for p in itertools.permutations((3, 1, 1))] + [p for p in itertools.permutations((-3, -1, 1))]
        P = list(set(P))
        # even permutations with sign pattern: all permutations of (3,1,1) and (3,-1,-1)?  use standard: (±3,±1,±1) with even number of minus signs
        P = [(a * 3, b, c) for a in (-1, 1) for b in (-1, 1) for c in (-1, 1) if a * b * c == 1]
        P = list(set(itertools.chain.from_iterable(itertools.permutations(p) for p in P)))
    elif name == 'truncated-cube':
        t = math.sqrt(2) - 1
        P = list(set(itertools.chain.from_iterable(itertools.permutations((s1 * t, s2, s3)) for s1 in (-1, 1) for s2 in (-1, 1) for s3 in (-1, 1))))
    elif name == 'truncated-octahedron':
        P = list(set(itertools.chain.from_iterable(itertools.permutations((0, s1, 2 * s2)) for s1 in (-1, 1) for s2 in (-1, 1))))
    elif name == 'rhombicuboctahedron':
        t = 1 + math.sqrt(2)
        P = list(set(itertools.chain.from_iterable(itertools.permutations((s1, s2, s3 * t)) for s1 in (-1, 1) for s2 in (-1, 1) for s3 in (-1, 1))))
    elif name.startswith('prism'):        # prism-m-h
        _, m, h = name.split('-'); m = int(m); h = float(h)
        P = [(math.cos(2 * math.pi * k / m), math.sin(2 * math.pi * k / m), h) for k in range(m)] + [(math.cos(2 * math.pi * k / m), math.sin(2 * math.pi * k / m), -h) for k in range(m)]
    elif name.startswith('antiprism'):    # antiprism-m-h
        _, m, h = name.split('-'); m = int(m); h = float(h)
        P = [(math.cos(2 * math.pi * k / m), math.sin(2 * math.pi * k / m), h) for k in range(m)] + [(math.cos((2 * k + 1) * math.pi / m), math.sin((2 * k + 1) * math.pi / m), -h) for k in range(m)]
    elif name.startswith('bipyramid'):    # bipyramid-m (m-gon on equator + poles)
        m = int(name.split('-')[1])
        P = [(math.cos(2 * math.pi * k / m), math.sin(2 * math.pi * k / m), 0) for k in range(m)] + [(0, 0, 1), (0, 0, -1)]
    else:
        raise ValueError(name)
    P = [tuple(float(c) for c in p) for p in P]
    P = [tuple(c / math.sqrt(sum(cc * cc for cc in p)) for c in p) for p in P]
    return name, P, [f"{name[:3]}{i}" for i in range(len(P))]


# ---------------------------------------------------------------- two circles + radical axis pencils
def two_circle_pencil_universe(d, r2, R_positions, angles, name=None):
    """Circle C1: centre (0,0) radius 1; circle C2: centre (d,0) radius r2.  Radical axis x = x0.
    Through each point R=(x0, y) draw lines with the given angles; add their intersections with C1 and
    C2 (all four points are pairwise concyclic in the pattern {chord of C1 through R} U {chord of C2
    through R}).  Also adds the points R, the centres and inf."""
    x0 = (1 - r2 * r2 + d * d) / (2 * d)
    pts, lab = [(0.0, 0.0), (d, 0.0), 'inf'], ['O1', 'O2', 'inf']
    for ri, y in enumerate(R_positions):
        R = (x0, y); pts.append(R); lab.append(f"R{ri}")
        for ai, th in enumerate(angles):
            ux, uy = math.cos(th), math.sin(th)
            for ci, (cx, cy, rr) in enumerate(((0.0, 0.0, 1.0), (d, 0.0, r2))):
                # |R + t u - c|^2 = rr^2
                fx, fy = R[0] - cx, R[1] - cy
                b = 2 * (fx * ux + fy * uy); c = fx * fx + fy * fy - rr * rr
                disc = b * b - 4 * c
                if disc < 0:
                    continue
                for s in (-1, 1):
                    t = (-b + s * math.sqrt(disc)) / 2
                    pts.append((R[0] + t * ux, R[1] + t * uy)); lab.append(f"C{ci}R{ri}a{ai}{'+' if s > 0 else '-'}")
    pts, lab = _dedupe(pts, lab)
    return (name or f"pencil-d{d}-r{r2}-R{len(R_positions)}-a{len(angles)}", pts, lab, None)


# ---------------------------------------------------------------- Klein-model regular polygon
def klein_polygon_universe(m2, a, name=None):
    """Regular m2-gon on the unit circle with its centre and the chord-direction points at infinity,
    mapped by the hyperbolic (Klein-model) isometry x -> ((x+a)/(1+ax), y sqrt(1-a^2)/(1+ax)), which
    preserves the unit circle and brings the points at infinity to a finite exterior line."""
    s = math.sqrt(1 - a * a)
    def g(x, y, w=1.0):
        # projective: (x, y, w) -> (x + a w, s y, w + a x)
        X, Y, W = x + a * w, s * y, w + a * x
        return (X / W, Y / W)
    pts, lab = [], []
    for k in range(m2):
        th = 2 * math.pi * k / m2
        pts.append(g(math.cos(th), math.sin(th))); lab.append(f"V{k}")
    pts.append(g(0.0, 0.0)); lab.append('O')
    # directions at infinity: chord directions of the m2-gon: angles j*pi/m2 (mod pi), j = 0..m2-1
    for j in range(m2):
        th = math.pi * j / m2 + math.pi / 2   # perpendicular to the axis at angle j*pi/m2 ... include all
        pts.append(g(math.cos(th), math.sin(th), 0.0)); lab.append(f"D{j}")
    pts.append('inf'); lab.append('inf')
    pts, lab = _dedupe(pts, lab)
    return (name or f"klein{m2}-a{a}", pts, lab, None)


if __name__ == "__main__":
    import sphere as SP
    from universe_opt import Universe
    for nm, pts, lab, ex in heronian_universes()[:3] + special_triangle_universes()[:1]:
        U = Universe(nm, SP.plane_to_sphere(pts), lab, ex)
        print(U.info())
    nm, pts, lab, ex = grid_universe(4, 4)
    U = Universe(nm, SP.plane_to_sphere(pts), lab, ex); print(U.info())
    for pn in ('cube', 'cube+octahedron', 'icosahedron', 'cuboctahedron'):
        nm, P, lab = polyhedron_universe(pn)
        U = Universe(nm, P, lab); print(U.info())
