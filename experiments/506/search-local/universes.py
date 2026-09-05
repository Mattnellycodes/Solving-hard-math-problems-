"""Universe builders for the local search (Erdős #506, search-local agent).

Each builder returns engine.Universe.  Exact coordinates (Fractions for rational universes, sympy
expressions for algebraic ones) are stored in U.exact for later certification.
"""
from fractions import Fraction as Fr
import math, itertools
import numpy as np
import sympy as sp
from engine import Universe


def _dedupe_exact(pts, labels):
    seen = {}
    out, lab = [], []
    for p, l in zip(pts, labels):
        if p not in seen:
            seen[p] = True; out.append(p); lab.append(l)
    return out, lab


def _from_exact(name, pts, labels, has_inf=True, extra_O=None, q=1e-7):
    P = np.array([[float(x), float(y)] for x, y in pts], dtype=np.float64)
    return Universe(name, P, has_inf=has_inf, exact=pts, labels=labels, q=q, extra_O=extra_O)


# ---------------------------------------------------------------- lattices
def grid(k, half=False, name=None, has_inf=True):
    pts, lab = [], []
    step = Fr(1, 2) if half else Fr(1)
    m = 2 * k - 1 if half else k
    for i in range(m):
        for j in range(m):
            pts.append((step * i, step * j)); lab.append(f"({step*i},{step*j})")
    return _from_exact(name or f"grid{k}x{k}{'h' if half else ''}", pts, lab, has_inf)


def rect(a, b, name=None):
    pts = [(Fr(i), Fr(j)) for i in range(a) for j in range(b)]
    lab = [f"({i},{j})" for i in range(a) for j in range(b)]
    return _from_exact(name or f"grid{a}x{b}", pts, lab, True)


def invert_pt(p, c, r2):
    dx, dy = p[0] - c[0], p[1] - c[1]
    d2 = dx * dx + dy * dy
    if d2 == 0:
        return None
    return (c[0] + r2 * dx / d2, c[1] + r2 * dy / d2)


def grid_inversions(k, centres, r2s, name=None, keep_grid=True):
    """k x k grid together with its images under inversions in the given lattice centres/radii."""
    base = [(Fr(i), Fr(j)) for i in range(k) for j in range(k)]
    pts, lab = [], []
    if keep_grid:
        for p in base:
            pts.append(p); lab.append(f"g({p[0]},{p[1]})")
    for c in centres:
        c = (Fr(c[0]), Fr(c[1]))
        for r2 in r2s:
            r2 = Fr(r2)
            for p in base:
                q = invert_pt(p, c, r2)
                if q is not None:
                    pts.append(q); lab.append(f"inv[{c[0]},{c[1]};{r2}]({p[0]},{p[1]})")
    pts, lab = _dedupe_exact(pts, lab)
    return _from_exact(name or f"grid{k}inv", pts, lab, True)


# ---------------------------------------------------------------- polygons (algebraic)
def polygons(m, radii, rots, centre=True, has_inf=True, name=None, extra_pts=()):
    """regular m-gons at the given radii (sympy expressions) with rotation offsets rots (in units
    of pi/m, integers), plus centre and infinity."""
    pts, lab, ex = [], [], []
    for r, rot in zip(radii, rots):
        rr = sp.nsimplify(r)
        for k in range(m):
            ang = sp.pi * (2 * k + rot) / m
            ex.append((rr * sp.cos(ang), rr * sp.sin(ang)))
            pts.append((float(rr) * math.cos(math.pi * (2 * k + rot) / m), float(rr) * math.sin(math.pi * (2 * k + rot) / m)))
            lab.append(f"P{m}[r={r},rot={rot}]#{k}")
    for (x, y, l) in extra_pts:
        ex.append((sp.nsimplify(x), sp.nsimplify(y))); pts.append((float(x), float(y))); lab.append(l)
    if centre:
        ex.append((sp.Integer(0), sp.Integer(0))); pts.append((0.0, 0.0)); lab.append('C')
    # dedupe by floats
    out, lo, le = [], [], []
    for p, l, e in zip(pts, lab, ex):
        if all(abs(p[0] - q[0]) > 1e-9 or abs(p[1] - q[1]) > 1e-9 for q in out):
            out.append(p); lo.append(l); le.append(e)
    P = np.array(out, dtype=np.float64)
    return Universe(name or f"poly{m}x{len(radii)}", P, has_inf=has_inf, exact=le, labels=lo)


def tri_lattice(R, name=None, has_inf=True):
    """triangular lattice points a + b*omega with |p|^2 <= R^2 (exact in Q(sqrt3))."""
    pts, lab, ex = [], [], []
    s3 = sp.sqrt(3)
    for a in range(-2 * R, 2 * R + 1):
        for b in range(-2 * R, 2 * R + 1):
            x = a + b / 2.0; y = b * math.sqrt(3) / 2
            if x * x + y * y <= R * R + 1e-9:
                pts.append((x, y)); lab.append(f"T({a},{b})")
                ex.append((sp.Rational(2 * a + b, 2), sp.Rational(b, 2) * s3))
    P = np.array(pts, dtype=np.float64)
    return Universe(name or f"trilat{R}", P, has_inf=has_inf, exact=ex, labels=lab)


# ---------------------------------------------------------------- orthocentric / triangle closure
def _circumcentre(a, b, c):
    ax, ay = a; bx, by = b; cx, cy = c
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    if d == 0:
        return None
    a2 = ax * ax + ay * ay; b2 = bx * bx + by * by; c2 = cx * cx + cy * cy
    return ((a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d, (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d)


def _orthocentre(a, b, c):
    o = _circumcentre(a, b, c)
    if o is None:
        return None
    return (a[0] + b[0] + c[0] - 2 * o[0], a[1] + b[1] + c[1] - 2 * o[1])


def _foot(p, a, b):
    """foot of the perpendicular from p to line ab"""
    dx, dy = b[0] - a[0], b[1] - a[1]
    t = ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / (dx * dx + dy * dy)
    return (a[0] + t * dx, a[1] + t * dy)


def _reflect_line(p, a, b):
    f = _foot(p, a, b)
    return (2 * f[0] - p[0], 2 * f[1] - p[1])


def _mid(a, b):
    return ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)


def ortho_closure(seed, levels=1, ops=('H', 'O', 'G', 'M', 'F', 'RH', 'RM', 'N'), cap=400, name=None, has_inf=True):
    """Start from rational points `seed`; for every triangle of the current set add orthocentre (H),
    circumcentre (O), centroid (G), side midpoints (M), feet of altitudes (F), reflections of H in
    sides (RH) and in side midpoints (RM), nine-point centre (N).  Repeat `levels` times (cap points)."""
    pts = [(Fr(x), Fr(y)) for x, y in seed]
    lab = [f"s{i}" for i in range(len(pts))]
    cur = list(pts)
    for lev in range(levels):
        new = {}
        base = list(cur)
        for a, b, c in itertools.combinations(base, 3):
            o = _circumcentre(a, b, c)
            if o is None:
                continue
            h = _orthocentre(a, b, c)
            cand = []
            if 'H' in ops: cand.append((h, 'H'))
            if 'O' in ops: cand.append((o, 'O'))
            if 'G' in ops: cand.append((((a[0] + b[0] + c[0]) / 3, (a[1] + b[1] + c[1]) / 3), 'G'))
            if 'M' in ops: cand += [(_mid(a, b), 'M'), (_mid(b, c), 'M'), (_mid(a, c), 'M')]
            if 'F' in ops: cand += [(_foot(a, b, c), 'F'), (_foot(b, a, c), 'F'), (_foot(c, a, b), 'F')]
            if 'RH' in ops: cand += [(_reflect_line(h, b, c), 'RH'), (_reflect_line(h, a, c), 'RH'), (_reflect_line(h, a, b), 'RH')]
            if 'RM' in ops:
                for x, y in ((a, b), (b, c), (a, c)):
                    m = _mid(x, y); cand.append(((2 * m[0] - h[0], 2 * m[1] - h[1]), 'RM'))
            if 'N' in ops: cand.append((_mid(o, h), 'N'))
            for p, l in cand:
                if p not in new:
                    new[p] = l
        for p, l in new.items():
            if p not in pts:
                pts.append(p); lab.append(f"{l}{lev}")
                if len(pts) >= cap:
                    break
        cur = list(pts)
        if len(pts) >= cap:
            break
    # dedupe exact
    pts, lab = _dedupe_exact(pts, lab)
    return _from_exact(name or f"ortho{len(pts)}", pts, lab, has_inf)


# ---------------------------------------------------------------- sphere universes
def _stereo(V, pole=None):
    """stereographic projection of unit-sphere points V (list of sympy 3-vectors) from `pole`
    (a unit vector, sympy) — pole itself maps to infinity.  Returns exact (x,y) pairs and floats."""
    # rotate so that pole -> (0,0,1) then (x,y,z) -> (x,y)/(1-z)
    pole = sp.Matrix(pole)
    # build orthonormal frame: e3 = pole
    e3 = pole
    # pick a vector not parallel
    tmp = sp.Matrix([1, 0, 0]) if abs(float(e3[0])) < 0.9 else sp.Matrix([0, 1, 0])
    e1 = tmp - tmp.dot(e3) * e3
    e1 = e1 / sp.sqrt(e1.dot(e1))
    e2 = e3.cross(e1)
    ex, fl, inf_idx = [], [], []
    for i, v in enumerate(V):
        v = sp.Matrix(v)
        x = sp.simplify(v.dot(e1)); y = sp.simplify(v.dot(e2)); z = sp.simplify(v.dot(e3))
        if sp.simplify(1 - z) == 0:
            inf_idx.append(i); continue
        X = sp.simplify(x / (1 - z)); Y = sp.simplify(y / (1 - z))
        ex.append((X, Y)); fl.append((float(X), float(Y)))
    return ex, fl, inf_idx


def _unit(v):
    v = sp.Matrix(v)
    return v / sp.sqrt(v.dot(v))


def icosa_family(which=('I', 'D'), pole='vertex', name=None):
    """union of icosahedron (I), dodecahedron (D), icosidodecahedron (E) vertices on the unit sphere,
    stereographically projected from a vertex (-> infinity) or from a generic point."""
    phi = (1 + sp.sqrt(5)) / 2
    V = []
    lab = []
    if 'I' in which:
        for s1 in (1, -1):
            for s2 in (1, -1):
                for perm in range(3):
                    v = [0, s1, s2 * phi]
                    v = v[perm:] + v[:perm]
                    V.append(_unit(v)); lab.append('I')
    if 'D' in which:
        for s in itertools.product((1, -1), repeat=3):
            V.append(_unit([s[0], s[1], s[2]])); lab.append('D')
        for s1 in (1, -1):
            for s2 in (1, -1):
                for perm in range(3):
                    v = [0, s1 / phi, s2 * phi]
                    v = v[perm:] + v[:perm]
                    V.append(_unit(v)); lab.append('D')
    if 'E' in which:
        for s in (1, -1):
            for perm in range(3):
                v = [s, 0, 0]
                v = v[perm:] + v[:perm]
                V.append(_unit(v)); lab.append('E')
        for s in itertools.product((1, -1), repeat=3):
            for perm in range(3):
                v = [s[0] / 2, s[1] * phi / 2, s[2] * phi * phi / 2]
                v = v[perm:] + v[:perm]
                V.append(_unit(v)); lab.append('E')
    # dedupe by float
    out, lo = [], []
    for v, l in zip(V, lab):
        fv = np.array([float(c) for c in v])
        if all(np.linalg.norm(fv - np.array([float(c) for c in w])) > 1e-9 for w in out):
            out.append(v); lo.append(l)
    if pole == 'vertex':
        p = out[0]
    elif pole == 'face':
        p = _unit([1, 1, 1])                       # icosahedron face centre / dodecahedron vertex direction
    else:
        p = _unit([3, 5, 7])                       # generic
    ex, fl, inf_idx = _stereo(out, p)
    labels = [lo[i] for i in range(len(out)) if i not in inf_idx]
    has_inf = len(inf_idx) > 0
    P = np.array(fl, dtype=np.float64)
    return Universe(name or f"icosa{''.join(which)}-{pole}", P, has_inf=has_inf, exact=ex, labels=labels, q=1e-7)


def octa_shells(norms=(1, 2, 3), pole='vertex', name=None):
    """integer points of the given norms scaled to the unit sphere (octahedral symmetry), projected."""
    V, lab = [], []
    for Nn in norms:
        r = int(math.isqrt(3 * Nn)) + 1
        for a in range(-r, r + 1):
            for b in range(-r, r + 1):
                for c in range(-r, r + 1):
                    if a * a + b * b + c * c == Nn:
                        V.append(sp.Matrix([a, b, c]) / sp.sqrt(Nn)); lab.append(f"({a},{b},{c})/√{Nn}")
    out, lo = [], []
    for v, l in zip(V, lab):
        fv = np.array([float(c) for c in v])
        if all(np.linalg.norm(fv - np.array([float(c) for c in w])) > 1e-9 for w in out):
            out.append(v); lo.append(l)
    if pole == 'vertex':
        p = sp.Matrix([0, 0, 1])
    elif pole == 'edge':
        p = _unit([1, 1, 0])
    elif pole == 'face':
        p = _unit([1, 1, 1])
    else:
        p = _unit([3, 5, 7])
    ex, fl, inf_idx = _stereo(out, p)
    labels = [lo[i] for i in range(len(out)) if i not in inf_idx]
    P = np.array(fl, dtype=np.float64)
    return Universe(name or f"octa{'-'.join(map(str,norms))}-{pole}", P, has_inf=len(inf_idx) > 0, exact=ex, labels=labels)


# ---------------------------------------------------------------- coincidence closure
def multiplicity_points(U, min_mult=3, max_pairs=3_000_000, tol=1e-7):
    """points (not in U) where >= min_mult blocks of U meet; returns (K,2) float array."""
    from engine import _intersections
    geoms = U.geom
    nb = len(geoms)
    scale = max(1.0, float(np.abs(U.P).max()))
    cands = {}
    sums = {}
    pairs = 0
    for i in range(nb):
        gi = geoms[i]
        for j in range(i + 1, nb):
            pairs += 1
            if pairs > max_pairs:
                break
            for (x, y) in _intersections(gi, geoms[j]):
                key = (round(x / (1e-6 * scale)), round(y / (1e-6 * scale)))
                s = cands.get(key)
                if s is None:
                    cands[key] = {i, j}; sums[key] = [x, y, 1]
                else:
                    s.add(i); s.add(j); sm = sums[key]; sm[0] += x; sm[1] += y; sm[2] += 1
        if pairs > max_pairs:
            break
    # merge clusters whose keys are neighbours (rounding-boundary effects), union-find
    parent = {k: k for k in cands}

    def find(k):
        while parent[k] != k:
            parent[k] = parent[parent[k]]; k = parent[k]
        return k
    for k in list(cands.keys()):
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                if dx == 0 and dy == 0:
                    continue
                k2 = (k[0] + dx, k[1] + dy)
                if k2 in cands:
                    a, b = find(k), find(k2)
                    if a != b:
                        parent[a] = b
    merged = {}
    msums = {}
    for k, s in cands.items():
        r = find(k)
        if r not in merged:
            merged[r] = set(); msums[r] = [0.0, 0.0, 0]
        merged[r] |= s
        sm = sums[k]; ms = msums[r]; ms[0] += sm[0]; ms[1] += sm[1]; ms[2] += sm[2]
    out = []
    for key, s in merged.items():
        if len(s) >= min_mult:
            ms = msums[key]
            x, y = ms[0] / ms[2], ms[1] / ms[2]        # mean of the actual intersection coordinates
            if any(abs(x - U.P[k, 0]) < 1e-5 * scale and abs(y - U.P[k, 1]) < 1e-5 * scale for k in range(U.N)):
                continue
            out.append((x, y, len(s)))
    return out


# ---------------------------------------------------------------- conics and polygon-diagonal closures
def ellipse(m, a, b, extras=('C',), name=None, has_inf=True, theta0=0):
    """points (a cos t, b sin t) at t = theta0 + 2*pi*k/m (4 points concyclic iff sum of t == 0 mod 2pi);
    extras: 'C' centre, 'F' foci, 'V' vertices at infinity? (no).  Exact in cyclotomic radicals."""
    pts, lab, ex = [], [], []
    A, B = sp.nsimplify(a), sp.nsimplify(b)
    for k in range(m):
        t = sp.pi * (2 * k) / m + sp.nsimplify(theta0)
        ex.append((A * sp.cos(t), B * sp.sin(t)))
        pts.append((float(A) * math.cos(float(t)), float(B) * math.sin(float(t)))); lab.append(f"E#{k}")
    if 'C' in extras:
        ex.append((sp.Integer(0), sp.Integer(0))); pts.append((0.0, 0.0)); lab.append('C')
    if 'F' in extras:
        c = sp.sqrt(A * A - B * B)
        for s in (1, -1):
            ex.append((s * c, sp.Integer(0))); pts.append((s * float(c), 0.0)); lab.append('F')
    out, lo, le = [], [], []
    for p, l, e in zip(pts, lab, ex):
        if all(abs(p[0] - q[0]) > 1e-9 or abs(p[1] - q[1]) > 1e-9 for q in out):
            out.append(p); lo.append(l); le.append(e)
    return Universe(name or f"ellipse{m}-{a}-{b}", np.array(out), has_inf=has_inf, exact=le, labels=lo)


def hyperbola_exact(ts, name=None, extras=('C',), has_inf=True):
    """rectangular hyperbola xy = 1 at parameters ts (4 points concyclic iff product == 1)."""
    pts, lab = [], []
    for t in ts:
        t = Fr(t)
        pts.append((t, 1 / t)); lab.append(f"H(t={t})")
    if 'C' in extras:
        pts.append((Fr(0), Fr(0))); lab.append('C')
    pts, lab = _dedupe_exact(pts, lab)
    return _from_exact(name or f"hyperbola{len(ts)}", pts, lab, has_inf)


def polygon_diag_closure(m, min_mult=3, max_pts=200, centre=True, name=None, has_inf=True, radius2=None):
    """regular m-gon vertices + centre + the intersection points of >= min_mult diagonals (the
    classical concurrency points of the regular polygon), as selectable points; optional second
    concentric m-gon of radius radius2 (same orientation) whose diagonals are included too."""
    V = [(math.cos(2 * math.pi * k / m), math.sin(2 * math.pi * k / m)) for k in range(m)]
    labs = [f"V#{k}" for k in range(m)]
    exact = [(sp.cos(2 * sp.pi * k / m), sp.sin(2 * sp.pi * k / m)) for k in range(m)]
    if radius2 is not None:
        r2 = sp.nsimplify(radius2)
        V += [(float(r2) * math.cos(2 * math.pi * k / m), float(r2) * math.sin(2 * math.pi * k / m)) for k in range(m)]
        labs += [f"W#{k}" for k in range(m)]
        exact += [(r2 * sp.cos(2 * sp.pi * k / m), r2 * sp.sin(2 * sp.pi * k / m)) for k in range(m)]
    # all lines through 2 vertices
    lines = []
    for i, j in itertools.combinations(range(len(V)), 2):
        (ax, ay), (bx, by) = V[i], V[j]
        A = by - ay; B = ax - bx; C = A * ax + B * ay
        nrm = math.hypot(A, B); A /= nrm; B /= nrm; C /= nrm
        if A < -1e-12 or (abs(A) <= 1e-12 and B < 0):
            A, B, C = -A, -B, -C
        key = (round(A / 1e-9), round(B / 1e-9), round(C / 1e-9))
        for L in lines:
            if L[0] == key:
                L[2].add(i); L[2].add(j); break
        else:
            lines.append([key, (A, B, C), {i, j}, (i, j)])
    # pairwise intersections
    cands = {}
    for (k1, g1, s1, p1), (k2, g2, s2, p2) in itertools.combinations(lines, 2):
        A1, B1, C1 = g1; A2, B2, C2 = g2
        det = A1 * B2 - A2 * B1
        if abs(det) < 1e-12:
            continue
        x = (C1 * B2 - C2 * B1) / det; y = (A1 * C2 - A2 * C1) / det
        key = (round(x / 1e-7), round(y / 1e-7))
        c = cands.get(key)
        if c is None:
            cands[key] = [x, y, {k1, k2}, (p1, p2)]
        else:
            c[2].add(k1); c[2].add(k2)
    pts = list(V); lab = list(labs); ex = list(exact)
    if centre:
        pts.append((0.0, 0.0)); lab.append('C'); ex.append((sp.Integer(0), sp.Integer(0)))
    extra = [(c[0], c[1], len(c[2]), c[3]) for c in cands.values() if len(c[2]) >= min_mult]
    extra = [e for e in extra if all(abs(e[0] - p[0]) > 1e-7 or abs(e[1] - p[1]) > 1e-7 for p in pts)]
    extra.sort(key=lambda e: -e[2])
    extra = extra[:max_pts]
    x, y = sp.symbols('x y')
    for (px, py, mult, (p1, p2)) in extra:
        pts.append((px, py)); lab.append(f"D{mult}")
        # exact: intersection of lines through exact vertex pairs p1, p2 (lazy: store the defining pairs)
        ex.append(('diag', p1, p2))
    U = Universe(name or f"pdiag{m}{'x2' if radius2 else ''}", np.array(pts), has_inf=has_inf, exact=None, labels=lab)
    U.exact_lazy = ex
    U.exact_vertices = exact
    return U


def resolve_lazy_exact(U, idx):
    """exact sympy coordinates of point idx for polygon_diag_closure universes."""
    e = U.exact_lazy[idx]
    if isinstance(e, tuple) and len(e) == 3 and e[0] == 'diag':
        (i1, j1), (i2, j2) = e[1], e[2]
        Vx = U.exact_vertices
        x, y = sp.symbols('x y')
        def lineq(P, Q):
            A = Q[1] - P[1]; B = P[0] - Q[0]; C = A * P[0] + B * P[1]
            return A * x + B * y - C
        sol = sp.solve([lineq(Vx[i1], Vx[j1]), lineq(Vx[i2], Vx[j2])], [x, y], dict=True)[0]
        return (sp.nsimplify(sp.simplify(sol[x])), sp.nsimplify(sp.simplify(sol[y])))
    return e
