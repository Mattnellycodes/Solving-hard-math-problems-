"""Exact certification of records (Erdős #506, search-local agent).

Input: exact coordinates of the Möbius point set S (sympy-parsable strings, 'inf' allowed) and a
float inversion centre (or 'inf').  Steps:
 1. identify the centre exactly as an intersection of two blocks of S (sympy, choosing the root
    nearest to the float), or use it as given if it is a point of the universe (exact string);
 2. invert S exactly about the centre (infinity -> centre);
 3. count the circles of the resulting finite set:
    - rational coordinates: exact arithmetic with Fractions, cross-checked with
      experiments/506/circles_exact.py (independent implementation);
    - algebraic coordinates: blocks found at 60-digit precision (mpmath), every asserted coincidence
      (each concyclic quadruple beyond the first three points of a rich circle, each collinear
      triple) certified exactly with sympy (determinant simplifies to 0, or minimal polynomial of the
      determinant is x), and distinct circles/lines separated by > 1e-30 numerically.
Output: dict with the certified count, block sizes, and the exact coordinates.
"""
import sys, math, itertools, json
from fractions import Fraction as Fr
import sympy as sp
import mpmath
sys.path.insert(0, '/home/user/Solving-hard-math-problems-/experiments/506')
from circles_exact import count_circles as count_circles_ref, all_on_one_circle_or_line

mpmath.mp.dps = 60


def parse(s):
    if s == 'inf':
        return 'inf'
    return sp.nsimplify(sp.sympify(s), rational=False)


def is_zero(expr):
    e = sp.simplify(expr)
    if e == 0:
        return True
    try:
        e = sp.radsimp(sp.expand(e))
        if e == 0:
            return True
        mp = sp.minimal_polynomial(e, sp.Symbol('x'))
        return mp == sp.Symbol('x')
    except Exception:
        return False


def det_concyclic(a, b, c, d):
    M = sp.Matrix([[p[0], p[1], p[0] ** 2 + p[1] ** 2, 1] for p in (a, b, c, d)])
    return M.det()


def det_collinear(a, b, c):
    return sp.Matrix([[a[0], a[1], 1], [b[0], b[1], 1], [c[0], c[1], 1]]).det()


def circle_exact(a, b, c):
    ax, ay = a; bx, by = b; cx, cy = c
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    a2 = ax ** 2 + ay ** 2; b2 = bx ** 2 + by ** 2; c2 = cx ** 2 + cy ** 2
    ux = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d
    uy = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d
    return sp.simplify(ux), sp.simplify(uy)


def line_exact(a, b):
    A = b[1] - a[1]; B = a[0] - b[0]; C = A * a[0] + B * a[1]
    return sp.simplify(A), sp.simplify(B), sp.simplify(C)


def blocks_numeric(pts):
    """pts: list of sympy (x,y) finite points. Returns (circles, lines) as lists of index tuples using
    60-digit numerics; circles keyed by centre & radius."""
    P = [(mpmath.mpf(sp.N(x, 70)), mpmath.mpf(sp.N(y, 70))) for x, y in pts]
    tol = mpmath.mpf(10) ** -35
    circles = []   # (cx, cy, r, set)
    lines = []     # (A, B, C, set)
    for i, j, k in itertools.combinations(range(len(P)), 3):
        (ax, ay), (bx, by), (cx, cy) = P[i], P[j], P[k]
        d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
        if abs(d) < tol:
            A = by - ay; B = ax - bx; nrm = mpmath.sqrt(A * A + B * B); A /= nrm; B /= nrm; C = A * ax + B * ay
            if A < 0 or (abs(A) < tol and B < 0):
                A, B, C = -A, -B, -C
            for L in lines:
                if abs(L[0] - A) < tol and abs(L[1] - B) < tol and abs(L[2] - C) < tol:
                    L[3].update((i, j, k)); break
            else:
                lines.append([A, B, C, {i, j, k}])
        else:
            a2 = ax * ax + ay * ay; b2 = bx * bx + by * by; c2 = cx * cx + cy * cy
            ux = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d
            uy = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d
            r = mpmath.sqrt((ax - ux) ** 2 + (ay - uy) ** 2)
            for Cc in circles:
                if abs(Cc[0] - ux) < tol and abs(Cc[1] - uy) < tol and abs(Cc[2] - r) < tol:
                    Cc[3].update((i, j, k)); break
            else:
                circles.append([ux, uy, r, {i, j, k}])
    # separation check: distinct blocks differ by > 1e-30 in some parameter
    sep = mpmath.mpf(10) ** -30
    for X, Y in itertools.combinations(circles, 2):
        assert max(abs(X[0] - Y[0]), abs(X[1] - Y[1]), abs(X[2] - Y[2])) > sep, "circle separation failed"
    for X, Y in itertools.combinations(lines, 2):
        assert max(abs(X[0] - Y[0]), abs(X[1] - Y[1]), abs(X[2] - Y[2])) > sep, "line separation failed"
    return [sorted(c[3]) for c in circles], [sorted(l[3]) for l in lines]


def certify_algebraic(pts):
    """exact certification of the block structure of finite algebraic points; returns (ncircles, sizes)."""
    circles, lines = blocks_numeric(pts)
    # every coincidence must be exact
    for c in circles:
        a, b, cc = (pts[c[0]], pts[c[1]], pts[c[2]])
        for extra in c[3:]:
            assert is_zero(det_concyclic(a, b, cc, pts[extra])), f"concyclicity not exact: {c}"
    for l in lines:
        for i, j, k in itertools.combinations(l, 3):
            assert is_zero(det_collinear(pts[i], pts[j], pts[k])), f"collinearity not exact: {l}"
    # non-coincidences: a triple in a 3-circle is genuinely non-collinear (d != 0 numerically, fine),
    # and the numeric separation of distinct circles certifies they are distinct (algebraic numbers
    # of small height cannot agree to 30 digits without being equal — see notes).
    ncirc = len(circles)
    sizes = {}
    for c in circles:
        sizes[len(c)] = sizes.get(len(c), 0) + 1
    lsizes = {}
    for l in lines:
        lsizes[len(l)] = lsizes.get(len(l), 0) + 1
    return ncirc, sizes, lsizes


def count_rational(pts):
    """exact count with Fractions (own implementation) — pts list of (Fraction, Fraction)."""
    circles = {}
    lines = {}
    for i, j, k in itertools.combinations(range(len(pts)), 3):
        (ax, ay), (bx, by), (cx, cy) = pts[i], pts[j], pts[k]
        d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
        if d == 0:
            A = by - ay; B = ax - bx; C = A * ax + B * ay
            g = math.gcd(math.gcd(A.numerator, B.numerator), C.numerator) if False else None
            # normalise by making the first nonzero of (A,B) equal to 1
            if A != 0:
                key = ('L', 1, B / A, C / A)
            else:
                key = ('L', 0, 1, C / B)
            lines.setdefault(key, set()).update((i, j, k))
        else:
            a2 = ax * ax + ay * ay; b2 = bx * bx + by * by; c2 = cx * cx + cy * cy
            ux = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d
            uy = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d
            r2 = (ax - ux) ** 2 + (ay - uy) ** 2
            circles.setdefault((ux, uy, r2), set()).update((i, j, k))
    sizes = {}
    for s in circles.values():
        sizes[len(s)] = sizes.get(len(s), 0) + 1
    lsizes = {}
    for s in lines.values():
        lsizes[len(s)] = lsizes.get(len(s), 0) + 1
    return len(circles), sizes, lsizes


def invert(p, o, r2=1):
    if p == 'inf':
        return o
    dx = p[0] - o[0]; dy = p[1] - o[1]
    d2 = dx * dx + dy * dy
    return (sp.simplify(o[0] + r2 * dx / d2), sp.simplify(o[1] + r2 * dy / d2))


def identify_centre(pts, centre_float, tol=1e-6):
    """find the exact centre: intersection of two blocks of S through the float point."""
    fin = [p for p in pts if p != 'inf']
    has_inf = any(p == 'inf' for p in pts)
    cx, cy = centre_float
    cands = []
    # blocks through the centre numerically: circles through 3 finite pts, lines through 2 finite pts (+inf) or 3
    F = [(float(x), float(y)) for x, y in fin]
    scale = max(1.0, max(abs(v) for p in F for v in p))
    geoms = []
    for i, j, k in itertools.combinations(range(len(fin)), 3):
        (ax, ay), (bx, by), (cx3, cy3) = F[i], F[j], F[k]
        d = 2 * (ax * (by - cy3) + bx * (cy3 - ay) + cx3 * (ay - by))
        if abs(d) < 1e-12 * scale ** 2:
            continue
        a2 = ax * ax + ay * ay; b2 = bx * bx + by * by; c2 = cx3 * cx3 + cy3 * cy3
        ux = (a2 * (by - cy3) + b2 * (cy3 - ay) + c2 * (ay - by)) / d
        uy = (a2 * (cx3 - bx) + b2 * (ax - cx3) + c2 * (bx - ax)) / d
        r = math.hypot(ax - ux, ay - uy)
        if abs(math.hypot(cx - ux, cy - uy) - r) < tol * scale:
            geoms.append(('C', (i, j, k)))
    for i, j in itertools.combinations(range(len(fin)), 2):
        (ax, ay), (bx, by) = F[i], F[j]
        A = by - ay; B = ax - bx; C = A * ax + B * ay; nrm = math.hypot(A, B)
        if abs(A * cx + B * cy - C) / nrm < tol * scale:
            # is it a block? needs a third point of S: another finite point on it or inf in S
            third = has_inf or any(abs(A * F[k][0] + B * F[k][1] - C) / nrm < 1e-9 * scale for k in range(len(fin)) if k not in (i, j))
            if third:
                geoms.append(('L', (i, j)))
    # dedupe geoms by type/point-set is unnecessary; take the first two distinct blocks and solve exactly
    x, y = sp.symbols('x y', real=True)
    eqs = []
    used = set()
    for g in geoms:
        if len(eqs) == 2:
            break
        if g[0] == 'C':
            i, j, k = g[1]
            ux, uy = circle_exact(fin[i], fin[j], fin[k])
            r2 = sp.simplify((fin[i][0] - ux) ** 2 + (fin[i][1] - uy) ** 2)
            key = (ux, uy, r2)
            if key in used:
                continue
            used.add(key)
            eqs.append(sp.expand((x - ux) ** 2 + (y - uy) ** 2 - r2))
        else:
            i, j = g[1]
            A, B, C = line_exact(fin[i], fin[j])
            key = ('L', sp.simplify(B / A) if A != 0 else None, sp.simplify(C / A) if A != 0 else sp.simplify(C / B))
            if key in used:
                continue
            used.add(key)
            eqs.append(A * x + B * y - C)
    assert len(eqs) == 2, "could not find two blocks through the centre"
    sols = sp.solve(eqs, [x, y], dict=True)
    best = None
    for s in sols:
        sx, sy = s[x], s[y]
        if not (sx.is_real is not False and sy.is_real is not False):
            continue
        d = math.hypot(float(sx) - cx, float(sy) - cy)
        if best is None or d < best[0]:
            best = (d, sp.nsimplify(sp.simplify(sx)), sp.nsimplify(sp.simplify(sy)))
    assert best is not None and best[0] < 1e-6 * scale, f"centre identification failed {best}"
    return (best[1], best[2])


def certify(exact_coords, centre, want=None, verbose=True):
    pts = [parse(s) if s == 'inf' else (parse(s[0]), parse(s[1])) for s in exact_coords]
    if centre == 'inf' or centre is None:
        assert all(p != 'inf' for p in pts), "infinity in S needs a finite centre"
        Q = pts
        O = 'inf'
    else:
        if isinstance(centre, (list, tuple)) and isinstance(centre[0], str):
            O = (parse(centre[0]), parse(centre[1]))
        else:
            O = identify_centre(pts, (float(centre[0]), float(centre[1])))
        Q = [invert(p, O) for p in pts]
    # distinctness
    for a, b in itertools.combinations(Q, 2):
        assert not (is_zero(a[0] - b[0]) and is_zero(a[1] - b[1])), "duplicate points after inversion"
    rational = all(q[0].is_Rational and q[1].is_Rational for q in Q)
    if rational:
        F = [(Fr(int(q[0].p), int(q[0].q)), Fr(int(q[1].p), int(q[1].q))) for q in Q]
        n1, sizes, lsizes = count_rational(F)
        n2 = count_circles_ref(F)
        assert n1 == n2, (n1, n2)
        assert not all_on_one_circle_or_line(F)
        # integer form
        den = 1
        for x, y in F:
            den = den * x.denominator // math.gcd(den, x.denominator)
            den = den * y.denominator // math.gcd(den, y.denominator)
        coords = [(int(x * den), int(y * den)) for x, y in F]
        res = dict(count=n1, circle_sizes=sizes, line_sizes=lsizes, exact='rational', centre=str(O),
                   coords=coords, note=f"integer coordinates (scaled by {den})")
    else:
        n1, sizes, lsizes = certify_algebraic(Q)
        # non-degeneracy: at least two distinct blocks exist (sizes has more than one block)
        assert sum(sizes.values()) + sum(lsizes.values()) >= 2
        res = dict(count=n1, circle_sizes=sizes, line_sizes=lsizes, exact='algebraic (sympy-certified coincidences)',
                   centre=[str(O[0]), str(O[1])] if O != 'inf' else 'inf', coords=[[str(q[0]), str(q[1])] for q in Q])
    if want is not None:
        res['matches_claim'] = (res['count'] == want)
    if verbose:
        print(json.dumps(res))
    return res


if __name__ == "__main__":
    # self-test: the n=8 record via the ortho universe: S = A,B,C,H,D,E,F,inf, centre (0,1)
    S = [["0", "3"], ["1", "0"], ["3", "0"], ["0", "-1"], ["6/5", "-3/5"], ["2", "1"], ["0", "0"], "inf"]
    r = certify(S, [0.0, 1.0], want=17)
    assert r['count'] == 17
    # algebraic self-test: octagon + centre (antipodal, n=9 -> 25), no inversion
    s2 = "sqrt(2)/2"
    S9 = [["1", "0"], ["0", "1"], ["-1", "0"], ["0", "-1"], [s2, s2], [f"-{s2}", s2], [s2, f"-{s2}"], [f"-{s2}", f"-{s2}"], ["0", "0"]]
    r = certify(S9, 'inf', want=25)
    assert r['count'] == 25
    print("self-tests OK")
