"""Exact verification of symmetric-orbit records (Erdos #506).

Given (m, orbits, float params, inversion centre O given by two of its blocks), this module
 1. recognises each parameter exactly (nsimplify against a list of constants, or as an algebraic
    root of the coincidence condition that vanishes there),
 2. builds the points in an algebraic number field K = Q(cos 2pi/m, sin 2pi/m, params...),
 3. computes the inversion centre O exactly (intersection of two blocks; may need a square root,
    handled by a quadratic extension K(sqrt d)),
 4. inverts exactly and counts circles with exact field arithmetic (certificate),
 5. returns exact sympy coordinates.
"""
import math, itertools, sys, json, os
import sympy as sp
from sympy import QQ
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import symcore as S

X = sp.Symbol('x')


# ------------------------------------------------------------------ quadratic extension K(sqrt d)
class QExt:
    """a + b*sqrt(d) with a,b in a sympy domain K (ANP elements); d in K, not a square in K."""
    __slots__ = ('a', 'b', 'F')
    def __init__(self, a, b, F): self.a = a; self.b = b; self.F = F
    def _c(self, o):
        if isinstance(o, QExt): return o
        return QExt(self.F.K.convert(o) if not hasattr(o, 'rep') else o, self.F.K.zero, self.F)
    def __add__(self, o): o = self._c(o); return QExt(self.a + o.a, self.b + o.b, self.F)
    __radd__ = __add__
    def __sub__(self, o): o = self._c(o); return QExt(self.a - o.a, self.b - o.b, self.F)
    def __rsub__(self, o): return self._c(o).__sub__(self)
    def __neg__(self): return QExt(-self.a, -self.b, self.F)
    def __mul__(self, o):
        o = self._c(o); d = self.F.d
        return QExt(self.a * o.a + self.b * o.b * d, self.a * o.b + self.b * o.a, self.F)
    __rmul__ = __mul__
    def inv(self):
        d = self.F.d; nrm = self.a * self.a - self.b * self.b * d
        if not nrm: raise ZeroDivisionError
        return QExt(self.a / nrm, -self.b / nrm, self.F)
    def __truediv__(self, o): return self * self._c(o).inv()
    def __rtruediv__(self, o): return self._c(o) * self.inv()
    def __eq__(self, o):
        o = self._c(o); return self.a == o.a and self.b == o.b
    def __hash__(self): return hash((self.a, self.b))
    def __bool__(self): return bool(self.a) or bool(self.b)
    def to_sympy(self):
        K = self.F.K
        return sp.radsimp(K.to_sympy(self.a) + K.to_sympy(self.b) * sp.sqrt(self.F.d_sym))


class QField:
    def __init__(self, K, d, d_sym): self.K = K; self.d = d; self.d_sym = d_sym
    def elem(self, a, b=None): return QExt(a, self.K.zero if b is None else b, self)


# ------------------------------------------------------------------ exact circle counting
def circumcircle(a, b, c):
    ax, ay = a; bx, by = b; cx, cy = c
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    if not d: return None
    a2 = ax * ax + ay * ay; b2 = bx * bx + by * by; c2 = cx * cx + cy * cy
    ux = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d
    uy = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d
    r2 = (ax - ux) * (ax - ux) + (ay - uy) * (ay - uy)
    return (ux, uy, r2)


def analyze_exact(pts):
    """pts: list of (x,y) with exact field elements. Returns circles dict, collinear-triple list,
    lines dict."""
    n = len(pts)
    assert all(pts[i] != pts[j] for i, j in itertools.combinations(range(n), 2)), "duplicate points"
    circles = defaultdict(set); lines = defaultdict(set); coll = 0
    for i, j, k in itertools.combinations(range(n), 3):
        cc = circumcircle(pts[i], pts[j], pts[k])
        if cc is None:
            coll += 1
            (ax, ay), (bx, by) = pts[i], pts[j]
            # line key: normalised (A,B,C) with A x + B y = C ; scale so that first nonzero of (A,B) is 1
            A = by - ay; B = ax - bx; C = A * ax + B * ay
            if A: B = B / A; C = C / A; A = A / A
            else: C = C / B; B = B / B
            lines[(A, B, C)].update((i, j, k))
        else:
            circles[cc].update((i, j, k))
    return circles, coll, lines


def all_concyclic_or_collinear(pts, circles, lines):
    n = len(pts)
    return any(len(s) == n for s in circles.values()) or any(len(s) == n for s in lines.values())


# ------------------------------------------------------------------ parameter recognition
def consts(m):
    L = [sp.sqrt(2), sp.sqrt(3), sp.sqrt(5), sp.sqrt(6), sp.sqrt(7), sp.pi, sp.GoldenRatio]
    for k in range(1, 2 * m + 1):
        L += [sp.cos(sp.pi * k / (2 * m)), sp.sin(sp.pi * k / (2 * m))]
    return L


def _simple(e):
    """accept only rationals combined with sqrt of small integers, pi, GoldenRatio, cos/sin(pi*q)."""
    if len(str(e)) > 45: return False
    for a in sp.preorder_traversal(e):
        if isinstance(a, sp.Pow):
            if a.exp not in (sp.Rational(1, 2), -sp.Rational(1, 2), -1): return False
            if a.exp.q == 2 and not (a.base.is_Integer and 0 < a.base < 100): return False
        elif isinstance(a, (sp.cos, sp.sin)):
            if not (a.args[0] / sp.pi).is_Rational: return False
        elif isinstance(a, sp.Function):
            return False
    return True


def recognise(val, m, kind):
    q = sp.nsimplify(val, tolerance=1e-10, rational=True)
    if q.is_Rational and q.q <= 1000 and abs(float(q) - val) < 1e-9:
        return q
    cands = []
    for cl in ([sp.sqrt(2), sp.sqrt(3), sp.sqrt(5), sp.GoldenRatio], [sp.sqrt(6), sp.sqrt(7), sp.sqrt(15), sp.sqrt(21)], consts(m)):
        for tol in (1e-12, 1e-10):
            try:
                e = sp.nsimplify(val, cl, tolerance=tol, rational=False)
            except Exception:
                continue
            if e.free_symbols or e.atoms(sp.Float) or not _simple(e): continue
            if abs(float(e.evalf(30)) - val) < 5e-9:
                cands.append(e)
    if kind in ('a', 'f'):
        q = sp.nsimplify(val / math.pi, tolerance=1e-9, rational=True)
        if q.is_Rational and q.q <= 48 * m and abs(float(q * sp.pi) - val) < 5e-9:
            cands.append(q * sp.pi)
    if cands:
        cands.sort(key=lambda e: len(str(e)))
        return cands[0]
    if kind == 'r':
        # algebraic number of small degree: guess the minimal polynomial from the double value
        # (the exact verification later certifies or refutes the guess)
        import mpmath
        mpmath.mp.dps = 15
        for deg in range(2, 7):
            try:
                P = mpmath.findpoly(val, deg, maxcoeff=3000, maxsteps=200)
            except Exception:
                P = None
            if P is None: continue
            poly = sp.Poly(list(P), X)
            for fac, _ in sp.factor_list(poly.as_expr(), X)[1]:
                fp = sp.Poly(fac, X)
                if fp.degree() < 1: continue
                for r in fp.real_roots():
                    if abs(float(r.evalf(30)) - val) < 1e-8:
                        rr = r
                        if isinstance(r, sp.CRootOf):
                            try:
                                for k in sp.roots(fp):
                                    if abs(complex(k.evalf(30)) - val) < 1e-8 and k.is_real is not False and _simple(k):
                                        rr = sp.radsimp(k); break
                            except Exception:
                                pass
                        return rr
    return None


def exact_points_sym(m, orbits, params_exact):
    pts = []; k = 0
    for ob in orbits:
        t = ob[0]
        if t == 'A':
            r = params_exact[k]; k += 1
            for j in range(m):
                a = 2 * sp.pi * j / m
                pts.append((r * sp.cos(a), r * sp.sin(a)))
        elif t == 'B':
            r = params_exact[k]; k += 1
            for j in range(m):
                a = (2 * j + 1) * sp.pi / m
                pts.append((r * sp.cos(a), r * sp.sin(a)))
        elif t == 'G':
            r, th = params_exact[k], params_exact[k + 1]; k += 2
            for j in range(m):
                for s in (1, -1):
                    a = s * th + 2 * sp.pi * j / m
                    pts.append((r * sp.cos(a), r * sp.sin(a)))
        elif t == 'R':
            r, th = params_exact[k], params_exact[k + 1]; k += 2
            for j in range(m):
                a = th + 2 * sp.pi * j / m
                pts.append((r * sp.cos(a), r * sp.sin(a)))
        elif t == 'C':
            pts.append((sp.Integer(0), sp.Integer(0)))
    return pts


def solve_param_algebraic(fam, params_exact_known, i, val):
    """Parameter i is a root of a coincidence condition; find it exactly.  Returns a sympy
    algebraic number (radical expression or CRootOf) or None."""
    t = sp.Symbol('t', real=True)
    pe = list(params_exact_known); pe[i] = t
    pts = exact_points_sym(fam.m, fam.orbits, pe)
    if not fam.has_c: pts.append((sp.Integer(0), sp.Integer(0)))
    # find vanishing conditions numerically
    pb = np.array([[float(x) if not isinstance(x, sp.Symbol) else val for x in pe]])
    vals = fam.conds_batch(pb)[0]
    order = np.argsort(np.abs(vals))
    for c in order[:6]:
        if abs(vals[c]) > 1e-8: break
        if c < len(fam.quads):
            idx = fam.quads[c]
            rows = [[pts[q][0], pts[q][1], pts[q][0] ** 2 + pts[q][1] ** 2, 1] for q in idx]
            expr = sp.Matrix(rows).det()
        else:
            idx = fam.trips[c - len(fam.quads)]
            (ax, ay), (bx, by), (cx, cy) = (pts[q] for q in idx)
            expr = (bx - ax) * (cy - ay) - (cx - ax) * (by - ay)
        if fam.kinds[i] in ('a', 'f'):
            # angle parameter: substitute u = tan(t/2)?  simpler: use s = sin t, c = cos t with c^2+s^2=1
            cs, sn = sp.symbols('c s', real=True)
            e2 = sp.expand(sp.expand_trig(expr).subs({sp.cos(t): cs, sp.sin(t): sn}))
            e2 = sp.expand(e2.subs(cs, sp.sqrt(1 - sn ** 2)))
            root_expr = _roots_near(e2, sn, math.sin(val))
            if root_expr is None: continue
            th = sp.asin(root_expr)
            if abs(float(th.evalf(30)) - val) > 1e-9: th = sp.pi - th
            if abs(float(th.evalf(30)) - val) > 1e-9: continue
            return th
        root_expr = _roots_near(sp.expand(expr), t, val)
        if root_expr is not None:
            return root_expr
    return None


def _roots_near(expr, t, val):
    """exact real root of expr(t)=0 nearest to val (expr has algebraic constant coefficients)."""
    expr = sp.nsimplify(expr)
    # bring to a polynomial in t over an algebraic field, then to a polynomial over Q via minimal polynomial
    try:
        P = sp.Poly(expr, t)
    except sp.PolynomialError:
        expr = sp.together(expr)
        P = sp.Poly(sp.numer(expr), t)
    if P.degree() < 1: return None
    coeffs = P.all_coeffs()
    gens = set()
    for c in coeffs:
        gens |= {a for a in c.atoms(sp.Function, sp.Pow) if not a.free_symbols}
    gens = [g for g in gens if not g.is_Rational]
    try:
        K = QQ.algebraic_field(*gens) if gens else QQ
        Pk = sp.Poly(expr, t, domain=K)
    except Exception:
        K = None; Pk = None
    if K is not None and K is not QQ:
        # norm down to QQ: resultant with the minimal polynomial of the primitive element
        z = sp.Symbol('z')
        alpha = K.ext
        mp = sp.Poly(sp.minimal_polynomial(alpha, z), z)
        # express coefficients as polynomials in z
        cz = [sp.Poly(K.to_sympy(cc), z) if False else sp.Poly(sp.Poly(cc.rep.to_list() if hasattr(cc, 'rep') else [cc], z).as_expr(), z) for cc in Pk.all_coeffs()]
        f = sum(c.as_expr() * t ** (P.degree() - i) for i, c in enumerate(cz))
        N = sp.Poly(sp.resultant(f, mp.as_expr(), z), t)
        Nq = N
    else:
        Nq = sp.Poly(expr, t, domain=QQ)
    facs = sp.factor_list(Nq.as_expr(), t)[1]
    best = None
    for fac, mult in facs:
        fp = sp.Poly(fac, t)
        if fp.degree() < 1: continue
        for r in fp.real_roots():
            fr = float(r.evalf(30))
            if abs(fr - val) < 1e-7 and (best is None or abs(fr - val) < best[0]):
                best = (abs(fr - val), r)
    if best is None: return None
    r = best[1]
    # prefer radical form for low degree
    if isinstance(r, sp.CRootOf):
        try:
            rad = sp.roots(sp.Poly(r.expr, r.poly.gen))
            for k in rad:
                if k.is_real is not False and abs(float(k.evalf(30)) - val) < 1e-9:
                    return sp.radsimp(k)
        except Exception:
            pass
    return r


# ------------------------------------------------------------------ main verification
def build_field(exprs):
    gens = set()
    for e in exprs:
        for a in sp.preorder_traversal(e):
            if isinstance(a, (sp.Pow, sp.Function, sp.CRootOf)) and not a.free_symbols and not a.is_Rational:
                if isinstance(a, sp.Pow) and not a.exp.is_Rational: continue
                gens.add(a)
    gens = [g for g in gens if g.is_real is not False]
    # replace nested things by their atoms: sympy handles cos(pi/7) etc. as generators
    return QQ.algebraic_field(*gens) if gens else QQ


def to_field(K, e):
    return K.from_sympy(sp.nsimplify(e))


def verify_entry(entry, log=print):
    from search import Family
    fam = Family(entry['m'], [tuple(o) for o in entry['orbits']])
    params = entry['params']
    pe = [None] * fam.p
    for i, v in enumerate(params):
        pe[i] = recognise(v, fam.m, fam.kinds[i])
    for i, v in enumerate(params):
        if pe[i] is None:
            known = [x if x is not None else sp.nsimplify(params[j]) for j, x in enumerate(pe)]
            pe[i] = solve_param_algebraic(fam, known, i, v)
            if pe[i] is None:
                log(f"  could not recognise parameter {i} = {v}"); return None
    log(f"  exact params: {pe}")
    pts_sym = exact_points_sym(fam.m, fam.orbits, pe)
    # numeric sanity: match floats
    P0 = fam.points(params)
    for (xs, ys), (xf, yf) in zip(pts_sym, P0):
        assert abs(float(xs.evalf(30)) - xf) < 1e-8 and abs(float(ys.evalf(30)) - yf) < 1e-8, "exact params do not match floats"
    # field
    gens_src = list(pe) + [sp.cos(2 * sp.pi / fam.m), sp.sin(2 * sp.pi / fam.m), sp.cos(sp.pi / fam.m), sp.sin(sp.pi / fam.m)]
    flat = []
    for e in gens_src:
        flat += [e]
    K = build_field([sp.expand_trig(sp.nsimplify(sp.expand(x))) for x in flat] + [sp.nsimplify(sp.expand_trig(sp.expand(c))) for p in pts_sym for c in p])
    log(f"  field: {K} (degree {K.mod.degree() if K is not QQ else 1})")
    pts = [(to_field(K, sp.expand_trig(x)), to_field(K, sp.expand_trig(y))) for x, y in pts_sym]
    circles, coll, lines = analyze_exact(pts)
    has_inf = fam.has_inf
    nb = len(circles) + len(lines) + (0 if has_inf else 0)
    if has_inf:
        # lines through 2 points count too; blocks = circles + all lines through >=2 finite points
        nl2 = 0
        seen = set()
        for i, j in itertools.combinations(range(len(pts)), 2):
            (ax, ay), (bx, by) = pts[i], pts[j]
            A = by - ay; B = ax - bx; C = A * ax + B * ay
            if A: B = B / A; C = C / A; A = A / A
            else: C = C / B; B = B / B
            seen.add((A, B, C))
        nb = len(circles) + len(seen)
    log(f"  Moebius blocks (exact): circles {len(circles)}, lines(>=3) {len(lines)}, total {nb}; float said {entry['nblocks']}")
    # inversion centre
    Ob = entry.get('O_blocks')
    result = {'exact_params': [str(x) for x in pe], 'nblocks_exact': nb}
    if entry.get('O_is_inf'):
        E = pts; F = None
        assert not has_inf
    else:
        assert Ob and len(Ob) >= 2
        # take two blocks; get exact geometry from members
        def geom(mem):
            fin = [q for q in mem if q >= 0]
            if len(fin) >= 3:
                cc = circumcircle(pts[fin[0]], pts[fin[1]], pts[fin[2]])
                if cc is not None: return ('C', cc)
            (ax, ay), (bx, by) = pts[fin[0]], pts[fin[1]]
            return ('L', (by - ay, ax - bx, (by - ay) * ax + (ax - bx) * ay))
        g1, g2 = geom(Ob[0]), geom(Ob[1])
        Of = entry['O']
        # solve symbolically in sympy using field->sympy conversions, then back to K or K(sqrt d)
        xs, ys = sp.symbols('xs ys', real=True)
        def eqn(g):
            if g[0] == 'L':
                A, B, C = (K.to_sympy(v) for v in g[1]); return A * xs + B * ys - C
            ux, uy, r2 = (K.to_sympy(v) for v in g[1]); return (xs - ux) ** 2 + (ys - uy) ** 2 - r2
        e1, e2 = eqn(g1), eqn(g2)
        if g1[0] == 'C' and g2[0] == 'C': e2 = sp.expand(e1 - e2)   # radical axis (linear)
        # linear one: solve for one variable
        lin = e2 if g2[0] == 'C' or g1[0] == 'C' else e1
        quad = e1 if lin is e2 else e2
        if g1[0] == 'L' and g2[0] == 'L':
            sol = sp.solve([e1, e2], [xs, ys], dict=True)
            cand = [(s[xs], s[ys]) for s in sol]
        else:
            lin = sp.expand(lin)
            if lin.coeff(ys) != 0:
                ysol = sp.solve(lin, ys)[0]
                q = sp.expand(quad.subs(ys, ysol))
                a2, a1, a0 = [sp.nsimplify(q.coeff(xs, k)) for k in (2, 1, 0)]
                disc = sp.expand(a1 ** 2 - 4 * a2 * a0)
                cand = []
                for sg in (1, -1):
                    xv = (-a1 + sg * sp.sqrt(disc)) / (2 * a2)
                    cand.append((xv, ysol.subs(xs, xv)))
            else:
                xsol = sp.solve(lin, xs)[0]
                q = sp.expand(quad.subs(xs, xsol))
                a2, a1, a0 = [sp.nsimplify(q.coeff(ys, k)) for k in (2, 1, 0)]
                disc = sp.expand(a1 ** 2 - 4 * a2 * a0)
                cand = []
                for sg in (1, -1):
                    yv = (-a1 + sg * sp.sqrt(disc)) / (2 * a2)
                    cand.append((xsol.subs(ys, yv), yv))
        best = None
        for cx, cy in cand:
            dd = abs(complex(cx.evalf(30)) - Of[0]) + abs(complex(cy.evalf(30)) - Of[1])
            if best is None or dd < best[0]: best = (dd, cx, cy)
        assert best[0] < 1e-6, f"could not match O exactly ({best[0]})"
        Ox, Oy = best[1], best[2]
        # try to express in K; else K(sqrt d)
        F = None
        try:
            OxK = to_field(K, Ox); OyK = to_field(K, Oy)
            E_O = (OxK, OyK); FK = K
            lift = lambda v: v
        except Exception:
            # need sqrt(disc): disc in K
            dK = to_field(K, disc)
            # is disc a square in K?  test via sympy factor over the field
            F = QField(K, dK, sp.nsimplify(disc))
            def lift(v):
                # v is a sympy expression a + b sqrt(disc): split
                v = sp.expand(v)
                b = sp.nsimplify(v.coeff(sp.sqrt(disc)))
                a = sp.nsimplify(sp.expand(v - b * sp.sqrt(disc)))
                return F.elem(to_field(K, a), to_field(K, b))
            # make sure sqrt(disc) truly is not in K by checking x^2 - disc irreducible over K
            fl = sp.factor_list(X ** 2 - sp.nsimplify(disc), X, extension=K.ext if K is not QQ else None)
            assert all(sp.Poly(f, X).degree() == 2 for f, _ in fl[1]), "disc is a square in K; use K"
            E_O = (lift(sp.radsimp(Ox)), lift(sp.radsimp(Oy)))
            FK = F
        # invert
        def conv(v): return v if F is None else F.elem(v)
        Ox_, Oy_ = E_O
        E = []
        for (x, y) in pts:
            dx = conv(x) - Ox_; dy = conv(y) - Oy_
            nrm = dx * dx + dy * dy
            assert nrm, "O coincides with a point of P"
            E.append((Ox_ + dx / nrm, Oy_ + dy / nrm))
        if has_inf:
            E.append((Ox_, Oy_))
    circles2, coll2, lines2 = analyze_exact(E)
    degenerate = all_concyclic_or_collinear(E, circles2, lines2)
    count = len(circles2)
    sizes = sorted((len(s) for s in circles2.values()), reverse=True)
    result.update({'circles_exact': count, 'collinear_triples': coll2, 'degenerate': degenerate,
                   'circle_sizes': sizes, 'line_sizes': sorted((len(s) for s in lines2.values()), reverse=True)})
    # exact coordinates as sympy strings
    def tos(v):
        if F is None: return sp.radsimp(K.to_sympy(v))
        return v.to_sympy()
    coords = [(tos(x), tos(y)) for x, y in E]
    result['coords_sympy'] = [(str(x), str(y)) for x, y in coords]
    result['coords_float'] = [[float(x.evalf(20)), float(y.evalf(20))] for x, y in coords]
    log(f"  EXACT: circles={count} collinear triples={coll2} degenerate={degenerate} circle sizes={sizes} line sizes={result['line_sizes']}")
    return result


if __name__ == '__main__':
    path = sys.argv[1]
    d = json.load(open(path))
    for n in sorted(d, key=int):
        for e in d[n][:1]:
            print(f"n={n} float best={e['best']} formula={e['formula']} {e['family']} params={e['params']}")
            try:
                r = verify_entry(e)
            except Exception as ex:
                import traceback; traceback.print_exc()
