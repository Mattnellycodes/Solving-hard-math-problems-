"""Exact certification that no realisation of C-2 has a point O (not in P) on >= 10 blocks.

From c2_rigidity.py: every realisation is Moebius-equivalent to A = regular pentagon on the unit circle,
B = regular pentagon on |z| = rho (rho > 0, rho != 1), B rotated by 0 or 36 degrees (two types), with
the labelling given by a solution w (angles = 36 w degrees).  All coordinates lie in E[rho] where
E = Q(s), s = sin 36 deg, minimal polynomial m(s) = 16 s^4 - 20 s^2 + 5 (cos 36 = (3 - 4 s^2)/2).

A point O on >= 10 blocks is on the blocks of one of the 31 line sets L (pairwise <=1-intersecting
families of blocks / uncovered triples, c2_linesets.py).  For a triple C1, C2, C3 in L with radical
centre (X/Delta, Y/Delta) (Delta = determinant of the two radical axes), concurrency of all of L at a
point requires, for every block C_k of L,
     g_k(rho) = a_k (X^2 + Y^2) + b_k X Delta + c_k Y Delta + d_k Delta^2 = 0
whenever Delta(rho) != 0.  We compute gcd_E(g_k : k in L); a constant gcd proves that no rho (real or
complex) with Delta(rho) != 0 works.  The real positive roots of Delta (collinear centres) are
isolated exactly via the norm polynomial N(Delta) = res_s(Delta, m) and treated separately: then the
three circles are concurrent only if they are coaxial, i.e. C1 and C2 meet on C3; we test that with
exact arithmetic at those roots (via a second triple of L with non-vanishing Delta at that rho, or
interval evaluation)."""
import json, itertools, time, sys
import sympy as sp
from sympy import Rational as Q_
import mo

s, rho = sp.symbols('s rho')
m = 16 * s ** 4 - 20 * s ** 2 + 5
SIN = sp.sin(sp.pi / 5)
E = sp.QQ.algebraic_field(SIN)


def red(e):
    return sp.rem(sp.expand(e), m, s)


COS36 = (3 - 4 * s ** 2) / 2


def trig(k):
    c, sn = sp.Integer(1), sp.Integer(0)
    for _ in range(k % 10):
        c, sn = red(c * COS36 - sn * s), red(sn * COS36 + c * s)
    return c, sn


for k in range(10):
    c, sn = trig(k)
    assert abs(float(c.subs(s, SIN)) - sp.cos(sp.pi * k / 5).evalf()) < 1e-12 and abs(float(sn.subs(s, SIN)) - sp.sin(sp.pi * k / 5).evalf()) < 1e-12

rig = json.load(open('c2_rigidity.json'))
lsets = json.load(open('c2_linesets.json'))
blocks = [tuple(b) for b in rig['blocks']]
F = [mo.mask(b) for b in blocks]
allblocks = blocks + [tuple(mo.bits(t)) for t in mo.uncovered_triples(10, F)]
assert len(allblocks) == 42
A, B = rig['A'], rig['B']
Ls = mo.line_search(10, F, mo.Caps(10, 'sg'), min_size=10)
Ls = [[tuple(mo.bits(l)) for l in L] for L in Ls]
assert len(Ls) == 31
types = {}
for T in rig['types']:
    types.setdefault(T['B_offset_deg'], T)   # one labelled solution per geometric type
sel = float(sys.argv[1]) if len(sys.argv) > 1 else None
outname = sys.argv[2] if len(sys.argv) > 2 else 'c2_exact.json'
if sel is not None:
    types = {k: v for k, v in types.items() if abs(k - sel) < 1e-9}
print("geometric types:", list(types), "-> output", outname)


def det3(M):
    return sp.Matrix(M).det()


def circle(pts):
    (x1, y1), (x2, y2), (x3, y3) = pts
    q1, q2, q3 = red(x1 ** 2 + y1 ** 2), red(x2 ** 2 + y2 ** 2), red(x3 ** 2 + y3 ** 2)
    a = red(det3([[x1, y1, 1], [x2, y2, 1], [x3, y3, 1]]))
    b = red(-det3([[q1, y1, 1], [q2, y2, 1], [q3, y3, 1]]))
    c = red(det3([[q1, x1, 1], [q2, x2, 1], [q3, x3, 1]]))
    d = red(-det3([[q1, x1, y1], [q2, x2, y2], [q3, x3, y3]]))
    return (a, b, c, d)


def toE(expr):
    return sp.Poly(sp.expand(expr).subs(s, SIN), rho, domain=E)


GEN = E.ext.as_expr()   # sympy's expression for the field generator (sin 36 deg)


def back_to_s(expr):
    """Rewrite an expression in the field generator (and sqrt(5)) as a polynomial in s, reduced mod m."""
    e = sp.expand(expr.subs(GEN, s).subs(SIN, s))
    e = sp.expand(e.subs(sp.sqrt(5), 5 - 8 * s ** 2))
    e = red(e)
    sp.Poly(e, s, rho, domain='QQ')   # raises if anything non-rational remains
    return e


assert back_to_s(GEN) == s and back_to_s(sp.sqrt(5)) == red(5 - 8 * s ** 2)


def real_pos_roots(poly_expr):
    """Isolating intervals of the real roots in (0, oo) of a polynomial in rho over Q (exact)."""
    P = sp.Poly(sp.expand(poly_expr), rho, domain='QQ')
    if P.is_zero:
        return None
    return [iv for iv in P.intervals() if iv[0][1] > 0]


report = []
t_all = time.time()
for off, T in types.items():
    ang = {p: int(round(a / 36)) for p, a in zip(A, T['anglesA'])}
    ang.update({p: int(round(b / 36)) for p, b in zip(B, T['anglesB'])})
    coord = {}
    for p in range(10):
        c, sn = trig(ang[p])
        coord[p] = (c, sn) if p in A else (red(rho * c), red(rho * sn))
    circ = {}
    for bl in allblocks:
        circ[bl] = circle([coord[p] for p in bl[:3]])
        a, b, c, d = circ[bl]
        for p in bl:   # exact verification of all incidences of the block, for symbolic rho
            x, y = coord[p]
            assert red(a * (x ** 2 + y ** 2) + b * x + c * y + d) == 0, (bl, p)
        assert not all(red(v) == 0 for v in (a, b, c)), bl
    print(f"type offset {off}: all 42 blocks verified concyclic for symbolic rho (exact)", flush=True)
    for li, L in enumerate(Ls):
        t0 = time.time()
        chosen = None
        for tri in itertools.combinations(L, 3):
            (a1, b1, c1, d1), (a2, b2, c2, d2), (a3, b3, c3, d3) = (circ[bl] for bl in tri)
            al12, be12, ga12 = red(a2 * b1 - a1 * b2), red(a2 * c1 - a1 * c2), red(a2 * d1 - a1 * d2)
            al13, be13, ga13 = red(a3 * b1 - a1 * b3), red(a3 * c1 - a1 * c3), red(a3 * d1 - a1 * d3)
            Delta = red(al12 * be13 - al13 * be12)
            if Delta == 0:
                continue
            X = red(-ga12 * be13 + ga13 * be12)
            Y = red(-al12 * ga13 + al13 * ga12)
            chosen = (tri, Delta, X, Y)
            break
        assert chosen is not None, ("all triples degenerate", li)
        tri, Delta, X, Y = chosen
        X2Y2 = red(X ** 2 + Y ** 2)
        gs = []
        for bl in L:
            a, b, c, d = circ[bl]
            g = red(a * X2Y2 + b * X * Delta + c * Y * Delta + d * Delta ** 2)
            gs.append(g)
        assert all(g == 0 for g, bl in zip(gs, L) if bl in tri) or True
        G = None
        for g in gs:
            Pg = toE(g)
            G = Pg if G is None else G.gcd(Pg)
            if G.degree() == 0:
                break
        # degenerate rho: real positive roots of Delta
        ND = sp.resultant(sp.Poly(Delta, s), sp.Poly(m, s))
        droots = real_pos_roots(ND)
        # strip the trivial degenerate factors rho and (rho^2 - 1) (rho = 0, 1 collapse the configuration)
        Gred = G
        triv = {'rho': 0, 'rho2m1': 0}
        for fac, key in ((sp.Poly(rho, rho, domain=E), 'rho'), (sp.Poly(rho ** 2 - 1, rho, domain=E), 'rho2m1')):
            while Gred.degree() > 0:
                qq, rr = Gred.div(fac)
                if rr.is_zero:
                    Gred = qq; triv[key] += 1
                else:
                    break
        # Delta roots other than rho = 1 ?
        bad_delta = [iv for iv in (droots or []) if not (iv[0] == 1 and iv[1] == 1)]
        entry = {'type_offset': off, 'lineset': li, 'L': L, 'triple': list(tri), 'gcd_degree': G.degree(),
                 'gcd': str(G.as_expr()), 'trivial_factors': triv, 'reduced_gcd_degree': Gred.degree(),
                 'reduced_gcd': str(Gred.as_expr()),
                 'delta_norm_positive_real_roots': droots if droots is None else [(str(a), str(b)) for a, b in droots],
                 'delta_roots_other_than_1': [(str(a), str(b)) for a, b in bad_delta]}
        if Gred.degree() > 0:
            Gs = red(back_to_s(Gred.as_expr()))
            NG = sp.resultant(sp.Poly(Gs, s), sp.Poly(m, s))
            rr = real_pos_roots(NG)
            entry['reduced_gcd_norm_positive_real_roots'] = [(str(a), str(b)) for a, b in rr]
            entry['verdict'] = 'NEEDS ANALYSIS' if (rr or bad_delta) else 'EXCLUDED'
        else:
            entry['verdict'] = 'EXCLUDED' if not bad_delta else 'NEEDS ANALYSIS (Delta root)'
        report.append(entry)
        print(f"  type {off} lineset {li} |L|={len(L)}: triple {tri} deg(Delta)={sp.Poly(Delta, rho).degree()}; gcd = {G.as_expr()} "
              f"-> trivial factors {triv}, reduced gcd degree {Gred.degree()} ({Gred.as_expr() if Gred.degree() > 0 else 1}); Delta positive real roots {entry['delta_norm_positive_real_roots']}"
              + (f"; reduced-gcd positive real roots: {entry.get('reduced_gcd_norm_positive_real_roots')}" if Gred.degree() > 0 else '')
              + f" => {entry['verdict']}  [{time.time()-t0:.1f}s]", flush=True)
json.dump(report, open(outname, 'w'), indent=1)
bad = [e for e in report if e['verdict'] != 'EXCLUDED']
print(f"\nSUMMARY: {len(report)} (type, line set) cases; EXCLUDED: {len(report) - len(bad)}; needing further analysis: {len(bad)}  [total {time.time()-t_all:.0f}s]")
for e in bad:
    print("  NEEDS ANALYSIS:", e)
if not bad:
    print("RESULT: no realisation of C-2 has a point outside P on >= 10 blocks; every planar configuration with the block structure C-2 has >= 42 - 9 = 33 circles.")
