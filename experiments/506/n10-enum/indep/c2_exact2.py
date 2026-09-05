"""Second pass for the (type, line set) cases left open by c2_exact.py: for each such case search the
triples of L for a *clean* certificate: Delta_T has no positive real root other than rho = 1 and the
gcd of the g_k (k in L) computed from T, after removing the trivial factors rho and rho^2 - 1, has no
positive real root other than rho = 1.  (Same mathematics as c2_exact.py, different triple.)
usage: python3 c2_exact2.py <type offset> <out.json>"""
import json, itertools, time, sys
import sympy as sp
import mo

s, rho = sp.symbols('s rho')
m = 16 * s ** 4 - 20 * s ** 2 + 5
SIN = sp.sin(sp.pi / 5)
E = sp.QQ.algebraic_field(SIN)
GEN = E.ext.as_expr()


def red(e):
    return sp.rem(sp.expand(e), m, s)


COS36 = (3 - 4 * s ** 2) / 2


def trig(k):
    c, sn = sp.Integer(1), sp.Integer(0)
    for _ in range(k % 10):
        c, sn = red(c * COS36 - sn * s), red(sn * COS36 + c * s)
    return c, sn


def back_to_s(expr):
    e = sp.expand(expr.subs(GEN, s).subs(SIN, s))
    e = sp.expand(e.subs(sp.sqrt(5), 5 - 8 * s ** 2))
    e = red(e)
    sp.Poly(e, s, rho, domain='QQ')
    return e


def toE(expr):
    return sp.Poly(sp.expand(expr).subs(s, SIN), rho, domain=E)


def pos_roots_not_1(poly_expr_in_s_rho):
    """Isolating intervals of positive real roots other than exactly rho = 1 of the norm polynomial."""
    N = sp.resultant(sp.Poly(poly_expr_in_s_rho, s), sp.Poly(m, s))
    P = sp.Poly(sp.expand(N), rho, domain='QQ')
    if P.is_zero:
        return None
    out = []
    for (a, b), mult in P.intervals():
        if b > 0 and not (a == 1 and b == 1):
            out.append((a, b))
    return out


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


sel = float(sys.argv[1]); outname = sys.argv[2]
open_cases = [int(x) for x in sys.argv[3].split(',')] if len(sys.argv) > 3 else [27, 28, 29, 30]
rig = json.load(open('c2_rigidity.json'))
blocks = [tuple(b) for b in rig['blocks']]
F = [mo.mask(b) for b in blocks]
allblocks = blocks + [tuple(mo.bits(t)) for t in mo.uncovered_triples(10, F)]
A, B = rig['A'], rig['B']
Ls = [[tuple(mo.bits(l)) for l in L] for L in mo.line_search(10, F, mo.Caps(10, 'sg'), min_size=10)]
T = next(T for T in rig['types'] if abs(T['B_offset_deg'] - sel) < 1e-9)
ang = {p: int(round(a / 36)) for p, a in zip(A, T['anglesA'])}
ang.update({p: int(round(b / 36)) for p, b in zip(B, T['anglesB'])})
coord = {}
for p in range(10):
    c, sn = trig(ang[p])
    coord[p] = (c, sn) if p in A else (red(rho * c), red(rho * sn))
circ = {bl: circle([coord[p] for p in bl[:3]]) for bl in allblocks}
report = []
for li in open_cases:
    L = Ls[li]
    t0 = time.time()
    found = None
    tried = 0
    for tri in itertools.combinations(L, 3):
        (a1, b1, c1, d1), (a2, b2, c2, d2), (a3, b3, c3, d3) = (circ[bl] for bl in tri)
        al12, be12, ga12 = red(a2 * b1 - a1 * b2), red(a2 * c1 - a1 * c2), red(a2 * d1 - a1 * d2)
        al13, be13, ga13 = red(a3 * b1 - a1 * b3), red(a3 * c1 - a1 * c3), red(a3 * d1 - a1 * d3)
        Delta = red(al12 * be13 - al13 * be12)
        if Delta == 0:
            continue
        dr = pos_roots_not_1(Delta)
        if dr is None or dr:
            continue          # Delta has positive real roots other than 1: not a clean triple
        tried += 1
        X = red(-ga12 * be13 + ga13 * be12)
        Y = red(-al12 * ga13 + al13 * ga12)
        X2Y2 = red(X ** 2 + Y ** 2)
        G = None
        for bl in L:
            a, b, c, d = circ[bl]
            g = red(a * X2Y2 + b * X * Delta + c * Y * Delta + d * Delta ** 2)
            Pg = toE(g)
            G = Pg if G is None else G.gcd(Pg)
            if G.degree() == 0:
                break
        Gred = G
        for fac in (sp.Poly(rho, rho, domain=E), sp.Poly(rho ** 2 - 1, rho, domain=E)):
            while Gred.degree() > 0:
                qq, rr = Gred.div(fac)
                if rr.is_zero:
                    Gred = qq
                else:
                    break
        gr = [] if Gred.degree() == 0 else pos_roots_not_1(back_to_s(Gred.as_expr()))
        print(f"  type {sel} lineset {li}: triple {tri}: Delta clean; reduced gcd {Gred.as_expr()}; its positive real roots != 1: {gr}", flush=True)
        if not gr:
            found = {'triple': list(tri), 'reduced_gcd': str(Gred.as_expr()), 'delta_positive_roots_not_1': dr}
            break
    entry = {'type_offset': sel, 'lineset': li, 'L': L, 'clean_certificate': found, 'triples_with_clean_delta_tried': tried}
    report.append(entry)
    print(f"type {sel} lineset {li} |L|={len(L)}: {'EXCLUDED by triple ' + str(found['triple']) if found else 'NO CLEAN TRIPLE FOUND'}  [{time.time()-t0:.0f}s]", flush=True)
json.dump(report, open(outname, 'w'), indent=1)
print("SUMMARY:", sum(1 for e in report if e['clean_certificate']), "of", len(report), "open cases excluded")
