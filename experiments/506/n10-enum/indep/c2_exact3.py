"""Third pass for the (type, line set) cases left open by c2_exact.py (line sets 27-30).
Mathematics: for a triple T of blocks of L with Delta_T != 0 identically, concurrency of all of L at a
point (for the radius ratio rho) forces  rho in Z_T := {positive real roots != 1 of Delta_T} u
{positive real roots != 1 of G_T}, where G_T is the gcd of the radical-centre polynomials g_k with the
trivial factors rho, rho^2-1 removed (c2_exact.py).  Hence rho lies in the intersection of Z_T over
several triples, which is contained in the set of positive real roots != 1 of
H = gcd_E(Delta_T1 * G_T1, Delta_T2 * G_T2, ...).  Real roots are now determined in the TRUE embedding
s = sin 36 deg (not via the norm alone): squarefree part over E, isolating intervals of the rational
norm polynomial, and exact sign changes of the E-polynomial at the rational endpoints (interval
arithmetic with an enclosure of sin 36 deg, exact zero test in E).
usage: python3 c2_exact3.py <type offset> <out.json> [linesets] [max_triples]"""
import json, itertools, time, sys
from fractions import Fraction
import sympy as sp
from mpmath import iv, mp
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


def sign_in_E(expr_in_s):
    """Exact sign of an element of E given as a polynomial in s with rational coefficients."""
    e = red(sp.expand(expr_in_s))
    if e == 0:
        return 0
    P = sp.Poly(e, s, domain='QQ')
    coeffs = [Fraction(int(c.p), int(c.q)) for c in P.all_coeffs()]   # highest degree first
    for dps in (60, 150, 400, 1000):
        iv.dps = dps
        sv = iv.sin(iv.pi / 5)
        val = iv.mpf(0)
        for c in coeffs:
            val = val * sv + iv.mpf(c.numerator) / iv.mpf(c.denominator)
        if val.a > 0:
            return 1
        if val.b < 0:
            return -1
    raise RuntimeError("sign undetermined for " + str(e))


def true_pos_roots_not_1(PE):
    """PE: Poly in rho over E.  Return isolating intervals (a, b) of its positive real roots != 1 in
    the true embedding s = sin 36 deg (exact).  Method: squarefree part over E; isolating intervals of
    the rational norm polynomial (each contains exactly one distinct root of the norm, hence at most
    one root of the squarefree part); an open interval is refined until neither endpoint is a root
    and it does not contain 0 or 1 in its interior; then the root belongs to PE iff the exact signs
    of the squarefree part at the two endpoints differ."""
    if PE.degree() <= 0:
        return []
    Psf = PE.quo(PE.gcd(PE.diff(rho)))
    Ps = back_to_s(Psf.as_expr())
    N = sp.resultant(sp.Poly(Ps, s), sp.Poly(m, s))
    NP = sp.Poly(sp.expand(N), rho, domain='QQ')
    assert not NP.is_zero
    NP = NP.sqf_part()          # same real roots; needed for isolation/refinement
    out = []
    for (a, b), mult in NP.intervals():
        if b <= 0:
            continue
        if a == b:
            if a > 0 and a != 1 and sign_in_E(Ps.subs(rho, a)) == 0:
                out.append((a, b))
            continue
        # open interval (a, b) isolating exactly one root of the norm, which is not rational
        for _ in range(300):
            if a == b or b <= 0:
                break
            if a > 0 and not (a < 1 < b) and sign_in_E(Ps.subs(rho, a)) != 0 and sign_in_E(Ps.subs(rho, b)) != 0:
                break
            a, b = NP.refine_root(a, b, steps=1)
        else:
            raise RuntimeError("could not refine interval")
        if a == b:
            if a > 0 and a != 1 and sign_in_E(Ps.subs(rho, a)) == 0:
                out.append((a, b))
            continue
        if b <= 0:
            continue            # negative root
        # now a > 0, 1 is not strictly inside (a, b), and neither endpoint is a root of Psf:
        # the (unique, irrational) root of the norm in (a, b) is a root of Psf iff the signs differ
        if sign_in_E(Ps.subs(rho, a)) * sign_in_E(Ps.subs(rho, b)) < 0:
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
max_triples = int(sys.argv[4]) if len(sys.argv) > 4 else 6
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
# self-test of the embedding-aware root finder: rho^2 - 2 rho cos36 + 1 has no real root; (rho - 2 s)(rho + 1) has root 2 sin36
assert true_pos_roots_not_1(toE(rho ** 2 - 2 * COS36 * rho + 1)) == []
tst = true_pos_roots_not_1(toE((rho - 2 * s) * (rho + 1)))
assert len(tst) == 1 and float(tst[0][0]) <= 2 * 0.5877852522924731 <= float(tst[0][1]), tst
tst2 = true_pos_roots_not_1(toE((rho - (3 - 4 * s ** 2) / 2 + 1 / sp.Integer(2)) ** 2 * (rho - 1) ** 2 * (rho + 1)))   # double root at cos36 - 1/2 = 0.309 (=cos72)
assert len(tst2) == 1 and float(tst2[0][0]) <= 0.30901699437 <= float(tst2[0][1]), tst2
print("root-finder self-test passed", flush=True)
report = []
for li in open_cases:
    L = Ls[li]
    t0 = time.time()
    H = None; used = []
    always_circle = {bl: not true_pos_roots_not_1(toE(circ[bl][0])) and toE(circ[bl][0]).degree() >= 0 and circ[bl][0] != 0 for bl in L}
    print(f"  type {sel} lineset {li}: blocks that are circles for all rho>0 (rho != 1): {[bl for bl in L if always_circle[bl]]}", flush=True)
    order = sorted(itertools.combinations(L, 3), key=lambda t: -sum(always_circle[bl] for bl in t))
    for tri in order:
        if len(used) >= max_triples:
            break
        (a1, b1, c1, d1), (a2, b2, c2, d2), (a3, b3, c3, d3) = (circ[bl] for bl in tri)
        al12, be12, ga12 = red(a2 * b1 - a1 * b2), red(a2 * c1 - a1 * c2), red(a2 * d1 - a1 * d2)
        al13, be13, ga13 = red(a3 * b1 - a1 * b3), red(a3 * c1 - a1 * c3), red(a3 * d1 - a1 * d3)
        Delta = red(al12 * be13 - al13 * be12)
        if Delta == 0:
            continue
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
        DE = toE(Delta)
        prod = DE * Gred
        H = prod if H is None else H.gcd(prod)
        droots = true_pos_roots_not_1(DE)
        groots = true_pos_roots_not_1(Gred)
        used.append({'triple': list(tri), 'delta_true_pos_roots_not_1': [(str(a), str(b)) for a, b in droots],
                     'reduced_gcd': str(Gred.as_expr()), 'reduced_gcd_true_pos_roots_not_1': [(str(a), str(b)) for a, b in groots]})
        print(f"  type {sel} lineset {li} triple {tri}: Delta true roots {droots}; reduced gcd {Gred.as_expr()} true roots {groots}; H degree now {H.degree()}", flush=True)
        if not droots and not groots:
            break          # clean triple: done
        if H.degree() == 0 or not true_pos_roots_not_1(H):
            break          # accumulated gcd has no admissible root: done
    Hroots = true_pos_roots_not_1(H) if H.degree() > 0 else []
    verdict = 'EXCLUDED' if not Hroots else 'NEEDS ANALYSIS'
    report.append({'type_offset': sel, 'lineset': li, 'L': L, 'triples': used, 'H': str(H.as_expr()),
                   'H_true_pos_roots_not_1': [(str(a), str(b)) for a, b in Hroots], 'verdict': verdict})
    print(f"type {sel} lineset {li} |L|={len(L)}: H = {H.as_expr()} ; positive real roots != 1 (true embedding): {Hroots} => {verdict}  [{time.time()-t0:.0f}s]", flush=True)
json.dump(report, open(outname, 'w'), indent=1)
print("SUMMARY:", sum(1 for e in report if e['verdict'] == 'EXCLUDED'), "of", len(report), "open cases excluded")
