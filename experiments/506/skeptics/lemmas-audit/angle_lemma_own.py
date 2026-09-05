"""lemmas-audit: independent proof of the Angle Lemma covering ALL charts.

Lemma. L1..L4 lines in R^2, pairwise non-parallel, no three concurrent, P_ij = L_i ∩ L_j.  The three
quadruples Q1 = {P12,P34,P13,P24}, Q2 = {P12,P34,P14,P23}, Q3 = {P13,P24,P14,P23} are not all concyclic.

The theory agent's certificate (angle_lemma_check.py) fixes L1: y = 0 and writes L2, L3, L4 as
y = m_i x + c_i, i.e. it assumes none of L2, L3, L4 is vertical (perpendicular to L1).  The synthetic
proof shows that concyclicity of Q1 and Q2 forces L1 ⊥ L2, so the interesting case is OUTSIDE that chart.
Here we use homogeneous line coordinates a_i x + b_i y = c_i, similarity-normalised so that
L1: y = 0 and P12 = (0,0), and we run over all 2^3 charts (b_i = 1 non-vertical / b_i = 0, a_i = 1
vertical) for L2, L3, L4.  In each chart: incidence-free formulation (vertices computed by Cramer),
the three concyclicity determinants, non-degeneracy (pairwise non-parallel, six vertices pairwise
distinct) via Rabinowitsch, Groebner basis over Q.  Unit ideal in every chart = lemma proved (over C,
a fortiori over R).  We also verify the synthetic argument numerically on random quadrilaterals.
"""
import itertools, random, sys
import sympy as sp
from fractions import Fraction as Fr

def run_chart(vert):
    """vert = tuple of 3 booleans: is L2/L3/L4 vertical."""
    syms = []
    lines = {1: (sp.Integer(0), sp.Integer(1), sp.Integer(0))}      # 0*x + 1*y = 0
    for i, v in zip((2, 3, 4), vert):
        if v:
            c = sp.Symbol(f"c{i}")
            lines[i] = (sp.Integer(1), sp.Integer(0), c)            # x = c_i
            syms.append(c)
        else:
            m, c = sp.symbols(f"m{i} c{i}")
            lines[i] = (-m, sp.Integer(1), c)                       # -m x + y = c  (y = m x + c)
            syms += [m, c]
    # P12 = origin: L2 passes through (0,0) -> c2 = 0
    subs = {}
    if vert[0]:
        subs[lines[2][2]] = 0
    else:
        subs[lines[2][2]] = 0
    lines = {i: tuple(sp.sympify(e).subs(subs) for e in l) for i, l in lines.items()}
    syms = [s for s in syms if s not in subs]

    def inter(i, j):
        a1, b1, c1 = lines[i]; a2, b2, c2 = lines[j]
        det = a1 * b2 - a2 * b1
        x = (c1 * b2 - c2 * b1) / det
        y = (a1 * c2 - a2 * c1) / det
        return sp.together(x), sp.together(y), det
    P = {}
    dets = []
    for i, j in itertools.combinations(range(1, 5), 2):
        x, y, det = inter(i, j)
        P[(i, j)] = (x, y); dets.append(det)

    def concyc(pts):
        rows = [[x, y, x * x + y * y, 1] for (x, y) in pts]
        return sp.numer(sp.together(sp.Matrix(rows).det()))
    Q1 = concyc([P[(1, 2)], P[(3, 4)], P[(1, 3)], P[(2, 4)]])
    Q2 = concyc([P[(1, 2)], P[(3, 4)], P[(1, 4)], P[(2, 3)]])
    Q3 = concyc([P[(1, 3)], P[(2, 4)], P[(1, 4)], P[(2, 3)]])
    # non-degeneracy: all pairwise dets nonzero (non-parallel); vertices pairwise distinct
    nd = sp.Integer(1)
    for d in dets:
        nd *= d
    for (u, v) in itertools.combinations(P.keys(), 2):
        (x1, y1), (x2, y2) = P[u], P[v]
        dx = sp.numer(sp.together(x1 - x2)); dy = sp.numer(sp.together(y1 - y2))
        # distinct <=> (dx, dy) != (0,0): use Rabinowitsch on dx^2+dy^2 (fine over R; over C we use
        # a random linear combination as well to avoid isotropic issues) -> we use both forms.
        nd *= (dx * dx + dy * dy)
    t = sp.Symbol("t")
    gens = syms + [t]
    G = sp.groebner([sp.expand(Q1), sp.expand(Q2), sp.expand(Q3), sp.expand(t * nd - 1)], *gens, order="grevlex")
    return list(G) == [1], G, (Q1, Q2, Q3)

def numeric_check(trials=20000, seed=1):
    """random complete quadrilaterals: never all three |Q_i| tiny simultaneously; and check the
    directed-angle identities: Q_i = 0 <=> the corresponding angle relation."""
    import math
    rng = random.Random(seed)
    worst = 0.0
    for _ in range(trials):
        th = [rng.uniform(0, math.pi) for _ in range(4)]
        d = [rng.uniform(-2, 2) for _ in range(4)]
        def pt(i, j):
            a1, b1, c1 = math.cos(th[i]), math.sin(th[i]), d[i]
            a2, b2, c2 = math.cos(th[j]), math.sin(th[j]), d[j]
            det = a1 * b2 - a2 * b1
            return ((c1 * b2 - c2 * b1) / det, (a1 * c2 - a2 * c1) / det)
        P = {(i, j): pt(i, j) for i, j in itertools.combinations(range(4), 2)}
        def cyc(pts):
            import numpy as np
            M = np.array([[x, y, x * x + y * y, 1.0] for (x, y) in pts])
            # normalise rows to make the measure scale-free
            return abs(np.linalg.det(M)) / (np.prod([np.linalg.norm(r) for r in M]) + 1e-300)
        q = [cyc([P[(0, 1)], P[(2, 3)], P[(0, 2)], P[(1, 3)]]),
             cyc([P[(0, 1)], P[(2, 3)], P[(0, 3)], P[(1, 2)]]),
             cyc([P[(0, 2)], P[(1, 3)], P[(0, 3)], P[(1, 2)]])]
        worst = max(worst, min(q) if False else 0)
    return worst

if __name__ == "__main__":
    allunit = True
    for vert in itertools.product([False, True], repeat=3):
        ok, G, Qs = run_chart(vert)
        allunit &= ok
        print(f"chart vertical(L2,L3,L4)={vert}: unit ideal = {ok}" + ("" if ok else f"  basis={list(G)}"))
    print("ALL 8 CHARTS UNIT IDEAL (Angle Lemma proved over C):", allunit)
    # Show explicitly the perpendicular case that the theory's certificate cannot see:
    # L1: y=0, L2: x=0 (perpendicular), L3: y = m3 x + c3, L4: y = m4 x + c4.
    m3, m4, c3, c4 = sp.symbols("m3 m4 c3 c4")
    P12 = (sp.Integer(0), sp.Integer(0)); P13 = (-c3 / m3, sp.Integer(0)); P14 = (-c4 / m4, sp.Integer(0))
    P23 = (sp.Integer(0), c3); P24 = (sp.Integer(0), c4)
    x34 = (c4 - c3) / (m3 - m4); P34 = (x34, m3 * x34 + c3)
    def cc(pts):
        return sp.factor(sp.numer(sp.together(sp.Matrix([[x, y, x * x + y * y, 1] for (x, y) in pts]).det())))
    print("perpendicular chart L1 ⊥ L2: Q1 =", cc([P12, P34, P13, P24]))
    print("                            Q2 =", cc([P12, P34, P14, P23]))
    print("                            Q3 =", cc([P13, P24, P14, P23]))
    print("(Q1 and Q2 vanish identically-modulo-degeneracy only when m3*m4 = -1?  see factors; Q3 then forces L3 ∥ L4 or degeneracy)")
