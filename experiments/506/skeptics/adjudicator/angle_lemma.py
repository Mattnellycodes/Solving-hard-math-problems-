"""Adjudicator's own algebraic certificate of the Angle Lemma, in ONE chart that covers every
complete quadrilateral (this closes the chart gap of the theory's angle_lemma_check.py, which
cannot represent a line perpendicular to L1).

Normalisation (all WLOG, exact): rotate so that L1 is y = 0.  Every line NOT parallel to L1 can be
written as x = m*y + c (vertical lines included: m = 0).  Since the four lines are pairwise
non-parallel, L2, L3, L4 are all of this form.  Translating along L1 and scaling fixes c2 = 0,
c3 = 1 (allowed because c2 != c3 is one of the non-degeneracy conditions: P12 != P13).
Unknowns: m2, m3, m4, c4.
Vertices: P1i = (c_i, 0);  P_ij (i,j >= 2): y = (c_j - c_i)/(m_i - m_j), x = m_i*y + c_i.
Non-degeneracy = pairwise non-parallel (m's distinct) + no three concurrent
  (c's distinct for triples containing L1; P23 not on L4 for {2,3,4}); this also makes the six
  vertices distinct.
Claim: the three concyclicity determinants Q1,Q2,Q3 (rows [x^2+y^2, x, y, 1]) cannot vanish
simultaneously under non-degeneracy: the ideal (Q1,Q2,Q3, t*nondeg - 1) is the unit ideal.
We also show the sharper facts: any TWO of the three conditions force a perpendicularity (as in
the synthetic proof), and two conditions ARE simultaneously realisable (the lemma is sharp).
"""
import sympy as sp
import random

m2, m3, m4, c4, t = sp.symbols('m2 m3 m4 c4 t')
c2, c3 = sp.Integer(0), sp.Integer(1)
M = {2: m2, 3: m3, 4: m4}
Cc = {2: c2, 3: c3, 4: c4}


def vertex(i, j):
    if i == 1:
        return (Cc[j], sp.Integer(0))
    y = (Cc[j] - Cc[i]) / (M[i] - M[j])
    return (sp.together(M[i] * y + Cc[i]), y)


P = {(1, 2): vertex(1, 2), (1, 3): vertex(1, 3), (1, 4): vertex(1, 4),
     (2, 3): vertex(2, 3), (2, 4): vertex(2, 4), (3, 4): vertex(3, 4)}


def concyc(*pts):
    rows = [[x * x + y * y, x, y, 1] for (x, y) in pts]
    return sp.factor(sp.numer(sp.together(sp.Matrix(rows).det())))


Q1 = concyc(P[(1, 2)], P[(3, 4)], P[(1, 3)], P[(2, 4)])
Q2 = concyc(P[(1, 2)], P[(3, 4)], P[(1, 4)], P[(2, 3)])
Q3 = concyc(P[(1, 3)], P[(2, 4)], P[(1, 4)], P[(2, 3)])

# concurrency of L2,L3,L4: P23 on L4  <=>  x23 - m4*y23 - c4 = 0
x23, y23 = P[(2, 3)]
conc234 = sp.factor(sp.numer(sp.together(x23 - m4 * y23 - c4)))
nondeg_factors = [m2 - m3, m2 - m4, m3 - m4, c4, c4 - 1, conc234]  # c2=0,c3=1 already distinct
print("non-degeneracy factors:", nondeg_factors)
for name, Q in (("Q1", Q1), ("Q2", Q2), ("Q3", Q3)):
    print(name, "=", Q)


def essential(Q):
    """strip factors that are (up to sign/power) non-degeneracy polynomials"""
    fs = sp.factor_list(Q)[1]
    ess = []
    for fac, mult in fs:
        if any(sp.simplify(fac - g) == 0 or sp.simplify(fac + g) == 0 for g in nondeg_factors):
            continue
        ess.append(fac)
    return ess


E = [essential(Q) for Q in (Q1, Q2, Q3)]
print("essential factors:", E)
nd = sp.prod(nondeg_factors)
gens = (m2, m3, m4, c4, t)
G = sp.groebner([sp.expand(Q1), sp.expand(Q2), sp.expand(Q3), sp.expand(t * nd - 1)], *gens, order='grevlex')
print("Groebner(Q1,Q2,Q3, t*nondeg-1) =", list(G))
print("ANGLE LEMMA CERTIFIED (unit ideal):", list(G) == [1])

# sharper: any two conditions -> perpendicularity;  In this chart L_i has direction (m_i, 1),
# L1 has direction (1,0): L1 perp Li <=> m_i = 0 ; Li perp Lj <=> m_i m_j + 1 = 0.
for (a, b, name) in ((Q1, Q2, "Q1,Q2"), (Q1, Q3, "Q1,Q3"), (Q2, Q3, "Q2,Q3")):
    G2 = sp.groebner([sp.expand(a), sp.expand(b), sp.expand(t * nd - 1)], *gens, order='lex')
    print(f"lex Groebner({name}, nondeg):", [sp.factor(g) for g in G2][:6])

# numeric sanity: random non-degenerate quadrilaterals never have all three near zero, and an
# explicit two-condition family exists (sharpness): L1 perp L2 (m2=0) then Q1,Q2 reduce to m3*m4+1=0.
random.seed(1)
worst = None
for _ in range(20000):
    vals = {m2: random.uniform(-3, 3), m3: random.uniform(-3, 3), m4: random.uniform(-3, 3), c4: random.uniform(-3, 3)}
    ndv = abs(float(nd.subs(vals)))
    if ndv < 1e-3:
        continue
    qs = [abs(float(Q.subs(vals))) for Q in (Q1, Q2, Q3)]
    mx = max(qs)
    if worst is None or mx < worst[0]:
        worst = (mx, vals)
print("numeric: min over random non-degenerate quadrilaterals of max|Q_k| =", worst[0])
fam = {m2: 0, m3: 2, m4: sp.Rational(-1, 2), c4: 3}
print("two-condition family (L1 perp L2, L3 perp L4): Q1,Q2,Q3 =", [Q.subs(fam) for Q in (Q1, Q2, Q3)],
      " nondeg =", nd.subs(fam))
