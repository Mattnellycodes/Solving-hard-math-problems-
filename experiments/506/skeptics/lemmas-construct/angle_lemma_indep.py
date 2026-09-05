"""Independent algebraic check of the Angle Lemma (my own parametrisation, not the report's).
Lines: L1: y=0, L2: y=m2 x, L3: y=m3 x+1 (scale), L4: y = m4 x + c.  Vertices P_ij = L_i ∩ L_j.
Concyclicity of a quadruple: det[x, y, x^2+y^2, 1] = 0.
Claim: the three quadruples {P12,P34,P13,P24}, {P12,P34,P14,P23}, {P13,P24,P14,P23} cannot all be
concyclic when the lines are pairwise non-parallel and no three concurrent.
Method: eliminate with resultants / Groebner basis, saturating by non-degeneracy.
"""
import sympy as sp
m2, m3, m4, c, t = sp.symbols('m2 m3 m4 c t')
def meet(l1, l2):
    (a1, b1), (a2, b2) = l1, l2   # y = a x + b
    x = (b2 - b1) / (a1 - a2); return sp.simplify(x), sp.simplify(a1 * x + b1)
L = {1: (0, 0), 2: (m2, 0), 3: (m3, 1), 4: (m4, c)}
P = {(i, j): meet(L[i], L[j]) for i in range(1, 5) for j in range(i + 1, 5)}
def conc(*pts):
    return sp.factor(sp.numer(sp.together(sp.Matrix([[x, y, x**2 + y**2, 1] for x, y in pts]).det())))
Q1 = conc(P[(1,2)], P[(3,4)], P[(1,3)], P[(2,4)])
Q2 = conc(P[(1,2)], P[(3,4)], P[(1,4)], P[(2,3)])
Q3 = conc(P[(1,3)], P[(2,4)], P[(1,4)], P[(2,3)])
print('Q1 =', Q1); print('Q2 =', Q2); print('Q3 =', Q3)
# non-degeneracy: slopes distinct and nonzero-difference (non-parallel): m2, m3, m4, m2-m3, m2-m4, m3-m4 != 0
# no three concurrent / vertices distinct: c != 0 (L4 != through origin... L4 through P12=(0,0) iff c=0), c != 1 (P13 != P14? P13=(-1/m3,0), P14=(-c/m4,0): equal iff m4 = c m3),
# P23 == P24 iff 1/(m2-m3) == c/(m2-m4); P34 on L1 iff ...; we include the vertex-distinctness conditions generically:
nd = m2 * m3 * m4 * (m2 - m3) * (m2 - m4) * (m3 - m4) * c * (m4 - c * m3) * ((m2 - m4) - c * (m2 - m3)) * (c - 1)
# also P34 != P12: P34 = ((c-1)/(m3-m4), ...) distinct from origin iff c != 1 (included).
G = sp.groebner([Q1, Q2, Q3, t * nd - 1], m2, m3, m4, c, t, order='grevlex')
print('Groebner basis with saturation:', list(G))
print('unit ideal ->', list(G) == [1])
# extra: show what Q1 and Q2 together force (without Q3)
G2 = sp.groebner([Q1, Q2, t * nd - 1], t, c, m4, m3, m2, order='lex')
print('Q1&Q2 saturated, lex basis:'); [print('  ', sp.factor(g)) for g in G2]
