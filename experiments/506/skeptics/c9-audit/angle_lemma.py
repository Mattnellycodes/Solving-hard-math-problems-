"""Own check of the Angle Lemma (Lemma 2.3) -- NOT load-bearing for n = 9, audited as requested.
Four lines y = m_i x + c_i (no line vertical after a rotation), pairwise non-parallel, no three
concurrent; P_ij = L_i ∩ L_j.  Claim: the three diagonal quadruples {P12,P34,P13,P24},
{P12,P34,P14,P23}, {P13,P24,P14,P23} are never all concyclic.
Method: (a) the directed-angle argument: concyclicity of {P12,P34,P13,P24} <=> angle(L1,L2) = angle(L3,L4)
as directed angles mod pi <=> (m2-m1)(1+m3 m4) = (m4-m3)(1+m1 m2).  We verify symbolically that the
concyclicity determinant equals this expression times non-degeneracy factors.  (b) the three relations
+ non-degeneracy generate the unit ideal (Groebner basis with Rabinowitsch).  (c) random numeric test."""
import sympy as sp, random
from fractions import Fraction as Fr
m1, m2, m3, m4, c1, c2, c3, c4, t = sp.symbols('m1 m2 m3 m4 c1 c2 c3 c4 t')
L = {1: (m1, c1), 2: (m2, c2), 3: (m3, c3), 4: (m4, c4)}
def P(i, j):
    mi, ci = L[i]; mj, cj = L[j]
    x = (cj - ci) / (mi - mj); return (x, mi * x + ci)
def concyc(a, b, c, d):
    return sp.Matrix([[x, y, x * x + y * y, 1] for (x, y) in (a, b, c, d)]).det()
def angle_rel(i, j, k, l):   # directed angle from L_i to L_j equals that from L_k to L_l (mod pi)
    return (L[j][0] - L[i][0]) * (1 + L[k][0] * L[l][0]) - (L[l][0] - L[k][0]) * (1 + L[i][0] * L[j][0])
Q = [concyc(P(1, 2), P(3, 4), P(1, 3), P(2, 4)), concyc(P(1, 2), P(3, 4), P(1, 4), P(2, 3)), concyc(P(1, 3), P(2, 4), P(1, 4), P(2, 3))]
A = [angle_rel(1, 2, 3, 4), angle_rel(1, 2, 4, 3), angle_rel(1, 3, 4, 2)]
# (concyclic {P12,P34,P13,P24}: inscribed angles at P13 and P24 subtending P12P34: angle(L1,L3)... we let sympy decide
#  which angle relation divides which determinant, trying all natural candidates.)
cands = {"A(12|34)": angle_rel(1, 2, 3, 4), "A(12|43)": angle_rel(1, 2, 4, 3), "A(13|24)": angle_rel(1, 3, 2, 4),
         "A(13|42)": angle_rel(1, 3, 4, 2), "A(14|23)": angle_rel(1, 4, 2, 3), "A(14|32)": angle_rel(1, 4, 3, 2)}
gens = (m1, m2, m3, m4, c1, c2, c3, c4)
factors_used = []
for idx, q in enumerate(Q):
    num = sp.Poly(sp.expand(sp.numer(sp.together(q))), *gens)
    print(f"Q{idx+1} numerator factorisation: {sp.factor(num.as_expr())}")
    for name, a in cands.items():
        quo, rem = sp.div(num, sp.Poly(sp.expand(a), *gens))
        if rem.is_zero:
            print(f"   divisible by {name}; cofactor = {sp.factor(quo.as_expr())}")
            factors_used.append(sp.expand(a))
            break
# the three angle relations that occur, saturated by pairwise non-parallel: unit ideal?
nd = 1
for i, j in [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]:
    nd *= (L[i][0] - L[j][0])
G = sp.groebner(factors_used + [t * nd - 1], m1, m2, m3, m4, t, order='grevlex')
print("Groebner basis of the three angle relations + 'pairwise non-parallel':", list(G), "=> unit ideal:", list(G) == [1])
# numeric random test: count how many of the three quadruples can be concyclic simultaneously
random.seed(3)
best = 0
for _ in range(20000):
    ms = [Fr(random.randint(-9, 9), random.randint(1, 5)) for _ in range(4)]
    cs = [Fr(random.randint(-9, 9), random.randint(1, 5)) for _ in range(4)]
    if len(set(ms)) < 4: continue
    Lx = {i + 1: (ms[i], cs[i]) for i in range(4)}
    def Pt(i, j):
        mi, ci = Lx[i]; mj, cj = Lx[j]; x = (cj - ci) / (mi - mj); return (x, mi * x + ci)
    pts = {(i, j): Pt(i, j) for i, j in [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]}
    if len(set(pts.values())) < 6: continue
    def cc(a, b, c, d):
        M = [[x, y, x * x + y * y, 1] for (x, y) in (a, b, c, d)]
        return sp.Matrix(M).det()
    k = sum(1 for (a, b, c, d) in [((1, 2), (3, 4), (1, 3), (2, 4)), ((1, 2), (3, 4), (1, 4), (2, 3)), ((1, 3), (2, 4), (1, 4), (2, 3))]
            if cc(pts[a], pts[b], pts[c], pts[d]) == 0)
    best = max(best, k)
print("random rational quadrilaterals: max number of simultaneously concyclic diagonal quadruples found:", best)
