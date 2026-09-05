import sympy as sp
m2, m3, m4, c, t = sp.symbols('m2 m3 m4 c t')
def meet(l1, l2):
    (a1, b1), (a2, b2) = l1, l2
    x = (b2 - b1) / (a1 - a2); return sp.cancel(x), sp.cancel(a1 * x + b1)
L = {1: (0, 0), 2: (m2, 0), 3: (m3, 1), 4: (m4, c)}
P = {(i, j): meet(L[i], L[j]) for i in range(1, 5) for j in range(i + 1, 5)}
def conc(*pts):
    return sp.factor(sp.numer(sp.cancel(sp.Matrix([[x, y, x**2 + y**2, 1] for x, y in pts]).det())))
Q1 = conc(P[(1,2)], P[(3,4)], P[(1,3)], P[(2,4)])
Q2 = conc(P[(1,2)], P[(3,4)], P[(1,4)], P[(2,3)])
Q3 = conc(P[(1,3)], P[(2,4)], P[(1,4)], P[(2,3)])
for name, Q in (('Q1', Q1), ('Q2', Q2), ('Q3', Q3)):
    print(name, '=', Q)
# tangent form of theta_i + theta_j = theta_k + theta_l (mod pi), with m1 = 0:
def tanform(mi, mj, mk, ml):   # (mi+mj)(1-mk ml) - (mk+ml)(1-mi mj)
    return sp.expand((mi + mj) * (1 - mk * ml) - (mk + ml) * (1 - mi * mj))
A1 = tanform(0, m4, m2, m3); A2 = tanform(0, m3, m2, m4); A3 = tanform(0, m2, m3, m4)
print('angle-form factors:', A1, '|', A2, '|', A3)
# check each Qi is (degeneracy factors) * Ai up to sign: divide
for Q, A in ((Q1, A1), (Q2, A2), (Q3, A3)):
    q, r = sp.div(sp.Poly(Q, m2, m3, m4, c), sp.Poly(A, m2, m3, m4, c))
    print('  remainder zero:', r.is_zero, ' cofactor:', sp.factor(q.as_expr()))
# now the three angle-form conditions alone, saturated by pairwise non-parallel (m2,m3,m4 nonzero & distinct)
nd = m2 * m3 * m4 * (m2 - m3) * (m2 - m4) * (m3 - m4)
G = sp.groebner([A1, A2, A3, t * nd - 1], m2, m3, m4, t, order='grevlex')
print('Groebner of angle conditions + non-parallel saturation:', list(G), ' unit ideal:', list(G) == [1])
# and the linear-in-angles argument, numerically over the unit circle: A1=A2=0 => ?
print('A1 - A2 =', sp.factor(A1 - A2), '   A1 + A3 =', sp.factor(A1 + A3), '  A2 + A3 =', sp.factor(A2 + A3))
