"""Independent algebraic certificate of the Angle Lemma (factorisation route, no Groebner).
Lines L_i: y = m_i x + c_i with m_1 = 0, c_1 = 0 (WLOG by rotation+translation).
For each diagonal quadruple the concyclicity determinant factors into non-degeneracy factors times a
linear-in-each-variable 'essential' factor A_i(m2,m3,m4); the three essential factors cannot vanish
simultaneously under m3 != m4 etc.
"""
import sympy as sp, itertools
m2, m3, m4, c2, c3, c4 = sp.symbols('m2 m3 m4 c2 c3 c4')
L = {1: (sp.Integer(0), sp.Integer(0)), 2: (m2, c2), 3: (m3, c3), 4: (m4, c4)}
def P(i, j):
    mi, ci = L[i]; mj, cj = L[j]; x = (cj - ci) / (mi - mj); return (x, mi * x + ci)
def concyc_num(pts):
    M = sp.Matrix([[x**2 + y**2, x, y, 1] for x, y in pts]); return sp.numer(sp.together(M.det()))
Q = {1: [P(1,2), P(3,4), P(1,3), P(2,4)], 2: [P(1,2), P(3,4), P(1,4), P(2,3)], 3: [P(1,3), P(2,4), P(1,4), P(2,3)]}
ess = {}
for i, q in Q.items():
    f = sp.factor_list(concyc_num(q))
    print(f"det{i}: constant {f[0]}, factors:")
    for g, e in f[1]:
        print(f"    ({g})^{e}")
    ess[i] = [g for g, e in f[1] if g.has(m2) and g.has(m3) and g.has(m4)]
print("essential factors:", ess)
A = {i: ess[i][0] for i in ess}
print("A1 - A2 =", sp.expand(A[1] - A[2]))
print("A1 - A3 =", sp.expand(A[1] - A[3]))
print("A2 - A3 =", sp.expand(A[2] - A[3]))
