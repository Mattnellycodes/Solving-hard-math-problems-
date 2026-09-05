"""Own re-proof of the Angle Lemma 2.3 (not load-bearing for n=9, but audited as requested).
Four lines L_i: y = m_i x + c_i (after rotation none vertical), pairwise non-parallel, no three concurrent.
Claim: the three diagonal quadruples cannot all be concyclic.
Method 1 (synthetic, checked symbolically): the concyclicity condition for {P12,P34,P13,P24} is
equivalent to tan-angle identity  (m1,m4) ~ (m2,m3):  the 'sum of directions' relation
theta1+theta4 = theta2+theta3 (mod pi)  <=>  (m1+m4)(1-m2 m3) = (m2+m3)(1-m1 m4).
We verify that identity symbolically (the concyclicity determinant factors into that expression times
non-degeneracy factors) and then that the three relations force two lines parallel.
Method 2: Groebner basis with Rabinowitsch saturation."""
import sympy as sp
m1, m2, m3, m4, c1, c2, c3, c4, t = sp.symbols('m1 m2 m3 m4 c1 c2 c3 c4 t')
L = {1: (m1, c1), 2: (m2, c2), 3: (m3, c3), 4: (m4, c4)}
def inter(i, j):
    mi, ci = L[i]; mj, cj = L[j]
    x = (cj - ci)/(mi - mj); return (x, mi*x + ci)
Pt = {(i, j): inter(i, j) for i, j in [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]}
def concyc(a, b, c, d):
    return sp.Matrix([[x, y, x*x+y*y, 1] for (x, y) in (a, b, c, d)]).det()
Q = {1: concyc(Pt[(1, 2)], Pt[(3, 4)], Pt[(1, 3)], Pt[(2, 4)]),
     2: concyc(Pt[(1, 2)], Pt[(3, 4)], Pt[(1, 4)], Pt[(2, 3)]),
     3: concyc(Pt[(1, 3)], Pt[(2, 4)], Pt[(1, 4)], Pt[(2, 3)])}
def angle_rel(i, j, k, l):   # theta_i + theta_j = theta_k + theta_l (mod pi) in slope form
    return (L[i][0] + L[j][0])*(1 - L[k][0]*L[l][0]) - (L[k][0] + L[l][0])*(1 - L[i][0]*L[j][0])
A = {1: angle_rel(1, 4, 2, 3), 2: angle_rel(1, 3, 2, 4), 3: angle_rel(1, 2, 3, 4)}
for k in (1, 2, 3):
    num = sp.factor(sp.numer(sp.together(Q[k])))
    print(f"Q{k} numerator factors:", num)
    quo, rem = sp.div(sp.Poly(sp.expand(num), m1, m2, m3, m4, c1, c2, c3, c4), sp.Poly(sp.expand(A[k]), m1, m2, m3, m4, c1, c2, c3, c4))
    print(f"   divisible by angle relation A{k}: remainder zero ->", rem.is_zero, "; cofactor:", sp.factor(quo.as_expr()))
# The cofactors must be products of non-degeneracy factors (slope differences, intercept differences ...).
# Now: A1, A2, A3 together with pairwise non-parallel => contradiction.
nd = 1
for i, j in [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]:
    nd *= (L[i][0] - L[j][0])
G = sp.groebner([sp.expand(A[1]), sp.expand(A[2]), sp.expand(A[3]), t*nd - 1], m1, m2, m3, m4, t, order='grevlex')
print("Groebner basis of (A1,A2,A3) saturated by 'pairwise non-parallel':", list(G))
# Show how: A1 + A2 (in angle form) forces L1 ⊥ L2 ; do it with the polynomial form
G12 = sp.groebner([sp.expand(A[1]), sp.expand(A[2]), t*nd - 1], m1, m2, m3, m4, t, order='lex')
print("From A1, A2 alone (lex basis):", [sp.factor(g) for g in G12])
