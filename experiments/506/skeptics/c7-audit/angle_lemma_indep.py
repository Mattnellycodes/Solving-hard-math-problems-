"""Skeptic c7-audit: independent algebraic verification of the Angle Lemma (REPORT Lemma 2.3).

Parametrisation (different from theory/angle_lemma_check.py, which uses slopes): a complete
quadrilateral is fixed up to similarity by its vertices
   P12 = (0,0), P13 = (1,0)            -> L1 = x-axis
   P23 = (u, v), v != 0                 -> L2 = line(P12,P23), L3 = line(P13,P23)
   P14 = (s, 0), P24 = t*(u, v)         -> L4 = line(P14,P24); P34 = L3 ∩ L4.
Non-degeneracy (pairwise non-parallel, no three concurrent, six distinct vertices) is exactly
   v * s * (s-1) * t * (t-1) * (t-s) != 0     (derived by hand, see comments).
The three diagonal quadruples are
   Q1 = {P12,P34,P13,P24}, Q2 = {P12,P34,P14,P23}, Q3 = {P13,P24,P14,P23}.
We show: the ideal (N1, N2, N3) saturated by the non-degeneracy product is the unit ideal, so no
complex (a fortiori no real) non-degenerate quadrilateral has all three quadruples concyclic.
We also show the sharper structure: N1 = N2 = 0 forces u = 0 (L1 ⊥ L2) and N1 = N3 = 0 forces
u = 1 (L1 ⊥ L3), which is the hand proof's mechanism, and give a numeric two-of-three example.
"""
import sympy as sp
u, v, s, t, z = sp.symbols('u v s t z')
P12 = (sp.Integer(0), sp.Integer(0)); P13 = (sp.Integer(1), sp.Integer(0)); P23 = (u, v)
P14 = (s, sp.Integer(0)); P24 = (t * u, t * v)
mu = (1 - s) / (t - s)                   # L4 parameter of P34 (solved by hand: lambda = mu*t)
P34 = (s + mu * (t * u - s), mu * t * v)
# sanity: P34 on L3 (through (1,0) and (u,v)) and on L4
def collin(a, b, c):
    return sp.simplify((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))
assert collin(P13, P23, P34) == 0 and collin(P14, P24, P34) == 0
def concyc(A, B, C, D):
    M = sp.Matrix([[x, y, x * x + y * y, 1] for (x, y) in (A, B, C, D)])
    return sp.factor(sp.numer(sp.together(M.det())))
N1 = concyc(P12, P34, P13, P24)
N2 = concyc(P12, P34, P14, P23)
N3 = concyc(P13, P24, P14, P23)
print("N1 =", N1); print("N2 =", N2); print("N3 =", N3)
nd = v * s * (s - 1) * t * (t - 1) * (t - s)
G = sp.groebner([sp.expand(N1), sp.expand(N2), sp.expand(N3), sp.expand(1 - z * nd)], u, v, s, t, z, order='grevlex')
print("Groebner basis of (N1,N2,N3) : nd^inf  =", list(G))
print("ANGLE LEMMA (unit ideal => no non-degenerate solution over C):", list(G) == [1])
G12 = sp.groebner([sp.expand(N1), sp.expand(N2), sp.expand(1 - z * nd)], z, v, s, t, u, order='lex')
print("N1=N2=0 saturated, lex basis:", [sp.factor(g) for g in G12])
G13 = sp.groebner([sp.expand(N1), sp.expand(N3), sp.expand(1 - z * nd)], z, v, s, t, u, order='lex')
print("N1=N3=0 saturated, lex basis:", [sp.factor(g) for g in G13])
# numeric two-of-three example: L1 ⊥ L2 (u=0) and L3 ⊥ L4: with u=0, L3 has direction (-1, v),
# L4 direction (t*0 - s, t*v) = (-s, t v); perpendicular iff s + t v^2 = 0.
ex = {u: 0, v: 2, s: -8, t: 1}  # t=1 is degenerate (L4 through P23) -> pick t=2: s = -t v^2 = -8
ex = {u: 0, v: 2, t: 2, s: -8}
print("example u=0,v=2,t=2,s=-8: nd =", nd.subs(ex), " N1,N2,N3 =", [sp.simplify(N.subs(ex)) for N in (N1, N2, N3)])
