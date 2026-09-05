"""Correct non-realisability proof of the Möbius–Kantor (8_3) configuration over R (2026-09-05).
Lines {i, i+1, i+3} mod 8. Frame p0=(1,0,0), p1=(0,1,0), p2=(0,0,1), p5=(1,1,1): valid, since no line
contains three of {0,1,2,5}. Forced coordinates:
  p3 on line p0p1 = {z=0}: p3 = (1, s, 0), s != 0 (s=0 would be p0)            [line {0,1,3}]
  p4 on line p1p2 = {x=0}: p4 = (0, 1, t), t != 0? (t=0 is p1, excluded)       [line {1,2,4}]
  p6 = (p5 x p0) x (p3 x p4)   [lines {5,6,0} and {3,4,6}]
  p7 = (p4 x p5) x (p0 x p2)   [lines {4,5,7} and {7,0,2}]
Remaining incidences: {2,3,5}: det(p2,p3,p5)=0 ; {6,7,1}: det(p6,p7,p1)=0.
Non-degeneracy: all 8 points distinct (as projective points), p6, p7 nonzero.
"""
import sympy as sp, itertools
s, t = sp.symbols('s t')
p = {0: sp.Matrix([1,0,0]), 1: sp.Matrix([0,1,0]), 2: sp.Matrix([0,0,1]), 5: sp.Matrix([1,1,1]),
     3: sp.Matrix([1, s, 0]), 4: sp.Matrix([0, 1, t])}
p[6] = (p[5].cross(p[0])).cross(p[3].cross(p[4]))
p[7] = (p[4].cross(p[5])).cross(p[0].cross(p[2]))
p[6] = sp.simplify(p[6]); p[7] = sp.simplify(p[7])
print("p6 =", p[6].T, " p7 =", p[7].T)
lines = [((i) % 8, (i + 1) % 8, (i + 3) % 8) for i in range(8)]
eqs = []
for (i, j, k) in lines:
    d = sp.expand(sp.Matrix.hstack(p[i], p[j], p[k]).det())
    print(f"line {i,j,k}: det = {sp.factor(d)}")
    if d != 0: eqs.append(d)
print("equations:", eqs)
sols = sp.solve(eqs, [s, t], dict=True)
print("all solutions over C:", sols)
for sol in sols:
    pts = {i: sp.simplify(v.subs(sol)) for i, v in p.items()}
    # distinctness check (projective): rank of pair matrix == 2
    ok = all(sp.Matrix.hstack(pts[i], pts[j]).rank() == 2 for i, j in itertools.combinations(range(8), 2))
    nonzero = all(any(c != 0 for c in pts[i]) for i in range(8))
    real = all(sp.im(v) == 0 for v in sol.values())
    print(f"  solution {sol}: points distinct={ok}, nonzero={nonzero}, real={real}")
# also cover the boundary charts: s=0 or t=0 mean p3=p0 or p4=p1 (excluded by distinctness);
# p6 or p7 = zero vector means two of the defining lines coincide -> then check separately:
print("p6 zero-vector condition:", sp.factor(sp.gcd(list(p[6]))), "| p7 zero-vector condition:", sp.factor(sp.gcd(list(p[7]))))
