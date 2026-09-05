"""lemmas-audit: own proof that the Moebius-Kantor configuration (8_3) has no realisation by 8
distinct points of RP^2 in which {p0,p1,p2,p5} is a projective frame (no three collinear).
In the application (derived structure of a degree-8 point: 8 four-blocks through p, no other block of
size >= 4 through p) the 8 derived lines are ALL the rich lines of the inverted 8-point set, so any four
points with no three on an MK line (e.g. 0,1,2,5) form a frame.  Hence the derived structure of a
point lying in 8 four-blocks cannot exist.

Case analysis (each step forced, checked with sympy):
  lines: {0,1,3},{1,2,4},{2,3,5},{3,4,6},{4,5,7},{5,6,0},{6,7,1},{7,0,2}
  p3 = line(0,1) ∩ line(2,5)                       (two distinct lines through frame points)
  p4 ∈ line(1,2) = {x=0}: p4 = (0,c,1), c != 0    ((0,1,0) = p1 and (0,0,1) = p2 are excluded)
  p6 ∈ line(5,0) = {y=z}: p6 = (a,1,1)            ((1,0,0) = p0 excluded)
  p7 ∈ line(0,2) = {y=0}: p7 = (b,0,1)            ((1,0,0) = p0 excluded)
  remaining incidences {3,4,6},{4,5,7},{6,7,1} -> polynomial system in a,b,c -> c^2 - c + 1 = 0.
Also: numerical random search for near-realisations (sanity, not proof).
"""
import itertools, random
import sympy as sp

MK = [frozenset({i, (i + 1) % 8, (i + 3) % 8}) for i in range(8)]
frame = (0, 1, 2, 5)
assert not any(len(l & set(frame)) >= 3 for l in MK), "frame has three points on an MK line"

a, b, c = sp.symbols("a b c")
p = {0: sp.Matrix([1, 0, 0]), 1: sp.Matrix([0, 1, 0]), 2: sp.Matrix([0, 0, 1]), 5: sp.Matrix([1, 1, 1])}
# p3 = line(0,1) x line(2,5)
l01 = p[0].cross(p[1]); l25 = p[2].cross(p[5])
p[3] = l01.cross(l25)
print("p3 =", list(p[3]))
# p4 on line(1,2): general point (x,y,z) with det(p1,p2,p4)=0 -> x = 0.  Charts: z=1 -> (0,c,1); z=0 -> (0,1,0)=p1.
# We record the chart argument by checking the excluded alternatives are frame points:
assert list(sp.Matrix([0, 1, 0])) == list(p[1])
p[4] = sp.Matrix([0, c, 1])
# p6 on line(5,0): det(p5,p0,p6)=0 -> y = z. Charts: z=1 -> (a,1,1); z=0 -> y=0 -> (1,0,0)=p0.
p[6] = sp.Matrix([a, 1, 1])
# p7 on line(0,2): y=0. Charts: z=1 -> (b,0,1); z=0 -> (1,0,0)=p0.
p[7] = sp.Matrix([b, 0, 1])
eqs = {}
for L in MK:
    i, j, k = sorted(L)
    e = sp.expand(sp.Matrix.hstack(p[i], p[j], p[k]).det())
    eqs[(i, j, k)] = e
    print(f"line {(i,j,k)}: incidence polynomial = {e}")
E = [e for e in eqs.values() if e != 0]
G = sp.groebner(E, a, b, c, order="lex")
print("Groebner basis:", list(G))
uni = [g for g in G if g.free_symbols <= {c}]
print("univariate in c:", uni, " real roots:", [sp.real_roots(sp.Poly(g, c)) for g in uni])
# distinctness sanity: c=0 would give p4=p2 (excluded anyway), check no other coincidences needed.
sols = sp.solve(E, [a, b, c], dict=True)
print("all complex solutions:", sols)
for s in sols:
    pts = {i: [sp.simplify(v.subs(s)) for v in p[i]] for i in p}
    print("  solution", s, "-> points", pts)
    for i, j in itertools.combinations(range(8), 2):
        M = sp.Matrix([pts[i], pts[j]])
        if M.rank() < 2:
            print("   coincident points", i, j)
print("CONCLUSION: every realisation has c^2 - c + 1 = 0 -> no real realisation with the frame {0,1,2,5}.")

# numeric sanity: random 8 real points near MK incidences?  Minimise sum of squared normalised
# determinants from random starts with scipy-free gradient-free perturbation; report best residual.
import numpy as np
def residual(X):
    r = 0.0
    for L in MK:
        i, j, k = sorted(L)
        M = np.array([[X[i][0], X[i][1], 1.0], [X[j][0], X[j][1], 1.0], [X[k][0], X[k][1], 1.0]])
        d = np.linalg.det(M)
        s = np.linalg.norm(X[i] - X[j]) * np.linalg.norm(X[j] - X[k]) * np.linalg.norm(X[i] - X[k]) + 1e-12
        r += (d / s) ** 2
    # penalise coincident points
    for i, j in itertools.combinations(range(8), 2):
        dd = np.linalg.norm(X[i] - X[j])
        if dd < 0.05:
            r += (0.05 - dd) * 100
    return r
rng = np.random.default_rng(0)
best = 1e9
for trial in range(30):
    X = rng.normal(size=(8, 2))
    step = 0.5
    val = residual(X)
    for it in range(3000):
        Y = X + rng.normal(size=(8, 2)) * step
        v = residual(Y)
        if v < val:
            X, val = Y, v
        else:
            step *= 0.999
    best = min(best, val)
print(f"numerical sanity: best normalised residual over 30 random descents = {best:.3e} (a realisation would give ~0)")
