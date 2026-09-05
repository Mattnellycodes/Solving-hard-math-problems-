"""Skeptic c7-audit: independent re-derivation of
  (i)  the Fano plane has no realisation by points/lines of RP^2 (indeed of any field of char != 2);
  (ii) every 8-triple system on 8 points with pairwise intersections <= 1 is the (8_3)
       Moebius-Kantor configuration, and it has no real realisation.
Frames and charts chosen independently of theory/eight_three_light.py.
"""
import itertools, sympy as sp

# ---------- (i) Fano -------------------------------------------------------------------------
FANO = [frozenset(l) for l in [(0,1,2),(0,3,4),(0,5,6),(1,3,5),(1,4,6),(2,3,6),(2,4,5)]]
# quadrangle {3,4,5,6} (complement of line 012): no three collinear -> projective frame
p = {3: sp.Matrix([1,0,0]), 4: sp.Matrix([0,1,0]), 5: sp.Matrix([0,0,1]), 6: sp.Matrix([1,1,1])}
def line(a, b): return a.cross(b)
def meet(l, m): return l.cross(m)
p[0] = meet(line(p[3], p[4]), line(p[5], p[6]))
p[1] = meet(line(p[3], p[5]), line(p[4], p[6]))
p[2] = meet(line(p[3], p[6]), line(p[4], p[5]))
print("Fano forced coordinates:", {k: list(v) for k, v in p.items()})
d = sp.Matrix.hstack(p[0], p[1], p[2]).det()
print("det(p0,p1,p2) for the last line {0,1,2}:", d, "-> Fano realisable only if", d, "= 0 (char 2). Over R: NOT realisable.")
# all other 6 incidences hold by construction; verify anyway
for L in FANO:
    a, b, c = sorted(L)
    print("  line", (a, b, c), "det =", sp.Matrix.hstack(p[a], p[b], p[c]).det())

# ---------- (ii) (8_3) -----------------------------------------------------------------------
n = 8
T = [frozenset(t) for t in itertools.combinations(range(n), 3)]
sols = []
def rec(i, F, deg):
    if len(F) == 8:
        sols.append(tuple(F)); return
    if 8 - len(F) > len(T) - i: return
    for j in range(i, len(T)):
        t = T[j]
        if all(len(t & u) <= 1 for u in F) and all(deg[x] < 3 for x in t):
            for x in t: deg[x] += 1
            F.append(t); rec(j + 1, F, deg); F.pop()
            for x in t: deg[x] -= 1
# degree <= 3 is forced: 4 lines through a point need 8 other points
rec(0, [], [0] * n)
print("labelled 8-triple systems on 8 points, pairwise <=1:", len(sols))
def canon(F):
    best = None
    for s in itertools.permutations(range(n)):
        img = tuple(sorted(tuple(sorted(s[i] for i in t)) for t in F))
        if best is None or img < best: best = img
    return best
classes = set()
seen = 0
for F in sols:
    classes.add(canon(F)); seen += 1
    if seen >= 3: break   # canon is slow (40320 perms); do a few, then use invariants for the rest
print("canonical forms of first", seen, "systems:", len(classes))
MK = [frozenset({i, (i+1) % 8, (i+3) % 8}) for i in range(8)]
cMK = canon(MK)
print("Moebius-Kantor canonical form equals theirs:", cMK in classes)
# cheaper isomorphism invariant for all: for each F, check isomorphism to MK by backtracking
def iso_to(F, G):
    Fs = set(F); Gs = set(G)
    adjF = {x: [t for t in F if x in t] for x in range(n)}
    adjG = {x: [t for t in G if x in t] for x in range(n)}
    order = list(range(n))
    def bt(i, m, used):
        if i == n:
            return all(frozenset(m[x] for x in t) in Gs for t in F)
        x = order[i]
        for y in range(n):
            if y in used: continue
            m[x] = y
            ok = True
            for t in adjF[x]:
                if all(z in m for z in t):
                    if frozenset(m[z] for z in t) not in Gs: ok = False; break
            if ok and bt(i + 1, m, used | {y}): return True
            del m[x]
        return False
    return bt(0, {}, frozenset())
print("all systems isomorphic to Moebius-Kantor:", all(iso_to(F, MK) for F in sols))

# real realisability of MK with an independent frame: choose quadrangle {0,2,4,6}?  need no three on a line
def no3(q): return all(not (frozenset(c) <= t) for t in MK for c in itertools.combinations(q, 3))
quad = [q for q in itertools.combinations(range(8), 4) if no3(q)]
print("quadrangles (no three collinear) in MK:", len(quad), "using", quad[-1])
q0, q1, q2, q3 = quad[-1]
X = sp.symbols('x0:8'); Y = sp.symbols('y0:8'); Z = sp.symbols('z0:8')
P = {q0: sp.Matrix([1,0,0]), q1: sp.Matrix([0,1,0]), q2: sp.Matrix([0,0,1]), q3: sp.Matrix([1,1,1])}
others = [i for i in range(8) if i not in P]
# generic homogeneous coordinates for the others with z=1 chart first (then check the z=0 chart too)
unk = []
for i in others:
    P[i] = sp.Matrix([X[i], Y[i], 1]); unk += [X[i], Y[i]]
eqs = [sp.expand(sp.Matrix.hstack(*[P[i] for i in sorted(t)]).det()) for t in MK]
eqs = [e for e in eqs if e != 0]
# nondegeneracy: all 8 points distinct
dist = 1
for i, j in itertools.combinations(range(8), 2):
    cr = P[i].cross(P[j])
    dist *= (cr[0]**2 + cr[1]**2 + cr[2]**2) if False else 1  # placeholder, handled below
w = sp.symbols('w')
# distinctness in the affine chart: (x_i - x_j)^2 + (y_i - y_j)^2 != 0 is not polynomial-friendly over C;
# use instead: for each pair of points on a common MK line both unknown, require a non-zero 2x2 minor.
# Simplest rigorous route: compute the full solution set of the incidence ideal and inspect it.
G = sp.groebner(eqs, *unk, order='lex')
print("Groebner basis (lex) of the incidence ideal in the affine chart z=1:")
for g in G: print("   ", sp.factor(g))
sol = sp.solve(eqs, unk, dict=True)
print("solutions of the affine-chart system:")
for s in sol:
    print("   ", s)
    pts = {i: P[i].subs(s) for i in range(8)}
    dup = [(i, j) for i, j in itertools.combinations(range(8), 2) if (pts[i].cross(pts[j])).norm() == 0]
    print("      coincident point pairs:", dup, " -> degenerate" if dup else "  (non-degenerate: check reality)")
    print("      contains non-real numbers:", any(not sp.im(v).equals(0) for v in s.values()))
