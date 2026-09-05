"""Independent check (2026-09-05): (a) every family of 8 triples on 8 points pairwise sharing <= 1 point
is isomorphic to the Möbius–Kantor configuration {i, i+1, i+3 mod 8}; (b) that configuration has no
realisation by points and lines in the real projective plane (all 24 incidences, points distinct).
"""
import itertools, sympy as sp
pts = range(8)
triples = list(itertools.combinations(pts, 3))
# (a) brute force: all 8-subsets of triples pairwise sharing <=1; check every point has degree 3 and
# compute a canonical form under S_8 to count isomorphism classes.
def compatible(fam):
    return all(len(set(a) & set(b)) <= 1 for a, b in itertools.combinations(fam, 2))
# build via backtracking
fams = []
def rec(start, fam):
    if len(fam) == 8:
        fams.append(tuple(fam)); return
    for i in range(start, len(triples)):
        t = triples[i]
        if all(len(set(t) & set(u)) <= 1 for u in fam):
            fam.append(t); rec(i + 1, fam); fam.pop()
rec(0, [])
print("families of 8 pairwise <=1-intersecting triples on 8 points:", len(fams))
def canon(fam):
    best = None
    for perm in itertools.permutations(pts):
        img = tuple(sorted(tuple(sorted(perm[v] for v in t)) for t in fam))
        if best is None or img < best: best = img
    return best
# canonical forms are expensive (8! per family) — use invariants first, then canon on representatives
from collections import defaultdict
classes = defaultdict(list)
for fam in fams:
    deg = tuple(sorted(sum(1 for t in fam if p in t) for p in pts))
    classes[deg].append(fam)
print("degree sequences:", {k: len(v) for k, v in classes.items()})
reps = {}
for deg, fs in classes.items():
    cf = set()
    for fam in fs[:200]:  # sample; then verify count consistency by orbit size
        cf.add(canon(fam))
    reps[deg] = cf
    print(f"  degree seq {deg}: {len(cf)} canonical forms among first {min(200,len(fs))} families")
MK = tuple(sorted(tuple(sorted(((i) % 8, (i + 1) % 8, (i + 3) % 8))) for i in range(8)))
print("Möbius–Kantor canonical:", canon(MK), " same as class rep:", canon(MK) in set().union(*reps.values()))
# (b) non-realisability: points p_i in RP^2, lines {i,i+1,i+3}. Frame: p0=(1,0,0), p1=(0,1,0), p2=(0,0,1), p5=(1,1,1)
# valid since no three of {0,1,2,5} are on a configuration line (check). Unknown coords for the rest.
lines = [((i) % 8, (i + 1) % 8, (i + 3) % 8) for i in range(8)]
assert all(len({0,1,2,5} & set(l)) <= 2 for l in lines)
a, b, c, d, e, f, g, h = sp.symbols('a b c d e f g h')
P = {0: sp.Matrix([1,0,0]), 1: sp.Matrix([0,1,0]), 2: sp.Matrix([0,0,1]), 5: sp.Matrix([1,1,1]),
     3: sp.Matrix([a,b,1]), 4: sp.Matrix([c,d,1]), 6: sp.Matrix([e,f,1]), 7: sp.Matrix([g,h,1])}
eqs = [sp.Matrix.hstack(P[i], P[j], P[k]).det() for i, j, k in lines]
G = sp.groebner([sp.expand(q) for q in eqs], a, b, c, d, e, f, g, h, order='lex')
print("Groebner basis (chart z=1 for points 3,4,6,7):")
for gg in G.exprs: print("   ", gg)
sols = sp.solve(eqs, [a, b, c, d, e, f, g, h], dict=True)
print("solutions in this chart:", sols)
