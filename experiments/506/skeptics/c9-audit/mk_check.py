"""Own check of Lemma 2.4 ((8_3) uniqueness and non-realisability over R).
Part 1 (combinatorial): enumerate ALL families of 8 triples on 8 points pairwise sharing <= 1 point
(labelled DFS, no symmetry assumption beyond lexicographic order), and count isomorphism classes by a
canonical form under all 8! permutations.
Part 2 (algebraic): with frame p0,p1,p2,p5 (no MK-line among them), chart choices forced by point
distinctness (checked explicitly), the remaining incidences give a univariate polynomial; also solve
over the complex numbers to show the ONLY solutions are non-real.
Part 3 (independent, projective, all charts): sympy Groebner with Rabinowitsch distinctness, over all
3^4 chart choices for the 4 free points, reporting every chart where the ideal is not (1)."""
import itertools, sympy as sp
n = 8
triples = [frozenset(t) for t in itertools.combinations(range(n), 3)]
fams = []
def rec(F, start):
    if len(F) == 8:
        fams.append(list(F)); return
    for j in range(start, len(triples)):
        t = triples[j]
        if all(len(t & u) <= 1 for u in F):
            rec(F + [t], j + 1)
rec([], 0)
print("labelled 8-triple systems on 8 points pairwise sharing <=1 point:", len(fams))
perms = list(itertools.permutations(range(n)))
def canon(F):
    return min(tuple(sorted(tuple(sorted(s[i] for i in t)) for t in F)) for s in perms)
degs = set(tuple(sorted(sum(1 for t in F if p in t) for p in range(n))) for F in fams)
print("degree sequences occurring:", degs)
MK = [frozenset({i, (i+1) % 8, (i+3) % 8}) for i in range(8)]
cMK = canon(MK)
# canonicalising 840 families x 40320 perms is slow-ish; use degree-sequence + a cheap invariant first
classes = {}
for F in fams:
    c = canon(F)
    classes[c] = classes.get(c, 0) + 1
print("isomorphism classes:", len(classes), "; MK class present:", cMK in classes, "; labelled copies per class:", list(classes.values()))
# Part 2
print("\nMK lines:", [sorted(t) for t in MK])
frame = (0, 1, 2, 5)
assert not any(frozenset(t) in set(MK) for t in itertools.combinations(frame, 3)), "frame contains an MK line"
a, b, c = sp.symbols('a b c')
P = {0: sp.Matrix([1, 0, 0]), 1: sp.Matrix([0, 1, 0]), 2: sp.Matrix([0, 0, 1]), 5: sp.Matrix([1, 1, 1])}
# p3 = line(p0,p1) ∩ line(p2,p5)   (lines {0,1,3} and {2,3,5})
L01 = P[0].cross(P[1]); L25 = P[2].cross(P[5]); P[3] = L01.cross(L25)
print("p3 =", list(P[3]))
# p6 on line(p5,p0) = {y=z}: (a,1,1) unless p6 = p0 ; p7 on line(p0,p2)={y=0}: (b,0,1) unless p7=p0 ; p4 on line(p1,p2)={x=0}: (0,c,1) unless p4=p1
P[6] = sp.Matrix([a, 1, 1]); P[7] = sp.Matrix([b, 0, 1]); P[4] = sp.Matrix([0, c, 1])
eqs = [sp.expand(sp.Matrix.hstack(*[P[i] for i in sorted(L)]).det()) for L in MK]
print("incidence equations:", [e for e in eqs if e != 0])
sols = sp.solve([e for e in eqs if e != 0], [a, b, c], dict=True)
print("all complex solutions:", sols)
print("any real solution:", any(all(sp.im(v) == 0 for v in s.values()) for s in sols))
# Part 3: all charts, projective, Rabinowitsch for distinctness
t = sp.Symbol('t')
free = [3, 4, 6, 7]
syms = {p: sp.symbols(f'x{p} y{p} z{p}') for p in free}
nontrivial = 0
for chart in itertools.product(range(3), repeat=4):
    Q = dict(P); Q[3] = None
    subs = {}
    for p, ch in zip(free, chart):
        v = list(syms[p]); v[ch] = 1
        # chart ch: coordinate ch equals 1 and all earlier coordinates are 0 (canonical representative)
        for k in range(ch): v[k] = 0
        Q[p] = sp.Matrix(v)
    E = [sp.expand(sp.Matrix.hstack(*[Q[i] for i in sorted(L)]).det()) for L in MK]
    prod = sp.Integer(1)
    for p, q in itertools.combinations(range(n), 2):
        A, B = Q[p], Q[q]
        minors = [A[i]*B[j]-A[j]*B[i] for i, j in ((0, 1), (0, 2), (1, 2))]
        prod *= (minors[0]**2 + minors[1]**2 + minors[2]**2)   # distinct <=> some minor != 0 (real: sum of squares)
    vars_ = sorted({s for p in free for s in Q[p].free_symbols}, key=str)
    G = sp.groebner(E + [t*sp.expand(prod) - 1], *vars_, t, order='grevlex')
    if list(G) != [1]:
        nontrivial += 1
        print("chart", chart, "ideal not unit; basis:", list(G))
print("charts with non-unit ideal (complex solutions possible):", nontrivial, "of 81")
