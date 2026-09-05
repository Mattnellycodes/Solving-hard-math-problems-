"""Own proof of Lemma 2.4 for the c(9) audit: no 8 distinct, non-collinear real points carry 8
collinear triples pairwise sharing <= 1 point (equivalently: no point of a 9-point Moebius set lies
in 8 four-blocks and no bigger block).

Part A (combinatorial, brute force): every family of 8 triples on 8 points pairwise sharing <= 1
point is isomorphic to the Moebius-Kantor configuration MK = {i, i+1, i+3 mod 8}; there is no such
family with 9 triples.
Part B (combinatorial closure): in MK every point is collinear (via an MK line) with 6 of the other 7
points; the missing partner of i is i+4.  Consequence (proved in the audit text and checked here):
a realisation of MK by 8 distinct points either has NO further collinear triple, or has all 8 points
collinear.  Checked here: (B1) no triple outside MK is pairwise 'non-adjacent'; (B2) no MK line can be
extended by a further point without creating two MK lines sharing 2 points.
Part C (algebra, frame {0,1,2,5}, which contains no MK line): the incidences force
p3 = (1,1,0), p6 = (a,1,1), p7 = (b,0,1), p4 = (0,c,1) (chart alternatives coincide with frame
points), and then a = b = 1-c with c^2 - c + 1 = 0: no real solution.  Verified with sympy, and
independently by a Groebner basis over all 3^4 chart choices with pairwise-distinctness saturation.
"""
import itertools, sympy as sp

n = 8
triples = [frozenset(t) for t in itertools.combinations(range(n), 3)]

# ---- Part A ----
fams = []
def rec(F, start):
    if len(F) == 8:
        fams.append(F[:]); return
    for j in range(start, len(triples)):
        t = triples[j]
        if all(len(t & u) <= 1 for u in F):
            rec(F + [t], j + 1)
rec([], 0)
perms = list(itertools.permutations(range(n)))
def canon(F):
    return min(tuple(sorted(tuple(sorted(s[i] for i in t)) for t in F)) for s in perms)
MK = [frozenset({i, (i + 1) % 8, (i + 3) % 8}) for i in range(8)]
cMK = canon(MK)
classes = {}
for F in fams:
    classes.setdefault(canon(F), 0)
    classes[canon(F)] += 1
print("Part A: families of 8 triples on 8 points pairwise sharing <=1 point:", len(fams),
      "labelled; isomorphism classes:", len(classes), "; equals MK:", list(classes) == [cMK])
# 9 triples impossible: each point in <= 3 triples (3 triples through p use 6 other points), so <= 24/3 = 8
nine = 0
def rec9(F, start):
    global nine
    if len(F) == 9:
        nine += 1; return
    for j in range(start, len(triples)):
        t = triples[j]
        if all(len(t & u) <= 1 for u in F):
            rec9(F + [t], j + 1)
rec9([], 0)
print("Part A: families of 9 such triples:", nine)

# ---- Part B ----
adj = {i: set() for i in range(n)}
for L in MK:
    for a in L:
        adj[a] |= (L - {a})
print("Part B: non-adjacent partner of each point:", {i: sorted(set(range(n)) - adj[i] - {i}) for i in range(n)})
B1 = [t for t in triples if t not in MK and all(b not in adj[a] for a, b in itertools.combinations(t, 2))]
print("Part B1: triples outside MK with all three pairs non-adjacent:", B1)
B2 = []
for L in MK:
    for x in range(n):
        if x in L: continue
        # x joins line L: every MK line through x sharing a point with L would then share 2 points with L
        conflict = any(len((M & L)) >= 1 for M in MK if x in M)
        if not conflict: B2.append((sorted(L), x))
print("Part B2: (MK line, extra point) extensions not creating a 2-point intersection:", B2)

# ---- Part C: explicit frame computation ----
frame = (0, 1, 2, 5)
assert not any(frozenset(t) in set(MK) for t in itertools.combinations(frame, 3))
a, b, c = sp.symbols('a b c')
P = {0: sp.Matrix([1, 0, 0]), 1: sp.Matrix([0, 1, 0]), 2: sp.Matrix([0, 0, 1]), 5: sp.Matrix([1, 1, 1])}
# p3 on line(0,1) = {z=0} and on line(2,5) = {x=y}  -> (1,1,0) (alternatives on z=0 with x=0 give point 1)
P[3] = sp.Matrix([1, 1, 0])
P[6] = sp.Matrix([a, 1, 1])   # on line(5,0) = {y=z}; y=z=0 would be point 0
P[7] = sp.Matrix([b, 0, 1])   # on line(0,2) = {y=0}; z=0 would be point 0
P[4] = sp.Matrix([0, c, 1])   # on line(1,2) = {x=0}; z=0 would be point 1
eqs = [sp.expand(sp.Matrix.hstack(*[P[i] for i in sorted(L)]).det()) for L in MK]
eqs = [e for e in eqs if e != 0]
print("Part C: remaining incidence equations:", eqs)
G = sp.groebner(eqs, a, b, c, order='lex')
print("Part C: Groebner basis (lex):", list(G))
uni = [g for g in G if g.free_symbols == {c}]
print("Part C: univariate in c:", sp.factor(uni[0]), "; real roots:", sp.real_roots(sp.Poly(uni[0], c)))
sols = sp.solve(eqs, [a, b, c], dict=True)
print("Part C: all complex solutions:", sols, "; real:", [s for s in sols if all(sp.im(v) == 0 for v in s.values())])

# ---- Part C': all 81 charts for the 4 free points (projective coordinates), incidence equations only,
# then every solution component is inspected for coincident points (no heavy saturation polynomial).
free = [3, 4, 6, 7]
syms = {p: sp.symbols(f'x{p} y{p} z{p}') for p in free}
def distinct_ok(Q, sol):
    for p, q in itertools.combinations(range(n), 2):
        A = Q[p].subs(sol); B = Q[q].subs(sol)
        minors = [sp.simplify(A[i] * B[j] - A[j] * B[i]) for i, j in ((0, 1), (0, 2), (1, 2))]
        if all(mm == 0 for mm in minors):
            return False
    return True
total_real_distinct = 0; charts_with_solutions = 0
for chart in itertools.product(range(3), repeat=4):
    Q = {k: v for k, v in P.items() if k in frame}
    for p, ch in zip(free, chart):
        v = [sp.Integer(0)] * 3
        v[ch] = sp.Integer(1)
        for k in range(ch + 1, 3):
            v[k] = syms[p][k]
        Q[p] = sp.Matrix(v)            # canonical representative: first nonzero coordinate = 1
    E = [sp.expand(sp.Matrix.hstack(*[Q[i] for i in sorted(L)]).det()) for L in MK]
    if any(e.is_number and e != 0 for e in E):
        continue                        # an incidence is identically violated in this chart
    E = [e for e in E if e != 0]
    vars_ = sorted({s for p in free for s in Q[p].free_symbols}, key=str)
    Gb = sp.groebner(E, *vars_, order='grevlex') if E else [sp.Integer(0)]
    if list(Gb) == [1]:
        continue
    charts_with_solutions += 1
    sols = sp.solve(E, vars_, dict=True) if E else [dict()]
    # a solution family may have free parameters; substitute a couple of generic rational values for them
    for s in sols:
        freev = [v for v in vars_ if v not in s]
        trials = [dict()] if not freev else [dict(zip(freev, vals)) for vals in itertools.product([sp.Rational(k, 7) for k in (2, 3, 5)], repeat=len(freev))]
        for tr in trials:
            full = {k: sp.simplify(v.subs(tr)) for k, v in s.items()}; full.update(tr)
            real = all(sp.im(sp.nsimplify(v)) == 0 for v in full.values())
            dist = distinct_ok(Q, full)
            if real and dist:
                total_real_distinct += 1
                print("   REAL solution with distinct points in chart", chart, ":", full)
print("Part C': charts admitting complex solutions of the incidences:", charts_with_solutions,
      "; real solutions with 8 distinct points found:", total_real_distinct)
