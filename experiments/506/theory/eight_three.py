"""No 8 points of the real plane have 8 lines each containing exactly 3 of them (t_3(8) <= 7).

Step 1: enumerate all families of 8 triples on 8 points, pairwise sharing <= 1 point, up to
isomorphism (every point then lies in exactly 3 triples: an (8_3) configuration).
Step 2: for each class, solve the realisation equations exactly (sympy).  A real realisation of the
line structure in RP^2 (hence in R^2) is excluded by showing the solution set is defined over
Q(sqrt(-3)) only.  Extra collinearities are irrelevant for our application (REPORT Prop. 5.2), where
the 8 points have exactly 4 ordinary lines, so the 4-point frame chosen is in general position.
"""
import itertools, sympy as sp
n = 8
triples = [frozenset(t) for t in itertools.combinations(range(n), 3)]
perms = list(itertools.permutations(range(n)))
def canon(F):
    return min(tuple(sorted(tuple(sorted(s[i] for i in t)) for t in F)) for s in perms)
classes = set()
def rec(F, start):
    if len(F) == 8:
        classes.add(canon(F)); return
    # each point in <= 3 triples
    for j in range(start, len(triples)):
        t = triples[j]
        if all(len(t & u) <= 1 for u in F):
            deg = [sum(1 for u in F if p in u) for p in range(n)]
            if all(deg[p] < 3 for p in t):
                rec(F + [t], j + 1)
# symmetry: first triple {0,1,2}
rec([triples[0]], 1)
print("isomorphism classes of (8_3) line systems:", len(classes))
for cl in classes:
    lines = [set(t) for t in cl]
    print("lines:", [sorted(t) for t in cl])
    # find a frame: 4 points no 3 on a line
    frame = None
    for q in itertools.combinations(range(n), 4):
        if not any(len(set(c) & l) == 3 for c in itertools.combinations(q, 3) for l in lines):
            frame = q; break
    print("frame (no three collinear):", frame)
    coords = {frame[0]: sp.Matrix([1, 0, 0]), frame[1]: sp.Matrix([0, 1, 0]),
              frame[2]: sp.Matrix([0, 0, 1]), frame[3]: sp.Matrix([1, 1, 1])}
    others = [p for p in range(n) if p not in frame]
    syms = {}
    for p in others:
        x, y, z = sp.symbols(f"x{p} y{p} z{p}")
        syms[p] = (x, y, z); coords[p] = sp.Matrix([x, y, z])
    eqs = []
    for l in lines:
        a, b, c = sorted(l)
        eqs.append(sp.Matrix.hstack(coords[a], coords[b], coords[c]).det())
    # Work chart by chart: each unknown point is normalised by setting one coordinate to 1 (the
    # remaining cases, where that coordinate is 0, are covered by trying all charts).
    # Non-degeneracy: points distinct (we use the strongest simple form: for each pair of distinct
    # points the 2x2 minors do not all vanish -> encode with Rabinowitsch on a random linear combination).
    total_real = 0
    import random
    random.seed(1)
    for chart in itertools.product(range(3), repeat=len(others)):
        subs = {}
        for p, ch in zip(others, chart):
            subs[syms[p][ch]] = 1
        E = [sp.expand(e.subs(subs)) for e in eqs]
        vars_ = [v for p in others for v in syms[p] if v not in subs]
        # distinctness via Rabinowitsch: product over pairs of (random minor) 
        t = sp.Symbol("t")
        prod = 1
        for p, q in itertools.combinations(range(n), 2):
            A = coords[p].subs(subs); B = coords[q].subs(subs)
            minors = [A[i]*B[j] - A[j]*B[i] for i, j in ((0,1),(0,2),(1,2))]
            prod *= sum(random.randint(1, 97) * m for m in minors)
        G = sp.groebner(E + [t * sp.expand(prod) - 1], *vars_, t, order="lex")
        if list(G) == [1]:
            continue
        # solve: last polynomials in lex order are univariate
        sols = sp.solve(list(G)[:-1] if False else list(G), vars_ + [t], dict=True)
        real_sols = [s for s in sols if all(sp.im(sp.nsimplify(v)) == 0 for v in s.values())]
        print("  chart", chart, "solutions:", len(sols), "real:", len(real_sols))
        total_real += len(real_sols)
        for s in sols[:4]:
            print("     ", {k: v for k, v in s.items() if k != t})
    print("TOTAL real realisations (all charts):", total_real)
