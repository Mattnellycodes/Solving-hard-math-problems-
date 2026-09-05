"""Step A.  Five points S on a conic, partition of the 10 chords into 5 perfect matchings M_v of S - v.
Skeleton-3 realisations (two-block lemma) need the five 'diagonal points' D_v = chord(M_v[0]) ∩ chord(M_v[1])
to be collinear (all on m = plane(S) ∩ plane(R)).  We solve this exactly: conic (1, t, t^2), S at t = 0, 1,
∞, a, b.  Output: all partitions, all real non-degenerate solutions (a, b), and whether each solution is
projectively the regular pentagon (the 5 points on the conic are then an orbit of a projectivity of order 5,
equivalently the pentagon's 'chord-direction' partition)."""
import itertools, sympy as sp
a, b = sp.symbols('a b')
pt = {0: sp.Matrix([1, 0, 0]), 1: sp.Matrix([1, 1, 1]), 2: sp.Matrix([0, 0, 1]), 3: sp.Matrix([1, a, a**2]), 4: sp.Matrix([1, b, b**2])}
V = range(5)
pairs = [frozenset(p) for p in itertools.combinations(V, 2)]
def pms(v):
    rest = [u for u in V if u != v]
    return [frozenset([frozenset([rest[0], rest[i]]), frozenset([u for u in rest[1:] if u != rest[i]])]) for i in (1, 2, 3)]
partitions = []
for choice in itertools.product(*[pms(v) for v in V]):
    used = [p for M in choice for p in M]
    if len(set(used)) == 10:
        partitions.append(choice)
print("partitions of E(K5) into perfect matchings of K5 - v:", len(partitions))
def line(p, q): return pt[p].cross(pt[q])
def dpoint(M):
    (c1, c2) = [sorted(x) for x in M]
    return line(*c1).cross(line(*c2))
nondeg = [a, b, a - 1, b - 1, a - b]
for k, part in enumerate(partitions):
    D = [dpoint(M) for M in part]
    Mx = sp.Matrix.hstack(*D).T   # 5 x 3
    eqs = [sp.expand(Mx.extract(list(rows), [0, 1, 2]).det()) for rows in itertools.combinations(range(5), 3)]
    eqs = [e for e in eqs if e != 0]
    G = sp.groebner(eqs, a, b, order='lex')
    sols = sp.solve(list(G), [a, b], dict=True)
    good = []
    for s in sols:
        if len(s) < 2:
            good.append(('positive-dimensional', s)); continue
        if all(sp.im(sp.N(s[x])) == 0 for x in (a, b)) and all(sp.N(nd.subs(s)) != 0 for nd in nondeg):
            good.append(s)
    print(f"partition {k}: M_v = {[sorted(sorted(x) for x in M) for M in part]}")
    print(f"   Groebner basis: {list(G)}")
    print(f"   real nondegenerate solutions: {good}")
    for s in good:
        if isinstance(s, dict):
            # cross-ratio test for the regular pentagon: for 5 points t0..t4 in cyclic order, the pentagon has
            # (t0,t1;t2,t3) = ... we test instead whether some projectivity of order 5 permutes the points cyclically
            ts = [0, 1, sp.oo, s[a], s[b]]
            found = False
            for perm in itertools.permutations(range(1, 5)):
                order = [0] + list(perm)
                # projectivity mapping t_order[0..2] -> t_order[1..3]; check it maps t_order[3] -> t_order[4] and t_order[4] -> t_order[0]
                def proj(x0, x1, x2, y0, y1, y2):
                    # find M with M(x_i) = y_i via cross-ratio: M(x) = y with (x0,x1;x2,x) = (y0,y1;y2,y)
                    x, y = sp.symbols('x y')
                    def cr(p, q, r, s_):
                        def sub(u, v):
                            return 1 if (u == sp.oo and v == sp.oo) else (1 if u == sp.oo else (-1 if v == sp.oo else u - v))
                        return sub(p, r) * sub(q, s_) / (sub(p, s_) * sub(q, r))
                    return sp.solve(sp.Eq(cr(x0, x1, x2, x), cr(y0, y1, y2, y)), y)
                t = [ts[i] for i in order]
                try:
                    img3 = proj(t[0], t[1], t[2], t[1], t[2], t[3])
                except Exception:
                    continue
                # simpler: check cross-ratio equality (t0,t1;t2,t3) == (t1,t2;t3,t4) == (t2,t3;t4,t0)
                def cr(p, q, r, s_):
                    def sub(u, v):
                        return 1 if (u == sp.oo and v == sp.oo) else (1 if u == sp.oo else (-1 if v == sp.oo else u - v))
                    return sp.nsimplify(sp.simplify(sub(p, r) * sub(q, s_) / (sub(p, s_) * sub(q, r))))
                c = [cr(t[i % 5], t[(i+1) % 5], t[(i+2) % 5], t[(i+3) % 5]) for i in range(5)]
                if all(sp.simplify(c[i] - c[0]) == 0 for i in range(5)):
                    found = True; break
            print(f"   solution {s}: cyclic-order-invariant (regular pentagon) = {found}; cross-ratio = {c[0] if found else None}")
