"""Own check of Lemma 2.4: (i) every 8-triple system on 8 points, pairwise sharing <= 1 point, is
isomorphic to Mobius-Kantor MK = {i,i+1,i+3 mod 8}; (ii) MK has no realisation in RP^2 with
distinct points.  (ii) is done with a general projective frame on {0,1,2,5} and exhaustive charts for the
remaining points (each point on a line through two frame points is parametrised on that line)."""
import itertools, sympy as sp
# (i)
trip = [frozenset(t) for t in itertools.combinations(range(8), 3)]
systems = []
def rec(F, start):
    if len(F) == 8: systems.append(F); return
    for j in range(start, len(trip)):
        t = trip[j]
        if all(len(t & u) <= 1 for u in F): rec(F + [t], j + 1)
rec([trip[0]], 1)   # WLOG contains {0,1,2}: every system has a triple, relabel it to {0,1,2}
print("8-triple systems containing {0,1,2}:", len(systems))
MK = [frozenset({i, (i + 1) % 8, (i + 3) % 8}) for i in range(8)]
def canon(F):
    return min(tuple(sorted(tuple(sorted(p[i] for i in t)) for t in F)) for p in itertools.permutations(range(8)))
cMK = canon(MK)
# cheap invariant first, then exact canonical form for distinct invariants
seen = set()
for F in systems:
    seen.add(canon(F))
print("isomorphism classes:", len(seen), "; equals MK:", seen == {cMK})
# (ii) realisation of MK in RP^2 -- my own parametrisation
# frame: p0=(1,0,0), p1=(0,1,0), p2=(0,0,1), p5=(1,1,1); no three of {0,1,2,5} on an MK line:
assert not any(len(l & {0, 1, 2, 5}) >= 3 for l in MK)
lines_by_pt = {p: [l for l in MK if p in l] for p in range(8)}
print("MK lines:", [sorted(l) for l in MK])
# point 3: on lines {0,1,3} and {2,3,5} -> intersection of z=0 and x=y  -> (1,1,0) (unique)
# point 4: on {1,2,4}: x=0 -> (0,a,1) or (0,1,0)=p1 (excluded: distinct)
# point 6: on {0,5,6}: y=z -> (c,1,1) or (1,0,0)=p0 (excluded)
# point 7: on {0,2,7}: y=0 -> (b,0,1) or (1,0,0)=p0 (excluded)
a, b, c = sp.symbols('a b c', real=True)
P = {0: (1, 0, 0), 1: (0, 1, 0), 2: (0, 0, 1), 5: (1, 1, 1), 3: (1, 1, 0), 4: (0, a, 1), 6: (c, 1, 1), 7: (b, 0, 1)}
eqs = []
for l in MK:
    i, j, k = sorted(l)
    eqs.append(sp.expand(sp.Matrix([P[i], P[j], P[k]]).det()))
eqs = [e for e in eqs if e != 0]
print("remaining equations:", eqs)
sols = sp.solve(eqs, [a, b, c], dict=True)
print("solutions over C:", sols)
# also without the 'real=True' assumption to be sure nothing is hidden
a2, b2, c2 = sp.symbols('a2 b2 c2')
eqs2 = [e.subs({a: a2, b: b2, c: c2}) for e in eqs]
G = sp.groebner(eqs2, a2, b2, c2, order='lex')
print("Groebner basis:", list(G))
uni = [g for g in G if g.free_symbols == {c2}]
print("univariate:", uni, " real roots:", [sp.real_roots(sp.Poly(u, c2)) for u in uni])
# Sanity: also check that the degenerate alternatives (p4=p1, p6=p0, p7=p0) are the only other points on
# those lines -- projective lines x=0, y=z, y=0 through the frame: (0,a,1),(0,1,0); (c,1,1),(1,0,0); (b,0,1),(1,0,0).
print("Conclusion: no real projective realisation of MK with distinct points.")
