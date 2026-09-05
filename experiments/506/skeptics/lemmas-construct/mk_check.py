"""Independent check of Lemma 2.4: (a) every 8-triple system on 8 points with pairwise
intersections <= 1 is isomorphic to Mobius-Kantor; (b) MK has no real realisation (my own frame)."""
import itertools, networkx as nx, sympy as sp
n = 8
T = [frozenset(t) for t in itertools.combinations(range(n), 3)]
sols = []
def rec(F, start):
    if len(F) == 8:
        sols.append(list(F)); return
    for j in range(start, len(T)):
        t = T[j]
        if all(len(t & u) <= 1 for u in F):
            rec(F + [t], j + 1)
rec([], 0)
print('labelled 8-triple systems (pairwise <=1):', len(sols))
def G_of(F):
    G = nx.Graph()
    for p in range(n): G.add_node(('p', p), k=0)
    for i, t in enumerate(F):
        G.add_node(('t', i), k=1)
        for p in t: G.add_edge(('t', i), ('p', p))
    return G
reps = []
nm = nx.algorithms.isomorphism.categorical_node_match('k', None)
for F in sols:
    G = G_of(F)
    if not any(nx.is_isomorphic(G, H, node_match=nm) for H in reps):
        reps.append(G)
print('isomorphism classes:', len(reps))
MK = [frozenset({i, (i+1) % 8, (i+3) % 8}) for i in range(8)]
print('MK is one of them:', any(nx.is_isomorphic(G_of(MK), H, node_match=nm) for H in reps))
# (b) realisability over R: points p_i = (x_i, y_i, 1) projective; MK lines: 8 determinant equations.
# frame: choose 4 points no three on an MK line: check {0,2,4,6}? lines containing pairs...
def ok_frame(fr): return not any(len(L & set(fr)) >= 3 for L in MK)
frame = next(fr for fr in itertools.combinations(range(8), 4) if ok_frame(fr))
print('frame (no three on a line):', frame)
xs = sp.symbols('x0:8'); ys = sp.symbols('y0:8')
P = {}
fixed = {frame[0]: (0, 0, 1), frame[1]: (1, 0, 1), frame[2]: (0, 1, 1), frame[3]: (1, 1, 1)}
for i in range(8):
    P[i] = sp.Matrix(fixed[i]) if i in fixed else sp.Matrix([xs[i], ys[i], 1])
eqs = [sp.expand(sp.Matrix.hstack(*[P[i] for i in sorted(L)]).det()) for L in MK]
vars_ = [v for i in range(8) if i not in fixed for v in (xs[i], ys[i])]
G = sp.groebner(eqs, *vars_, order='lex')
print('Groebner basis (lex):'); [print('  ', sp.factor(g)) for g in G]
uni = [g for g in G if len(g.free_symbols) == 1]
for g in uni:
    print('univariate:', sp.factor(g), ' real roots:', sp.Poly(g, list(g.free_symbols)[0]).real_roots())
