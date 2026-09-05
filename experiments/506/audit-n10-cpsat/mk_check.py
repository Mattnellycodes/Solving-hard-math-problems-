"""Auditor's own check of the (8_3) lemma used by the MK constraint:
 (a) every family of 8 triples on 8 points pairwise sharing <= 1 point is isomorphic to the
     Moebius-Kantor configuration {i, i+1, i+3 mod 8};
 (b) the Moebius-Kantor configuration has no realisation by 8 distinct points of the real
     projective plane (hence none in R^2): with the frame p0=(1,0,0), p1=(0,1,0), p2=(0,0,1),
     p5=(1,1,1) (no three of them on a configuration line) the incidences force
     p3=(1,1,0), p4=(0,c,1), p6=(a,1,1), p7=(b,0,1) and a=b=1-c, c^2-c+1=0."""
import itertools, sympy as sp
import networkx as nx
from networkx.algorithms import isomorphism

MK = [frozenset(((i) % 8, (i + 1) % 8, (i + 3) % 8)) for i in range(8)]
assert all(len(a & b) <= 1 for a, b in itertools.combinations(MK, 2))

def graph(trs):
    G = nx.Graph()
    for p in range(8): G.add_node(('p', p), c=0)
    for i, t in enumerate(trs):
        G.add_node(('t', i), c=1)
        for p in t: G.add_edge(('t', i), ('p', p))
    return G
GMK = graph(MK)
nm = isomorphism.categorical_node_match('c', -1)
# (a): enumerate all 8-packings by DFS over triples in lexicographic order
trip = [frozenset(t) for t in itertools.combinations(range(8), 3)]
count = 0; nonmk = 0
def dfs(start, chosen):
    global count, nonmk
    if len(chosen) == 8:
        count += 1
        if not isomorphism.GraphMatcher(graph(chosen), GMK, node_match=nm).is_isomorphic():
            nonmk += 1
        return
    for i in range(start, len(trip)):
        t = trip[i]
        if all(len(t & c) <= 1 for c in chosen):
            chosen.append(t); dfs(i + 1, chosen); chosen.pop()
dfs(0, [])
print('(a) 8-packings of triples on 8 points:', count, '; not isomorphic to MK:', nonmk)
# (b)
a, b, c = sp.symbols('a b c')
P = {0: (1, 0, 0), 1: (0, 1, 0), 2: (0, 0, 1), 5: (1, 1, 1), 3: (1, 1, 0), 4: (0, c, 1), 6: (a, 1, 1), 7: (b, 0, 1)}
# check the forced coordinates: p3 on lines {0,1,3},{2,3,5}; p4 on {1,2,4}; p6 on {5,6,0}; p7 on {7,0,2}
def coll(x, y, z): return sp.Matrix([P[x], P[y], P[z]]).det()
frame_lines = [(0, 1, 3), (2, 3, 5), (1, 2, 4), (5, 6, 0), (7, 0, 2)]
assert all(sp.simplify(coll(*l)) == 0 for l in frame_lines)
# generic form check: a point on line {0,1} is (x,y,0); on line {2,5} (x=y): (1,1,0) up to scale (x=y=0 impossible)
# a point on line {1,2} is (0,y,z) with z != 0 (z=0 gives p1): (0,c,1); on {5,0}: (x,y,y), y != 0: (a,1,1); on {0,2}: (x,0,z), z != 0: (b,0,1)
rest = [(3, 4, 6), (4, 5, 7), (6, 7, 1)]
eqs = [sp.expand(coll(*l)) for l in rest]
print('(b) remaining incidence equations:', eqs)
G = sp.groebner(eqs, a, b, c, order='lex')
print('    Groebner basis:', list(G.exprs))
sols = sp.solve(eqs, [a, b, c], dict=True)
print('    solutions:', sols, ' real:', [all(sp.im(v) == 0 for v in s.values()) for s in sols])
