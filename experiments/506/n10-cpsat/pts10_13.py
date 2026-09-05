"""All maximum partial triple systems on 10 points (13 triples pairwise sharing <= 1 point) up to
isomorphism, and whether each contains a Fano plane (7 points, 7 triples) or an (8_3) (8 points,
8 triples).  Counting: 13 triples cover 39 pairs, leave 6 edges; every point has odd leave-degree,
so nine points have leave-degree 1 and one point has leave-degree 3: the leave is K_{1,3} plus a
perfect matching of the remaining 6 points.  WLOG leave = {01,02,03} ∪ {45,67,89}."""
import itertools, json
from ortools.sat.python import cp_model
import networkx as nx
P = range(10)
T = [frozenset(t) for t in itertools.combinations(P, 3)]
leave = {frozenset(e) for e in [(0,1),(0,2),(0,3),(4,5),(6,7),(8,9)]}
m = cp_model.CpModel(); x = {t: m.NewBoolVar('') for t in T}
for e in itertools.combinations(P, 2):
    e = frozenset(e); s = sum(x[t] for t in T if e <= t)
    m.Add(s == (0 if e in leave else 1))
sols = []
class CB(cp_model.CpSolverSolutionCallback):
    def on_solution_callback(s): sols.append([t for t in T if s.Value(x[t])])
sv = cp_model.CpSolver(); sv.parameters.enumerate_all_solutions = True; sv.parameters.num_workers = 1
st = sv.Solve(m, CB()); print("status", sv.StatusName(st), "labelled systems with this leave:", len(sols))
def G(F):
    g = nx.Graph()
    for p in P: g.add_node(('p', p), c=0)
    for i, t in enumerate(F):
        g.add_node(('t', i), c=1)
        for p in t: g.add_edge(('t', i), ('p', p))
    return g
nm = nx.algorithms.isomorphism.categorical_node_match('c', 0)
reps = []
for F in sols:
    g = G(F)
    if not any(nx.is_isomorphic(g, r[1], node_match=nm) for r in reps): reps.append((F, g))
print("isomorphism classes:", len(reps))
for F, g in reps:
    fano = any(sum(1 for t in F if t <= frozenset(S)) >= 7 for S in itertools.combinations(P, 7))
    mk = any(sum(1 for t in F if t <= frozenset(S)) >= 8 for S in itertools.combinations(P, 8))
    print(" system", sorted(sorted(t) for t in F), "contains Fano:", fano, " contains (8_3):", mk)
