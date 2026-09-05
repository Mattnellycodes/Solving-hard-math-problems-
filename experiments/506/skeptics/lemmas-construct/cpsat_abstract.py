"""Independent CP-SAT model of the abstract Mobius block relaxation (written from scratch).
Variables: x_S (S subset, 4<=|S|<=n-1) = S is a block (circle or line with >=4 points);
           y_T (T subset, 3<=|T|<=n-1) = T is a line (block through infinity).
Constraints: every triple in at most one block of size>=4, and if the triple is itself a 3-line it
is in no bigger block; every pair in at most one line; a line of size>=4 is a block.
Optional: derived Sylvester-Gallai caps (--sg) at every point and at infinity.
Objective: maximise D + ell; circles = C(n,3) - (D + ell).
Modes: 'opt' (optimum), 'enum T' (enumerate all solutions with D+ell >= T, collect iso classes via
incidence-graph isomorphism).
"""
import sys, itertools, json
from math import comb
from ortools.sat.python import cp_model
import networkx as nx

OL = {2: 0, 3: 3, 4: 3, 5: 4, 6: 3, 7: 3, 8: 4, 9: 6, 10: 5}

def build(n, sg, extra_maxdeg4=None):
    m = cp_model.CpModel()
    pts = range(n)
    X = {}
    for k in range(4, n):
        for S in itertools.combinations(pts, k):
            X[S] = m.NewBoolVar('x' + ''.join(map(str, S)))
    Y = {}
    for k in range(3, n):
        for T in itertools.combinations(pts, k):
            Y[T] = m.NewBoolVar('y' + ''.join(map(str, T)))
    # per triple
    for t in itertools.combinations(pts, 3):
        big = [X[S] for S in X if set(t) <= set(S)]
        m.Add(sum(big) + Y[t] <= 1)
    # per pair: at most one line
    for p in itertools.combinations(pts, 2):
        m.Add(sum(Y[T] for T in Y if set(p) <= set(T)) <= 1)
    for T in Y:
        if len(T) >= 4:
            m.AddImplication(Y[T], X[T])
    if sg:
        for p in pts:
            m.Add(sum(X[S] * comb(len(S) - 1, 2) for S in X if p in S) <= comb(n - 1, 2) - OL[n - 1])
        m.Add(sum(Y[T] * comb(len(T), 2) for T in Y) <= comb(n, 2) - OL[n])
    obj = sum(X[S] * (comb(len(S), 3) - 1) for S in X) + sum(Y[T] for T in Y)
    return m, X, Y, obj

def incidence_graph(n, blocks, lines):
    G = nx.Graph()
    for p in range(n):
        G.add_node(('p', p), kind=0)
    for i, B in enumerate(blocks):
        G.add_node(('b', i), kind=1 + len(B))
        for p in B:
            G.add_edge(('b', i), ('p', p))
    G.add_node('inf', kind=99)
    for i, L in enumerate(lines):
        # a line is a block through infinity; if it's also a size>=4 block, connect that block node to inf
        idx = None
        for j, B in enumerate(blocks):
            if set(B) == set(L):
                idx = j
        if idx is None:
            G.add_node(('l', i), kind=1 + len(L))
            for p in L:
                G.add_edge(('l', i), ('p', p))
            G.add_edge(('l', i), 'inf')
        else:
            G.add_edge(('b', idx), 'inf')
    return G

class Collector(cp_model.CpSolverSolutionCallback):
    def __init__(self, n, X, Y, obj):
        super().__init__()
        self.n, self.X, self.Y, self.obj = n, X, Y, obj
        self.count = 0
        self.classes = {}  # hash -> list of (G, record)
    def on_solution_callback(self):
        self.count += 1
        blocks = [S for S in self.X if self.Value(self.X[S])]
        lines = [T for T in self.Y if self.Value(self.Y[T])]
        G = incidence_graph(self.n, blocks, lines)
        h = nx.weisfeiler_lehman_graph_hash(G, node_attr='kind', iterations=4)
        lst = self.classes.setdefault(h, [])
        nm = nx.algorithms.isomorphism.categorical_node_match('kind', None)
        for G2, rec in lst:
            if nx.is_isomorphic(G, G2, node_match=nm):
                rec['mult'] += 1
                return
        val = self.Value(self.obj)
        lst.append((G, dict(blocks=[list(S) for S in blocks], lines=[list(T) for T in lines],
                            D_plus_l=val, circles=comb(self.n, 3) - val, mult=1)))

if __name__ == '__main__':
    n = int(sys.argv[1]); mode = sys.argv[2]; sg = '--sg' in sys.argv
    workers = 2
    m, X, Y, obj = build(n, sg)
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = workers
    solver.parameters.max_time_in_seconds = float([a for a in sys.argv if a.startswith('--time=')][0][7:]) if any(a.startswith('--time=') for a in sys.argv) else 600
    if mode == 'opt':
        m.Maximize(obj)
        st = solver.Solve(m)
        print('n', n, 'sg', sg, 'status', solver.StatusName(st), 'obj', solver.ObjectiveValue(), 'bound', solver.BestObjectiveBound(),
              'circles', comb(n, 3) - solver.ObjectiveValue(), 'time', solver.WallTime())
        blocks = [S for S in X if solver.Value(X[S])]; lines = [T for T in Y if solver.Value(Y[T])]
        print('blocks', blocks); print('lines', lines)
    else:
        T = int(sys.argv[3])
        m.Add(obj >= T)
        solver.parameters.enumerate_all_solutions = True
        solver.parameters.num_search_workers = 1
        col = Collector(n, X, Y, obj)
        st = solver.Solve(m, col)
        print('n', n, 'sg', sg, 'threshold', T, 'status', solver.StatusName(st), 'labelled solutions', col.count, 'time', solver.WallTime())
        recs = [rec for lst in col.classes.values() for G, rec in lst]
        recs.sort(key=lambda r: -r['D_plus_l'])
        print('iso classes:', len(recs))
        for r in recs:
            print(' D+l', r['D_plus_l'], 'circles', r['circles'], 'mult', r['mult'], 'sizes', sorted(len(b) for b in r['blocks']), 'nlines', len(r['lines']), 'linesizes', sorted(len(l) for l in r['lines']))
        out = [a for a in sys.argv if a.startswith('--out=')]
        if out:
            json.dump(recs, open(out[0][6:], 'w'))
