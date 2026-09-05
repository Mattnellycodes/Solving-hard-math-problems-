"""Independent CP-SAT model (written from scratch, 2026-09-05) for the combinatorial relaxation
used in the lower-bound arguments for Erdős #506.

Möbius setting: P has n points; blocks = circles-or-lines through >= 3 points. Every triple lies in
exactly one block. Two blocks share <= 2 points. Lines (blocks through infinity) pairwise share <= 1.
circles(P) = C(n,3) - D - l,  D = sum_{blocks of size k>=4} (C(k,3) - 1),  l = number of lines of size >= 3.
Necessary conditions from Sylvester-Gallai applied to the derived structure at each point p
(inversion about p turns blocks through p into lines of an (n-1)-point non-collinear set, which has
at least o(n-1) ordinary lines): rich blocks through p cover <= C(n-1,2) - o(n-1) pairs of P\{p}.
Lines cover <= C(n,2) - o(n) pairs.  o(m) = min number of ordinary lines of m non-collinear points:
Sylvester-Gallai gives o(m) >= 1; the known exact table gives o(3..14) = 3,3,4,3,3,4,6,5,6,6,6,7.

Usage: python3 cpsat_relaxation.py n [sg|table] [threshold]
  Maximises D + l and prints the lower bound circles >= C(n,3) - opt.  With a threshold T it instead
  enumerates all labelled structures with D + l >= T and prints isomorphism-invariant summaries.
"""
import sys, itertools, math
from collections import Counter
from ortools.sat.python import cp_model

O_TABLE = {3: 3, 4: 3, 5: 4, 6: 3, 7: 3, 8: 4, 9: 6, 10: 5, 11: 6, 12: 6, 13: 6, 14: 7}

def build(n, mode):
    o = (lambda m: 1) if mode == 'sg' else (lambda m: O_TABLE[m])
    pts = range(n)
    rich = [frozenset(c) for k in range(4, n) for c in itertools.combinations(pts, k)]
    lines = [frozenset(c) for k in range(3, n) for c in itertools.combinations(pts, k)]
    m = cp_model.CpModel()
    x = {B: m.NewBoolVar(f"x{sorted(B)}") for B in rich}
    y = {S: m.NewBoolVar(f"y{sorted(S)}") for S in lines}
    richset = set(rich)
    # two rich blocks share <= 2 points
    for B1, B2 in itertools.combinations(rich, 2):
        if len(B1 & B2) >= 3:
            m.AddBoolOr([x[B1].Not(), x[B2].Not()])
    # lines: size>=4 line must be a rich block; size-3 line must be an uncovered triple
    for S in lines:
        if len(S) >= 4:
            m.AddImplication(y[S], x[S])
        else:
            for B in rich:
                if S <= B:
                    m.AddBoolOr([y[S].Not(), x[B].Not()])
    # two lines share <= 1 point
    for S1, S2 in itertools.combinations(lines, 2):
        if len(S1 & S2) >= 2:
            m.AddBoolOr([y[S1].Not(), y[S2].Not()])
    # derived Sylvester-Gallai caps
    for p in pts:
        m.Add(sum(math.comb(len(B) - 1, 2) * x[B] for B in rich if p in B) <= math.comb(n - 1, 2) - o(n - 1))
    # line cap
    m.Add(sum(math.comb(len(S), 2) * y[S] for S in lines) <= math.comb(n, 2) - o(n))
    obj = sum((math.comb(len(B), 3) - 1) * x[B] for B in rich) + sum(y[S] for S in lines)
    return m, x, y, obj

class Collector(cp_model.CpSolverSolutionCallback):
    def __init__(self, x, y, n):
        super().__init__(); self.x, self.y, self.n = x, y, n; self.count = 0; self.inv = Counter(); self.examples = {}
    def on_solution_callback(self):
        self.count += 1
        F = [B for B, v in self.x.items() if self.Value(v)]
        L = [S for S, v in self.y.items() if self.Value(v)]
        deg = Counter(); 
        for B in F:
            for p in B: deg[p] += 1
        inv = (tuple(sorted(len(B) for B in F)), tuple(sorted(deg[p] for p in range(self.n))), tuple(sorted(len(S) for S in L)))
        self.inv[inv] += 1
        self.examples.setdefault(inv, (F, L))

if __name__ == "__main__":
    n = int(sys.argv[1]); mode = sys.argv[2] if len(sys.argv) > 2 else 'sg'
    thr = int(sys.argv[3]) if len(sys.argv) > 3 else None
    workers = int(sys.argv[4]) if len(sys.argv) > 4 else 2
    m, x, y, obj = build(n, mode)
    solver = cp_model.CpSolver(); solver.parameters.num_search_workers = workers
    solver.parameters.max_time_in_seconds = 3600
    if thr is None:
        m.Maximize(obj)
        st = solver.Solve(m)
        print(f"n={n} mode={mode} status={solver.StatusName(st)} opt(D+l)={solver.ObjectiveValue():.0f} bound={solver.BestObjectiveBound():.0f}  =>  circles >= {math.comb(n,3) - round(solver.ObjectiveValue())}  (formula {math.comb(n-1,2)+1-(n-1)//2})")
        F = [sorted(B) for B, v in x.items() if solver.Value(v)]; L = [sorted(S) for S, v in y.items() if solver.Value(v)]
        print("  an optimal structure: blocks", F, " lines", L)
    else:
        m.Add(obj >= thr)
        solver.parameters.enumerate_all_solutions = True
        solver.parameters.num_search_workers = 1
        col = Collector(x, y, n)
        st = solver.Solve(m, col)
        print(f"n={n} mode={mode} D+l>={thr}: status={solver.StatusName(st)} labelled solutions={col.count}")
        for inv, c in sorted(col.inv.items(), key=lambda t: -t[1]):
            F, L = col.examples[inv]
            print(f"  invariant blocksizes={inv[0]} degrees={inv[1]} linesizes={inv[2]} : {c} labelled copies; count = {math.comb(n,3) - sum(math.comb(len(B),3)-1 for B in F) - len(L)}")
            print("     example blocks:", [sorted(B) for B in F], "lines:", [sorted(S) for S in L])
