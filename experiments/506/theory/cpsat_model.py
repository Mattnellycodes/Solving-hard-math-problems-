"""Independent formulation of the combinatorial relaxation as a CP-SAT model (OR-tools).

Variables: x_B for every subset B of size 4..n-1 (B is a block);  y_S for every subset S of size
>= 3 (S is a LINE, i.e. a block through infinity).  Constraints:
  (C1) x_B + x_B' <= 1 if |B ∩ B'| >= 3;  y_S + y_S' <= 1 if |S ∩ S'| >= 2;
       y_S <= x_S for |S| >= 4;  y_t + x_B <= 1 for every triple t ⊂ B (a triple is a block iff it is
       in no bigger block);
  (C2) for every point p: sum_{B ∋ p} C(|B|-1,2) x_B <= C(n-1,2) - o(n-1);
  (C3) sum_S C(|S|,2) y_S <= C(n,2) - o(n).
Objective: maximise  D + ell = sum_B (C(|B|,3)-1) x_B + sum_S y_S ;  circles = C(n,3) - (D + ell).
Usage: python3 cpsat_model.py n [time_limit_s] [--all target]  (--all enumerates all solutions with
D+ell >= C(n,3)-target, no symmetry reduction; use only when the count is small)
"""
import sys, itertools, time
from math import comb
from ortools.sat.python import cp_model
from mobius_enum import o_lower

def build(n, o_func=o_lower, max_block=None):
    max_block = max_block or n - 1
    m = cp_model.CpModel()
    subsets = {k: [frozenset(c) for c in itertools.combinations(range(n), k)] for k in range(3, max_block + 1)}
    x = {B: m.NewBoolVar(f"x{sorted(B)}") for k in range(4, max_block + 1) for B in subsets[k]}
    y = {S: m.NewBoolVar(f"y{sorted(S)}") for k in range(3, max_block + 1) for S in subsets[k]}
    blocks = list(x)
    # (C1) block conflicts
    for i in range(len(blocks)):
        for j in range(i + 1, len(blocks)):
            if len(blocks[i] & blocks[j]) >= 3:
                m.Add(x[blocks[i]] + x[blocks[j]] <= 1)
    lines = list(y)
    for i in range(len(lines)):
        for j in range(i + 1, len(lines)):
            if len(lines[i] & lines[j]) >= 2:
                m.Add(y[lines[i]] + y[lines[j]] <= 1)
    for S in lines:
        if len(S) >= 4:
            m.Add(y[S] <= x[S])
        else:
            for B in blocks:
                if S <= B:
                    m.Add(y[S] + x[B] <= 1)
    # (C2)
    capd = comb(n - 1, 2) - o_func(n - 1)
    for p in range(n):
        m.Add(sum(comb(len(B) - 1, 2) * x[B] for B in blocks if p in B) <= capd)
    # (C3)
    m.Add(sum(comb(len(S), 2) * y[S] for S in lines) <= comb(n, 2) - o_func(n))
    obj = sum((comb(len(B), 3) - 1) * x[B] for B in blocks) + sum(y[S] for S in lines)
    return m, x, y, obj

class Collector(cp_model.CpSolverSolutionCallback):
    def __init__(self, x, y, limit=200):
        super().__init__(); self.x = x; self.y = y; self.sols = []; self.limit = limit
    def on_solution_callback(self):
        F = sorted(sorted(B) for B, v in self.x.items() if self.Value(v))
        L = sorted(sorted(S) for S, v in self.y.items() if self.Value(v))
        self.sols.append((F, L))
        if len(self.sols) >= self.limit:
            self.StopSearch()

def main():
    n = int(sys.argv[1]); tl = float(sys.argv[2]) if len(sys.argv) > 2 and not sys.argv[2].startswith("--") else 600.0
    m, x, y, obj = build(n)
    N3 = comb(n, 3)
    if "--all" in sys.argv:
        target = int(sys.argv[sys.argv.index("--all") + 1])
        m.Add(obj >= N3 - target)
        solver = cp_model.CpSolver(); solver.parameters.max_time_in_seconds = tl
        solver.parameters.enumerate_all_solutions = True; solver.parameters.num_workers = 1
        col = Collector(x, y)
        st = solver.Solve(m, col)
        print(f"n={n}: status {solver.StatusName(st)}, {len(col.sols)} structures with circles <= {target} (labelled, no symmetry reduction)")
        for F, L in col.sols[:10]:
            print("  blocks:", F, " lines:", L)
        return
    m.Maximize(obj)
    solver = cp_model.CpSolver(); solver.parameters.max_time_in_seconds = tl; solver.parameters.num_workers = 2
    t0 = time.time(); st = solver.Solve(m)
    print(f"n={n}: status={solver.StatusName(st)} best D+ell={solver.ObjectiveValue():.0f} bound={solver.BestObjectiveBound():.0f} "
          f"=> combinatorial min circles = {N3 - solver.ObjectiveValue():.0f} (lower bound {N3 - solver.BestObjectiveBound():.0f}) in {time.time()-t0:.1f}s")
    F = sorted(sorted(B) for B, v in x.items() if solver.Value(v)); L = sorted(sorted(S) for S, v in y.items() if solver.Value(v))
    print("  optimal blocks:", F); print("  lines:", L)

if __name__ == "__main__":
    main()
