"""Independent CP-SAT model of the combinatorial relaxation for Erdos #506, n = 8.
Formulated differently from the theory agent's model: triple-assignment form.
  x_B  (|B| in 4..n-1): B is a block (circle or line with >= 4 points)
  z_t  (|t| = 3):       t is a 3-block (its circle/line contains no 4th point)
  y_S  (|S| in 3..n-1): S is a line (a block through infinity)
  every triple lies in exactly one block:   z_t + sum_{B ⊃ t} x_B = 1
  blocks pairwise share <= 2 points:        x_B + x_B' <= 1 if |B∩B'| >= 3
  lines are blocks:                         y_S <= x_S (|S|>=4),  y_t <= z_t (|t|=3)
  lines pairwise share <= 1 point:          y_S + y_S' <= 1 if |S∩S'| >= 2
  derived Sylvester-Gallai at p:            sum_{B∋p} C(|B|-1,2) x_B <= C(n-1,2) - o(n-1)
  Sylvester-Gallai at infinity:             sum_S C(|S|,2) y_S <= C(n,2) - o(n)
  circles = (#blocks) - (#lines) = sum x_B + sum z_t - sum y_S   (minimised).
Usage: python3 cpsat_audit.py MODE [--enum TARGET [LIMIT]]   MODE in strong|weak|none
"""
import sys, itertools, time
from math import comb
from ortools.sat.python import cp_model

N = 8
OVALS = {"strong": {7: 3, 8: 4},   # Kelly-Moser 3m/7: o(7)>=3, o(8)>=4
         "weak":   {7: 1, 8: 1},   # Sylvester-Gallai only
         "none":   {7: 0, 8: 0}}   # no SG at all (for comparison only; NOT a valid relaxation)

def build(mode):
    o = OVALS[mode]
    capd, capl = comb(N - 1, 2) - o[N - 1], comb(N, 2) - o[N]
    m = cp_model.CpModel()
    subs = lambda k: [frozenset(c) for c in itertools.combinations(range(N), k)]
    big = [B for k in range(4, N) for B in subs(k)]
    tri = subs(3)
    x = {B: m.NewBoolVar("x" + "".join(map(str, sorted(B)))) for B in big}
    z = {t: m.NewBoolVar("z" + "".join(map(str, sorted(t)))) for t in tri}
    y = {S: m.NewBoolVar("y" + "".join(map(str, sorted(S)))) for S in big + tri}
    for t in tri:
        m.Add(z[t] + sum(x[B] for B in big if t <= B) == 1)
    for B1, B2 in itertools.combinations(big, 2):
        if len(B1 & B2) >= 3:
            m.AddBoolOr([x[B1].Not(), x[B2].Not()])
    for S in big:
        m.AddImplication(y[S], x[S])
    for t in tri:
        m.AddImplication(y[t], z[t])
    for S1, S2 in itertools.combinations(big + tri, 2):
        if len(S1 & S2) >= 2:
            m.AddBoolOr([y[S1].Not(), y[S2].Not()])
    for p in range(N):
        m.Add(sum(comb(len(B) - 1, 2) * x[B] for B in big if p in B) <= capd)
    m.Add(sum(comb(len(S), 2) * y[S] for S in big + tri) <= capl)
    circles = sum(x.values()) + sum(z.values()) - sum(y.values())
    return m, x, z, y, circles, capd, capl

class Collector(cp_model.CpSolverSolutionCallback):
    def __init__(self, x, y, limit):
        super().__init__(); self.x, self.y, self.limit = x, y, limit; self.sols = []
    def on_solution_callback(self):
        F = tuple(sorted(sum(1 << i for i in B) for B, v in self.x.items() if self.Value(v)))
        L = tuple(sorted(sum(1 << i for i in S) for S, v in self.y.items() if self.Value(v)))
        self.sols.append((F, L))
        if len(self.sols) >= self.limit:
            self.StopSearch()

def main():
    mode = sys.argv[1]
    m, x, z, y, circles, capd, capl = build(mode)
    print(f"mode={mode}: cap_derived={capd} cap_lines={capl}", flush=True)
    if "--enum" in sys.argv:
        i = sys.argv.index("--enum"); target = int(sys.argv[i + 1])
        limit = int(sys.argv[i + 2]) if len(sys.argv) > i + 2 else 100000
        m.Add(circles <= target)
        s = cp_model.CpSolver(); s.parameters.num_workers = 1
        s.parameters.enumerate_all_solutions = True; s.parameters.max_time_in_seconds = 1500
        col = Collector(x, y, limit)
        t0 = time.time(); st = s.Solve(m, col)
        print(f"enumeration circles<={target}: status={s.StatusName(st)} labelled solutions={len(col.sols)} "
              f"time={time.time()-t0:.1f}s", flush=True)
        import pickle
        with open(f"cpsat_{mode}_enum{target}.pkl", "wb") as f:
            pickle.dump(col.sols, f)
        return
    m.Minimize(circles)
    s = cp_model.CpSolver(); s.parameters.num_workers = 2; s.parameters.max_time_in_seconds = 1500
    t0 = time.time(); st = s.Solve(m)
    print(f"status={s.StatusName(st)} min circles={s.ObjectiveValue():.0f} bound={s.BestObjectiveBound():.0f} "
          f"time={time.time()-t0:.1f}s")
    F = sorted(sorted(B) for B, v in x.items() if s.Value(v)); L = sorted(sorted(S) for S, v in y.items() if s.Value(v))
    print("  blocks:", F); print("  lines:", L)
    print("  D =", sum(comb(len(B), 3) - 1 for B in F), " ell =", len(L))

if __name__ == "__main__":
    main()
