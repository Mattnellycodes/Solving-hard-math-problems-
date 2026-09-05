"""lemmas-audit: own CP-SAT model of the combinatorial relaxation (independent second check of the
optimum D + ell for n = 6, 7, 8 under weak (SG-only) and strong caps).
usage: python3 cpsat_own.py n capmode [time_limit] [maxsize]
Model: x_B (blocks, sizes 4..maxsize), y_S (lines, sizes 3..maxsize).
  blocks pairwise share <= 2 points;  lines pairwise share <= 1 point;
  y_S -> x_S (|S| >= 4);  y_t -> no block contains t (|t| = 3);
  derived cap at every point; line cap.  maximise D + ell.
"""
import sys, itertools, time
from math import comb
from ortools.sat.python import cp_model

def olower(m, capmode):
    if m <= 2: return 0
    if capmode == 0: return 0
    if capmode == 1: return 1
    tab = {3: 3, 4: 3, 5: 4, 6: 3, 7: 3, 8: 4, 9: 6, 10: 5}
    return tab[m]

def main():
    n = int(sys.argv[1]); capmode = int(sys.argv[2])
    tl = float(sys.argv[3]) if len(sys.argv) > 3 else 600
    maxsize = int(sys.argv[4]) if len(sys.argv) > 4 else n - 1
    m = cp_model.CpModel()
    subsets = {k: [frozenset(c) for c in itertools.combinations(range(n), k)] for k in range(3, maxsize + 1)}
    x = {B: m.NewBoolVar("x" + "".join(map(str, sorted(B)))) for k in range(4, maxsize + 1) for B in subsets[k]}
    y = {S: m.NewBoolVar("y" + "".join(map(str, sorted(S)))) for k in range(3, maxsize + 1) for S in subsets[k]}
    xs = list(x); ys = list(y)
    for A, B in itertools.combinations(xs, 2):
        if len(A & B) >= 3:
            m.AddBoolOr([x[A].Not(), x[B].Not()])
    for A, B in itertools.combinations(ys, 2):
        if len(A & B) >= 2:
            m.AddBoolOr([y[A].Not(), y[B].Not()])
    for S in ys:
        if len(S) >= 4:
            m.AddImplication(y[S], x[S])
        else:
            for B in xs:
                if S <= B:
                    m.AddBoolOr([y[S].Not(), x[B].Not()])
    capd = comb(n - 1, 2) - olower(n - 1, capmode)
    for p in range(n):
        m.Add(sum(comb(len(B) - 1, 2) * x[B] for B in xs if p in B) <= capd)
    capl = comb(n, 2) - olower(n, capmode)
    m.Add(sum(comb(len(S), 2) * y[S] for S in ys) <= capl)
    obj = sum((comb(len(B), 3) - 1) * x[B] for B in xs) + sum(y[S] for S in ys)
    m.Maximize(obj)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = tl
    solver.parameters.num_workers = 2
    t0 = time.time()
    st = solver.Solve(m)
    N3 = comb(n, 3)
    print(f"n={n} capmode={capmode} maxsize={maxsize}: status={solver.StatusName(st)} best D+ell={solver.ObjectiveValue():.0f} "
          f"bound={solver.BestObjectiveBound():.0f} => min circles (relaxation) = {N3 - solver.ObjectiveValue():.0f}, "
          f"proven lower bound {N3 - solver.BestObjectiveBound():.0f}; time {time.time() - t0:.1f}s", flush=True)
    F = [sorted(B) for B in xs if solver.Value(x[B])]; L = [sorted(S) for S in ys if solver.Value(y[S])]
    print("  blocks:", F)
    print("  lines:", L)

if __name__ == "__main__":
    main()
