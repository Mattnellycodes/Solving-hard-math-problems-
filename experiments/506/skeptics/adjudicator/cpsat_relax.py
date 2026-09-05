"""Adjudicator's own CP-SAT model of the combinatorial relaxation (written from scratch).

Variables: x_B for every subset B of size 4..n-1 (block of size >= 4), y_S for every subset S of
size 3..n-1 (line with >= 3 points).  Constraints:
  (C1) x_B + x_B' <= 1 if |B ∩ B'| >= 3
  (C1') y_S <= x_S for |S| >= 4 (a rich line is a block);  y_t + x_B <= 1 for |t| = 3, t ⊂ B
  (C3) y_S + y_S' <= 1 if |S ∩ S'| >= 2
  caps (regime sg, Sylvester-Gallai with o = 1): for every p, sum_{B ∋ p} C(|B|-1,2) x_B <= C(n-1,2)-1;
       sum_S C(|S|,2) y_S <= C(n,2)-1.   regime none: no caps.
Objective: maximise D + ell = sum (C(|B|,3)-1) x_B + sum y_S ; min circles = C(n,3) - optimum.
Usage: python3 cpsat_relax.py n regime [time_limit_s] [forced block as comma list]
"""
import sys, time
from itertools import combinations
from math import comb
from ortools.sat.python import cp_model

n = int(sys.argv[1]); regime = sys.argv[2]
tl = float(sys.argv[3]) if len(sys.argv) > 3 else 600.0
forced = [int(x) for x in sys.argv[4].split(',')] if len(sys.argv) > 4 else None
maxsize = int(sys.argv[5]) if len(sys.argv) > 5 else n - 1
pts = range(n)
blocks = [frozenset(c) for k in range(4, maxsize + 1) for c in combinations(pts, k)]
lsets = [frozenset(c) for k in range(3, maxsize + 1) for c in combinations(pts, k)]
m = cp_model.CpModel()
x = {B: m.NewBoolVar(f"x{sorted(B)}") for B in blocks}
y = {S: m.NewBoolVar(f"y{sorted(S)}") for S in lsets}
for A, B in combinations(blocks, 2):
    if len(A & B) >= 3:
        m.AddBoolOr([x[A].Not(), x[B].Not()])
for S in lsets:
    if len(S) >= 4:
        m.AddImplication(y[S], x[S])
    else:
        for B in blocks:
            if S <= B:
                m.AddBoolOr([y[S].Not(), x[B].Not()])
for S, T in combinations(lsets, 2):
    if len(S & T) >= 2:
        m.AddBoolOr([y[S].Not(), y[T].Not()])
if regime == 'sg':
    for p in pts:
        m.Add(sum(comb(len(B) - 1, 2) * x[B] for B in blocks if p in B) <= comb(n - 1, 2) - 1)
    m.Add(sum(comb(len(S), 2) * y[S] for S in lsets) <= comb(n, 2) - 1)
if forced is not None:
    m.Add(x[frozenset(forced)] == 1)
obj = sum((comb(len(B), 3) - 1) * x[B] for B in blocks) + sum(y[S] for S in lsets)
m.Maximize(obj)
solver = cp_model.CpSolver()
solver.parameters.num_workers = 2
solver.parameters.max_time_in_seconds = tl
t0 = time.time()
st = solver.Solve(m)
name = solver.StatusName(st)
val = solver.ObjectiveValue(); bd = solver.BestObjectiveBound()
print(f"n={n} regime={regime} forced={forced} maxsize={maxsize}: status={name} max D+ell={val:.0f} bound={bd:.0f} "
      f"=> min circles={comb(n,3)-val:.0f} (proved lower bound {comb(n,3)-bd:.0f}) time={time.time()-t0:.1f}s")
if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    F = [sorted(B) for B in blocks if solver.Value(x[B])]
    L = [sorted(S) for S in lsets if solver.Value(y[S])]
    print("   blocks:", F)
    print("   lines:", L)
    print("   degrees:", [sum(1 for B in F if p in B) for p in pts])
