"""Skeptic c7-audit: own CP-SAT model of the combinatorial relaxation for n = 7.
Maximise D + ell under (C1) and a cap regime; also enumerate all labelled solutions with
D + ell >= 25 (circles <= 10) under (C1) only, to cross-check enum_c7.py's labelled counts.
"""
import itertools, sys
from math import comb
from ortools.sat.python import cp_model
n = 7
subs = {k: [frozenset(c) for c in itertools.combinations(range(n), k)] for k in range(3, 7)}
blocks = subs[4] + subs[5] + subs[6]
lines = subs[3] + subs[4] + subs[5] + subs[6]
def build(dcap, lcap):
    m = cp_model.CpModel()
    x = {B: m.NewBoolVar(str(sorted(B))) for B in blocks}
    y = {S: m.NewBoolVar("y" + str(sorted(S))) for S in lines}
    for A, B in itertools.combinations(blocks, 2):
        if len(A & B) >= 3: m.AddBoolOr([x[A].Not(), x[B].Not()])
    for A, B in itertools.combinations(lines, 2):
        if len(A & B) >= 2: m.AddBoolOr([y[A].Not(), y[B].Not()])
    for S in lines:
        if len(S) >= 4: m.AddImplication(y[S], x[S])
        else:
            for B in blocks:
                if S < B: m.AddBoolOr([y[S].Not(), x[B].Not()])
    for p in range(n):
        m.Add(sum(comb(len(B) - 1, 2) * x[B] for B in blocks if p in B) <= dcap)
    m.Add(sum(comb(len(S), 2) * y[S] for S in lines) <= lcap)
    obj = sum((comb(len(B), 3) - 1) * x[B] for B in blocks) + sum(y[S] for S in lines)
    return m, x, y, obj
for name, (dcap, lcap) in {"none": (1000, 1000), "weak(SG only)": (14, 20), "strong(o6=o7=3)": (12, 18)}.items():
    m, x, y, obj = build(dcap, lcap); m.Maximize(obj)
    sol = cp_model.CpSolver(); sol.parameters.num_workers = 2; sol.parameters.max_time_in_seconds = 300
    st = sol.Solve(m)
    Fb = sorted(sorted(B) for B in blocks if sol.Value(x[B])); Lb = sorted(sorted(S) for S in lines if sol.Value(y[S]))
    print(f"regime {name}: {sol.StatusName(st)} max D+ell = {sol.ObjectiveValue():.0f} (bound {sol.BestObjectiveBound():.0f}) "
          f"=> min circles = {35 - sol.ObjectiveValue():.0f}; blocks={Fb} lines={Lb}")
# enumerate all labelled solutions with D+ell >= 25 under (C1) only
m, x, y, obj = build(1000, 1000); m.Add(obj >= 25)
class Coll(cp_model.CpSolverSolutionCallback):
    def __init__(s): super().__init__(); s.cnt = {}
    def on_solution_callback(s):
        D = sum((comb(len(B), 3) - 1) for B in blocks if s.Value(x[B])); ell = sum(1 for S in lines if s.Value(y[S]))
        nb = sum(1 for B in blocks if s.Value(x[B]))
        s.cnt[(nb, D, ell)] = s.cnt.get((nb, D, ell), 0) + 1
sol = cp_model.CpSolver(); sol.parameters.enumerate_all_solutions = True; sol.parameters.num_workers = 1
c = Coll(); st = sol.Solve(m, c)
print("all labelled (F,L) with D+ell>=25 under (C1) only:", sol.StatusName(st), "total", sum(c.cnt.values()))
for k in sorted(c.cnt): print("  (#blocks, D, ell) =", k, "circles =", 35 - k[1] - k[2], ":", c.cnt[k], "labelled")
