"""Independent CP-SAT feasibility check for the all-4-block branch b4 = 25 with every point in <= 10 four-blocks
(this is the only all-4-block branch in which no point has degree 11): is there a family of >= 25 four-subsets of
10 points, every triple in at most one of them, all point degrees <= 10?"""
import itertools
from ortools.sat.python import cp_model
n = 10
subs = list(itertools.combinations(range(n), 4))
m = cp_model.CpModel()
x = [m.NewBoolVar(f"x{i}") for i in range(len(subs))]
for t in itertools.combinations(range(n), 3):
    m.Add(sum(x[i] for i, s in enumerate(subs) if set(t) <= set(s)) <= 1)
for p in range(n):
    m.Add(sum(x[i] for i, s in enumerate(subs) if p in s) <= 10)
m.Add(sum(x) >= 25)
s = cp_model.CpSolver(); s.parameters.num_workers = 2; s.parameters.max_time_in_seconds = 1200
st = s.Solve(m)
print("status:", s.StatusName(st))
if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    print([subs[i] for i in range(len(subs)) if s.Value(x[i])])
