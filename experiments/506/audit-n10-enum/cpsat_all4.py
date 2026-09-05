"""Independent CP-SAT checks for the all-4-block case n = 10 (SG-only caps: degree <= 11, line pairs <= 44).
Variables: x_c (4-subset c is a block), y_t (uncovered triple t is a 3-line), z_c (block c is a 4-line).
Constraints: each triple in <= 1 block; y_t + sum_{c ⊇ t} x_c <= 1; z_c <= x_c; each pair in <= 1 line;
3 sum y + 6 sum z <= 44; degrees <= dmax.
Check 1: dmax = 10, sum x >= 25, lines >= 13  (the branch with no degree-11 point; need D + ell >= 88).
Check 2: dmax = 11, 3 sum x + lines >= 89     (would give count <= 31)."""
import itertools, sys, time
from ortools.sat.python import cp_model
n = 10
subs = list(itertools.combinations(range(n), 4))
trips = list(itertools.combinations(range(n), 3))
def build(dmax):
    m = cp_model.CpModel()
    x = [m.NewBoolVar(f"x{i}") for i in range(len(subs))]
    y = [m.NewBoolVar(f"y{i}") for i in range(len(trips))]
    z = [m.NewBoolVar(f"z{i}") for i in range(len(subs))]
    for j, t in enumerate(trips):
        sup = [x[i] for i, s in enumerate(subs) if set(t) <= set(s)]
        m.Add(sum(sup) <= 1)
        m.Add(y[j] + sum(sup) <= 1)
    for i in range(len(subs)):
        m.Add(z[i] <= x[i])
    for p, q in itertools.combinations(range(n), 2):
        m.Add(sum(y[j] for j, t in enumerate(trips) if p in t and q in t) + sum(z[i] for i, s in enumerate(subs) if p in s and q in s) <= 1)
    for p in range(n):
        m.Add(sum(x[i] for i, s in enumerate(subs) if p in s) <= dmax)
    m.Add(3 * sum(y) + 6 * sum(z) <= 44)
    return m, x, y, z
def solve(m, label):
    s = cp_model.CpSolver(); s.parameters.num_workers = 2; s.parameters.max_time_in_seconds = 1500
    t0 = time.time(); st = s.Solve(m)
    print(f"{label}: status {s.StatusName(st)} [{time.time()-t0:.1f}s]", flush=True)
    return s, st
m, x, y, z = build(10)
m.Add(sum(x) >= 25); m.Add(sum(y) + sum(z) >= 13)
s, st = solve(m, "check 1 (degrees <= 10, >= 25 blocks, >= 13 lines)")
if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    print("  blocks", [subs[i] for i in range(len(subs)) if s.Value(x[i])]); print("  lines", [trips[j] for j in range(len(trips)) if s.Value(y[j])] + [subs[i] for i in range(len(subs)) if s.Value(z[i])])
m, x, y, z = build(11)
m.Add(3 * sum(x) + sum(y) + sum(z) >= 89)
s, st = solve(m, "check 2 (degrees <= 11, 3*b4 + ell >= 89, i.e. count <= 31)")
if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    print("  blocks", [subs[i] for i in range(len(subs)) if s.Value(x[i])]); print("  lines", [trips[j] for j in range(len(trips)) if s.Value(y[j])] + [subs[i] for i in range(len(subs)) if s.Value(z[i])])
