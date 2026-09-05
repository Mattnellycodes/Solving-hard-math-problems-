"""Confirm the counting lemma for case A by CP-SAT: with a 6-block {0..5} and D + l >= 88 (SG-only caps),
every solution has exactly 18 four-blocks with 2 points in {0..5} and 2 in R = {6..9}, and R is a block.
Two feasibility runs on the complement: (a) R not a block; (b) R a block and at most 17 (2,2)-blocks."""
import sys, time
from ortools.sat.python import cp_model
sys.path.insert(0, '/home/user/Solving-hard-math-problems-/experiments/506/n10-cpsat')
from n10model import Model
R = frozenset({6, 7, 8, 9})
for tag in ('a', 'b'):
    M = Model('A', 'sg', threshold=88, orchard9=True)
    n22 = sum(M.x[B] for B in M.blocks if len(B) == 4 and len(B & R) == 2)
    if tag == 'a':
        M.m.Add(M.x[R] == 0)
    else:
        M.m.Add(M.x[R] == 1); M.m.Add(n22 <= 17)
    s = cp_model.CpSolver(); s.parameters.num_workers = 2; s.parameters.max_time_in_seconds = 1200
    t0 = time.time(); st = s.Solve(M.m)
    print(f"case A complement ({tag}): status={s.StatusName(st)} in {time.time()-t0:.0f}s", flush=True)
