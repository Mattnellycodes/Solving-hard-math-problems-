"""Extra (outside the assigned task): largest block 4.  Feasibility of D + l >= 88 with only 4-blocks
under SG-only caps + (F7)/(M8) + orchard10 (t3(10) <= 12, proved in t3_10_proof.py).  Counting: 25 four-blocks
(every point in 10) and 13 lines = 12 three-lines + one 4-line are needed; with o(10) >= 5 (Kelly–Moser)
this is already excluded (42 > 40 pairs)."""
import sys, time
from ortools.sat.python import cp_model
sys.path.insert(0, '/home/user/Solving-hard-math-problems-/experiments/506/n10-cpsat')
from n10model import Model
for omode in ('sg', 'table'):
    M = Model('C', omode, orchard10=True, threshold=88)
    s = cp_model.CpSolver(); s.parameters.num_workers = 1; s.parameters.max_time_in_seconds = 1500
    t0 = time.time(); st = s.Solve(M.m)
    print(f"case C (largest block 4) omode={omode} orchard10: status={s.StatusName(st)} in {time.time()-t0:.0f}s", flush=True)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        F, L = M.extract(lambda v: s.Value(v)); print("  blocks", sorted(sorted(B) for B in F)); print("  lines", sorted(sorted(S) for S in L))
