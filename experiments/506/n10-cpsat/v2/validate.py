"""Validation of model.py on n = 8 and n = 9 against known values:
 n = 8 (SG-only caps, no skeleton): optimum D + l = 39 (17 circles = c(8)) [theory REPORT / verify_independent].
 n = 9 (SG-only caps): with the MK-derived constraints the structures with D + l >= 60 (<= 24 circles) should all be
       'two 5-blocks through a point + 12 four-blocks' (count 24) [verify_independent/n9_maxdeg7_sg.out]."""
import sys, time, math
sys.path.insert(0, '/home/user/Solving-hard-math-problems-/experiments/506/n10-cpsat/v2')
from model import Relaxation, ClassStore, count
from ortools.sat.python import cp_model

for n, sizes, thr in ((8, (4, 5, 6, 7), None), (9, (4, 5, 6, 7, 8), 60)):
    t0 = time.time()
    M = Relaxation(n=n, sizes=sizes, o='sg', threshold=thr, e11=True)
    print(f"n={n}: built in {time.time()-t0:.1f}s; constraint counts {M.ncons}", flush=True)
    s = cp_model.CpSolver(); s.parameters.num_workers = 2; s.parameters.max_time_in_seconds = 900
    if thr is None:
        M.m.Maximize(M.obj)
        st = s.Solve(M.m)
        F, L = M.extract(lambda v: s.Value(v))
        print(f"   status={s.StatusName(st)} opt D+l={s.ObjectiveValue():.0f} bound={s.BestObjectiveBound():.0f} -> circles >= {math.comb(n,3)-round(s.ObjectiveValue())}; "
              f"example blocks {sorted(sorted(B) for B in F)} lines {sorted(sorted(S) for S in L)}  [{time.time()-t0:.0f}s]", flush=True)
    else:
        s.parameters.num_workers = 1; s.parameters.enumerate_all_solutions = True
        store = ClassStore(n)
        class CB(cp_model.CpSolverSolutionCallback):
            def on_solution_callback(self):
                F, L = M.extract(lambda v: self.Value(v)); store.add(F, L)
        st = s.Solve(M.m, CB())
        print(f"   status={s.StatusName(st)} labelled={store.total} classes={len(store.classes)}  [{time.time()-t0:.0f}s]")
        for r in store.records():
            print("   class:", r['blocks'], r['lines'], "copies", r['copies'], "count", r['count'])
