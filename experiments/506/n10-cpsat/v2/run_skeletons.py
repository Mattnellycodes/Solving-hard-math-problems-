"""Per-skeleton feasibility of D + l >= thr (default 88, i.e. <= 32 circles) for n = 10.
Usage: python3 run_skeletons.py omode [--time T] [--only i,j] [--orchard10] [--thr K] [--workers W] [--out file]
For every skeleton (skeletons.json) the big blocks are fixed exactly (all other 5/6-blocks forbidden); block sizes
are (4,5,6) in case A and (4,5) in case B.  INFEASIBLE is a proof (modulo the solver) that no structure with these
big blocks reaches the threshold.
"""
import sys, json, time, argparse
sys.path.insert(0, '/home/user/Solving-hard-math-problems-/experiments/506/n10-cpsat/v2')
from model import Relaxation, count
from ortools.sat.python import cp_model

ap = argparse.ArgumentParser()
ap.add_argument('omode'); ap.add_argument('--time', type=float, default=1500); ap.add_argument('--only', default=None)
ap.add_argument('--orchard10', action='store_true'); ap.add_argument('--thr', type=int, default=88)
ap.add_argument('--workers', type=int, default=2); ap.add_argument('--out', default=None)
ap.add_argument('--no-e11', action='store_true')
args = ap.parse_args()
sk = json.load(open('/home/user/Solving-hard-math-problems-/experiments/506/n10-cpsat/v2/skeletons.json'))
idx = list(range(len(sk))) if args.only is None else [int(t) for t in args.only.split(',')]
results = []
for i in idx:
    rec = sk[i]
    sizes = (4, 5, 6) if rec['case'] == 'A' else (4, 5)
    t0 = time.time()
    M = Relaxation(n=10, sizes=sizes, o=args.omode, threshold=args.thr, skeleton=rec['blocks'], orchard10=args.orchard10,
                   e11=not args.no_e11)
    tb = time.time() - t0
    s = cp_model.CpSolver(); s.parameters.num_workers = args.workers; s.parameters.max_time_in_seconds = args.time
    st = s.Solve(M.m)
    line = f"skeleton {i} case {rec['case']} k={rec['k']} blocks={rec['blocks']}: status={s.StatusName(st)} (build {tb:.1f}s, solve {s.WallTime():.0f}s)"
    out = {'i': i, 'case': rec['case'], 'k': rec['k'], 'blocks': rec['blocks'], 'status': s.StatusName(st), 'time': s.WallTime()}
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        F, L = M.extract(lambda v: s.Value(v))
        out['example'] = {'blocks': sorted(sorted(B) for B in F), 'lines': sorted(sorted(S) for S in L), 'count': count(10, F, L)}
        line += f"\n    example: count={out['example']['count']} blocks={out['example']['blocks']} lines={out['example']['lines']}"
    print(line, flush=True)
    results.append(out)
    if args.out:
        json.dump(results, open(args.out, 'w'), indent=0)
print("done")
