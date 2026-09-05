"""Case B by skeleton: for each family of 5-blocks (skeletons.json) fix exactly those 5-blocks,
forbid all other 5-blocks, and (1) maximise D + l with 2 workers; (2) if the proven optimum is >= thr,
enumerate all labelled solutions with D + l >= thr and store isomorphism classes.
Usage: python3 extend.py omode [--thr 88] [--time T] [--only i,j,...] [--orchard10] [--no-mk] [--no-fano]
       [--out survivors_<omode>.json]
"""
import sys, json, time, math, argparse, itertools
from ortools.sat.python import cp_model
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from n10model import Model, ClassStore, structure_count, PTS

ap = argparse.ArgumentParser()
ap.add_argument('omode'); ap.add_argument('--thr', type=int, default=88); ap.add_argument('--time', type=float, default=1200)
ap.add_argument('--only', default=None); ap.add_argument('--orchard10', action='store_true')
ap.add_argument('--no-mk', action='store_true'); ap.add_argument('--no-fano', action='store_true')
ap.add_argument('--no-orchard9', action='store_true')
ap.add_argument('--out', default=None); ap.add_argument('--workers', type=int, default=2)
ap.add_argument('--skip-opt', action='store_true'); ap.add_argument('--bundle', action='store_true')
args = ap.parse_args()
sk = json.load(open(__file__.rsplit('/', 1)[0] + '/skeletons.json'))
idx = list(range(len(sk))) if args.only is None else [int(t) for t in args.only.split(',')]
FIVES = [frozenset(c) for c in itertools.combinations(PTS, 5)]
results = []
survivors = []
T0 = time.time()
for i in idx:
    rec = sk[i]
    fam = [frozenset(B) for B in rec['blocks']]
    fixed = {tuple(sorted(B)): (1 if B in fam else 0) for B in FIVES}
    t0 = time.time()
    M = Model('B', args.omode, fano=not args.no_fano, mk=not args.no_mk, orchard9=not args.no_orchard9,
              orchard10=args.orchard10, degree_order=(i in (0, 3)), threshold=None, extra_fixed=fixed, bundle=args.bundle)
    line = f"skeleton {i} (k={rec['k']}) blocks={rec['blocks']}"
    opt_ok = True
    if not args.skip_opt:
        # phase 1: feasibility of D + l >= thr with several workers (INFEASIBLE is a proof)
        M.m.Add(M.obj >= args.thr)
        solver = cp_model.CpSolver(); solver.parameters.num_workers = args.workers
        solver.parameters.max_time_in_seconds = args.time
        st = solver.Solve(M.m)
        line += f" | feas(D+l>={args.thr}) status={solver.StatusName(st)} ({time.time()-t0:.0f}s)"
        if st == cp_model.INFEASIBLE:
            opt_ok = False
        results.append({'i': i, 'k': rec['k'], 'status': solver.StatusName(st), 'time': time.time() - t0})
    print(line, flush=True)
    if not opt_ok:
        continue
    # enumeration
    M = Model('B', args.omode, fano=not args.no_fano, mk=not args.no_mk, orchard9=not args.no_orchard9,
              orchard10=args.orchard10, degree_order=(i in (0, 3)), threshold=args.thr, extra_fixed=fixed, bundle=args.bundle)
    solver = cp_model.CpSolver(); solver.parameters.num_workers = 1
    solver.parameters.enumerate_all_solutions = True; solver.parameters.max_time_in_seconds = args.time
    store = ClassStore()
    class CB(cp_model.CpSolverSolutionCallback):
        def on_solution_callback(s):
            F, L = M.extract(lambda v: s.Value(v))
            if store.add(F, L):
                print(f"    [{time.time()-t0:.0f}s] class #{len(store.classes)} (labelled {store.total}) sizes={store.classes[-1][2][0]} "
                      f"lines={store.classes[-1][2][1]} count={structure_count(F, L)}", flush=True)
    st = solver.Solve(M.m, CB())
    print(f"    enum status={solver.StatusName(st)} labelled={store.total} classes={len(store.classes)} ({time.time()-t0:.0f}s)", flush=True)
    for c in store.classes:
        survivors.append({'skeleton': i, 'blocks': sorted(sorted(B) for B in c[0]), 'lines': sorted(sorted(S) for S in c[1]),
                          'copies': c[3], 'count': structure_count(c[0], c[1]), 'enum_status': solver.StatusName(st)})
    if args.out:
        json.dump({'results': results, 'survivors': survivors}, open(args.out, 'w'))
print(f"done in {time.time()-T0:.0f}s; survivors={len(survivors)}")
if args.out:
    json.dump({'results': results, 'survivors': survivors}, open(args.out, 'w'))
