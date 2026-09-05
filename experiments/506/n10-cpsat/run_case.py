"""Runner.  python3 run_case.py A|B omode [--opt|--enum] [--time T] [--workers W] [--thr K]
   [--no-fano] [--no-mk] [--no-orchard9] [--orchard10] [--no-pairdeg] [--no-order] [--fix B=v ...]
   [--b5 K] [--out file.json]
--opt : maximise D + l (with the threshold constraint removed) and report the optimum / bound.
--enum: enumerate all labelled solutions with D + l >= thr (default 88); store isomorphism classes.
"""
import sys, json, time, math, argparse
from ortools.sat.python import cp_model
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from n10model import Model, ClassStore, structure_count, PTS

ap = argparse.ArgumentParser()
ap.add_argument('case'); ap.add_argument('omode')
ap.add_argument('--opt', action='store_true'); ap.add_argument('--enum', action='store_true')
ap.add_argument('--time', type=float, default=1800); ap.add_argument('--workers', type=int, default=2)
ap.add_argument('--thr', type=int, default=88)
ap.add_argument('--no-fano', action='store_true'); ap.add_argument('--no-mk', action='store_true')
ap.add_argument('--no-orchard9', action='store_true'); ap.add_argument('--orchard10', action='store_true')
ap.add_argument('--no-pairdeg', action='store_true'); ap.add_argument('--no-order', action='store_true')
ap.add_argument('--bundle', action='store_true'); ap.add_argument('--fix', nargs='*', default=[]); ap.add_argument('--b5', type=int, default=None)
ap.add_argument('--b5min', type=int, default=None); ap.add_argument('--b5max', type=int, default=None)
ap.add_argument('--out', default=None); ap.add_argument('--log', action='store_true')
args = ap.parse_args()

extra = {}
for f in args.fix:
    B, v = f.split('=')
    extra[tuple(int(c) for c in B)] = int(v)
t0 = time.time()
M = Model(args.case, args.omode, fano=not args.no_fano, mk=not args.no_mk, orchard9=not args.no_orchard9,
          orchard10=args.orchard10, pairdeg=not args.no_pairdeg, degree_order=not args.no_order,
          threshold=None if args.opt else args.thr, extra_fixed=extra, bundle=args.bundle)
m = M.m
b5 = sum(M.x[B] for B in M.blocks if len(B) == 5)
if args.b5 is not None:
    m.Add(b5 == args.b5)
if args.b5min is not None:
    m.Add(b5 >= args.b5min)
if args.b5max is not None:
    m.Add(b5 <= args.b5max)
print(f"case {args.case} omode {args.omode} o={M.o} subset-caps={M.n_sub} bundle={M.n_bundle} vars={len(M.blocks)}+{len(M.lines)} "
      f"built in {time.time()-t0:.1f}s; options={vars(args)}", flush=True)
solver = cp_model.CpSolver()
solver.parameters.max_time_in_seconds = args.time
solver.parameters.log_search_progress = args.log
if args.opt:
    m.Maximize(M.obj)
    solver.parameters.num_workers = args.workers
    st = solver.Solve(m)
    print(f"status={solver.StatusName(st)} obj={solver.ObjectiveValue():.0f} bound={solver.BestObjectiveBound():.0f} "
          f"time={solver.WallTime():.1f}s", flush=True)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        F, L = M.extract(lambda v: solver.Value(v))
        print("  blocks:", sorted(sorted(B) for B in F)); print("  lines:", sorted(sorted(S) for S in L))
        print("  circles =", structure_count(F, L))
else:
    solver.parameters.num_workers = 1
    solver.parameters.enumerate_all_solutions = True
    store = ClassStore()
    class CB(cp_model.CpSolverSolutionCallback):
        def __init__(s):
            super().__init__(); s.t = time.time()
        def on_solution_callback(s):
            F, L = M.extract(lambda v: s.Value(v))
            new = store.add(F, L)
            if new:
                print(f"  [{time.time()-t0:.0f}s] new class #{len(store.classes)} (labelled so far {store.total}): "
                      f"sizes={store.classes[-1][2][0]} lines={store.classes[-1][2][1]} count={structure_count(F, L)}", flush=True)
    cb = CB()
    st = solver.Solve(m, cb)
    print(f"status={solver.StatusName(st)} labelled={store.total} classes={len(store.classes)} time={solver.WallTime():.1f}s", flush=True)
    if args.out:
        json.dump([{'blocks': sorted(sorted(B) for B in c[0]), 'lines': sorted(sorted(S) for S in c[1]),
                    'copies': c[3], 'count': structure_count(c[0], c[1])} for c in store.classes], open(args.out, 'w'))
    for c in store.classes:
        print(" class: blocks", sorted(sorted(B) for B in c[0]), "lines", sorted(sorted(S) for S in c[1]), "copies", c[3])
