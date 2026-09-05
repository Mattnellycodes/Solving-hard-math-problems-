"""Complete enumeration, up to isomorphism, of all structures (F, L) with D + l >= thr whose big blocks are exactly
a given skeleton — by CP-SAT with ORBIT NOGOODS:
  repeat: solve the model; if INFEASIBLE stop; else take the 4-block family F of the solution, compute its orbit
  under the automorphism group G of the skeleton (all point permutations preserving the skeleton), and add for
  every image F' the clause 'the set of 4-blocks is not exactly F''.
Since every labelled structure with these big blocks has its 4-block family in the G-orbit of some found F, and a
clause only excludes exact copies, the loop terminates with a complete list of 4-block-family classes.
Then for each class all line sets L with |L| >= thr - D(F) are enumerated exhaustively (small CP-SAT over the y's
with the same line constraints), and (F, L) pairs are reduced to isomorphism classes.
Usage: python3 enum_skeleton.py omode i [--time T] [--orchard10] [--thr K] [--workers W] [--out file]
"""
import sys, json, time, argparse, itertools, math
sys.path.insert(0, '/home/user/Solving-hard-math-problems-/experiments/506/n10-cpsat/v2')
from model import Relaxation, ClassStore, count, automorphisms
from ortools.sat.python import cp_model

ap = argparse.ArgumentParser()
ap.add_argument('omode'); ap.add_argument('i', type=int); ap.add_argument('--time', type=float, default=2400)
ap.add_argument('--orchard10', action='store_true'); ap.add_argument('--thr', type=int, default=88)
ap.add_argument('--workers', type=int, default=2); ap.add_argument('--out', default=None)
ap.add_argument('--no-e11', action='store_true')
args = ap.parse_args()
sk = json.load(open('/home/user/Solving-hard-math-problems-/experiments/506/n10-cpsat/v2/skeletons.json'))
rec = sk[args.i]
sizes = (4, 5, 6) if rec['case'] == 'A' else (4, 5)
skel = [frozenset(B) for B in rec['blocks']]
T0 = time.time()
G = automorphisms(10, skel, [])
print(f"skeleton {args.i} {rec['blocks']}: |Aut| = {len(G)}  [{time.time()-T0:.0f}s]", flush=True)
M = Relaxation(n=10, sizes=sizes, o=args.omode, threshold=args.thr, skeleton=rec['blocks'], orchard10=args.orchard10,
               e11=not args.no_e11)
fours = [B for B in M.blocks if len(B) == 4]
fstore = ClassStore(10)
found = []          # list of (F class record)
n_nogood = 0
status_final = None
while True:
    s = cp_model.CpSolver(); s.parameters.num_workers = args.workers; s.parameters.max_time_in_seconds = args.time
    st = s.Solve(M.m)
    if st == cp_model.INFEASIBLE:
        status_final = 'COMPLETE'; print(f"  no further 4-block family: INFEASIBLE  [{time.time()-T0:.0f}s]", flush=True); break
    if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        status_final = 'INCOMPLETE:' + s.StatusName(st); print(f"  solver stopped with {s.StatusName(st)} -> enumeration INCOMPLETE", flush=True); break
    F, L = M.extract(lambda v: s.Value(v))
    F4 = frozenset(B for B in F if len(B) == 4)
    idx, new = fstore.add(list(F), [])
    assert new, "found a 4-block family already excluded?!"
    orbit = {frozenset(frozenset(g[p] for p in B) for B in F4) for g in G}
    for img in orbit:
        M.m.AddBoolOr([M.x[B].Not() if B in img else M.x[B] for B in fours])
    n_nogood += len(orbit)
    D = sum(math.comb(len(B), 3) - 1 for B in F)
    print(f"  F-class #{len(fstore.classes)}: b4={len(F4)} D={D} (example l={len(L)}, count={count(10, F, L)}) orbit size {len(orbit)}; "
          f"nogoods so far {n_nogood}  [{time.time()-T0:.0f}s]", flush=True)
    found.append({'F': F, 'D': D})
print(f"4-block-family classes: {len(found)}; enumeration status {status_final}", flush=True)

# ---- line sets for each F class
store = ClassStore(10)
per_class = []
for k, fc in enumerate(found):
    F = fc['F']
    Mx = Relaxation(n=10, sizes=sizes, o=args.omode, threshold=args.thr, skeleton=rec['blocks'], orchard10=args.orchard10,
                    e11=not args.no_e11, degree_order=False)
    for B in Mx.blocks:
        if len(B) == 4:
            Mx.m.Add(Mx.x[B] == (1 if B in F else 0))
    s = cp_model.CpSolver(); s.parameters.num_workers = 1; s.parameters.enumerate_all_solutions = True
    s.parameters.max_time_in_seconds = args.time
    sols = []
    class CB(cp_model.CpSolverSolutionCallback):
        def on_solution_callback(self):
            sols.append([S for S in Mx.lines if self.Value(Mx.y[S])])
    st = s.Solve(Mx.m, CB())
    before = len(store.classes)
    for L in sols:
        store.add(F, L)
    per_class.append({'F': sorted(sorted(B) for B in F), 'D': fc['D'], 'line_sets': len(sols), 'status': s.StatusName(st),
                      'new_FL_classes': len(store.classes) - before})
    print(f"  F-class {k}: D={fc['D']}: {len(sols)} labelled line sets (status {s.StatusName(st)}), (F,L) classes so far {len(store.classes)}  [{time.time()-T0:.0f}s]", flush=True)
recs = store.records()
print(f"TOTAL (F, L) classes: {len(recs)}; F classes {len(found)}; enumeration status {status_final}; time {time.time()-T0:.0f}s")
for r in recs:
    print("   count", r['count'], "copies", r['copies'], "blocks", r['blocks'], "lines", r['lines'])
if args.out:
    json.dump({'skeleton': args.i, 'omode': args.omode, 'orchard10': args.orchard10, 'thr': args.thr, 'status': status_final,
               'F_classes': per_class, 'structures': recs}, open(args.out, 'w'), indent=0)
