"""Driver: two-phase orderly enumeration of all abstract structures with count <= target.
usage: python3 run_enum.py n mode target [--sizes 5,6,7,8,9] [--empty auto|skip|force] [--out file.json]
Outputs every family F (big blocks + 4-blocks) with D(F) + ellmax(F) >= C(n,3) - target, one per
isomorphism class (for the big-block part; see mo.py), together with ellmax(F) and an example line set.
"""
import sys, math, time, json, argparse
import mo

ap = argparse.ArgumentParser()
ap.add_argument('n', type=int); ap.add_argument('mode'); ap.add_argument('target', type=int)
ap.add_argument('--sizes', default=None)
ap.add_argument('--empty', default='auto', help="empty big-block family: auto (run iff root bound allows and |S_n| <= 50000), skip, force")
ap.add_argument('--out', default=None)
args = ap.parse_args()
n, mode, target = args.n, args.mode, args.target
sizes = [int(s) for s in args.sizes.split(',')] if args.sizes else list(range(5, n))
caps = mo.Caps(n, mode)
need = math.comb(n, 3) - target
print(f"n={n} mode={mode} target={target} need=D+ell>={need} {caps} big sizes={sizes}", flush=True)
t0 = time.time()
reps = mo.phase1(n, caps, sizes, verbose=True)
print(f"phase 1: {len(reps)} classes of big-block families (incl. empty) [{time.time()-t0:.1f}s]", flush=True)
allres = []
tot_nodes = 0
for ri, R in enumerate(reps):
    t1 = time.time()
    D0 = mo.deficit(R)
    bound, M = mo.root_bound(n, caps, R, need)
    desc = f"[{ri}] sizes={[mo.popcount(b) for b in R]} blocks={[mo.bits(b) for b in R]} D0={D0} cands={M} rootD<={bound} (need D>={need-caps.ellmax})"
    if bound + caps.ellmax < need:
        print(desc + " -> pruned at root", flush=True); continue
    if not R:
        if args.empty == 'skip' or (args.empty == 'auto' and math.factorial(n) > 50000):
            print(desc + " -> EMPTY FAMILY SKIPPED (all-4-block case not enumerated here)", flush=True); continue
    aut = mo.automorphisms(n, R)
    res, st = mo.phase2(n, caps, R, aut, need)
    tot_nodes += st['nodes']
    for r in res:
        r['rep'] = ri; r['sizes'] = sorted((len(b) for b in r['blocks']), reverse=True)
        r['degrees'] = [sum(1 for b in r['blocks'] if p in b) for p in range(n)]
    allres.extend(res)
    print(desc + f" |Aut|={len(aut)} nodes={st['nodes']} evals={st['evals']} candidates={len(res)} [{time.time()-t1:.1f}s]", flush=True)
print(f"TOTAL phase-2 nodes={tot_nodes}; candidate structures (count <= {target}): {len(allres)} [{time.time()-t0:.1f}s]")
for i, r in enumerate(allres):
    print(f"  cand {i}: sizes={r['sizes']} D={r['D']} ellmax={r['ellmax']} count={r['count']} degrees={r['degrees']}")
    print(f"     blocks={r['blocks']}")
if args.out:
    json.dump({'n': n, 'mode': mode, 'target': target, 'need': need, 'caps': repr(caps), 'sizes': sizes,
               'phase1_reps': [[mo.bits(b) for b in R] for R in reps], 'candidates': allres}, open(args.out, 'w'))
