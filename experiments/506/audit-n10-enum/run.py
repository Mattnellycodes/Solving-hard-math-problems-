"""Usage: python3 run.py n target mode [--all4] [--sizes 5,6,...]
Enumerates all classes of (F,L) with circles <= target under caps `mode` ('sg' or 'table').
--all4: also run the empty big family with G = S_n (all-4-block case) -- feasible for n <= 9."""
import sys, json, time
from math import comb
from mob import *
n = int(sys.argv[1]); target = int(sys.argv[2]); mode = sys.argv[3]
all4 = '--all4' in sys.argv
sizes = None
for a in sys.argv:
    if a.startswith('--sizes'):
        sizes = [int(x) for x in a.split('=')[1].split(',')]
capD, capL = caps(n, mode)
need = comb(n, 3) - target
print(f"n={n} target={target} mode={mode} capD={capD} capL={capL} ellmax={capL//3} need D+ell>={need}", flush=True)
PT = PermTable(n)
res = run_big(n, target, mode, sizes=sizes, PT=PT)
if all4:
    t0 = time.time()
    S = Stage2(n, [], PT.P, capD, capL, need, verbose=True)
    c4 = S.run()
    print(f"[all-4-block] nodes={S.nodes} evals={S.evals} candidates={len(c4)} [{time.time()-t0:.1f}s]", flush=True)
    for c in c4: c['big'] = []
    res.extend(c4)
for i, c in enumerate(res):
    print(f"cand {i}: sizes={[len(b) for b in c['blocks']]} D={c['D']} ell_max={c['ell_max']} count={c['count_min']} degrees={c['degrees']} #line_sets={len(c['line_sets'])}")
    print(f"   blocks={c['blocks']}")
json.dump(dict(n=n, target=target, mode=mode, capD=capD, capL=capL, results=res), open(f"out_n{n}_t{target}_{mode}{'_all4' if all4 else ''}.json", 'w'))
