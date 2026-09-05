"""Apply the necessary realisability conditions (common.full_report) to every (F, L) class of the given
enumeration files.  Tier T = theorems only (hereditary SG on all derived spaces, orchard t3(7..9) proved and
t3(10) <= 13 by parity, Miquel closure); tier C adds the cited tables (t3(10) = 12; o(9) = 6, o(10) = 5).
Usage: python3 postfilter.py runs/enum_*.json --out runs/post.json"""
import sys, json, argparse, math, time
sys.path.insert(0, '/home/user/Solving-hard-math-problems-/experiments/506/n10-continue')
from common import full_report, ClassStore, count, deficit, fs
ap = argparse.ArgumentParser(); ap.add_argument('files', nargs='+'); ap.add_argument('--out', default=None)
a = ap.parse_args()
store = ClassStore(); src = {}
for fn in a.files:
    d = json.load(open(fn))
    for r in d['structures']:
        idx, new = store.add(r['blocks'], r['lines'])
        src.setdefault(idx, []).append((fn, d['stage1_status']))
print(f"{len(store.classes)} distinct (F, L) classes from {len(a.files)} files")
out = []; surv_T = []; surv_C = []
t0 = time.time()
for k, c in enumerate(store.classes):
    F, L = c['F'], c['L']
    D = sum(deficit(len(B)) for B in F)
    repT = full_report(F, L, cited=False)
    repC = full_report(F, L, cited=True, miquel=False)
    line = f"structure {k}: sizes={sorted(len(B) for B in F)} lines={sorted(len(S) for S in L)} D={D} l={len(L)} count={count(F, L)} -> theorems: {repT['kills'] or 'SURVIVES'}; +cited: {repC['kills'] or 'SURVIVES'}"
    if repT['hsg']: line += f"\n     HSG e.g. {list(repT['hsg'].items())[0]}"
    if repT['orchard']: line += f"\n     orchard e.g. {list(repT['orchard'].items())[0]}"
    if repT['miquel']: line += f"\n     Miquel cube faces {repT['miquel']}"
    if repC['ordinary']: line += f"\n     ordinary-lines (cited) e.g. {list(repC['ordinary'].items())[0]}"
    print(line, flush=True)
    if not repT['kills']: surv_T.append(k)
    if not repT['kills'] and not repC['kills']: surv_C.append(k)
    out.append({'index': k, 'blocks': sorted(sorted(B) for B in F), 'lines': sorted(sorted(S) for S in L), 'D': D, 'l': len(L),
                'count': count(F, L), 'theorem_kills': repT['kills'], 'cited_kills': repC['kills'], 'sources': src[k]})
print(f"survivors (theorems only): {surv_T}\nsurvivors (theorems + cited tables): {surv_C}   [{time.time()-t0:.0f}s]")
if a.out:
    json.dump(out, open(a.out, 'w'), indent=0)
