"""Apply the necessary conditions (filters.py) and the numerical realiser (realise.py) to the structures found by
enum_skeleton.py.  Usage: python3 postprocess.py enum_file.json [...] [--tries N] [--cited] [--out report.json]"""
import sys, json, argparse, math
sys.path.insert(0, '/home/user/Solving-hard-math-problems-/experiments/506/n10-cpsat/v2')
from filters import full_report
from realise import try_realise
from model import ClassStore

ap = argparse.ArgumentParser()
ap.add_argument('files', nargs='+'); ap.add_argument('--tries', type=int, default=200); ap.add_argument('--cited', action='store_true')
ap.add_argument('--out', default=None); ap.add_argument('--no-realise', action='store_true')
args = ap.parse_args()
store = ClassStore(10)
src = {}
for fn in args.files:
    data = json.load(open(fn))
    for rec in data['structures']:
        idx, new = store.add(rec['blocks'], rec['lines'])
        src.setdefault(idx, []).append((fn, data['status']))
print(f"{len(store.classes)} distinct (F, L) isomorphism classes from {len(args.files)} file(s)")
out = []
alive = []
for k, c in enumerate(store.classes):
    F, L = c['F'], c['L']
    D = sum(math.comb(len(B), 3) - 1 for B in F)
    rep = full_report(F, L, cited=args.cited)
    line = f"structure {k}: sizes={sorted(len(B) for B in F)} lines={sorted(len(S) for S in L)} D={D} l={len(L)} count={120-D-len(L)} -> kills={rep['kills']}"
    if rep['hereditary_sg']:
        line += f"\n    hereditary SG violations at {list(rep['hereditary_sg'].keys())}: e.g. {list(rep['hereditary_sg'].items())[0]}"
    if rep['orchard']:
        line += f"\n    orchard violations at {list(rep['orchard'].keys())}: e.g. {list(rep['orchard'].items())[0]}"
    if rep['n_miquel']:
        line += f"\n    Miquel violations: {rep['n_miquel']}, e.g. cube faces {rep['miquel_example'][0]}"
    print(line, flush=True)
    r = None
    if not rep['kills'] and not args.no_realise:
        r = try_realise(F, L, tries=args.tries, verbose=True)
        print(f"    SURVIVES all filters; numerical realisation: realised={r['realised']}" + ("" if r['realised'] else f" best residual {r['best_residual']:.2e} (min dist {r['min_dist']:.2e})"), flush=True)
        alive.append(k)
    out.append({'index': k, 'blocks': sorted(sorted(B) for B in F), 'lines': sorted(sorted(S) for S in L), 'D': D, 'l': len(L),
                'count': 120 - D - len(L), 'report': rep, 'realise': r, 'sources': src[k]})
print("survivors of all filters:", alive)
if args.out:
    json.dump(out, open(args.out, 'w'), default=str, indent=0)
