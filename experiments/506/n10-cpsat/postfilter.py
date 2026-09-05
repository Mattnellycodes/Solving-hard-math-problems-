"""Apply the necessary conditions (filters.py) and the numerical realiser to a JSON list of structures.
Usage: python3 postfilter.py file.json [--key survivors] [--tries N] [--out report.json]
"""
import sys, json, argparse, itertools, math
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from filters import full_report, hereditary_sg, three_line_excess, miquel_violations, bundle_violations, fregier_violations
from realise import try_realise
from n10model import O_MODES

ap = argparse.ArgumentParser()
ap.add_argument('file'); ap.add_argument('--key', default=None); ap.add_argument('--tries', type=int, default=100)
ap.add_argument('--out', default=None); ap.add_argument('--no-realise', action='store_true')
args = ap.parse_args()
data = json.load(open(args.file))
if args.key:
    data = data[args.key]
out = []
for i, rec in enumerate(data):
    F = [frozenset(B) for B in rec['blocks']]; L = [frozenset(S) for S in rec['lines']]
    D = sum(math.comb(len(B), 3) - 1 for B in F)
    # which o-modes admit this structure (SG caps)
    modes = []
    for name, o in O_MODES.items():
        okp = all(sum(math.comb(len(B) - 1, 2) for B in F if p in B) <= 36 - o[9] for p in range(10))
        okl = sum(math.comb(len(S), 2) for S in L) <= 45 - o[10]
        if okp and okl:
            modes.append(name)
    print(f"structure {i}: sizes={sorted(len(B) for B in F)} lines={sorted(len(S) for S in L)} D={D} l={len(L)} count={120-D-len(L)} modes={modes}")
    rep = full_report(F, L)
    rep['modes'] = modes
    bv = bundle_violations(F, L); fv = fregier_violations(F, L)
    rep['n_bundle'] = len(bv); rep['n_fregier'] = len(fv); rep['fregier_example'] = fv[:1]
    kills = [k for k, c in (('hereditarySG', bool(rep['hereditary_sg'])), ('Fano/8_3', bool(rep['three_line_excess'])),
             ('Miquel', rep['n_miquel'] > 0), ('bundle', len(bv) > 0), ('Fregier', len(fv) > 0)) if c]
    rep['kills'] = kills
    print(f"   bundle violations: {len(bv)}  Fregier-component violations: {len(fv)}  => killed by {kills}")
    dead = bool(kills)
    rep['dead'] = dead
    if not dead and not args.no_realise:
        r = try_realise(F, L, tries=args.tries, verbose=True)
        rep['realise'] = r
        print(f"   SURVIVES all combinatorial filters; numerical realisation: {r}")
    out.append({'index': i, 'blocks': rec['blocks'], 'lines': rec['lines'], 'report': rep})
if args.out:
    json.dump(out, open(args.out, 'w'), default=str)
print("dead:", sum(1 for o in out if o['report']['dead']), "alive:", sum(1 for o in out if not o['report']['dead']))
