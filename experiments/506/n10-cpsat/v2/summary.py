"""Collect the results of the pipeline into one overview (prints a table; robust to missing files)."""
import json, glob, os, sys
V = '/home/user/Solving-hard-math-problems-/experiments/506/n10-cpsat/v2'
def load(fn):
    try: return json.load(open(fn))
    except Exception: return None
sk = load(V + '/skeletons.json')
res = {}
for fn in (V + '/skel_sg_part1.json', V + '/skel_sg_part2.json'):
    d = load(fn)
    if d:
        for r in d: res[r['i']] = r
print("feasibility of D + l >= 88 per skeleton (mode sg, o = 1; E11, HSG, MK, O9, LPK; no O10):")
for i, s in enumerate(sk):
    r = res.get(i)
    print(f"  {i:2d} case {s['case']} k={s['k']} {s['blocks']}: {r['status'] if r else 'not run'}" + (f" ({r['time']:.0f}s)" if r else ''))
for fn in sorted(glob.glob(V + '/enum_sg_*.json')):
    d = load(fn)
    if d:
        print(f"enumeration {os.path.basename(fn)}: skeleton {d['skeleton']} status {d['status']}: {len(d['F_classes'])} 4-block-family classes, {len(d['structures'])} (F,L) classes; counts {sorted(set(s['count'] for s in d['structures']))}")
p = load(V + '/post_sg.json')
if p:
    alive = [r['index'] for r in p if not r['report']['kills']]
    print(f"post-processing: {len(p)} (F,L) classes; killed by filters: {len(p) - len(alive)}; survivors {alive}")
    for r in p:
        print(f"   #{r['index']}: count {r['count']} sizes {sorted(len(B) for B in r['blocks'])} lines {sorted(len(S) for S in r['lines'])} kills {r['report']['kills']}" + (f" realised={r['realise']['realised']}" if r.get('realise') else ''))
