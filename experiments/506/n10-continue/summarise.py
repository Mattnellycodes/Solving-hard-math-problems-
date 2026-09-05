"""Per-skeleton summary table from runs/enum_*.json and runs/post_all.json."""
import json, glob, os, sys, math
sys.path.insert(0, '/home/user/Solving-hard-math-problems-/experiments/506/n10-continue')
from common import ClassStore, automorphisms, fs
sk = json.load(open('skeletons_own.json'))
post = json.load(open('runs/post_all.json'))
pst = ClassStore()
pidx = {}
for r in post:
    i, _ = pst.add(r['blocks'], r['lines']); pidx[i] = r
print("| own idx | v2 idx | big blocks | D_skel | b4_min | Aut | stage 1 | labelled | F-classes | (F,L) classes | kills of the (F,L) classes (theorems only) |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
v2 = json.load(open('/home/user/Solving-hard-math-problems-/experiments/506/n10-cpsat/v2/skeletons.json'))
v2st = ClassStore(); v2map = {}
for j, r in enumerate(v2):
    i, _ = v2st.add(r['blocks'], []); v2map[i] = j
for i, rec in enumerate(sk):
    fn = f'runs/enum_{i}.json'
    vi, new = v2st.add(rec['blocks'], []); v2i = '-' if new else str(v2map[vi])
    if not os.path.exists(fn):
        print(f"| {i} | {v2i} | {rec['sizes']} | | | | excluded (block >= 7: Lemma A/B/C) | | | | |"); continue
    d = json.load(open(fn))
    D = sum(math.comb(len(B), 3) - 1 for B in rec['blocks']); b4min = max(0, math.ceil((88 - 14 - D) / 3))
    aut = len(automorphisms([fs(B) for B in rec['blocks']]))
    lab = sum(c['copies'] for c in d['F_classes'])
    kills = []
    for s in d['structures']:
        j, _ = pst.add(s['blocks'], s['lines'])
        kills.append('+'.join(pidx[j]['theorem_kills']) if pidx[j]['theorem_kills'] else 'SURVIVES->exact')
    from collections import Counter
    kc = ', '.join(f"{k} x{v}" for k, v in sorted(Counter(kills).items()))
    print(f"| {i} | {v2i} | {rec['sizes']} | {D} | {b4min} | {aut} | {d['stage1_status']} | {lab} | {len(d['F_classes'])} | {len(d['structures'])} | {kc} |")
