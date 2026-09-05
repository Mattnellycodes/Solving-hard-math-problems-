"""Stage 2 (auditor's own): for every F-class of a stage-1 output, enumerate ALL line sets L with
|L| >= 88 - D(F):
  * a line of size >= 4 is a member of F; a 3-line is a triple contained in no member of F;
  * lines pairwise share <= 1 point;  * sum_S C(|S|,2) <= C(10,2) - 1 = 44 (Sylvester-Gallai).
Exhaustive DFS over candidate lines in a fixed order with the bound (remaining candidates).
Then (F, L) reduced to isomorphism classes and written to fl_<tag>.json.

usage: python3 lines.py f4_<tag>.json
"""
import sys, json, itertools, time
from math import comb
from common import N, ClassCollector, D_of, check_structure, count_of

src = sys.argv[1]
data = json.load(open(src))
t0 = time.time()
coll = ClassCollector()
tot_labelled = 0
for ci, cl in enumerate(data['classes']):
    F = [frozenset(B) for B in cl['F']]
    D = D_of(F)
    need = 88 - D
    cands = [B for B in F]
    for T in itertools.combinations(range(N), 3):
        Ts = frozenset(T)
        if not any(Ts <= B for B in F):
            cands.append(Ts)
    cands.sort(key=lambda s: (-len(s), sorted(s)))
    # conflict graph
    nc = len(cands)
    conflict = [[len(cands[i] & cands[j]) > 1 for j in range(nc)] for i in range(nc)]
    w = [comb(len(s), 2) for s in cands]
    sols = []

    def dfs(i, chosen, pairs):
        if len(chosen) + (nc - i) < need:
            return
        if i == nc:
            if len(chosen) >= need:
                sols.append(list(chosen))
            return
        # take i
        if pairs + w[i] <= 44 and all(not conflict[i][j] for j in chosen):
            chosen.append(i)
            dfs(i + 1, chosen, pairs + w[i])
            chosen.pop()
        dfs(i + 1, chosen, pairs)

    dfs(0, [], 0)
    tot_labelled += len(sols)
    for s in sols:
        L = [sorted(cands[i]) for i in s]
        Fl = [sorted(B) for B in F]
        check_structure(Fl, L)
        coll.add(Fl, L)
    print('F-class %d: D=%d need l>=%d: %d line sets  (%.0fs)' % (ci, D, need, len(sols), time.time() - t0), flush=True)
cl = coll.classes()
print('total labelled (F,L):', tot_labelled, ' (F,L)-classes:', len(cl))
out = []
for F, L, cnt in cl:
    print('  D=%d l=%d count=%d  L=%s' % (D_of(F), len(L), count_of(F, L), L))
    out.append({'F': F, 'L': L, 'D': D_of(F), 'l': len(L), 'count': count_of(F, L), 'labelled': cnt})
tag = src[3:-5]
json.dump({'skeleton': data['skeleton'], 'classes': out}, open('fl_%s.json' % tag, 'w'))
