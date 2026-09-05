"""Independent (SAT-free) enumeration for the skeleton of two disjoint 5-blocks S = {0..4},
R = {5..9}.  Every 4-block meets S and R in 2 points each (3 points in S would share 3 points with
the block S).  Two 4-blocks with the same S-pair must have disjoint R-pairs (else a common triple),
and vice versa.  So a 4-block family = set of (S-pair, R-pair) incidences with every S-pair in <= 2
incidences with disjoint R-pairs and every R-pair in <= 2 incidences with disjoint S-pairs.
D = 18 + 3 b4 >= 75 needs b4 >= 19, so b4 in {19, 20}.
Symmetry breaking (valid): some S-pair is in 2 incidences; by S5 x S5 we may take it to be {0,1}
with R-pairs {5,6}, {7,8}.  Then canonical dedupe.  Filters applied afterwards: SG cap (automatic:
d_p = 6 + 3*8 = 30 <= 35), (MK) t3(8) <= 7 on every derived 8-subset, (CAP) <= 10 four-blocks/point."""
import json
import time
from itertools import combinations
from lib import mask, bits, popcount, canon, mk_violation, mobius_blocks, derived

S = list(range(5))
R = list(range(5, 10))
Sp = [frozenset(c) for c in combinations(S, 2)]
Rp = [frozenset(c) for c in combinations(R, 2)]
# options for a row: 0, 1, or 2 pairwise disjoint R-pairs
row_opts = [()] + [(i,) for i in range(10)] + [(i, j) for i, j in combinations(range(10), 2) if not (Rp[i] & Rp[j])]
t0 = time.time()
sols = []
col_use = [[] for _ in range(10)]  # list of S-pair indices per R-pair


def rec(r, chosen, total):
    if total + 2 * (10 - r) < 19:
        return
    if r == 10:
        if total >= 19:
            sols.append(list(chosen))
        return
    opts = row_opts
    if r == 0:
        opts = [(Rp.index(frozenset((5, 6))), Rp.index(frozenset((7, 8))))]
    for opt in opts:
        ok = True
        for j in opt:
            if len(col_use[j]) >= 2 or any(Sp[i] & Sp[r] for i in col_use[j]):
                ok = False
                break
        if not ok:
            continue
        for j in opt:
            col_use[j].append(r)
        chosen.append(opt)
        rec(r + 1, chosen, total + len(opt))
        chosen.pop()
        for j in opt:
            col_use[j].pop()


rec(0, [], 0)
print(f"labelled families (with the symmetry-breaking normalisation): {len(sols)}  ({time.time() - t0:.1f}s)")
skel = [mask(S), mask(R)]
classes = {}
for sol in sols:
    F4 = [mask(Sp[i] | Rp[j]) for i, opt in enumerate(sol) for j in opt]
    F = skel + F4
    cf = canon(F, 10)[0]
    if cf in classes:
        continue
    # filters
    cap_ok = all(sum(1 for b in F4 if b >> p & 1) <= 10 for p in range(10))
    mk_ok = True
    for p in range(10):
        lines = derived(F, p)
        if mk_violation(lines, ((1 << 10) - 1) & ~(1 << p)) is not None:
            mk_ok = False
            break
    classes[cf] = {"F": [bits(b) for b in F], "b4": len(F4), "cap_ok": cap_ok, "mk_ok": mk_ok}
n_all = len(classes)
n_mk = sum(1 for c in classes.values() if c["mk_ok"] and c["cap_ok"])
print(f"F-classes (T only): {n_all}; by b4: { {k: sum(1 for c in classes.values() if c['b4'] == k) for k in (19, 20)} }")
print(f"F-classes surviving (MK)+(CAP): {n_mk}; by b4: { {k: sum(1 for c in classes.values() if c['b4'] == k and c['mk_ok'] and c['cap_ok']) for k in (19, 20)} }")
print(f"time {time.time() - t0:.1f}s")
json.dump({"classes": [c for c in classes.values() if c["mk_ok"] and c["cap_ok"]],
           "n_all": n_all}, open("runs/two5_dfs.json", "w"))
