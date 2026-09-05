"""Stage 1 (own): for a labelled skeleton (blocks of size >= 5), enumerate ALL families of
4-blocks F4 such that F = skel u F4 satisfies
  (T)   every triple in <= 1 block,
  (SG9) per point p: sum_{B ∋ p} C(|B|-1,2) <= 35,
  (MK)  per point p and 8-subset E of P-p: #{B ∋ p : |(B-p) ∩ E| = 3} <= 7   [t3(8) <= 7, from (8_3)]
  (CAP, optional) per point: #4-blocks through p <= 10                       [t3(9) <= 10]
  (LB)  D(F) >= LB (default 75: l <= 13 because every point is on <= 4 lines)
up to Aut(skel), using a SAT solver (pysat / CaDiCaL) with cardinality encodings, enumerating
b4 from its maximum downwards, and blocking every Aut(skel)-image of each solution (superset
blocking within a fixed b4 round).  Output: one representative per isomorphism class of F.

usage: python3 stage1_sat.py <skeleton index in runs/skeletons.json | 'empty'> [--nocap] [--nomk] [--lb=74] [--tag=x]
"""
import json
import sys
import time
from itertools import combinations
from math import comb
from pysat.solvers import Solver
from pysat.card import CardEnc, EncType
from lib import mask, bits, popcount, automorphisms, canon

N = 10
args = sys.argv[1:]
idx = args[0]
use_cap = "--nocap" not in args
use_mk = "--nomk" not in args
LB = 75
tag = ""
solver_name = "cadical153"
for a in args:
    if a.startswith("--solver="):
        solver_name = a[9:]
    if a.startswith("--lb="):
        LB = int(a[5:])
    if a.startswith("--tag="):
        tag = a[6:]

if idx == "empty":
    skel = []
else:
    skel = [mask(b) for b in json.load(open("runs/skeletons.json"))[int(idx)]["blocks"]]
D_skel = sum(comb(popcount(b), 3) - 1 for b in skel)
b4min = max(0, -(-(LB - D_skel) // 3))
t0 = time.time()

# candidates
cands = [mask(c) for c in combinations(range(N), 4) if all(popcount(mask(c) & s) <= 2 for s in skel)]
var = {c: i + 1 for i, c in enumerate(cands)}
top = len(cands)
base = []
# (T)
for a, b in combinations(cands, 2):
    if popcount(a & b) == 3:
        base.append([-var[a], -var[b]])


def atmost(lits, k):
    global top
    if k < 0:
        return [[]]  # unsatisfiable
    if not lits or k >= len(lits):
        return []
    if k == 0:
        return [[-l] for l in lits]
    enc = CardEnc.atmost(lits=lits, bound=k, top_id=top, encoding=EncType.seqcounter)
    top = max(top, enc.nv)
    return enc.clauses


def atleast(lits, k):
    global top
    if k <= 0:
        return []
    if k > len(lits):
        return [[]]
    enc = CardEnc.atleast(lits=lits, bound=k, top_id=top, encoding=EncType.seqcounter)
    top = max(top, enc.nv)
    return enc.clauses


for p in range(N):
    through = [var[c] for c in cands if c >> p & 1]
    d = sum(comb(popcount(s) - 1, 2) for s in skel if s >> p & 1)
    base += atmost(through, (35 - d) // 3)
    if use_cap:
        base += atmost(through, 10)
    if use_mk:
        others = [q for q in range(N) if q != p]
        for E8 in combinations(others, 8):
            E = mask(E8)
            const = sum(1 for s in skel if s >> p & 1 and popcount((s & ~(1 << p)) & E) == 3)
            lits = [var[c] for c in cands if c >> p & 1 and popcount((c & ~(1 << p)) & E) == 3]
            base += atmost(lits, 7 - const)

b4max = len(cands)
# cheap upper bound on b4 from the caps: sum over points of cap / 4
capsum = 0
for p in range(N):
    d = sum(comb(popcount(s) - 1, 2) for s in skel if s >> p & 1)
    c = (35 - d) // 3
    if use_cap:
        c = min(c, 10)
    capsum += c
b4max = min(b4max, capsum // 4, (120 - sum(comb(popcount(s), 3) for s in skel)) // 4)

aut = automorphisms(skel, N) if skel else None
print(f"skeleton {idx}: sizes {sorted(popcount(s) for s in skel)} D_skel={D_skel} b4 in [{b4min},{b4max}] "
      f"cands={len(cands)} base clauses={len(base)} |Aut|={len(aut) if aut else 'S10'} cap={use_cap} mk={use_mk} LB={LB} solver={solver_name}",
      flush=True)
# permutation action on candidate indices
if aut:
    cidx = {c: i for i, c in enumerate(cands)}
    perm_maps = []
    for g in aut:
        m = [0] * len(cands)
        for i, c in enumerate(cands):
            m[i] = cidx[mask(g[p] for p in bits(c))]
        perm_maps.append(m)

classes = []
total_labelled = 0
for k in range(b4max, b4min - 1, -1):
    solver = Solver(name=solver_name, bootstrap_with=base)
    saved_top = top
    lits = [var[c] for c in cands]
    for cl in atmost(lits, k) + atleast(lits, k):
        solver.add_clause(cl)
    nsol = 0
    nlab = 0
    while solver.solve():
        model = solver.get_model()
        sol = [i for i, c in enumerate(cands) if model[var[c] - 1] > 0]
        assert len(sol) == k
        nsol += 1
        F = skel + [cands[i] for i in sol]
        cf, _, na = canon(F, N)
        classes.append({"b4": k, "F": [bits(b) for b in F], "canon": [bits(b) for b in cf], "aut": na})
        if aut:
            images = {tuple(sorted(m[i] for i in sol)) for m in perm_maps}
        else:
            images = {tuple(sorted(sol))}
        nlab += len(images)
        for img in images:
            solver.add_clause([-var[cands[i]] for i in img])
    solver.delete()
    top = saved_top
    total_labelled += nlab
    print(f"  b4={k}: {nsol} classes, {nlab} labelled solutions  ({time.time() - t0:.1f}s)", flush=True)

# sanity: class representatives pairwise non-isomorphic
cfs = [tuple(mask(b) for b in c["canon"]) for c in classes]
assert len(set(cfs)) == len(cfs), "duplicate classes (orbit blocking failed?)"
print(f"skeleton {idx}: {len(classes)} F-classes, {total_labelled} labelled solutions, {time.time() - t0:.1f}s")
name = f"runs/stage1_{idx}{tag}.json"
json.dump({"skeleton": [bits(b) for b in skel], "D_skel": D_skel, "b4min": b4min, "cap": use_cap, "mk": use_mk,
           "LB": LB, "classes": classes, "labelled": total_labelled}, open(name, "w"))
