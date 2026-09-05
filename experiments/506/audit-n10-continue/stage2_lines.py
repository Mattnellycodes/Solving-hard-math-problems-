"""Stage 2 (own DFS): for every F-class of a stage-1 file, enumerate all line sets L with
|L| >= 88 - D(F) (circles <= 32):
  rich lines = blocks of F, 3-lines = triples in no block of F, lines pairwise share <= 1 point,
  (SG10)  lines of P cover <= 44 pairs,
  (SGQ)   for each p: d_p(F) + sum_{l ∋ p} w(l) <= 44 with w(3-line) = 3, w(rich line of size k) = k-1
          (pairs of Q_p = (P-p) u {inf'} covered by its rich lines).
(F, L) classes are deduped by the canonical form of the Möbius structure with inf marked.
usage: python3 stage2_lines.py runs/stage1_<idx>.json [--lb=88]
"""
import json
import sys
import time
from itertools import combinations
from math import comb
from lib import mask, bits, popcount, canon, mobius_blocks, circles_count, check_structure

N = 10
fname = sys.argv[1]
TARGET = 88
for a in sys.argv[2:]:
    if a.startswith("--lb="):
        TARGET = int(a[5:])
data = json.load(open(fname))


def line_sets(F):
    D = sum(comb(popcount(b), 3) - 1 for b in F)
    need = TARGET - D
    covered = set()
    for b in F:
        for t in combinations(bits(b), 3):
            covered.add(mask(t))
    cands = list(F) + [mask(t) for t in combinations(range(N), 3) if mask(t) not in covered]
    w = [popcount(c) - 1 if popcount(c) >= 4 else 3 for c in cands]
    pairs = [comb(popcount(c), 2) for c in cands]
    dF = [sum(comb(popcount(b) - 1, 2) for b in F if b >> p & 1) for p in range(N)]
    compat = [[popcount(a & b) <= 1 for b in cands] for a in cands]
    out = []
    m = len(cands)

    def rec(start, chosen, tot, q):
        if len(chosen) + (m - start) < need:
            return
        if len(chosen) >= need:
            out.append(list(chosen))
        for i in range(start, m):
            if any(not compat[i][j] for j in chosen):
                continue
            if tot + pairs[i] > 44:
                continue
            ok = True
            for p in bits(cands[i]):
                if q[p] + w[i] > 44:
                    ok = False
                    break
            if not ok:
                continue
            for p in bits(cands[i]):
                q[p] += w[i]
            chosen.append(i)
            rec(i + 1, chosen, tot + pairs[i], q)
            chosen.pop()
            for p in bits(cands[i]):
                q[p] -= w[i]

    rec(0, [], 0, list(dF))
    return [[cands[i] for i in sol] for sol in out], need


t0 = time.time()
results = []
seen = {}
for ci, cl in enumerate(data["classes"]):
    F = [mask(b) for b in cl["F"]]
    sols, need = line_sets(F)
    ncls = 0
    for L in sols:
        check_structure(F, L)
        M = mobius_blocks(F, L, N) + [1 << N]
        cf = canon(M, N + 1)[0]
        if cf in seen:
            continue
        seen[cf] = True
        ncls += 1
        results.append({"F": [bits(b) for b in F], "L": [bits(l) for l in L], "D": sum(comb(popcount(b), 3) - 1 for b in F),
                        "l": len(L), "circles": circles_count(F, L), "F_class": ci})
    print(f"F-class {ci} (b4={cl['b4']}, |Aut|={cl.get('aut', '?')}): need l >= {need}: {len(sols)} labelled line sets, {ncls} new (F,L) classes",
          flush=True)
print(f"{fname}: {len(results)} (F,L) classes in total, {time.time() - t0:.1f}s")
json.dump(results, open(fname.replace("stage1_", "stage2_"), "w"))
