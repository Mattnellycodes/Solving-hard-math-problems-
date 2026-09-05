"""Exhaustive 10-subsets of the 5x5 grid (pure mode): stage 1 (|B(S)|, L, Euclidean count) and exact
stage 2 for all subsets with |B(S)| <= THRESH + 15  (deg <= 15 for 10 points: sum C(k_i,2) <= 45)."""
import sys, time, json
from itertools import combinations
from math import comb
from collections import Counter
import numpy as np
sys.path.insert(0, "/home/user/Solving-hard-math-problems-/experiments/506/audit-search-local")
from grid5_exhaustive import grid_blocks, GRID, N
from stage2_exact import analyse_subset

THRESH = int(sys.argv[1]) if len(sys.argv) > 1 else 44
t0 = time.time()
blocks = grid_blocks(False)
keys = list(blocks)
masks = np.array([blocks[k] for k in keys], dtype=np.uint32)
is_line = np.array([k[0] == 0 for k in keys])
m = 10
S = np.zeros(comb(N, m), dtype=np.uint32)
pos = 0
CH = 500000
buf = []
for c in combinations(range(N), m):
    buf.append(sum(1 << t for t in c))
    if len(buf) == CH:
        S[pos:pos + CH] = buf; pos += CH; buf = []
S[pos:pos + len(buf)] = buf
print(f"subsets {len(S)} t={time.time()-t0:.0f}s", flush=True)
hit3 = np.array([1 if k >= 3 else 0 for k in range(N + 1)], dtype=np.int32)
NB = np.zeros(len(S), dtype=np.int32); L = np.zeros(len(S), dtype=np.int32); E = np.zeros(len(S), dtype=np.int32)
for bi in range(len(keys)):
    k = np.bitwise_count(S & masks[bi]).astype(np.int32)
    h = hit3[k]
    NB += h
    if is_line[bi]: L += h
    else: E += h
print(f"stage 1 done t={time.time()-t0:.0f}s", flush=True)
hist = lambda a: sorted({int(v): int(c) for v, c in zip(*np.unique(a, return_counts=True))}.items())
print("Euclidean hist head:", hist(E)[:12])
print("NB hist head:", hist(NB)[:15])
print("min Euclid", E.min(), "count", int((E == E.min()).sum()), " min NB", NB.min())
cand = np.nonzero(NB <= THRESH + 15)[0]
print(f"candidates with NB <= {THRESH+15}: {len(cand)}", flush=True)
rows = []
for ci, idx in enumerate(cand.tolist()):
    pts = [GRID[t] for t in range(N) if (int(S[idx]) >> t) & 1]
    nb, l, md, P, bl = analyse_subset(pts, False)
    assert nb == NB[idx] and l == L[idx]
    best = max(md, l)
    rows.append((pts, nb, l, best, nb - best, sorted((len(v) for v in bl.values()), reverse=True)))
    if ci % 2000 == 0:
        print(f"  {ci}/{len(cand)} t={time.time()-t0:.0f}s", flush=True)
h = Counter(r[4] for r in rows)
print("histogram of per-subset min over all centres (complete for values <= %d):" % THRESH, sorted(h.items()))
print("(NB,best_deg) hist:", sorted(Counter((r[1], r[3]) for r in rows).items()))
best = min(r[4] for r in rows)
for r in rows:
    if r[4] == best:
        print("minimiser:", r[0], "NB", r[1], "L", r[2], "best_deg", r[3], "sizes", Counter(r[5]))
json.dump([{"pts": r[0], "NB": r[1], "L": r[2], "best_deg": r[3], "min_circles": r[4]} for r in rows if r[4] <= THRESH],
          open("stage2_n10_pure.json", "w"))
print(f"total t={time.time()-t0:.0f}s")
