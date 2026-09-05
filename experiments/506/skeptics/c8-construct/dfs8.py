"""Solver-free DFS enumeration (third method) of block families F on 8 points with
  (I1) pairwise |B∩B'| <= 2, sizes 4..7,
  (SG) for every point p: sum_{B∋p} C(|B|-1,2) <= 20   (Sylvester-Gallai at p, o>=1),
and D(F) >= 31, reporting D + ell_max >= 39 (circles <= 17).
Symmetry: any such family contains a 4-block (families of >=5-blocks only have D <= 19),
so WLOG {0,1,2,3} (code 15, the smallest code of any subset of size >= 4) is in F and is the first
block in code order.  Blocks are added in increasing code order.
Pruning: D + 0.34 * sum_p (20 - cov_p) < 31  => cut  (deficit per unit of point-coverage is
3/12, 9/30, 19/60, 34/105 for sizes 4,5,6,7, all <= 0.324).
ell_max: max family of pairwise <=1-intersecting members of F ∪ {uncovered triples}."""
import itertools, sys, time
from math import comb
n = 8
def pc(x): return bin(x).count("1")
cands = sorted([c for c in range(1 << n) if 4 <= pc(c) <= 7])
size = {c: pc(c) for c in cands}
deficit = {c: comb(pc(c), 3) - 1 for c in cands}
cov = {c: comb(pc(c) - 1, 2) for c in cands}
triples = [sum(1 << i for i in t) for t in itertools.combinations(range(n), 3)]
def ell_max(F):
    unc = [t for t in triples if not any(t & b == t for b in F)]
    items = sorted(F) + unc
    best = 0
    def rec(i, chosen):
        nonlocal best
        if len(chosen) > best: best = len(chosen)
        if len(chosen) + (len(items) - i) <= best: return
        for j in range(i, len(items)):
            c = items[j]
            if all(pc(c & d) <= 1 for d in chosen): rec(j + 1, chosen + [c])
    rec(0, [])
    return best
found = {}
nodes = 0
t0 = time.time()
def dfs(F, last, D, covp):
    global nodes
    nodes += 1
    if D >= 31:
        e = ell_max(F)
        if D + e >= 39:
            key = tuple(F)
            found[key] = (D, e)
            print(f"  found D={D} ell={e} circles={56-D-e} sizes={sorted(size[c] for c in F)} F={[[i for i in range(n) if c>>i&1] for c in F]}", flush=True)
    # potential bound
    slack = sum(20 - cp for cp in covp)
    if D + 0.34 * slack < 31: return
    for c in cands:
        if c <= last: continue
        if any(pc(c & b) >= 3 for b in F): continue
        ok = True
        newcov = list(covp)
        for p in range(n):
            if c >> p & 1:
                newcov[p] += cov[c]
                if newcov[p] > 20: ok = False; break
        if not ok: continue
        dfs(F + [c], c, D + deficit[c], newcov)
start = 15
dfs([start], start, deficit[start], [cov[start] if start >> p & 1 else 0 for p in range(n)])
print(f"nodes {nodes}, time {time.time()-t0:.1f}s, families with circles<=17: {len(found)}")
# isomorphism classes
perms = list(itertools.permutations(range(n)))
def img(c, s): return sum(1 << s[i] for i in range(n) if c >> i & 1)
def canon(F): return min(tuple(sorted(img(c, s) for c in F)) for s in perms)
classes = {}
for F, (D, e) in found.items():
    classes.setdefault(canon(F), (D, e, F))
print("iso classes:", len(classes))
for cf, (D, e, F) in classes.items():
    print(f"  D={D} ell={e} circles={56-D-e} blocks={[[i for i in range(n) if c>>i&1] for c in F]}")
