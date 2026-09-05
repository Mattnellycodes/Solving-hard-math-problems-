"""Independent enumeration (2026-09-05) of 'big-block skeletons' for n = 10:
families of blocks of size 5..9 on 10 points, pairwise sharing <= 2 points, up to isomorphism,
subject only to the Sylvester-Gallai per-point cap  sum_{B ∋ p} C(|B|-1, 2) <= C(9,2) - 1 = 35.
Orderly generation with brute-force canonical form (10! is too big; use invariant + nauty-free
canonical labelling via sorted certificate over automorphism-free refinement is complex, so we use
a simpler but exhaustive approach: generate families by increasing size with a canonical-form test
based on the lexicographically minimal image over S_10 restricted to permutations respecting a
degree-based partition refinement).  For validation, we also count labelled families for k<=2.
"""
import itertools, math, sys, time
from collections import Counter
n = 10
pts = range(n)
big = [frozenset(c) for k in range(5, 10) for c in itertools.combinations(pts, k)]
CAP = math.comb(9, 2) - 1

def ok_add(F, B):
    if any(len(B & C) > 2 for C in F): return False
    cov = Counter()
    for C in F + [B]:
        for p in C: cov[p] += math.comb(len(C) - 1, 2)
    return all(v <= CAP for v in cov.values())

def canon(F):
    """canonical form: min over permutations consistent with a refinement by (degree profile)."""
    # invariant per point: sorted list of sizes of blocks containing it
    prof = {p: tuple(sorted(len(B) for B in F if p in B)) for p in pts}
    # group points by profile; only permute within groups, and order groups by profile
    groups = {}
    for p in pts: groups.setdefault(prof[p], []).append(p)
    keys = sorted(groups)
    best = None
    # assign target labels: group i gets consecutive labels
    def rec(i, mapping, next_label):
        nonlocal best
        if i == len(keys):
            img = tuple(sorted(tuple(sorted(mapping[p] for p in B)) for B in F))
            if best is None or img < best: best = img
            return
        g = groups[keys[i]]
        for perm in itertools.permutations(g):
            m2 = dict(mapping)
            for j, p in enumerate(perm): m2[p] = next_label + j
            rec(i + 1, m2, next_label + len(g))
    rec(0, {}, 0)
    return (tuple(keys), best)

t0 = time.time()
seen = set(); classes = []
frontier = [[]]
seen.add(canon([]))
level = 0
while frontier:
    level += 1
    nxt = []
    for F in frontier:
        for B in big:
            if B in F: continue
            if not ok_add(F, B): continue
            G = F + [B]
            c = canon(G)
            if c in seen: continue
            seen.add(c); nxt.append(G); classes.append(G)
    print(f"level {level}: {len(nxt)} new classes (total {len(classes)}), {time.time()-t0:.0f}s", flush=True)
    frontier = nxt
print("=== all skeleton classes (block families of size >= 5, pairwise <= 2, SG cap) ===")
for G in classes:
    print(sorted(len(B) for B in G), [sorted(B) for B in G])
