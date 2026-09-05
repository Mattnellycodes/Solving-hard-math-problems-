"""Own enumeration of skeletons: families of blocks of size >= 5 on 10 points, pairwise sharing
<= 2 points, under the per-point SG cap sum_{B ∋ p} C(|B|-1, 2) <= 35, up to isomorphism.
BFS over the number of blocks with canonical-form dedupe (lib.canon)."""
import json
from itertools import combinations
from math import comb
from lib import canon, mask, bits, popcount, automorphisms

N = 10
CAP = 35
cands = [mask(c) for k in range(5, 10) for c in combinations(range(N), k)]


def dcap_ok(blocks):
    for p in range(N):
        if sum(comb(popcount(b) - 1, 2) for b in blocks if b >> p & 1) > CAP:
            return False
    return True


level = {(): None}
classes = []
k = 0
while level:
    for fam in level:
        classes.append(fam)
    nxt = {}
    for fam in level:
        for c in cands:
            if c in fam:
                continue
            if any(popcount(c & b) > 2 for b in fam):
                continue
            new = list(fam) + [c]
            if not dcap_ok(new):
                continue
            cf, _, _ = canon(new, N)
            nxt[cf] = None
    level = nxt
    k += 1

out = []
for fam in classes:
    if not fam:
        continue
    sizes = sorted(popcount(b) for b in fam)
    aut = len(automorphisms(fam, N))
    D = sum(comb(popcount(b), 3) - 1 for b in fam)
    out.append({"sizes": sizes, "blocks": [bits(b) for b in fam], "D": D, "aut": aut,
                "intersections": sorted(popcount(a & b) for a, b in combinations(fam, 2))})
out.sort(key=lambda r: (len(r["sizes"]), r["sizes"], r["intersections"]))
for i, r in enumerate(out):
    print(i, r["sizes"], "D =", r["D"], "|Aut| =", r["aut"], "inters", r["intersections"], r["blocks"])
print("total non-empty skeleton classes:", len(out))
json.dump(out, open("runs/skeletons.json", "w"), indent=1)
