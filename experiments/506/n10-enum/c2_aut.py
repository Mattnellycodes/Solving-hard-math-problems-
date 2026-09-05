#!/usr/bin/env python3
"""c2_aut.py -- automorphism group of the abstract structure C-2 (all block-preserving permutations of the 10 points)."""
import json, itertools, numpy as np
data = json.load(open('n10_table_C.json')); rec = data['results'][2]
blocks = [frozenset(b) for b in rec['blocks']]
big = [b for b in blocks if len(b) == 5]; fours = [b for b in blocks if len(b) == 4]
bset = set(blocks)
auts = []
# permutations preserving {A, B}: A -> A or A -> B
A, B = sorted(big[0]), sorted(big[1])
for swap in (False, True):
    tgtA, tgtB = (A, B) if not swap else (B, A)
    for pa in itertools.permutations(tgtA):
        for pb in itertools.permutations(tgtB):
            perm = {}
            for s, t in zip(A, pa): perm[s] = t
            for s, t in zip(B, pb): perm[s] = t
            if all(frozenset(perm[x] for x in b) in bset for b in fours):
                auts.append(tuple(perm[i] for i in range(10)))
print("|Aut(C-2)| =", len(auts))
# orbit structure on line sets of size >= 10
ls = [frozenset(frozenset(l) for l in L) for L in rec['line_sets']]
print("line sets of size >= 10:", len(ls), " sizes:", sorted(len(L) for L in ls))
orbits = []
seen = set()
for i, L in enumerate(ls):
    if i in seen: continue
    orb = set()
    for a in auts:
        img = frozenset(frozenset(a[x] for x in l) for l in L)
        for j, L2 in enumerate(ls):
            if L2 == img: orb.add(j)
    seen |= orb; orbits.append(sorted(orb))
print("orbits of Aut on the line sets:", orbits)
for orb in orbits:
    L = rec['line_sets'][orb[0]]
    print(f"   representative (size {len(L)}): {L}")
json.dump({"auts": auts, "orbits": orbits}, open('c2_aut.json', 'w'))
