#!/usr/bin/env python3
"""Automorphism groups of the surviving structures and the orbits of their surviving line sets."""
import json, itertools, numpy as np
from math import comb
from oe import list_to_mask, mask_to_list, images, rows_less_eq, all_perms
from hcheck import check_structure, cumulative

d = json.load(open('n10_sg_t32.json')); n = d['n']
perms = all_perms(n)
def aut_group(fam):
    fam = sorted(fam); auts = []
    for s in range(0, len(perms), 400000):
        chunk = perms[s:s+400000]
        rows = np.sort(images(chunk, fam), axis=1)
        _, eq = rows_less_eq(rows, fam)
        auts.append(chunk[eq])
    return np.concatenate(auts)
for i in (4, 5, 6, 7, 8):
    r = d['results'][i]
    F = [list_to_mask(b) for b in r['blocks']]
    A = aut_group(F)
    rep = check_structure(n, r['blocks'], r['line_sets']); cum = cumulative(rep)
    print(f"structure {i}: |Aut(F)| = {len(A)}")
    for tier in ('B', 'C'):
        Ls = [frozenset(list_to_mask(l) for l in r['line_sets'][j]) for j in cum[tier]]
        # orbits under Aut(F)
        seen = set(); orbits = 0
        for L in Ls:
            if L in seen: continue
            orbits += 1
            for sigma in A:
                img = frozenset(int(x) for x in images(sigma[None, :], list(L))[0])
                seen.add(img)
        print(f"   tier {'A+B' if tier=='B' else 'A+B+C'}: {len(Ls)} surviving line sets, {orbits} orbits under Aut(F)")
