#!/usr/bin/env python3
"""iso_check.py -- match the candidate structures of two enumeration runs up to isomorphism (brute-force canonical
form over the permutations mapping the big blocks onto each other, via orderly2.Phase1)."""
import sys, json, itertools, numpy as np
sys.path.insert(0, '.')
from orderly2 import Params, Phase1, popcount
def canon(rec, ph1):
    n = rec['n']
    masks = [sum(1 << i for i in b) for b in rec['blocks']]
    big = [m for m in masks if popcount(m) >= 5]
    # candidates: all permutations mapping some smallest big block to {0..k-1}, then compare full sorted images
    S = ph1.relevant_perms(big)
    pw = np.int32(1) << S.astype(np.int32)
    M = np.stack([pw[:, [i for i in range(n) if m >> i & 1]].sum(axis=1) for m in masks], axis=1)
    M.sort(axis=1)
    idx = np.lexsort(M.T[::-1])
    return tuple(int(x) for x in M[idx[0]])
P = Params(10, 32, 'table'); ph1 = Phase1(P)
runs = {}
for name in sys.argv[1:]:
    d = json.load(open(name))
    runs[name] = [(i, canon(r, ph1), r['count_min'], r['ell_max']) for i, r in enumerate(d['results'])]
names = list(runs)
for i, c, cnt, ell in runs[names[0]]:
    matches = [(j, cnt2, ell2) for j, c2, cnt2, ell2 in runs[names[1]] if c2 == c]
    print(f"{names[0]}[{i}] (count {cnt}, ell {ell})  ==  {names[1]}{[m[0] for m in matches]} (count {[m[1] for m in matches]}, ell {[m[2] for m in matches]})")
