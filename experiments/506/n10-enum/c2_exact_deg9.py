#!/usr/bin/env python3
"""c2_exact_deg9.py -- same exact test as c2_exact_linesets.py but for line sets of size >= 9 (i.e. whether some
realisation of C-2 has a point on 9 blocks, which would give a planar realisation with exactly 33 circles)."""
import sys, json, time, itertools
sys.path.insert(0, '.')
import sympy as sp
from orderly2 import Params, line_sets, popcount
from c2_torsion import torsion_solutions
from c2_exact_linesets import test_line_set, test_line_set_infinity

data = json.load(open('n10_table_C.json')); rec = data['results'][2]
P = Params(10, 32, 'table')
masks = [sum(1 << i for i in b) for b in rec['blocks']]
ell, sets = line_sets(masks, P, 9)
sets = [frozenset(frozenset(i for i in range(10) if m >> i & 1) for m in L) for L in sets]
sets = sorted(set(sets), key=lambda L: (len(L), sorted(sorted(l) for l in L)))
print(f"line sets of size >= 9: {len(sets)} (ell_max = {ell})")
auts = [tuple(a) for a in json.load(open('c2_aut.json'))['auts']]
orbits = []; seen = set()
for L in sets:
    if L in seen: continue
    orb = {frozenset(frozenset(a[x] for x in l) for l in L) for a in auts}
    seen |= orb; orbits.append((L, len(orb)))
print(f"orbits under Aut(C-2): {len(orbits)}")
M, A, B, n2, n5, n10, good = torsion_solutions(rec)
A, B = list(A), list(B)
allok = True; t0 = time.time()
for k, (L, osz) in enumerate(orbits):
    Ll = [sorted(l) for l in L]
    res = []
    for sol in good:
        s_ang, t_ang = sol[:5], sol[5:]
        res.append(test_line_set(Ll, s_ang, t_ang, A, B, None, verbose=False))
        res.append(test_line_set_infinity(Ll, s_ang, t_ang, A, B, verbose=False))
    ok = all(res); allok &= ok
    print(f"orbit {k} (size {osz}, |L| = {len(L)}): {Ll} -> {'impossible for all labellings' if ok else 'NOT DECIDED'}  [{time.time() - t0:.0f}s]", flush=True)
print("RESULT:", "no realisation of C-2 has a point on >= 9 blocks (every planar realisation has >= 34 circles)" if allok else "SOME CASE NOT DECIDED")
