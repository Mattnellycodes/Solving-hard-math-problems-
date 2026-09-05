#!/usr/bin/env python3
"""c2_exact_nocaps.py -- assumption-free version of c2_exact_deg9.py: the candidate sets of blocks through a
common point are ALL sets of >= 9 blocks of C-2 pairwise sharing <= 1 point (no Sylvester-Gallai or orchard caps
at all), reduced modulo Aut(C-2), and each orbit representative is tested exactly (finite O and O = infinity) for
every torsion labelling."""
import sys, json, time
sys.path.insert(0, '.')
from orderly2 import Params, line_sets, INF
from c2_torsion import torsion_solutions
from c2_exact_linesets import test_line_set, test_line_set_infinity

need = int(sys.argv[1]) if len(sys.argv) > 1 else 9
data = json.load(open('n10_table_C.json')); rec = data['results'][2]
P = Params(10, 32, 'none')            # no caps
masks = [sum(1 << i for i in b) for b in rec['blocks']]
ell, sets = line_sets(masks, P, need, max_sets=10 ** 6)
sets = [frozenset(frozenset(i for i in range(10) if m >> i & 1) for m in L) for L in sets]
sets = sorted(set(sets), key=lambda L: (len(L), sorted(sorted(l) for l in L)))
print(f"sets of >= {need} pairwise <=1-intersecting blocks (no caps): {len(sets)} (maximum size {ell})", flush=True)
auts = [tuple(a) for a in json.load(open('c2_aut.json'))['auts']]
orbits = []; seen = set()
for L in sets:
    if L in seen: continue
    orb = {frozenset(frozenset(a[x] for x in l) for l in L) for a in auts}
    seen |= orb; orbits.append((L, len(orb)))
print(f"orbits under Aut(C-2): {len(orbits)}", flush=True)
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
    print(f"orbit {k} (size {osz}, |L| = {len(L)}): {'impossible for all labellings' if ok else 'NOT DECIDED: ' + str(Ll)}  [{time.time() - t0:.0f}s]", flush=True)
print("RESULT:", f"no realisation of C-2 has a point on >= {need} blocks (no caps assumed)" if allok else "SOME CASE NOT DECIDED")
