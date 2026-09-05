"""All-4-block case for n = 10 (SG-only caps), own method:
D + ell >= 88 with ell <= 14 forces b4 >= 25, so sum of degrees >= 100 and some point has degree d0 in {10, 11}
(degree <= 11 by the derived cap 35).  Take point 0 of maximum degree; its derived structure is a partial triple
system on the other 9 points with d0 triples (pairwise sharing <= 1 point).  Step 1 enumerates all such PSTS(9, d0)
up to isomorphism (orderly generation, brute-force S_9 canonicity).  Step 2 extends each (blocks through 0 fixed,
symmetry group = Aut(PSTS) fixing 0, every other point of degree <= d0) by 4-blocks avoiding 0, exactly as in the
big-block case, and evaluates line sets exactly.  Results are deduplicated by brute-force S_10 canonical form."""
import sys, json, time, itertools
import numpy as np
from math import comb
from mob import *

def psts9(kmin, PT9):
    n = 9
    trip = sorted(mask_of(t) for t in itertools.combinations(range(n), 3))
    tests = [0]
    res = {}
    def bound(F, deg, last):
        cnt = [0] * n
        for m in trip:
            if m > last and all(popcount(m & b) <= 1 for b in F) and all(deg[p] < 4 for p in bits(m)):
                for p in bits(m): cnt[p] += 1
        return sum(min(cnt[p], 4 - deg[p]) for p in range(n)) // 3
    def rec(F, deg, aut):
        k = len(F)
        if k >= kmin:
            res.setdefault(k, []).append((list(F), aut))
        last = F[-1] if F else -1
        if k + bound(F, deg, last) < kmin:
            return
        for m in trip:
            if m <= last: continue
            if any(popcount(m & b) > 1 for b in F): continue
            if any(deg[p] >= 4 for p in bits(m)): continue
            F2 = F + [m]
            if F2[0] != 7: continue
            tests[0] += 1
            imgs = PT9.image_masks(F2)
            less, equal = rows_less_or_equal(imgs, F2)
            if less: continue
            deg2 = deg[:]
            for p in bits(m): deg2[p] += 1
            rec(F2, deg2, PT9.P[equal])
    rec([], [0] * n, PT9.P)
    print(f"PSTS(9,k) classes with k >= {kmin}: " + ", ".join(f"k={k}: {len(v)}" for k, v in sorted(res.items())) + f"  ({tests[0]} brute-force tests)", flush=True)
    return res

if __name__ == '__main__':
    n = 10; mode = sys.argv[1] if len(sys.argv) > 1 else 'sg'; target = int(sys.argv[2]) if len(sys.argv) > 2 else 32
    capD, capL = caps(n, mode); need = comb(n, 3) - target
    print(f"all-4-block case n={n} target={target} mode={mode} capD={capD} capL={capL} need={need}", flush=True)
    t0 = time.time()
    PT9 = PermTable(9)
    P = psts9(10, PT9)
    print(f"step 1 done [{time.time()-t0:.1f}s]", flush=True)
    allc = []
    tot = 0
    for d0 in sorted(P):
        if d0 > capD // 3:
            print(f"  d0={d0}: skipped (exceeds derived cap {capD}: a point of degree {d0} covers {3*d0} > {capD} pairs)"); continue
        for T, aut in P[d0]:
            fixed = [(m << 1) | 1 for m in T]  # relabel 0..8 -> 1..9, add point 0
            G = np.zeros((aut.shape[0], n), dtype=np.int8)
            G[:, 1:] = aut + 1
            S = Stage2(n, fixed, G, capD, capL, need, degcap=[d0] * n, verbose=False)
            cands = S.run(); tot += S.nodes
            print(f"  d0={d0} PSTS={[bits(m) for m in T]} |Aut|={aut.shape[0]} #cands={S.K} nodes={S.nodes} evals={S.evals} candidates={len(cands)}", flush=True)
            allc.extend(cands)
    print(f"step 2 done: total nodes {tot}, {len(allc)} candidate structures before dedupe [{time.time()-t0:.1f}s]", flush=True)
    PT = PermTable(10)
    classes = {}
    for c in allc:
        classes.setdefault(canon_full(PT, c['blocks']), c)
    res = list(classes.values())
    print(f"{len(res)} isomorphism classes with count <= {target}:")
    for i, c in enumerate(res):
        print(f"cand {i}: b4={len(c['blocks'])} D={c['D']} ell_max={c['ell_max']} count={c['count_min']} degrees={c['degrees']} #line_sets={len(c['line_sets'])}")
        print(f"   blocks={c['blocks']}")
    json.dump(dict(n=n, target=target, mode=mode, capD=capD, capL=capL, results=res), open(f"out_n{n}_t{target}_{mode}_all4only.json", 'w'))
    print(f"total time {time.time()-t0:.1f}s")
