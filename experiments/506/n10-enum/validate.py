#!/usr/bin/env python3
"""validate.py -- brute-force validation of the canonical-form machinery of orderly2.py.

For n = 6, 7 (and optionally 8 with restricted sizes) we enumerate ALL labelled families F of rich blocks
(sizes 4..n-1, pairwise sharing <= 2 points; optionally under the derived-SG caps of a mode), compute for
each its brute-force canonical form (lexicographically least sorted mask tuple over all n! permutations,
computed with a full permutation table -- a completely different computation from the two-phase form)
and collect the set of isomorphism classes.  Then we run the two-phase orderly generator with no deficit
pruning (target = huge) including the empty big family, map every visited (F_big, F4) pair to its
brute-force canonical form, and check that (i) the generator visits each class exactly once and (ii) the
sets of classes coincide.

Usage: python3 validate.py [n ...]   (default 5 6 7)
"""
import sys, itertools, time
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from orderly2 import Params, Phase1, Phase2, popcount, INF


def brute_classes(n, P):
    cands = [m for m in range(1 << n) if 4 <= popcount(m) <= n - 1]
    perms = np.array(list(itertools.permutations(range(n))), dtype=np.int64)
    pw = np.int64(1) << perms                                                # (n!, n)
    pts = {m: [i for i in range(n) if m >> i & 1] for m in cands}
    tab = np.zeros((len(perms), 1 << n), dtype=np.int64)
    for m in cands:
        tab[:, m] = pw[:, pts[m]].sum(axis=1)
    from math import comb
    inc = {m: comb(popcount(m) - 1, 2) for m in cands}

    def canon(F):
        M = tab[:, F]
        M.sort(axis=1)
        # lexicographic minimum row
        idx = np.lexsort(M.T[::-1])
        return tuple(int(x) for x in M[idx[0]])

    classes = {}
    count_labelled = [0]

    def rec(F, last, cov):
        count_labelled[0] += 1
        key = canon(F) if F else ()
        classes[key] = classes.get(key, 0) + 1
        for b in cands:
            if b <= last:
                continue
            if any(popcount(a & b) > 2 for a in F):
                continue
            if any(cov[p] + inc[b] > P.cap_derived for p in pts[b]):
                continue
            cov2 = cov[:]
            for p in pts[b]:
                cov2[p] += inc[b]
            rec(F + [b], b, cov2)

    rec([], -1, [0] * n)
    return classes, count_labelled[0], canon


def generator_classes(n, P, canon):
    ph1 = Phase1(P)
    fams = ph1.run()
    seen = {}
    nodes = 0
    for Fbig, cov, D0 in fams:
        aut = ph1.automorphisms(Fbig)
        ph2 = Phase2(P, Fbig, cov, D0, aut)
        res = ph2.run()
        nodes += ph2.nodes
        for r in res:
            F = list(r["Fbig"]) + list(r["F4"])
            key = canon(F) if F else ()
            seen[key] = seen.get(key, 0) + 1
    return seen, nodes, len(fams)


if __name__ == "__main__":
    ns = [int(x) for x in sys.argv[1:]] or [5, 6, 7]
    for n in ns:
        for mode in ("none", "table"):
            P = Params(n, 10 ** 6, mode)            # target huge -> no deficit pruning
            P.ell_max = INF; P.D_min = -INF        # make sure no pruning at all in phase 2
            if mode == "none":
                P.cap_derived = INF; P.max4 = INF
            t0 = time.time()
            classes, labelled, canon = brute_classes(n, P)
            t1 = time.time()
            seen, nodes, nfam = generator_classes(n, P, canon)
            t2 = time.time()
            dup = [k for k, v in seen.items() if v != 1]
            ok = (set(seen) == set(classes)) and not dup
            print(f"n={n} mode={mode} (cap_derived={P.cap_derived}, max4={P.max4}): labelled families={labelled}, "
                  f"brute-force classes={len(classes)} [{t1 - t0:.1f}s]; generator: phase-1 families={nfam}, "
                  f"phase-2 nodes={nodes}, distinct classes={len(seen)}, duplicates={len(dup)} [{t2 - t1:.1f}s]  "
                  f"=> {'MATCH' if ok else 'MISMATCH'}", flush=True)
            if not ok:
                missing = set(classes) - set(seen)
                extra = set(seen) - set(classes)
                print("   missing:", list(missing)[:5], " extra:", list(extra)[:5], " dup:", dup[:5])
