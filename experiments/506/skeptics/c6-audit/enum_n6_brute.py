"""Hostile-referee brute force for the combinatorial relaxation at n = 6 (and n = 5 as a control).

Written from scratch; no orderly generation, no pruning: we enumerate ALL labelled families F of
"rich blocks" (subsets of size 4..n-1, pairwise sharing <= 2 points) and, for each, ALL labelled
line sets L (subsets of size >= 3 that are blocks, i.e. members of F or triples contained in no
member of F, pairwise sharing <= 1 point).  count = C(n,3) - D(F) - |L|.

Three cap regimes are reported:
   none   : no Sylvester-Gallai information at all (o(m) = 0)
   weak   : Sylvester-Gallai only, o(m) >= 1 for m >= 3
   strong : o(5) = 4, o(6) = 3, o(7) = 3, o(8) = 4 (Kelly-Moser / Csima-Sawyer tables)
Caps used: for every p, sum_{B in F, p in B} C(|B|-1, 2) <= C(n-1,2) - o(n-1);
           sum_{S in L} C(|S|,2) <= C(n,2) - o(n).
We list every labelled (F, L) with count <= target, and the isomorphism classes under S_n.
"""
import sys
from itertools import combinations, permutations
from math import comb

O_STRONG = {3: 3, 4: 3, 5: 4, 6: 3, 7: 3, 8: 4, 9: 6}


def o_val(regime, m):
    if regime == "none":
        return 0
    if regime == "weak":
        return 1 if m >= 3 else 0
    return O_STRONG[m]


def enumerate_structures(n, target, regime):
    pts = range(n)
    rich = [frozenset(c) for k in range(4, n) for c in combinations(pts, k)]
    triples = [frozenset(c) for c in combinations(pts, 3)]
    cap_d = comb(n - 1, 2) - o_val(regime, n - 1)
    cap_l = comb(n, 2) - o_val(regime, n)
    N3 = comb(n, 3)
    results = []
    families = 0

    def ok_family(F):
        for p in pts:
            if sum(comb(len(B) - 1, 2) for B in F if p in B) > cap_d:
                return False
        return True

    def rec_F(F, start):
        nonlocal families
        families += 1
        if ok_family(F):
            handle(F)
        for j in range(start, len(rich)):
            B = rich[j]
            if all(len(B & A) <= 2 for A in F):
                rec_F(F + [B], j + 1)

    def handle(F):
        D = sum(comb(len(B), 3) - 1 for B in F)
        cands = list(F) + [t for t in triples if not any(t <= B for B in F)]
        # enumerate all line sets (pairwise <= 1 common point, pair budget cap_l)
        def rec_L(L, start, pairs):
            cnt = N3 - D - len(L)
            if cnt <= target:
                results.append((tuple(sorted(tuple(sorted(B)) for B in F)),
                                tuple(sorted(tuple(sorted(S)) for S in L)), D, len(L), cnt))
            for j in range(start, len(cands)):
                S = cands[j]
                if pairs + comb(len(S), 2) > cap_l:
                    continue
                if all(len(S & T) <= 1 for T in L):
                    rec_L(L + [S], j + 1, pairs + comb(len(S), 2))
        rec_L([], 0, 0)

    rec_F([], 0)
    return families, results


def canon(F, L, n):
    best = None
    for s in permutations(range(n)):
        img = (tuple(sorted(tuple(sorted(s[i] for i in B)) for B in F)),
               tuple(sorted(tuple(sorted(s[i] for i in S)) for S in L)))
        if best is None or img < best:
            best = img
    return best


if __name__ == "__main__":
    n = int(sys.argv[1]); target = int(sys.argv[2])
    for regime in ("none", "weak", "strong"):
        fams, res = enumerate_structures(n, target, regime)
        classes = {}
        for F, L, D, ell, cnt in res:
            classes.setdefault(canon(F, L, n), []).append((D, ell, cnt))
        print(f"n={n} target={target} regime={regime}: families of rich blocks scanned={fams}, "
              f"labelled (F,L) with count<={target}: {len(res)}, isomorphism classes: {len(classes)}")
        for key, lst in sorted(classes.items()):
            F, L = key
            print(f"   class: blocks={[list(b) for b in F]} lines={[list(s) for s in L]} "
                  f"D={lst[0][0]} ell={lst[0][1]} count={lst[0][2]} (labelled copies: {len(lst)})")
        # minimum count over everything
        allres = enumerate_structures(n, 10**9, regime)[1]
        mn = min(r[4] for r in allres)
        print(f"   minimum combinatorial count under regime {regime}: {mn}")
