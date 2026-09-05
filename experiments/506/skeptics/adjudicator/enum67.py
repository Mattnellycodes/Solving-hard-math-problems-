"""Adjudicator's own exhaustive labelled enumeration of abstract Moebius structures for small n.

Structure = (F, L): F = family of subsets of size 4..n-1 (blocks of size >= 4) pairwise sharing
<= 2 points; L = set of "lines" chosen among F and the triples not inside any block of F, pairwise
sharing <= 1 point.  count = C(n,3) - D(F) - |L|,  D(F) = sum (C(|B|,3)-1).
Regimes: 'none' = only the intersection axioms;  'sg' = plus Sylvester-Gallai caps with o = 1
(at every point the blocks of size>=4 through it cover <= C(n-1,2)-1 pairs; lines cover
<= C(n,2)-1 pairs).  No orderly generation: plain DFS over all families, canonical forms by brute
force over S_n.  Usage: python3 enum67.py n target
"""
import sys
from itertools import combinations, permutations
from math import comb


def run(n, target, regime):
    pts = range(n)
    cands = [frozenset(c) for k in range(4, n) for c in combinations(pts, k)]
    triples = [frozenset(c) for c in combinations(pts, 3)]
    cap_pt = comb(n - 1, 2) - (1 if regime == 'sg' else 0)
    cap_ln = comb(n, 2) - (1 if regime == 'sg' else 0)
    results = []
    nfam = [0]

    def cov_ok(F):
        for p in pts:
            if sum(comb(len(B) - 1, 2) for B in F if p in B) > cap_pt:
                return False
        return True

    def ell_max(F):
        lc = list(F) + [t for t in triples if not any(t <= B for B in F)]
        lc.sort(key=len)
        best = [0, None]

        def rec(i, chosen, pairs):
            if len(chosen) > best[0]:
                best[0] = len(chosen); best[1] = list(chosen)
            for j in range(i, len(lc)):
                if len(chosen) + (len(lc) - j) <= best[0]:
                    return
                S = lc[j]
                if pairs + comb(len(S), 2) > cap_ln:
                    continue
                if all(len(S & T) <= 1 for T in chosen):
                    chosen.append(S); rec(j + 1, chosen, pairs + comb(len(S), 2)); chosen.pop()
        rec(0, [], 0)
        return best

    def dfs(start, F, D):
        nfam[0] += 1
        e, L = ell_max(F)
        cnt = comb(n, 3) - D - e
        if cnt <= target:
            results.append((cnt, D, e, sorted(sorted(B) for B in F), sorted(sorted(S) for S in L)))
        for j in range(start, len(cands)):
            B = cands[j]
            if all(len(B & A) <= 2 for A in F):
                F2 = F + [B]
                if cov_ok(F2):
                    dfs(j + 1, F2, D + comb(len(B), 3) - 1)
    dfs(0, [], 0)

    # canonical forms under S_n
    perms = list(permutations(pts))

    def canon(F, L):
        best = None
        for s in perms:
            key = (tuple(sorted(tuple(sorted(s[i] for i in B)) for B in F)),
                   tuple(sorted(tuple(sorted(s[i] for i in S)) for S in L)))
            if best is None or key < best:
                best = key
        return best
    classes = {}
    for cnt, D, e, F, L in results:
        classes.setdefault(canon(F, L), (cnt, D, e, F, L))
    print(f"n={n} target={target} regime={regime}: families visited={nfam[0]}, labelled (F,L_max) with count<=target: {len(results)}, iso classes: {len(classes)}")
    for key, (cnt, D, e, F, L) in sorted(classes.items(), key=lambda kv: kv[1][0]):
        print(f"   count={cnt} D={D} ell_max={e} blocks={F} lines(example max set)={L}")
    minimum = min((comb(n, 3) - D - e) for cnt, D, e, F, L in results) if results else None
    return classes


if __name__ == '__main__':
    n = int(sys.argv[1]); target = int(sys.argv[2])
    for regime in ('none', 'sg'):
        run(n, target, regime)
