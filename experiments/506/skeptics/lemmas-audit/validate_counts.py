"""lemmas-audit: independent Python brute force of the number of labelled families of >=4-subsets
pairwise sharing <= 2 points (no other constraint), to validate enum_own's raw DFS (capmode 0,
target = C(n,3) records every family, including the empty one)."""
import itertools, sys
def count(n):
    cands = [frozenset(c) for k in range(4, n) for c in itertools.combinations(range(n), k)]
    cands.sort(key=lambda s: sum(1 << i for i in s))
    tot = 0
    def rec(F, start):
        nonlocal tot
        tot += 1
        for j in range(start, len(cands)):
            if all(len(cands[j] & b) <= 2 for b in F):
                rec(F + [cands[j]], j + 1)
    rec([], 0)
    return tot
for n in (5, 6, 7):
    print(f"n={n}: labelled families (incl. empty) = {count(n)}")
