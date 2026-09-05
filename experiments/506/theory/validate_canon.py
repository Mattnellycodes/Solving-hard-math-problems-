"""(V1) Validate orderly generation: for n=6,7 count isomorphism classes of families of >=4-blocks
(pairwise sharing <=2 points, no block of size n) by (a) brute-force labelled enumeration + canonical
form via all permutations, (b) the Enumerator's DFS with D_min=0 (no pruning)."""
import itertools, sys
from mobius_enum import Enumerator, popcount
def brute(n, maxk):
    cands=[c for c in range(1<<n) if 4<=popcount(c)<=n-1]
    perms=list(itertools.permutations(range(n)))
    def img(c,s): return sum(1<<s[i] for i in range(n) if c>>i&1)
    def canon(F): return min(tuple(sorted(img(c,s) for c in F)) for s in perms)
    classes=set()
    def rec(F,last):
        if len(F)>maxk: return
        classes.add(canon(F))
        for b in cands:
            if b<=last: continue
            if all(popcount(a&b)<=2 for a in F): rec(F+[b],b)
    rec([],-1)
    return classes
for n in (5,6,7):
    E=Enumerator(n, target=10**6, verbose=False)   # D_min very negative -> no pruning by deficit
    E.D_min=-10**9
    # override process to just record
    seen=[]
    E.process=lambda F,D: seen.append(tuple(sorted(F)))
    E.run()
    # also verify each recorded family is canonical under brute force and all distinct
    perms=list(itertools.permutations(range(n)))
    def img(c,s): return sum(1<<s[i] for i in range(n) if c>>i&1)
    def canon(F): return min(tuple(sorted(img(c,s) for c in F)) for s in perms)
    assert all(canon(F)==F for F in seen), "non-canonical family output"
    assert len(set(seen))==len(seen)
    B=brute(n, 100)
    print(f"n={n}: orderly generation classes={len(seen)} brute-force classes={len(B)} match={set(seen)==B}")
