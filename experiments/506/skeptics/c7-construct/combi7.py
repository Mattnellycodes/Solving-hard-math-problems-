"""Independent brute-force enumeration (labelled, no canonicity tricks) of all abstract Mobius block
structures on 7 points that could give <= 10 circles.
Structure: family F of blocks of size 4,5,6 pairwise sharing <=2 points (bitmasks).
count = 35 - D(F) - ell, D = sum (C(k,3)-1), ell = max number of 'lines' = blocks (from F or
uncovered triples) pairwise sharing <= 1 point.
Only the trivial necessary condition (pairwise intersections) is used; NO Sylvester-Gallai caps.
Then we list all F with 35 - D - ell_max <= 10."""
from itertools import combinations
from math import comb
n=7
pc=lambda x: bin(x).count('1')
cands=[m for m in range(1<<n) if 4<=pc(m)<=6]
cands.sort(key=lambda m:(pc(m),m))
deficit={m:comb(pc(m),3)-1 for m in cands}
ok={}
for a in cands:
    ok[a]=[b for b in cands if b>a and pc(a&b)<=2]
triples=[sum(1<<i for i in t) for t in combinations(range(n),3)]

def max_lines(F):
    # candidates: F blocks + triples not inside any block of F
    unc=[t for t in triples if not any(t&b==t for b in F)]
    L=list(F)+unc
    m=len(L)
    adj=[[pc(L[i]&L[j])<=1 for j in range(m)] for i in range(m)]
    best=0
    def rec(i,chosen,k):
        nonlocal best
        if k+ (m-i) <= best: return
        if i==m:
            best=max(best,k); return
        if all(adj[c][i] for c in chosen):
            chosen.append(i); rec(i+1,chosen,k+1); chosen.pop()
        rec(i+1,chosen,k)
    rec(0,[],0)
    return best

results=[]
nfam=0
def dfs(F,D,avail):
    global nfam
    nfam+=1
    if D>=35-10-7:   # ell <= 7 trivially (7 lines pairwise sharing <=1 on 7 points is the max, 21 pairs)
        ell=max_lines(F)
        cnt=35-D-ell
        if cnt<=10:
            results.append((sorted(F),D,ell,cnt))
    for idx,b in enumerate(avail):
        if D+sum(deficit[c] for c in avail[idx:])<35-10-7: break
        nav=[c for c in avail[idx+1:] if pc(b&c)<=2]
        dfs(F+[b],D+deficit[b],nav)
dfs([],0,cands)
print("families visited:",nfam)
print("structures with 35-D-ell_max <= 10:",len(results))
import collections
iso=collections.Counter()
for F,D,ell,cnt in results:
    sizes=tuple(sorted(pc(b) for b in F))
    iso[(sizes,D,ell,cnt)]+=1
for k,v in iso.items(): print(k,"labelled copies:",v)
# show one example of each
seen=set()
for F,D,ell,cnt in results:
    sizes=tuple(sorted(pc(b) for b in F))
    if sizes in seen: continue
    seen.add(sizes)
    print("example:",[ [i for i in range(n) if b>>i&1] for b in F],"D",D,"ell",ell,"count",cnt)
