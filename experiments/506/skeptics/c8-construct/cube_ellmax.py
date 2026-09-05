"""Brute-force: the cube structure on 8 points (12 four-blocks = 6 faces + 6 diagonal planes,
8 corner triples), what is the maximum number of blocks pairwise sharing <= 1 point?
Also: is the block set really the one from an actual cube (verify from coordinates)."""
import itertools
V = list(itertools.product([0,1],repeat=3))
idx = {v:i for i,v in enumerate(V)}
# all planes through >=3 vertices
import numpy as np
blocks=set()
for t in itertools.combinations(range(8),3):
    p=np.array([V[i] for i in t],float)
    nrm=np.cross(p[1]-p[0],p[2]-p[0])
    memb=frozenset(i for i in range(8) if abs(np.dot(np.array(V[i],float)-p[0],nrm))<1e-9)
    blocks.add(memb)
blocks=sorted(blocks,key=lambda b:(len(b),sorted(b)))
print(len(blocks),"blocks; sizes",sorted(len(b) for b in blocks))
four=[b for b in blocks if len(b)==4]; three=[b for b in blocks if len(b)==3]
assert all(len(a&b)<=2 for a,b in itertools.combinations(blocks,2))
best=0;bestsets=[]
for k in range(1,9):
    for S in itertools.combinations(blocks,k):
        if all(len(a&b)<=1 for a,b in itertools.combinations(S,2)):
            if k>best: best=k;bestsets=[]
            if k==best: bestsets.append([sorted(b) for b in S])
print("ell_max =",best,"; number of maximum line-sets:",len(bestsets))
print("circles = 56 - 36 - ell_max =",56-36-best)
# also: with only the 12 four-blocks fixed but the 8 triples free (any triple not inside a 4-block)?  In the
# cube structure the uncovered triples are exactly the 8 corner triples, so nothing else is possible.
unc=[frozenset(t) for t in itertools.combinations(range(8),3) if not any(set(t)<=b for b in four)]
print("uncovered triples:",len(unc),"== corner triples?",set(unc)==set(three))
