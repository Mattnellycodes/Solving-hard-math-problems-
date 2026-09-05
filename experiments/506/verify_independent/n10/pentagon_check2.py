import math, itertools, numpy as np, sys
from collections import defaultdict, Counter
sys.path.insert(0,'.')
from pentagon_check import circ, blocks, pent
FAM = {
 'continue': [[0,1,6,9],[0,1,7,8],[0,2,5,9],[0,2,6,8],[0,3,5,7],[0,3,8,9],[0,4,5,6],[0,4,7,9],[1,2,5,8],[1,2,7,9],[1,3,5,9],[1,3,6,7],[1,4,5,7],[1,4,6,8],[2,3,5,6],[2,3,7,8],[2,4,6,7],[2,4,8,9],[3,4,5,8],[3,4,6,9]],
 'enum5':    [[0,1,5,6],[0,1,7,8],[0,2,5,7],[0,2,6,9],[0,3,5,8],[0,3,7,9],[0,4,5,9],[0,4,6,8],[1,2,6,7],[1,2,8,9],[1,3,5,9],[1,3,6,8],[1,4,5,7],[1,4,6,9],[2,3,5,6],[2,3,7,8],[2,4,5,8],[2,4,7,9],[3,4,6,7],[3,4,8,9]],
}
for name,F in FAM.items():
    deg=Counter(x for b in F for x in b); pairs=Counter(frozenset(p) for b in F for p in itertools.combinations(b,2))
    print(name, "blocks", len(F), "degrees", sorted(deg.values()), "max pair multiplicity", max(pairs.values()),
          "triples covered twice?", any(v>1 for v in Counter(frozenset(t) for b in F for t in itertools.combinations(b,3)).values()))
for phi in (0.0, math.pi/5):
    pts=pent(0.61,phi); B=blocks(pts); four=sorted([sorted(b) for b in B if len(b)==4])
    print(f"\nphi={phi:.3f}: geometric 4-blocks (0-4 outer, 5-9 inner):", four)
    deg=Counter(x for b in four for x in b); print("   degrees", sorted(deg.values()))
    G=set(frozenset(b) for b in four)
    for name,F in FAM.items():
        Fs=set(frozenset(b) for b in F)
        for orient in ('outer->S','outer->R'):
            found=0
            for p in itertools.permutations(range(5)):
                for q in itertools.permutations(range(5)):
                    if orient=='outer->S': lab={k:p[k] for k in range(5)}; lab.update({5+k:5+q[k] for k in range(5)})
                    else: lab={k:5+p[k] for k in range(5)}; lab.update({5+k:q[k] for k in range(5)})
                    if set(frozenset(lab[x] for x in b) for b in G)==Fs: found+=1
            print(f"   {name} {orient}: {found} isomorphisms")
