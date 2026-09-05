"""Numerical scan of the two-concentric-pentagon family R = ±lambda*S: for each lambda, compute all
blocks (rich blocks and 3-point circles/lines) of the 10-point Möbius set and the maximum number of blocks
through a point Z not in P (candidates: pairwise intersections of blocks, and infinity).  A structure with
<= 32 circles needs a Z on >= 10 blocks."""
import math, itertools, sys
import numpy as np
sys.path.insert(0, '/home/user/Solving-hard-math-problems-/experiments/506')
from explore_families import evaluate
best = []
for sign in (+1, -1):
    for lam in np.concatenate([np.linspace(0.05, 0.99, 380), np.linspace(1.01, 6.0, 500)]):
        S = [(math.cos(math.pi/2 + 2*math.pi*k/5), math.sin(math.pi/2 + 2*math.pi*k/5)) for k in range(5)]
        R = [(sign*lam*x, sign*lam*y) for x, y in S]
        res = evaluate(S + R)
        if res is None: continue
        bestc, eucl, nb, O = res[0], res[1], res[2], res[3]
        best.append((bestc, sign, lam, nb, eucl))
best.sort()
print("best (count, sign, lambda, nblocks, euclid) over the family:")
for b in best[:12]: print("  ", b)
print("minimum count over scanned family:", best[0][0])
