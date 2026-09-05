import math, sys, itertools
import numpy as np
sys.path.insert(0, '/home/user/Solving-hard-math-problems-/experiments/506')
sys.path.insert(0, '.')
from explore_families import evaluate
from realise import structure_of
def config(sign, lam):
    S = [(math.cos(math.pi/2 + 2*math.pi*k/5), math.sin(math.pi/2 + 2*math.pi*k/5)) for k in range(5)]
    R = [(sign*lam*x, sign*lam*y) for x, y in S]
    return S + R
def invert(pts, O, r2=1.0):
    out = []
    for x, y in pts:
        dx, dy = x - O[0], y - O[1]; d2 = dx*dx + dy*dy
        out.append((O[0] + r2*dx/d2, O[1] + r2*dy/d2))
    return out
for sign, lam in [(-1, 0.10456464379947229), (-1, 0.9676781002638523), (-1, 0.3848284960422163), (+1, 0.5)]:
    pts = config(sign, lam)
    res = evaluate(pts)
    print(f"sign={sign} lam={lam:.5f}: evaluate -> best={res[0]} euclid={res[1]} nblocks={res[2]} O={res[3]}")
    O = res[3]
    if O is None or O == 'inf':
        continue
    inv = invert(pts, O)
    P = np.array(inv)
    for tol in (1e-7, 1e-9, 1e-11):
        F, L, co = structure_of(P, tol=tol)
        D = sum(math.comb(len(B), 3) - 1 for B in F)
        print(f"   inverted about O, tol={tol}: rich blocks {sorted(len(B) for B in F)} lines {len(L)} (sizes {sorted(len(l) for l in L)}) -> circles = 120 - {D} - {len(L)} = {120 - D - len(L)}; min pairwise dist {min(np.hypot(*(P[i]-P[j])) for i,j in itertools.combinations(range(10),2)):.3e}")
