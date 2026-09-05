import sys, time, numpy as np
sys.path.insert(0, '/home/user/Solving-hard-math-problems-/experiments/506/search-local')
from engine import *
import universes as UV

# 1. the 8-point record (17) as a subset of an orthocentric closure universe
U = UV.ortho_closure([(0, 3), (1, 0), (3, 0)], levels=1)
print(U.summary())
# the record set scaled: (0,0),(0,5),(0,10),(0,15),(3,6),(3,9),(5,5),(5,10) = inversion of {A,B,C,H,D,E,F,inf} about O=(0,1)
# In the universe: A=(0,3),B=(1,0),C=(3,0),H=(0,-1),D=(6/5,-3/5),E=(2,1),F=(0,0), inf
from fractions import Fraction as Fr
want = [(Fr(0),Fr(3)),(Fr(1),Fr(0)),(Fr(3),Fr(0)),(Fr(0),Fr(-1)),(Fr(6,5),Fr(-3,5)),(Fr(2),Fr(1)),(Fr(0),Fr(0))]
idx = [U.exact.index(w) for w in want] + [U.N]
S = np.array(idx, dtype=np.int64)
v, o, nb = k_eval(S, len(S), U.T, U.bptr, U.bmem, U.pptr, U.pblk)
print("proxy value", v, "centre", U.label(o), "nblocks", nb, "true:", true_count(U, S)[:4], describe(U, S))
# 2. grid 4x4 exhaustive n=6..8 (compare with subset-search agent: 8, 11, 18)
G = UV.grid(4)
print(G.summary())
for n in (6, 7, 8):
    t = time.time()
    total, best, hist, hist_e, good = k_exhaustive(n, G.M, G.T, G.bptr, G.bmem, G.pptr, G.pblk, best_thresh := 0)
    print(f"grid4 n={n}: subsets={total} best={best} time={time.time()-t:.1f}s hist={[(i,int(hist[i])) for i in range(400) if hist[i]][:6]}")
