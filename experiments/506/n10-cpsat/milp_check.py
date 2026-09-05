"""Independent cross-check of the CP-SAT infeasibility results with a different solver (HiGHS MILP via
scipy.optimize.milp).  Same combinatorial relaxation, case B with a fixed skeleton of 5-blocks, written
as a pure 0/1 integer program:  maximise D + l  subject to the linear constraints (T),(L),(SG),(SGL),
(PD),(F7),(M8),(O9),(O10 optional),(DD).  A proven optimum <= 87 confirms the CP-SAT INFEASIBLE result.
Usage: python3 milp_check.py omode skeleton_index [--orchard10] [--time T]
"""
import sys, json, math, itertools, time, argparse
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix, csr_matrix

ap = argparse.ArgumentParser()
ap.add_argument('omode'); ap.add_argument('idx', type=int); ap.add_argument('--orchard10', action='store_true')
ap.add_argument('--time', type=float, default=1200)
args = ap.parse_args()
O = {'table': {9: 6, 10: 5}, 'cs': {9: 5, 10: 5}, 'sg': {9: 1, 10: 1}}[args.omode]
D = '/home/user/Solving-hard-math-problems-/experiments/506/n10-cpsat'
sk = json.load(open(D + '/skeletons.json'))[args.idx]
fam = [frozenset(B) for B in sk['blocks']]
P = range(10)
fours = [frozenset(c) for c in itertools.combinations(P, 4)]
threes = [frozenset(c) for c in itertools.combinations(P, 3)]
# variables: x4[B] (210), y3[T] (120), y4[B] (210), y5[B] for B in fam
var = {}
for B in fours: var[('x', B)] = len(var)
for T in threes: var[('y', T)] = len(var)
for B in fours: var[('y', B)] = len(var)
for B in fam: var[('y', B)] = len(var)
nv = len(var)
rows = []; lo = []; hi = []
def add(coefs, l, h):
    rows.append(coefs); lo.append(l); hi.append(h)
blocks_all = fours + fam           # rich blocks (5-blocks fixed = fam)
def xcoef(B):                      # coefficient dict for indicator of rich block B (fixed 1 for fam)
    return ({var[('x', B)]: 1}, 0) if len(B) == 4 else ({}, 1)
lines_all = threes + fours + fam
# (T): every triple in <= 1 rich block, and y3[T] excluded if covered
for T in threes:
    c = {var[('y', T)]: 1}; const = 0
    for B in blocks_all:
        if T <= B:
            d, k = xcoef(B); c.update(d); const += k
    add(c, -np.inf, 1 - const)
# (L): y4[B] <= x4[B]; two lines share <= 1 point
for B in fours:
    add({var[('y', B)]: 1, var[('x', B)]: -1}, -np.inf, 0)
for Q in itertools.combinations(P, 2):
    Q = frozenset(Q)
    add({var[('y', S)]: 1 for S in lines_all if Q <= S}, -np.inf, 1)
# (SG) derived caps and (DD) local bounds
for p in P:
    c = {}; const = 0
    for B in blocks_all:
        if p in B:
            d, k = xcoef(B); w = math.comb(len(B) - 1, 2)
            for v in d: c[v] = c.get(v, 0) + w
            const += w * k
    add(c, -np.inf, 36 - O[9] - const)
    d5 = sum(1 for B in fam if p in B)
    d4 = {var[('x', B)]: 1 for B in fours if p in B}
    add(d4, -np.inf, 13 - 3 * d5)
    if d5 >= 1: add(d4, -np.inf, 8)
    add(d4, -np.inf, 10)                                   # (O9)
    # (PD) lines through p
    add({var[('y', S)]: (len(S) - 1) for S in lines_all if p in S}, -np.inf, 9)
# (SGL)
add({var[('y', S)]: math.comb(len(S), 2) for S in lines_all}, -np.inf, 45 - O[10])
if args.orchard10:
    add({var[('y', T)]: 1 for T in threes}, -np.inf, 12)
# (PD) pair degree
for Q in itertools.combinations(P, 2):
    Q = frozenset(Q); c = {}; const = 0
    for B in blocks_all:
        if Q <= B:
            d, k = xcoef(B)
            for v in d: c[v] = c.get(v, 0) + (len(B) - 2)
            const += (len(B) - 2) * k
    add(c, -np.inf, 8 - const)
# (F7),(M8)
for r, cap in ((7, 6), (8, 7)):
    for p in P:
        others = [q for q in P if q != p]
        for S in itertools.combinations(others, r):
            S = frozenset(S); c = {}; const = 0
            for B in blocks_all:
                if p in B and len(B & S) == 3:
                    d, k = xcoef(B); c.update(d); const += k
            if len(c) + const > cap: add(c, -np.inf, cap - const)
    for S in itertools.combinations(P, r):
        S = frozenset(S)
        c = {var[('y', L)]: 1 for L in lines_all if len(L & S) == 3}
        if len(c) > cap: add(c, -np.inf, cap)
A = lil_matrix((len(rows), nv))
for i, c in enumerate(rows):
    for v, w in c.items(): A[i, v] = w
A = csr_matrix(A)
obj = np.zeros(nv)
for B in fours: obj[var[('x', B)]] = 3
for S in lines_all: obj[var[('y', S)]] = 1
const_obj = 9 * len(fam)
t0 = time.time()
res = milp(-obj, constraints=LinearConstraint(A, lo, hi), integrality=np.ones(nv), bounds=Bounds(0, 1),
           options={'time_limit': args.time, 'disp': False})
val = -res.fun + const_obj if res.x is not None else None
bound = None
try:
    bound = -res.mip_dual_bound + const_obj
except Exception:
    pass
print(f"skeleton {args.idx} omode {args.omode} orchard10={args.orchard10}: status={res.status} ({res.message.strip()}) "
      f"opt D+l={val} dual bound={bound} gap={getattr(res,'mip_gap',None)} time={time.time()-t0:.0f}s vars={nv} rows={len(rows)}")
