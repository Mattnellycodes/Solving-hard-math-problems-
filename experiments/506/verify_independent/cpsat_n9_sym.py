"""n = 9, complete enumeration with symmetry breaking (2026-09-05).
Same model as cpsat_n9.py (SG-only or table caps, "<= 7 four-blocks through any point" by the (8_3) lemma,
D + l >= 60, i.e. <= 24 circles).
Symmetry breaking: with <= 7 four-blocks per point, b4 <= floor(9*7/4) = 15, so D <= 45 from 4-blocks, and
l <= floor((C(9,2) - 1)/3) = 11; D + l >= 60 therefore forces a block of size >= 5. Relabel so that it
contains {0,1,2,3,4}: add  OR_{B ⊇ {0,1,2,3,4}} x[B].
Enumerate ALL solutions (no time limit) and apply the hereditary Sylvester–Gallai check to each.
"""
import sys, itertools, math, time
from collections import Counter
from ortools.sat.python import cp_model
sys.path.insert(0, '.')
from cpsat_relaxation import build
from cpsat_n9 import hereditary_sg_ok
n = 9; thr = 60; mode = sys.argv[1] if len(sys.argv) > 1 else 'sg'
m, x, y, obj = build(n, mode)
for p in range(n):
    m.Add(sum(x[B] for B in x if p in B and len(B) == 4) <= 7)
m.Add(obj >= thr)
core = frozenset(range(5))
m.AddBoolOr([x[B] for B in x if core <= B])
solver = cp_model.CpSolver(); solver.parameters.enumerate_all_solutions = True
solver.parameters.num_search_workers = 1
class Col(cp_model.CpSolverSolutionCallback):
    def __init__(s): super().__init__(); s.count = 0; s.inv = Counter(); s.ex = {}; s.fails = 0; s.passes = []
    def on_solution_callback(s):
        s.count += 1
        F = [B for B, v in x.items() if s.Value(v)]; L = [S for S, v in y.items() if s.Value(v)]
        deg = Counter(q for B in F for q in B)
        inv = (tuple(sorted(len(B) for B in F)), tuple(sorted(deg[q] for q in range(n))), len(L))
        s.inv[inv] += 1
        if inv not in s.ex: s.ex[inv] = (F, L)
        ok, p, S = hereditary_sg_ok(n, F)
        if ok: s.passes.append((F, L))
        else: s.fails += 1
t0 = time.time(); col = Col(); st = solver.Solve(m, col)
print(f"n=9 mode={mode} D+l>={thr}, maxdeg4<=7, symmetry-broken: status={solver.StatusName(st)} labelled solutions={col.count} time={time.time()-t0:.0f}s", flush=True)
print(f"hereditary-SG: {col.fails} solutions FAIL, {len(col.passes)} PASS")
for inv, c in sorted(col.inv.items(), key=lambda t: -t[1]):
    F, L = col.ex[inv]
    cnt = math.comb(n, 3) - sum(math.comb(len(B), 3) - 1 for B in F) - len(L)
    print(f"  blocksizes={inv[0]} degrees={inv[1]} lines={inv[2]} copies={c} count={cnt}")
    print("     blocks:", [sorted(B) for B in F], "lines:", [sorted(S) for S in L])
for F, L in col.passes[:5]:
    print("PASSING structure (needs geometric analysis):", [sorted(B) for B in F], [sorted(S) for S in L])
