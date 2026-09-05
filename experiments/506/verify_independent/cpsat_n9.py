"""n = 9 independent verification (2026-09-05).
Relaxation as in cpsat_relaxation.py (SG-only caps) PLUS the provably valid constraint that no point
lies in >= 8 four-blocks (8 derived 3-point lines on 8 points pairwise sharing <= 1 point form the
Möbius–Kantor (8_3) configuration, which has no real realisation; verified separately in mk_unrealisable.py).
Enumerate all labelled structures with D + l >= 60 (i.e. <= 24 circles) and apply an independent
hereditary Sylvester–Gallai check to the derived structure at every point:
  for every subset S of the 8 derived points, |S| >= 3, not contained in one derived line, some pair of S
  must lie on a derived line containing no third point of S (an 'ordinary line' of S).
Structures failing this at some point are unrealisable. Report survivors.
"""
import sys, itertools, math
from collections import Counter
from ortools.sat.python import cp_model
sys.path.insert(0, '.')
from cpsat_relaxation import build

def hereditary_sg_ok(n, F):
    """True if every derived structure passes hereditary SG."""
    for p in range(n):
        others = [q for q in range(n) if q != p]
        dlines = [frozenset(B - {p}) for B in F if p in B]   # derived rich lines (size >= 3)
        for r in range(3, len(others) + 1):
            for S in itertools.combinations(others, r):
                S = frozenset(S)
                if any(S <= L for L in dlines):
                    continue  # collinear subset: SG does not apply
                ordinary = False
                for a, b in itertools.combinations(S, 2):
                    L = next((L for L in dlines if a in L and b in L), None)
                    if L is None or len(L & S) == 2:
                        ordinary = True; break
                if not ordinary:
                    return False, p, sorted(S)
    return True, None, None

if __name__ == "__main__":
    n = 9; thr = 60; mode = sys.argv[1] if len(sys.argv) > 1 else 'sg'
    m, x, y, obj = build(n, mode)
    for p in range(n):
        m.Add(sum(x[B] for B in x if p in B and len(B) == 4) <= 7)
    m.Add(obj >= thr)
    solver = cp_model.CpSolver(); solver.parameters.enumerate_all_solutions = True
    solver.parameters.num_search_workers = 1; solver.parameters.max_time_in_seconds = 7200
    class Col(cp_model.CpSolverSolutionCallback):
        def __init__(s): super().__init__(); s.count = 0; s.inv = Counter(); s.ex = {}
        def on_solution_callback(s):
            s.count += 1
            F = [B for B, v in x.items() if s.Value(v)]; L = [S for S, v in y.items() if s.Value(v)]
            deg = Counter(q for B in F for q in B)
            inv = (tuple(sorted(len(B) for B in F)), tuple(sorted(deg[q] for q in range(n))), len(L))
            s.inv[inv] += 1; s.ex.setdefault(inv, (F, L))
    col = Col(); st = solver.Solve(m, col)
    print(f"n=9 mode={mode} D+l>={thr} with maxdeg4<=7: status={solver.StatusName(st)} labelled solutions={col.count}", flush=True)
    for inv, c in sorted(col.inv.items(), key=lambda t: -t[1]):
        F, L = col.ex[inv]
        cnt = math.comb(n, 3) - sum(math.comb(len(B), 3) - 1 for B in F) - len(L)
        ok, p, S = hereditary_sg_ok(n, F)
        print(f"  blocksizes={inv[0]} degrees={inv[1]} lines={inv[2]} copies={c} count={cnt}  hereditarySG={'PASS' if ok else 'FAIL at point '+str(p)+' subset '+str(S)}")
        print("     blocks:", [sorted(B) for B in F], "lines:", [sorted(S) for S in L], flush=True)
