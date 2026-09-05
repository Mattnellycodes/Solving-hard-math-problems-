"""lemmas-audit: demonstrate that the theory agent's cpsat_hereditary.py encodes the hereditary
Sylvester-Gallai constraint WITHOUT the hypothesis 'S not contained in one block through p', and is
therefore unsound: it declares the (realisable!) antipodal 9-point structure (one 8-block + centre,
4 lines, 25 circles) infeasible.  (This tool was NOT used for the n <= 9 results in REPORT.md, but is
proposed there as the tool for n = 10.)
We import their model builder only to exhibit the bug; nothing here is used as evidence for the
lemmas.
"""
import sys, itertools, time
sys.path.insert(0, "/home/user/Solving-hard-math-problems-/experiments/506/theory")
from ortools.sat.python import cp_model
from cpsat_model import build
from cpsat_hereditary import add_hereditary

n = 9
m, x, y, obj = build(n)
t0 = time.time()
c = add_hereditary(m, x, y, n, use_mk=False)     # Fano constraints only (H1),(H1')
print(f"model with {c} Fano-type hereditary constraints built in {time.time()-t0:.1f}s")
# fix the antipodal structure: 8-block {1..8}, lines {0,1,2},{0,3,4},{0,5,6},{0,7,8}
big = frozenset(range(1, 9))
for B, v in x.items():
    m.Add(v == (1 if B == big else 0))
lines = [frozenset({0, 1, 2}), frozenset({0, 3, 4}), frozenset({0, 5, 6}), frozenset({0, 7, 8})]
for S, v in y.items():
    m.Add(v == (1 if S in lines else 0))
solver = cp_model.CpSolver(); solver.parameters.max_time_in_seconds = 300; solver.parameters.num_workers = 1
st = solver.Solve(m)
print("status for the antipodal structure under the theory's hereditary model:", solver.StatusName(st))
print("(INFEASIBLE means the model wrongly excludes a realisable configuration: 8 concyclic points + centre)")
