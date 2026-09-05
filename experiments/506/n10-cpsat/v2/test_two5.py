"""Test of analyse_two5.py on a known two-disjoint-5-block structure (the 'C-2' block family reported by the
n10-enum agent; only the block list is taken from their log as TEST INPUT, the analysis is ours).  A line set of
size >= 10 is generated here by a tiny CP-SAT packing model (lines = uncovered triples or blocks, pairwise <= 1)."""
import itertools, json, sys
from ortools.sat.python import cp_model
F = [[0,1,2,3,4],[5,6,7,8,9],[0,1,5,6],[0,1,7,8],[0,2,5,7],[0,2,6,9],[0,3,5,8],[0,3,7,9],[0,4,5,9],[0,4,6,8],[1,2,6,7],[1,2,8,9],[1,3,5,9],[1,3,6,8],[1,4,5,7],[1,4,6,9],[2,3,5,6],[2,3,7,8],[2,4,5,8],[2,4,7,9],[3,4,6,7],[3,4,8,9]]
Fs = [frozenset(B) for B in F]
covered = {frozenset(T) for B in Fs for T in itertools.combinations(sorted(B), 3)}
cand = Fs + [frozenset(T) for T in itertools.combinations(range(10), 3) if frozenset(T) not in covered]
m = cp_model.CpModel(); y = [m.NewBoolVar('') for _ in cand]
for i, j in itertools.combinations(range(len(cand)), 2):
    if len(cand[i] & cand[j]) >= 2: m.AddBoolOr([y[i].Not(), y[j].Not()])
m.Add(sum(y) >= 10); m.Add(sum(len(c) * (len(c) - 1) // 2 * y[i] for i, c in enumerate(cand)) <= 44)
s = cp_model.CpSolver(); s.parameters.num_workers = 1; st = s.Solve(m)
L = [sorted(cand[i]) for i in range(len(cand)) if s.Value(y[i])]
print("test line set:", L)
json.dump([{'index': 0, 'blocks': F, 'lines': L, 'count': 120 - (18 + 60) - len(L), 'report': {'kills': []}}], open('test_two5_input.json', 'w'))
