import itertools, sys
from ortools.sat.python import cp_model
S = frozenset(range(5)); R = frozenset(range(5, 10))
fours = [frozenset(c) for c in itertools.combinations(range(10), 4) if len(frozenset(c) & S) <= 2]
m = cp_model.CpModel()
x = {B: m.NewBoolVar('') for B in fours}
for T in itertools.combinations(range(10), 3):
    T = frozenset(T)
    m.AddAtMostOne([x[B] for B in fours if T <= B])
for p in range(10):
    bl = [B for B in fours if p in B]
    m.Add(3 * sum(x[B] for B in bl) + (6 if p in S else 0) <= 35)
    if p in R:
        m.Add(sum(x[B] for B in bl) <= 10)
m.Add(sum(x.values()) >= 22)
s = cp_model.CpSolver(); s.parameters.num_workers = 2
st = s.Solve(m)
F = sorted(sorted(B) for B in fours if s.Value(x[B]))
print(len(F), F)
from collections import Counter
print("d4:", [sum(1 for B in F if p in B) for p in range(10)])
print("types:", Counter(len(set(B) & S) for B in F))
