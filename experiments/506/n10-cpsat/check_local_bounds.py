"""Brute-force check of local packing bounds at a point p (derived structure on the 9 other points):
k four-lines and t three-lines on 9 points, all pairwise sharing <= 1 point.  Report max t for each k.
This is a pure packing computation (no geometry)."""
import itertools
from ortools.sat.python import cp_model
pts = range(9)
fours = [frozenset(c) for c in itertools.combinations(pts, 4)]
threes = [frozenset(c) for c in itertools.combinations(pts, 3)]
for k in range(0, 5):
    m = cp_model.CpModel()
    x = {B: m.NewBoolVar('') for B in fours}; y = {T: m.NewBoolVar('') for T in threes}
    for Q in itertools.combinations(pts, 2):
        Q = frozenset(Q)
        m.AddAtMostOne([x[B] for B in fours if Q <= B] + [y[T] for T in threes if Q <= T])
    m.Add(sum(x.values()) == k)
    # symmetry breaking: none needed, just maximise
    m.Maximize(sum(y.values()))
    s = cp_model.CpSolver(); s.parameters.num_workers = 2; s.parameters.max_time_in_seconds = 120
    st = s.Solve(m)
    print(f"k={k} four-lines: status={s.StatusName(st)} max three-lines={s.ObjectiveValue():.0f} (pure packing, no MK/SG)")
