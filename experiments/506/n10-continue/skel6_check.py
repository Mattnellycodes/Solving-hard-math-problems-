"""Computational cross-check of the single-5-block case (skeleton 6): with B5 = {0..4} the only block of size
>= 5, no family of 4-blocks pairwise sharing <= 2 points (and <= 2 with B5) has b4 >= 22 when every point off B5
lies in <= 10 four-blocks (t3(9) <= 10, from the (8_3) lemma).  Since D + l >= 88 with D = 9 + 3 b4 and l <= 14
(SG) needs b4 >= 22, this settles the skeleton.  Two CP-SAT models: (a) explicit cap d4(r) <= 10; (b) instead the
(8_3) constraint 'at most 7 derived 3-lines inside any 8-subset of the derived 9-set' (t3(8) <= 7).
"""
import itertools, math, sys, time
from ortools.sat.python import cp_model
S = frozenset(range(5)); R = frozenset(range(5, 10))
fours = [frozenset(c) for c in itertools.combinations(range(10), 4) if len(frozenset(c) & S) <= 2]
print("candidate 4-blocks:", len(fours))
def build(mode):
    m = cp_model.CpModel()
    x = {B: m.NewBoolVar('') for B in fours}
    for T in itertools.combinations(range(10), 3):
        T = frozenset(T)
        m.AddAtMostOne([x[B] for B in fours if T <= B])          # pairwise <= 2 (triples of S are in B5: no 4-block has 3 pts of S anyway)
    for p in range(10):
        bl = [B for B in fours if p in B]
        # SG cap on the derived 9-point set: 3 d4 + 6 [p in S] <= 35
        m.Add(3 * sum(x[B] for B in bl) + (6 if p in S else 0) <= 35)
        if mode == 'a' and p in R:
            m.Add(sum(x[B] for B in bl) <= 10)
        if mode == 'b':
            others = [q for q in range(10) if q != p]
            for E in itertools.combinations(others, 8):
                E = frozenset(E)
                terms = [x[B] for B in bl if (B - {p}) <= E]
                if len(terms) > 7:
                    m.Add(sum(terms) <= 7)
    m.Add(sum(x.values()) >= 22)
    return m
for mode in 'ab':
    t0 = time.time()
    m = build(mode)
    s = cp_model.CpSolver(); s.parameters.num_workers = 2; s.parameters.max_time_in_seconds = 1800
    st = s.Solve(m)
    print(f"mode {mode}: status {s.StatusName(st)}  [{time.time()-t0:.0f}s]", flush=True)
