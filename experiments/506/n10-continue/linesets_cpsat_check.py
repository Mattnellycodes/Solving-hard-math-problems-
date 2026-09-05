"""Independent cross-check of the stage-2 line-set enumeration (enum4.line_sets) by CP-SAT enumeration."""
import sys, json, itertools, math
sys.path.insert(0, '/home/user/Solving-hard-math-problems-/experiments/506/n10-continue')
from enum4 import line_sets
from common import fs, deficit
from ortools.sat.python import cp_model
def cpsat_count(F, need):
    F = [fs(B) for B in F]
    covered = {fs(T) for B in F for T in itertools.combinations(sorted(B), 3)}
    cand = list(F) + [fs(T) for T in itertools.combinations(range(10), 3) if fs(T) not in covered]
    m = cp_model.CpModel(); y = {S: m.NewBoolVar('') for S in cand}
    for Q in itertools.combinations(range(10), 2):
        m.AddAtMostOne([y[S] for S in cand if set(Q) <= S])
    m.Add(sum(math.comb(len(S), 2) * y[S] for S in cand) <= 44)
    m.Add(sum(y.values()) >= need)
    n = [0]
    class CB(cp_model.CpSolverSolutionCallback):
        def on_solution_callback(self): n[0] += 1
    s = cp_model.CpSolver(); s.parameters.enumerate_all_solutions = True; s.parameters.num_workers = 1
    st = s.Solve(m, CB()); assert st in (cp_model.OPTIMAL, cp_model.INFEASIBLE), s.StatusName(st)
    return n[0]
for fn in sys.argv[1:]:
    d = json.load(open(fn))
    for k, rec in enumerate(d['F_classes']):
        need = rec['need_l']
        a = len(line_sets(rec['blocks'], need)); b = cpsat_count(rec['blocks'], need)
        flag = "OK" if a == b else "MISMATCH"
        if a or b or flag != "OK":
            print(f"{fn} F-class {k}: need l>={need}: DFS {a}, CP-SAT {b}: {flag}")
    print(fn, "checked", len(d['F_classes']), "F-classes")
