"""Joint subset + centre optimisation (one CP-SAT model per (universe, n)).

Variables: x_p (p in S), z_c (c is the inversion centre; c ranges over the universe points and the
extra sphere points), y_b (block b of the universe counted as a circle of the picture).
For every block b:  sum_{p in b} x_p  <=  2 + (|b| - 2) (y_b + sum_{c on b} z_c)
i.e. a block with >= 3 chosen points is counted unless the centre lies on it.  Exactly one centre,
which is not a chosen point.  Minimise sum y_b.
"""
import time
import numpy as np
from ortools.sat.python import cp_model
import sphere as SP


def optimise_joint(U, n, time_limit=60.0, workers=2, forbid=None, log=False, hint=None):
    m = cp_model.CpModel()
    x = {p: m.NewBoolVar(f"x{p}") for p in range(U.N)}
    zc = {('point', i): m.NewBoolVar(f"zp{i}") for i in range(U.N)}
    for j in range(len(U.extras)):
        zc[('extra', j)] = m.NewBoolVar(f"ze{j}")
    m.Add(sum(x.values()) == n)
    m.AddExactlyOne(zc.values())
    for i in range(U.N):
        m.AddImplication(zc[('point', i)], x[i].Not())
    # blocks through each centre
    through = {}
    for i in range(U.N):
        for t in U.blocks_of_point[i]:
            through.setdefault(t, []).append(zc[('point', i)])
    for j, bl in enumerate(U.extra_blocks):
        for t in bl:
            through.setdefault(t, []).append(zc[('extra', j)])
    y = {}
    for t, b in enumerate(U.blocks):
        pts = sorted(b)
        if len(pts) >= n:
            m.Add(sum(x[p] for p in pts) <= n - 1)          # non-degeneracy
        y[t] = m.NewBoolVar(f"y{t}")
        m.Add(sum(x[p] for p in pts) <= 2 + (len(pts) - 2) * (y[t] + sum(through.get(t, []))))
    if forbid:
        for S0 in forbid:
            m.Add(sum(x[p] for p in S0) <= len(S0) - 1)
    if hint:
        S0, c0 = hint
        for p in range(U.N):
            m.AddHint(x[p], 1 if p in S0 else 0)
        for c, v in zc.items():
            m.AddHint(v, 1 if c == c0 else 0)
    m.Minimize(sum(y.values()))
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_workers = workers
    solver.parameters.log_search_progress = log
    st = solver.Solve(m)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        S = sorted(p for p in x if solver.Value(x[p]))
        c = [k for k, v in zc.items() if solver.Value(v)][0]
        return ('OPTIMAL' if st == cp_model.OPTIMAL else 'FEASIBLE', int(round(solver.ObjectiveValue())), S, c,
                int(round(solver.BestObjectiveBound())))
    return ('INFEASIBLE' if st == cp_model.INFEASIBLE else 'UNKNOWN', None, None, None, None)


if __name__ == "__main__":
    import universes as UV
    from universe_opt import Universe
    nm, P, lab = UV.polyhedron_universe('cube+octahedron')
    U = Universe(nm, np.array(P, float), lab, None)
    for n in (9, 10):
        t0 = time.time()
        r = optimise_joint(U, n, time_limit=60, workers=2)
        print(n, r[:2], r[3], "bound", r[4], f"{time.time()-t0:.1f}s")
