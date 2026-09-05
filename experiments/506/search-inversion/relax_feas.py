"""Feasibility of the combinatorial relaxation for beating f(n): is there an abstract Moebius
block structure (rich blocks of size 4..mmax pairwise sharing <= 2 points, lines pairwise sharing
<= 1, derived Sylvester-Gallai caps at every point, line cap) with D + l >= C(n-1,3) + floor((n-1)/2)?

Reuses ../verify_independent/cpsat_relaxation.build (independent model) but restricts block sizes to
<= mmax (largest-block Lemmas A, B, C of ../theory/REPORT.md give mmax = 6 for n = 10, 11, 12) and
solves the pure feasibility problem.  A proof of infeasibility shows that no planar n-point set beats
f(n), modulo the ordinary-line values used (mode 'sg': o = 1 only, i.e. Sylvester-Gallai alone).

Usage: python3 relax_feas.py n mmax mode[sg|table] seconds [workers]
"""
import sys, os, math, itertools, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'verify_independent'))
from ortools.sat.python import cp_model
from cpsat_relaxation import O_TABLE


def build(n, mmax, mode):
    o = (lambda m: 1) if mode == 'sg' else (lambda m: O_TABLE[m])
    pts = range(n)
    rich = [frozenset(c) for k in range(4, mmax + 1) for c in itertools.combinations(pts, k)]
    lines = [frozenset(c) for k in range(3, mmax + 1) for c in itertools.combinations(pts, k)]
    m = cp_model.CpModel()
    x = {B: m.NewBoolVar(f"x{sorted(B)}") for B in rich}
    y = {S: m.NewBoolVar(f"y{sorted(S)}") for S in lines}
    # two rich blocks share <= 2 points: for every triple T, at most one rich block contains T
    by_triple = {}
    for B in rich:
        for T in itertools.combinations(sorted(B), 3):
            by_triple.setdefault(T, []).append(x[B])
    for T, vs in by_triple.items():
        if len(vs) > 1:
            m.AddAtMostOne(vs)
    # lines: size>=4 line must be a rich block; size-3 line must be an uncovered triple
    for S in lines:
        if len(S) >= 4:
            m.AddImplication(y[S], x[S])
        else:
            T = tuple(sorted(S))
            for v in by_triple.get(T, []):
                m.AddBoolOr([y[S].Not(), v.Not()])
    # two lines share <= 1 point: for every pair, at most one line contains it
    by_pair = {}
    for S in lines:
        for pr in itertools.combinations(sorted(S), 2):
            by_pair.setdefault(pr, []).append(y[S])
    for pr, vs in by_pair.items():
        if len(vs) > 1:
            m.AddAtMostOne(vs)
    # derived SG caps and line cap
    for p in pts:
        m.Add(sum(math.comb(len(B) - 1, 2) * x[B] for B in rich if p in B) <= math.comb(n - 1, 2) - o(n - 1))
    m.Add(sum(math.comb(len(S), 2) * y[S] for S in lines) <= math.comb(n, 2) - o(n))
    D = sum((math.comb(len(B), 3) - 1) * x[B] for B in rich)
    L = sum(y.values())
    target = math.comb(n - 1, 3) + (n - 1) // 2
    m.Add(D + L >= target)
    # symmetry breaking: point 0 has the maximum rich degree (weighted coverage)  -- valid up to relabelling
    cov = [sum(math.comb(len(B) - 1, 2) * x[B] for B in rich if p in B) for p in pts]
    for p in range(1, n):
        m.Add(cov[0] >= cov[p])
    return m, x, y, target


if __name__ == "__main__":
    n = int(sys.argv[1]); mmax = int(sys.argv[2]); mode = sys.argv[3]
    secs = float(sys.argv[4]); workers = int(sys.argv[5]) if len(sys.argv) > 5 else 2
    t0 = time.time()
    m, x, y, target = build(n, mmax, mode)
    print(f"n={n} mmax={mmax} mode={mode}: {len(x)} block vars, {len(y)} line vars, target D+l>={target}, built in {time.time()-t0:.1f}s", flush=True)
    solver = cp_model.CpSolver()
    solver.parameters.num_workers = workers
    solver.parameters.max_time_in_seconds = secs
    solver.parameters.log_search_progress = True
    st = solver.Solve(m)
    print(f"STATUS {solver.StatusName(st)} after {time.time()-t0:.1f}s", flush=True)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        F = [sorted(B) for B, v in x.items() if solver.Value(v)]
        L = [sorted(S) for S, v in y.items() if solver.Value(v)]
        D = sum(math.comb(len(B), 3) - 1 for B in F)
        print("FEASIBLE structure: D=", D, "l=", len(L), "circles=", math.comb(n, 3) - D - len(L))
        print("blocks", F); print("lines", L)
