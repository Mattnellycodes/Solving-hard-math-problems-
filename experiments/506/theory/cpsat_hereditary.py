"""CP-SAT model of the combinatorial relaxation PLUS hereditary non-realisability constraints:
  (H1) no derived Fano plane: for every point p, every 7-subset S of P\{p} and every Fano plane
       structure on S, at most 6 of its 7 lines l can satisfy 'l ∪ {p} is contained in a block';
  (H2) no derived (8_3): for every p, 8-subset S, and Mobius-Kantor labelling on S, at most 7 of the
       8 lines are derived lines (t_3(8) <= 7, proved in eight_three_light.py);
  (H1'),(H2') the same for the line set at infinity (y variables).
Usage: python3 cpsat_hereditary.py n [time_limit] [--nomk] [--enumerate target limit]
"""
import sys, itertools, time
from math import comb
from ortools.sat.python import cp_model
from cpsat_model import build, Collector

FANO = [{0,1,2},{0,3,4},{0,5,6},{1,3,5},{1,4,6},{2,3,6},{2,4,5}]
MK = [{i % 8, (i + 1) % 8, (i + 3) % 8} for i in range(8)]

def labellings(structure, k):
    """distinct images of the structure under all bijections {0..k-1} -> {0..k-1}."""
    out = set()
    for s in itertools.permutations(range(k)):
        out.add(frozenset(frozenset(s[i] for i in l) for l in structure))
    return [sorted(tuple(sorted(l)) for l in st) for st in out]

def add_hereditary(m, x, y, n, use_mk=True):
    blocks = list(x); lines = list(y)
    sup = {}   # frozenset -> list of block vars containing it
    def supersets_of(S):
        if S not in sup:
            sup[S] = [x[B] for B in blocks if S <= B]
        return sup[S]
    def linesup(S):
        return [y[L] for L in lines if S <= L]
    fano_l = labellings(FANO, 7); mk_l = labellings(MK, 8)
    cnt = 0
    for p in range(n):
        others = [q for q in range(n) if q != p]
        for S in itertools.combinations(others, 7):
            for lab in fano_l:
                terms = []
                for l in lab:
                    terms += supersets_of(frozenset(S[i] for i in l) | {p})
                m.Add(sum(terms) <= 6); cnt += 1
        if use_mk and len(others) >= 8:
            for S in itertools.combinations(others, 8):
                for lab in mk_l:
                    terms = []
                    for l in lab:
                        terms += supersets_of(frozenset(S[i] for i in l) | {p})
                    m.Add(sum(terms) <= 7); cnt += 1
    # lines at infinity
    for S in itertools.combinations(range(n), 7):
        for lab in fano_l:
            terms = []
            for l in lab:
                terms += linesup(frozenset(S[i] for i in l))
            m.Add(sum(terms) <= 6); cnt += 1
    if use_mk and n >= 8:
        for S in itertools.combinations(range(n), 8):
            for lab in mk_l:
                terms = []
                for l in lab:
                    terms += linesup(frozenset(S[i] for i in l))
                m.Add(sum(terms) <= 7); cnt += 1
    return cnt

def main():
    n = int(sys.argv[1]); tl = float(sys.argv[2]) if len(sys.argv) > 2 else 600.0
    use_mk = "--nomk" not in sys.argv
    t0 = time.time()
    m, x, y, obj = build(n)
    c = add_hereditary(m, x, y, n, use_mk)
    print(f"n={n}: model built with {c} hereditary constraints in {time.time()-t0:.1f}s", flush=True)
    N3 = comb(n, 3)
    if "--enumerate" in sys.argv:
        i = sys.argv.index("--enumerate"); target = int(sys.argv[i + 1]); limit = int(sys.argv[i + 2])
        m.Add(obj >= N3 - target)
        solver = cp_model.CpSolver(); solver.parameters.max_time_in_seconds = tl
        solver.parameters.enumerate_all_solutions = True; solver.parameters.num_workers = 1
        col = Collector(x, y, limit); st = solver.Solve(m, col)
        print(f"status {solver.StatusName(st)}: {len(col.sols)} labelled structures with circles <= {target}")
        for F, L in col.sols[:5]: print("  blocks:", F, "lines:", L)
        return
    m.Maximize(obj)
    solver = cp_model.CpSolver(); solver.parameters.max_time_in_seconds = tl; solver.parameters.num_workers = 2
    st = solver.Solve(m)
    print(f"n={n}: status={solver.StatusName(st)} best D+ell={solver.ObjectiveValue():.0f} bound={solver.BestObjectiveBound():.0f} "
          f"=> min circles = {N3 - solver.ObjectiveValue():.0f} (lower bound {N3 - solver.BestObjectiveBound():.0f}); "
          f"formula f(n)={comb(n-1,2)+1-(n-1)//2}; time {time.time()-t0:.1f}s", flush=True)
    F = sorted(sorted(B) for B, v in x.items() if solver.Value(v)); L = sorted(sorted(S) for S, v in y.items() if solver.Value(v))
    print("  optimal blocks:", F); print("  lines:", L)

if __name__ == "__main__":
    main()
