"""From-scratch CP-SAT model, n=8.  Mode A: intersection conditions only, maximise D+ell.
Mode B: + 'no derived Fano' (every point in <= 6 four-blocks; the only consequence of
Sylvester-Gallai that can bite on 8 points), prove D+ell <= 39 (circles >= 17) and enumerate
all solutions with D+ell >= thr (collapsed to block families / iso classes).
Usage: python3 opt8_cpsat.py A|B [threshold]
"""
import sys, itertools, time
from math import comb
from ortools.sat.python import cp_model
n = 8
mode = sys.argv[1]; thr = int(sys.argv[2]) if len(sys.argv) > 2 else 39
def build(no_fano):
    subsets = {k: [frozenset(c) for c in itertools.combinations(range(n), k)] for k in range(3, n)}
    m = cp_model.CpModel()
    x = {B: m.NewBoolVar("x%s" % sorted(B)) for k in range(4, n) for B in subsets[k]}
    y = {S: m.NewBoolVar("y%s" % sorted(S)) for k in range(3, n) for S in subsets[k]}
    Bs = list(x); Ls = list(y)
    for A, B in itertools.combinations(Bs, 2):
        if len(A & B) >= 3: m.AddBoolOr([x[A].Not(), x[B].Not()])
    for A, B in itertools.combinations(Ls, 2):
        if len(A & B) >= 2: m.AddBoolOr([y[A].Not(), y[B].Not()])
    for S in Ls:
        if len(S) >= 4: m.AddImplication(y[S], x[S])
        else:
            for B in Bs:
                if S <= B: m.AddBoolOr([y[S].Not(), x[B].Not()])
    if no_fano:
        for p in range(n):
            m.Add(sum(x[B] for B in Bs if p in B and len(B) == 4) <= 6)
    obj = sum((comb(len(B), 3) - 1) * x[B] for B in Bs) + sum(y[S] for S in Ls)
    return m, x, y, Bs, Ls, obj
def show(F, L):
    def bits(c): return [i for i in range(n) if c >> i & 1]
    return [bits(b) for b in F], [bits(b) for b in L]
if mode == "A":
    m, x, y, Bs, Ls, obj = build(False)
    m.Maximize(obj)
    s = cp_model.CpSolver(); s.parameters.num_workers = 2; s.parameters.max_time_in_seconds = 600
    t = time.time(); st = s.Solve(m)
    print("A (intersection only): status", s.StatusName(st), "max D+ell =", s.ObjectiveValue(), "bound", s.BestObjectiveBound(), "=> min circles", 56 - s.ObjectiveValue(), f"({time.time()-t:.1f}s)")
    F = sorted(sum(1 << i for i in B) for B in Bs if s.Value(x[B])); L = sorted(sum(1 << i for i in S) for S in Ls if s.Value(y[S]))
    print("  blocks, lines:", show(F, L))
    print("  four-block degrees:", [sum(1 for B in Bs if s.Value(x[B]) and p in B and len(B) == 4) for p in range(n)])
else:
    m, x, y, Bs, Ls, obj = build(True)
    m.Maximize(obj)
    s = cp_model.CpSolver(); s.parameters.num_workers = 2; s.parameters.max_time_in_seconds = 900
    t = time.time(); st = s.Solve(m)
    print("B (+ no derived Fano): status", s.StatusName(st), "max D+ell =", s.ObjectiveValue(), "bound", s.BestObjectiveBound(), "=> min circles", 56 - s.ObjectiveValue(), f"({time.time()-t:.1f}s)", flush=True)
    m2, x2, y2, Bs2, Ls2, obj2 = build(True)
    m2.Add(obj2 >= thr)
    class Col(cp_model.CpSolverSolutionCallback):
        def __init__(s_): super().__init__(); s_.sols = {}
        def on_solution_callback(s_):
            F = tuple(sorted(sum(1 << i for i in B) for B in Bs2 if s_.Value(x2[B])))
            L = tuple(sorted(sum(1 << i for i in S) for S in Ls2 if s_.Value(y2[S])))
            D = sum(comb(bin(b).count("1"), 3) - 1 for b in F)
            if F not in s_.sols or s_.sols[F][0] < D + len(L): s_.sols[F] = (D + len(L), L)
    s2 = cp_model.CpSolver(); s2.parameters.enumerate_all_solutions = True; s2.parameters.num_workers = 1
    s2.parameters.max_time_in_seconds = 900
    col = Col(); t = time.time(); st = s2.Solve(m2, col)
    print("enumeration of D+ell >=", thr, ": status", s2.StatusName(st), "labelled block families:", len(col.sols), f"({time.time()-t:.1f}s)")
    perms = list(itertools.permutations(range(n)))
    def img(c, sg): return sum(1 << sg[i] for i in range(n) if c >> i & 1)
    def canon(F): return min(tuple(sorted(img(c, sg) for c in F)) for sg in perms)
    classes = {}
    for F, (val, L) in col.sols.items():
        cf = canon(F)
        if cf not in classes or classes[cf][0] < val: classes[cf] = (val, F, L)
    print("isomorphism classes:", len(classes))
    for cf, (val, F, L) in classes.items():
        print(f"  D+ell={val} circles={56-val} sizes={sorted(bin(b).count('1') for b in F)} ell={len(L)}")
        print("   ", show(F, L))
