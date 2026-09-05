"""From-scratch CP-SAT enumeration of abstract Moebius structures on n=8 points using ONLY
the trivially necessary conditions:
  (I1) blocks (subsets of size 4..7) pairwise share <= 2 points;
  (I2) lines (subsets of size 3..7) pairwise share <= 1 point; a line of size >=4 must be a block;
       a line of size 3 must not be contained in a block (its triple would then lie in that block).
Objective D + ell where D = sum_B (C(|B|,3)-1), ell = #lines. circles = 56 - D - ell.
No Sylvester-Gallai / orchard caps at all.  Enumerate ALL labelled solutions with D+ell >= T,
then collapse to (F, best ell) and to isomorphism classes (canonical form under S_8, brute force).
Usage: python3 enum8_cpsat.py T
"""
import sys, itertools
from math import comb
from ortools.sat.python import cp_model

n = 8
T = int(sys.argv[1])
subsets = {k: [frozenset(c) for c in itertools.combinations(range(n), k)] for k in range(3, n)}
m = cp_model.CpModel()
x = {B: m.NewBoolVar("x%s" % sorted(B)) for k in range(4, n) for B in subsets[k]}
y = {S: m.NewBoolVar("y%s" % sorted(S)) for k in range(3, n) for S in subsets[k]}
Bs = list(x); Ls = list(y)
for A, B in itertools.combinations(Bs, 2):
    if len(A & B) >= 3:
        m.AddBoolOr([x[A].Not(), x[B].Not()])
for A, B in itertools.combinations(Ls, 2):
    if len(A & B) >= 2:
        m.AddBoolOr([y[A].Not(), y[B].Not()])
for S in Ls:
    if len(S) >= 4:
        m.AddImplication(y[S], x[S])
    else:
        for B in Bs:
            if S <= B:
                m.AddBoolOr([y[S].Not(), x[B].Not()])
obj = sum((comb(len(B), 3) - 1) * x[B] for B in Bs) + sum(y[S] for S in Ls)
m.Add(obj >= T)

class Col(cp_model.CpSolverSolutionCallback):
    def __init__(s):
        super().__init__(); s.sols = {}
    def on_solution_callback(s):
        F = tuple(sorted(sum(1 << i for i in B) for B in Bs if s.Value(x[B])))
        L = tuple(sorted(sum(1 << i for i in S) for S in Ls if s.Value(y[S])))
        D = sum(comb(bin(b).count("1"), 3) - 1 for b in F)
        key = F
        if key not in s.sols or s.sols[key][0] < D + len(L):
            s.sols[key] = (D + len(L), L)

solver = cp_model.CpSolver()
solver.parameters.enumerate_all_solutions = True
solver.parameters.num_workers = 1
col = Col()
st = solver.Solve(m, col)
print("status", solver.StatusName(st), "labelled block-families with D+ell >=", T, ":", len(col.sols))

perms = list(itertools.permutations(range(n)))
def img(c, s): return sum(1 << s[i] for i in range(n) if c >> i & 1)
def canon(F): return min(tuple(sorted(img(c, s) for c in F)) for s in perms)
classes = {}
for F, (val, L) in col.sols.items():
    cf = canon(F)
    if cf not in classes or classes[cf][0] < val:
        classes[cf] = (val, F, L)
print("isomorphism classes:", len(classes))
def bits(c): return [i for i in range(n) if c >> i & 1]
for cf, (val, F, L) in sorted(classes.items(), key=lambda kv: -kv[1][0]):
    sizes = sorted(bin(b).count("1") for b in F)
    deg = [sum(1 for b in F if b >> p & 1) for p in range(n)]
    print(f"D+ell={val} circles={56-val} sizes={sizes} degrees={deg} ell={len(L)}")
    print("   blocks:", [bits(b) for b in F])
    print("   lines :", [bits(b) for b in L])
