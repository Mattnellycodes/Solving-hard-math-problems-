"""Independent cross-check of skeleton infeasibility with a different solver (z3, pseudo-Boolean encoding of the
same constraint set as model.py, written separately).  Usage: python3 z3_check.py omode i [--orchard10] [--time T]"""
import sys, json, time, argparse, itertools, math
import z3
sys.path.insert(0, '/home/user/Solving-hard-math-problems-/experiments/506/n10-cpsat/v2')
from model import O_MODES, T3_CITED

ap = argparse.ArgumentParser(); ap.add_argument('omode'); ap.add_argument('i', type=int)
ap.add_argument('--orchard10', action='store_true'); ap.add_argument('--time', type=float, default=1800); ap.add_argument('--thr', type=int, default=88)
args = ap.parse_args()
sk = json.load(open('/home/user/Solving-hard-math-problems-/experiments/506/n10-cpsat/v2/skeletons.json'))[args.i]
if sk['case'] == 'A': sizes, forbid = (4, 5, 6), (5, 6)
elif sk['case'] == 'B': sizes, forbid = (4, 5), (5,)
elif sk['case'] == 'C': sizes, forbid = (4,), ()
else:
    mm = max(len(B) for B in sk['blocks']); sizes, forbid = tuple(range(4, mm + 1)), (mm,)
o = O_MODES[args.omode]; n = 10; PTS = range(n)
skel = [frozenset(B) for B in sk['blocks']]
blocks = [frozenset(c) for k in sizes for c in itertools.combinations(PTS, k)]
lines = [frozenset(c) for k in (3,) + sizes for c in itertools.combinations(PTS, k)]
x = {B: z3.Bool('x' + ''.join(map(str, sorted(B)))) for B in blocks}
y = {S: z3.Bool('y' + ''.join(map(str, sorted(S)))) for S in lines}
s = z3.Solver(); s.set('timeout', int(args.time * 1000))
def le(terms, bound):   # sum c_i * lit_i <= bound
    s.add(z3.PbLe([(v, c) for c, v in terms], bound))
def ge(terms, bound):
    s.add(z3.PbGe([(v, c) for c, v in terms], bound))
for T in itertools.combinations(PTS, 3):
    T = frozenset(T); le([(1, y[T])] + [(1, x[B]) for B in blocks if T <= B], 1)
for S in lines:
    if len(S) >= 4: s.add(z3.Implies(y[S], x[S]))
for Q in itertools.combinations(PTS, 2):
    Q = frozenset(Q); le([(1, y[S]) for S in lines if Q <= S], 1)
    le([(len(B) - 2, x[B]) for B in blocks if Q <= B], n - 2)
for p in PTS:
    le([(len(S) - 1, y[S]) for S in lines if p in S], n - 1)
    le([(math.comb(len(B) - 1, 2), x[B]) for B in blocks if p in B], math.comb(n - 1, 2) - o(n - 1))
    d4 = [(1, x[B]) for B in blocks if p in B and len(B) == 4]; d5 = [(1, x[B]) for B in blocks if p in B and len(B) == 5]
    d6 = [(1, x[B]) for B in blocks if p in B and len(B) == 6]
    le([(2, v) for _, v in d4] + [(3, v) for _, v in d5] + [(8, v) for _, v in d6], 20)
    le([(1, v) for _, v in d5] + [(2, v) for _, v in d6], 3)
le([(math.comb(len(S), 2), y[S]) for S in lines], math.comb(n, 2) - o(n))
def subset_caps(point_set, linevars):
    for r in (7, 8, 9):
        for Ssub in itertools.combinations(sorted(point_set), r):
            Ssub = frozenset(Ssub)
            t3 = [(1, v) for (Ls, v) in linevars if len(Ls & Ssub) == 3]
            if r in (7, 8):
                terms = [(math.comb(len(Ls & Ssub), 2), v) for (Ls, v) in linevars if len(Ls & Ssub) >= 3]
                if sum(c for c, _ in terms) > math.comb(r, 2) - 1: le(terms, math.comb(r, 2) - 1)
            if r == 8 and len(t3) > 7: le(t3, 7)
            if r == 9 and len(t3) > 10: le(t3, 10)
for p in PTS:
    others = frozenset(q for q in PTS if q != p)
    rich = [B for B in blocks if p in B]
    subset_caps(others, [(B - {p}, x[B]) for B in rich])
    three = [S for S in lines if p in S and len(S) == 3]
    le([(math.comb(len(B) - 1, 2), x[B]) for B in rich] + [(len(B) - 1, y[B]) for B in rich] + [(3, y[S]) for S in three], math.comb(n, 2) - o(n))
    for r in (7, 8, 9):
        for S0 in itertools.combinations(sorted(others), r - 1):
            S0 = frozenset(S0); pair_terms = []; three_terms = []
            for B in rich:
                k = len((B - {p}) & S0)
                if k >= 3:
                    pair_terms += [(math.comb(k, 2), x[B]), (k, y[B])]
                    if k == 3: three_terms += [(1, x[B]), (-1, y[B])]
                elif k == 2:
                    pair_terms.append((3, y[B])); three_terms.append((1, y[B]))
            for S in three:
                if len((S - {p}) & S0) == 2: pair_terms.append((3, y[S])); three_terms.append((1, y[S]))
            if r in (7, 8) and sum(c for c, _ in pair_terms if c > 0) > math.comb(r, 2) - 1: le(pair_terms, math.comb(r, 2) - 1)
            if r == 8 and sum(c for c, _ in three_terms if c > 0) > 7: le(three_terms, 7)
            if r == 9 and sum(c for c, _ in three_terms if c > 0) > 10: le(three_terms, 10)
    if args.orchard10:
        le([(1, x[B]) for B in rich if len(B) == 4] + [(-1, y[B]) for B in rich if len(B) == 4] + [(1, y[S]) for S in three], T3_CITED[10])
subset_caps(frozenset(PTS), [(S, y[S]) for S in lines])
if args.orchard10: le([(1, y[S]) for S in lines if len(S) == 3], T3_CITED[10])
for B in skel: s.add(x[B])
for B in blocks:
    if len(B) in forbid and B not in skel: s.add(z3.Not(x[B]))
ge([(math.comb(len(B), 3) - 1, x[B]) for B in blocks] + [(1, y[S]) for S in lines], args.thr)
t0 = time.time(); res = s.check()
print(f"z3 skeleton {args.i} {sk['blocks']} omode={args.omode} orchard10={args.orchard10}: {res}  [{time.time()-t0:.0f}s]")
if res == z3.sat:
    mdl = s.model(); F = [sorted(B) for B in blocks if z3.is_true(mdl.eval(x[B], model_completion=True))]
    L = [sorted(S) for S in lines if z3.is_true(mdl.eval(y[S], model_completion=True))]
    print("   blocks", F, "lines", L)
