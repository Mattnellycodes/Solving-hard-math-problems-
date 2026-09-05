"""Stage 1 (auditor's own): for a fixed skeleton (the exact family of blocks of size >= 5), enumerate
ALL families of 4-blocks such that

  (T)   every triple lies in <= 1 block (4-blocks and skeleton blocks together)
  (SG)  for every p: sum_{B ∋ p} C(|B|-1, 2) <= C(9,2) - 1 = 35        (Sylvester-Gallai at p)
  (MK)  for every p and every q != p: the number of blocks B ∋ p whose derived line B - p meets
        E = P - {p, q} in exactly 3 points is <= 7                       ((8_3) is not real)
  (LB)  D_skel + 3 b4 >= 88 - 14                                         (l <= 14 by SG on lines)

up to the automorphism group of the skeleton, by CP-SAT all-solution enumeration with lex-leader
constraints x <=_lex g.x for a subset of Aut(skeleton) (sound for any subset: the lex-smallest
member of every orbit satisfies all of them), followed by exact isomorphism reduction.

usage: python3 enum_f4.py <skeleton index into skeletons_audit.json> [--nolex] [--nomk] [--t3cap] [--first]
"""
import sys, json, time, itertools, random
from math import comb
from ortools.sat.python import cp_model
from common import N, automorphisms, ClassCollector, D_of

args = sys.argv[1:]
idx = int(args[0])
NOLEX = '--nolex' in args
NOMK = '--nomk' in args
T3CAP = '--t3cap' in args
FIRST = '--first' in args
SEED = 1
for a in args:
    if a.startswith('--seed='):
        SEED = int(a[7:])
skel = [] if idx < 0 else [frozenset(B) for B in json.load(open('skeletons_audit.json'))[idx]]   # idx = -1: no block of size >= 5
assert all(len(B) <= 6 for B in skel)
Dsk = D_of(skel)
t0 = time.time()
print('skeleton', idx, [sorted(B) for B in skel], 'D_skel =', Dsk, flush=True)

cands = [frozenset(c) for c in itertools.combinations(range(N), 4)]
cands = [c for c in cands if all(len(c & B) <= 2 for B in skel)]
cid = {c: i for i, c in enumerate(cands)}
print('candidate 4-blocks:', len(cands))

m = cp_model.CpModel()
x = [m.NewBoolVar('x%d' % i) for i in range(len(cands))]

# (T)
for T in itertools.combinations(range(N), 3):
    Ts = frozenset(T)
    if any(Ts <= B for B in skel):
        continue
    lst = [x[cid[c]] for c in cands if Ts <= c]
    if len(lst) > 1:
        m.Add(sum(lst) <= 1)
# (SG)
for p in range(N):
    base = sum(comb(len(B) - 1, 2) for B in skel if p in B)
    m.Add(3 * sum(x[cid[c]] for c in cands if p in c) <= 35 - base)
    if T3CAP:
        m.Add(sum(x[cid[c]] for c in cands if p in c) <= 10)
# (MK)
if not NOMK:
    for p in range(N):
        for q in range(N):
            if q == p:
                continue
            E = frozenset(range(N)) - {p, q}
            base = sum(1 for B in skel if p in B and len((B - {p}) & E) == 3)
            lst = [x[cid[c]] for c in cands if p in c and len((c - {p}) & E) == 3]
            m.Add(sum(lst) + base <= 7)
# (LB)
m.Add(3 * sum(x) >= 88 - 14 - Dsk)

# symmetry breaking
random.seed(SEED)
if skel:
    aut = automorphisms(skel, [])
    print('|Aut(skeleton)| =', len(aut), flush=True)
else:
    # Aut = S_10: transpositions, 3-cycles, double transpositions and random permutations
    aut = set()
    for a, b in itertools.combinations(range(N), 2):
        g = list(range(N)); g[a], g[b] = b, a; aut.add(tuple(g))
    for a, b, c in itertools.permutations(range(N), 3):
        g = list(range(N)); g[a], g[b], g[c] = b, c, a; aut.add(tuple(g))
    for (a, b), (c, d) in itertools.combinations(list(itertools.combinations(range(N), 2)), 2):
        if len({a, b, c, d}) == 4:
            g = list(range(N)); g[a], g[b], g[c], g[d] = b, a, d, c; aut.add(tuple(g))
    for _ in range(3000):
        g = list(range(N)); random.shuffle(g); aut.add(tuple(g))
    aut = sorted(aut)
    print('empty skeleton: Aut = S_10, using', len(aut), 'elements')
if NOLEX:
    use = []
elif len(aut) <= 1500:
    use = [g for g in aut if any(g[i] != i for i in range(N))]
else:
    def nmoved(g):
        return sum(1 for i in range(N) if g[i] != i)
    small = [g for g in aut if 0 < nmoved(g) <= (3 if not skel else 4)]
    rest = [g for g in aut if nmoved(g) > (3 if not skel else 4)]
    use = small + random.sample(rest, min(400, len(rest)))
print('lex-leader constraints for', len(use), 'group elements', flush=True)


def apply(g, c):
    return frozenset(g[i] for i in c)


order = list(range(len(cands)))   # variable order for lex comparison
naux = 0
for g in use:
    ginv = [0] * N
    for i in range(N):
        ginv[g[i]] = i
    # (g.x)_B = x_{g^{-1} B}
    pos = [(k, cid[apply(ginv, cands[k])]) for k in order]
    pos = [(a, b) for a, b in pos if a != b]
    e = m.NewBoolVar('e'); m.Add(e == 1)
    for a, b in pos:
        v, w = x[a], x[b]
        # e -> v <= w
        m.AddBoolOr([e.Not(), v.Not(), w])
        eq = m.NewBoolVar('q')
        m.AddBoolOr([v, w, eq]); m.AddBoolOr([v.Not(), w.Not(), eq])
        m.AddBoolOr([v.Not(), w, eq.Not()]); m.AddBoolOr([v, w.Not(), eq.Not()])
        e2 = m.NewBoolVar('e')
        m.AddBoolOr([e.Not(), eq.Not(), e2]); m.AddImplication(e2, e); m.AddImplication(e2, eq)
        naux += 2
        e = e2
print('aux vars:', naux, 'build time %.1fs' % (time.time() - t0), flush=True)


class CB(cp_model.CpSolverSolutionCallback):
    def __init__(self):
        super().__init__()
        self.sols = []

    def on_solution_callback(self):
        self.sols.append([i for i in range(len(cands)) if self.Value(x[i])])
        if FIRST:
            self.StopSearch()
        if len(self.sols) % 500 == 0:
            print('  ...', len(self.sols), 'solutions  %.0fs' % (time.time() - t0), flush=True)


solver = cp_model.CpSolver()
solver.parameters.enumerate_all_solutions = True
solver.parameters.num_workers = 1
solver.parameters.max_time_in_seconds = 2400
cb = CB()
st = solver.Solve(m, cb)
print('status:', solver.StatusName(st), 'labelled solutions:', len(cb.sols), ' %.0fs' % (time.time() - t0), flush=True)
assert FIRST or st in (cp_model.OPTIMAL, cp_model.INFEASIBLE), 'INCOMPLETE'

coll = ClassCollector()
for s in cb.sols:
    F = [sorted(B) for B in skel] + [sorted(cands[i]) for i in s]
    coll.add(F, [])
cl = coll.classes()
print('F-classes:', len(cl))
out = []
for F, _, cnt in cl:
    deg = [sum(1 for B in F if p in B and len(B) == 4) for p in range(N)]
    b4 = sum(1 for B in F if len(B) == 4)
    print('  b4=%d D=%d deg4=%s labelled=%d' % (b4, D_of(F), deg, cnt))
    out.append({'F': F, 'b4': b4, 'D': D_of(F), 'labelled': cnt})
tag = ('empty' if idx < 0 else str(idx)) + ('_seed%d' % SEED if SEED != 1 else '') + ('_nolex' if NOLEX else '') + ('_nomk' if NOMK else '') + ('_t3cap' if T3CAP else '') + ('_first' if FIRST else '')
json.dump({'skeleton': [sorted(B) for B in skel], 'status': solver.StatusName(st),
           'labelled': len(cb.sols), 'classes': out}, open('f4_%s.json' % tag, 'w'))
print('done %.0fs' % (time.time() - t0))
