"""Exhaustive enumeration for skeleton 3 (two disjoint 5-blocks S = {0..4}, R = {5..9}).
Counting (valid in all modes, uses only packing + (M8)/(O9) + t3(10) <= 12):  D = 18 + 3 b4 and every point
lies on one 5-block, so d4(p) <= 8 and b4 <= 20; l <= 12 forces D >= 76, i.e. b4 = 20 and l >= 10.
Blocks of type (1,3),(3,1),(0,4),(4,0) share >= 3 points with S or R, so all 20 four-blocks have two points
in S and two in R.  For an S-pair {a,b} the R-pairs of its blocks are pairwise disjoint (else two blocks
share 3 points), so an S-pair is in <= 2 blocks; 20 blocks on 10 S-pairs gives exactly 2 each; likewise for
R-pairs.  So the block structure is a 2-regular bipartite graph between S-pairs and R-pairs in which the
two neighbours of any vertex are disjoint pairs.
Step 1: enumerate all such graphs (labelled DFS) up to isomorphism (S5 x S5 x swap).
Step 2: for each, enumerate all line systems (CP-SAT): lines from {S, R, the 20 blocks, the 20 uncovered
triples}, pairwise sharing <= 1 point, >= 10 lines, pair budget 45 - o(10) (o(10) = 1: SG only), <= 12
three-lines (t3(10) <= 12), <= 6 / <= 7 three-lines inside any 7- / 8-subset (Fano / (8_3)).
Step 3: isomorphism classes of (blocks, lines); apply all filters (filters.py) and the numerical realiser.
"""
import itertools, json, sys, time, math
from ortools.sat.python import cp_model
sys.path.insert(0, '/home/user/Solving-hard-math-problems-/experiments/506/n10-cpsat')
from n10model import ClassStore, structure_count
from filters import full_report, bundle_violations, fregier_violations
from realise import try_realise

S = list(range(5)); R = list(range(5, 10))
SP = [frozenset(t) for t in itertools.combinations(S, 2)]
RP = [frozenset(t) for t in itertools.combinations(R, 2)]
FIVE = [frozenset(S), frozenset(R)]

# ---- step 1
structs = []
def dfs(i, use, nb_r, chosen):
    if i == len(SP):
        if all(use[Q] == 2 for Q in RP):
            structs.append(list(chosen))
        return
    A = SP[i]
    for Q1, Q2 in itertools.combinations(RP, 2):
        if Q1 & Q2 or use[Q1] >= 2 or use[Q2] >= 2:
            continue
        # S-neighbours of Q1, Q2 must be disjoint from A
        if any(A & A2 for A2 in nb_r[Q1]) or any(A & A2 for A2 in nb_r[Q2]):
            continue
        use[Q1] += 1; use[Q2] += 1; nb_r[Q1].append(A); nb_r[Q2].append(A)
        chosen.append(A | Q1); chosen.append(A | Q2)
        dfs(i + 1, use, nb_r, chosen)
        chosen.pop(); chosen.pop()
        use[Q1] -= 1; use[Q2] -= 1; nb_r[Q1].pop(); nb_r[Q2].pop()

t0 = time.time()
dfs(0, {Q: 0 for Q in RP}, {Q: [] for Q in RP}, [])
print(f"labelled block structures: {len(structs)} ({time.time()-t0:.0f}s)", flush=True)
bstore = ClassStore()
for F4 in structs:
    bstore.add(FIVE + F4, [])
print(f"block-structure isomorphism classes: {len(bstore.classes)}", flush=True)

def fregier_ok(F4):
    """5 four-cycles (two-block lemma applied to S, R)."""
    adj = {}
    for B in F4:
        A = frozenset(B & set(S)); Q = frozenset(B & set(R))
        adj.setdefault(A, set()).add(Q); adj.setdefault(Q, set()).add(A)
    seen = set()
    for A in SP:
        if A in seen: continue
        comp, st = set(), [A]
        while st:
            u = st.pop()
            if u in comp: continue
            comp.add(u); st.extend(adj[u])
        seen |= comp
        ins = [u for u in comp if u in set(SP)]; outs = [u for u in comp if u in set(RP)]
        if any(a & b for a, b in itertools.combinations(ins, 2)) or any(a & b for a, b in itertools.combinations(outs, 2)):
            return False
    return True

# ---- step 2
O10 = int(sys.argv[1]) if len(sys.argv) > 1 else 1
store = ClassStore()
summary = []
for ci, c in enumerate(bstore.classes):
    F = c[0]; F4 = [B for B in F if len(B) == 4]
    covered = set()
    for B in F:
        for T in itertools.combinations(sorted(B), 3):
            covered.add(frozenset(T))
    unc = [frozenset(T) for T in itertools.combinations(range(10), 3) if frozenset(T) not in covered]
    cand = FIVE + F4 + unc
    m = cp_model.CpModel(); y = {i: m.NewBoolVar('') for i in range(len(cand))}
    for i, j in itertools.combinations(range(len(cand)), 2):
        if len(cand[i] & cand[j]) >= 2:
            m.AddBoolOr([y[i].Not(), y[j].Not()])
    m.Add(sum(y.values()) >= 10)
    m.Add(sum(math.comb(len(cand[i]), 2) * y[i] for i in range(len(cand))) <= 45 - O10)
    m.Add(sum(y[i] for i in range(len(cand)) if len(cand[i]) == 3) <= 12)
    for r, cap in ((7, 6), (8, 7)):
        for Sub in itertools.combinations(range(10), r):
            Sub = frozenset(Sub)
            terms = [y[i] for i in range(len(cand)) if len(cand[i] & Sub) == 3]
            if len(terms) > cap: m.Add(sum(terms) <= cap)
    for p in range(10):
        m.Add(sum((len(cand[i]) - 1) * y[i] for i in range(len(cand)) if p in cand[i]) <= 9)
    sols = []
    class CB(cp_model.CpSolverSolutionCallback):
        def on_solution_callback(s): sols.append([cand[i] for i in range(len(cand)) if s.Value(y[i])])
    sv = cp_model.CpSolver(); sv.parameters.enumerate_all_solutions = True; sv.parameters.num_workers = 1
    st = sv.Solve(m, CB())
    n_before = len(store.classes)
    for L in sols:
        store.add(F, L)
    summary.append({'block_class': ci, 'fregier_ok': fregier_ok(F4), 'line_solutions': len(sols), 'status': sv.StatusName(st), 'new_classes': len(store.classes) - n_before})
    print(f"block class {ci}: Fregier-4-cycle structure={fregier_ok(F4)}  line systems status={sv.StatusName(st)} labelled={len(sols)} -> total classes {len(store.classes)} ({time.time()-t0:.0f}s)", flush=True)

json.dump([{'blocks': sorted(sorted(B) for B in c[0]), 'lines': sorted(sorted(S_) for S_ in c[1]), 'copies': c[3],
            'count': structure_count(c[0], c[1])} for c in store.classes], open(f'sk3_classes_o{O10}.json', 'w'))
print(f"total (blocks, lines) classes with >= 10 lines: {len(store.classes)}", flush=True)
# ---- step 3
alive = []
for k, c in enumerate(store.classes):
    F, L = c[0], c[1]
    print(f"class {k}: lines={sorted(len(S_) for S_ in L)} count={structure_count(F, L)} copies={c[3]}")
    rep = full_report(F, L)
    bv = bundle_violations(F, L); fv = fregier_violations(F, L)
    kills = [kk for kk, cc in (('hereditarySG', bool(rep['hereditary_sg'])), ('Fano/8_3', bool(rep['three_line_excess'])),
             ('Miquel', rep['n_miquel'] > 0), ('bundle', len(bv) > 0), ('Fregier', len(fv) > 0)) if cc]
    print(f"   bundle={len(bv)} Fregier={len(fv)} => killed by {kills}", flush=True)
    if not kills:
        alive.append(k)
        r = try_realise(F, L, tries=100, verbose=True)
        print("   SURVIVOR of all filters; numerical realisation:", r, flush=True)
print("survivors of all filters:", alive)
