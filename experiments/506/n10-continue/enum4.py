"""Two-stage complete enumeration, for one skeleton (the exact family of blocks of size >= 5), of all abstract
structures (F, L) with D + l >= thr (default 88, i.e. <= 32 circles).

Stage 1 (CP-SAT, all solutions): 4-block families F4 such that
  (T)   every triple lies in <= 1 block (skeleton blocks + 4-blocks) -> blocks pairwise share <= 2 points;
  (SG)  at every p: sum_{B ∋ p} C(|B|-1, 2) <= 35 (derived 9-point set has an ordinary line);
  (MK)  at every p and every 8-subset E of the other 9 points: at most 7 blocks B ∋ p with |(B-p) ∩ E| = 3
        (the (8_3) configuration is not realisable);  and #4-blocks through p <= 10 (t3(9) <= 10, consequence);
  (LB)  D_skel + 3 b4 + 14 >= thr (l <= 14 because lines cover <= 44 pairs, each >= 3).
Symmetry breaking: lex-leader constraints x <=_lex g(x) for a generating set of Aut(skeleton) plus all
transpositions in Aut(skeleton) (valid: the lex-smallest member of every orbit satisfies them all).
Solutions are reduced to isomorphism classes (common.ClassStore).
Stage 2 (own DFS): for each F-class all line sets L with |L| >= thr - D: lines pairwise share <= 1 point, 3-lines
are uncovered triples, rich lines are blocks, sum C(|S|,2) <= 44.  (F, L) reduced to isomorphism classes.
Usage: python3 enum4.py <skeleton index into skeletons_own.json | 'blocks json'> [--thr 88] [--out f.json]
"""
import sys, json, time, argparse, itertools, math
sys.path.insert(0, '/home/user/Solving-hard-math-problems-/experiments/506/n10-continue')
from common import ClassStore, fs, N, count, automorphisms, apply_perm, deficit
from ortools.sat.python import cp_model

def generators(G):
    """small generating set of the permutation group G (list of tuples) by greedy closure."""
    ident = tuple(range(N))
    gens = []
    closure = {ident}
    def compose(a, b):  # a after b
        return tuple(a[b[i]] for i in range(N))
    for g in G:
        if g in closure:
            continue
        gens.append(g)
        # regenerate closure
        closure = {ident}
        frontier = [ident]
        while frontier:
            new = []
            for h in frontier:
                for s in gens:
                    k = compose(s, h)
                    if k not in closure:
                        closure.add(k); new.append(k)
            frontier = new
        if len(closure) == len(G):
            break
    assert len(closure) == len(G)
    return gens

def add_lex(m, xs, ys, tag):
    """x <=_lex y for lists of BoolVars (same length)."""
    n = len(xs)
    eq = m.NewBoolVar(tag + 'e0'); m.Add(eq == 1)
    for i in range(n):
        if xs[i].Index() == ys[i].Index():
            continue          # identical variable: equal, prefix status unchanged
        m.Add(xs[i] <= ys[i]).OnlyEnforceIf(eq)
        b = m.NewBoolVar(''); m.Add(xs[i] == ys[i]).OnlyEnforceIf(b); m.Add(xs[i] != ys[i]).OnlyEnforceIf(b.Not())
        eq2 = m.NewBoolVar('')
        m.AddImplication(eq2, eq); m.AddImplication(eq2, b); m.AddBoolOr([eq.Not(), b.Not(), eq2])
        eq = eq2

def stage1(skel, thr, verbose=True, time_limit=3600, lex=True, cap10=True):
    skel = [fs(B) for B in skel]
    Dsk = sum(deficit(len(B)) for B in skel)
    b4min = max(0, math.ceil((thr - 14 - Dsk) / 3))
    cands = [fs(c) for c in itertools.combinations(range(N), 4) if all(len(fs(c) & B) <= 2 for B in skel)]
    G = automorphisms(skel, [])
    gens = generators(G)
    transp = [g for g in G if sum(1 for i in range(N) if g[i] != i) == 2]
    if verbose:
        print(f"skeleton {sorted(sorted(B) for B in skel)}: D_skel={Dsk}, b4 >= {b4min}, {len(cands)} candidate 4-blocks, |Aut|={len(G)}, "
              f"{len(gens)} generators, {len(transp)} transpositions", flush=True)
    m = cp_model.CpModel()
    x = {B: m.NewBoolVar('x%d' % i) for i, B in enumerate(cands)}
    for T in itertools.combinations(range(N), 3):
        T = fs(T)
        if any(T <= B for B in skel):
            for B in cands:
                if T <= B:
                    m.Add(x[B] == 0)      # cannot happen (filtered), kept for safety
        else:
            m.AddAtMostOne([x[B] for B in cands if T <= B])
    for p in range(N):
        bl = [B for B in cands if p in B]
        base = sum(math.comb(len(B) - 1, 2) for B in skel if p in B)
        m.Add(3 * sum(x[B] for B in bl) + base <= 35)
        if cap10:
            m.Add(sum(x[B] for B in bl) <= 10)
        others = [q for q in range(N) if q != p]
        skel_p = [B - {p} for B in skel if p in B]
        for E in itertools.combinations(others, 8):
            E = fs(E)
            fixed = sum(1 for Ls in skel_p if len(Ls & E) == 3)
            terms = [x[B] for B in bl if (B - {p}) <= E]
            if fixed + len(terms) > 7:
                m.Add(sum(terms) <= 7 - fixed)
    m.Add(sum(x.values()) >= b4min)
    idx = {B: i for i, B in enumerate(cands)}
    xs = [x[B] for B in cands]
    for k, g in enumerate((gens + transp) if lex else []):
        ginv = [0] * N
        for i in range(N): ginv[g[i]] = i
        ys = [x[fs(ginv[p] for p in B)] for B in cands]
        add_lex(m, xs, ys, 'g%d' % k)
    sols = []
    class CB(cp_model.CpSolverSolutionCallback):
        def on_solution_callback(self):
            sols.append([B for B in cands if self.Value(x[B])])
            if len(sols) % 5000 == 0 and verbose:
                print(f"   ... {len(sols)} labelled solutions", flush=True)
    s = cp_model.CpSolver(); s.parameters.enumerate_all_solutions = True; s.parameters.num_workers = 1
    s.parameters.max_time_in_seconds = time_limit
    t0 = time.time()
    st = s.Solve(m, CB())
    status = s.StatusName(st)
    store = ClassStore(N)
    for F4 in sols:
        store.add(skel + F4, [])
    if verbose:
        print(f"   stage 1: status {status}, {len(sols)} labelled solutions, {len(store.classes)} classes  [{time.time()-t0:.0f}s]", flush=True)
    return status, store, Dsk

def line_sets(F, need, cap=44):
    """all sets L of lines with |L| >= need: rich lines from F, 3-lines = uncovered triples; pairwise <= 1 point."""
    F = [fs(B) for B in F]
    covered = set()
    for B in F:
        for T in itertools.combinations(sorted(B), 3):
            covered.add(fs(T))
    cand = sorted(F, key=lambda B: (-len(B), sorted(B))) + [fs(T) for T in itertools.combinations(range(N), 3) if fs(T) not in covered]
    pairs_of = [frozenset(fs(Q) for Q in itertools.combinations(sorted(S), 2)) for S in cand]
    out = []
    def rec(i, chosen, used, npairs):
        if len(chosen) + (len(cand) - i) < need:
            return
        if i == len(cand):
            if len(chosen) >= need:
                out.append(list(chosen))
            return
        # take cand[i]
        if not (pairs_of[i] & used) and npairs + len(pairs_of[i]) <= cap:
            chosen.append(cand[i]); rec(i + 1, chosen, used | pairs_of[i], npairs + len(pairs_of[i])); chosen.pop()
        rec(i + 1, chosen, used, npairs)
    rec(0, [], frozenset(), 0)
    return out

def run(skel, thr=88, out=None, verbose=True, lex=True, cap10=True):
    T0 = time.time()
    status, store, Dsk = stage1(skel, thr, verbose, lex=lex, cap10=cap10)
    result = {'skeleton': sorted(sorted(B) for B in skel), 'thr': thr, 'stage1_status': status, 'F_classes': [], 'structures': []}
    fl = ClassStore(N)
    for k, c in enumerate(store.classes):
        F = c['F']
        D = sum(deficit(len(B)) for B in F)
        need = thr - D
        Ls = line_sets(F, need)
        before = len(fl.classes)
        for L in Ls:
            fl.add(F, L)
        b4 = sum(1 for B in F if len(B) == 4)
        degs = [sum(1 for B in F if p in B) for p in range(N)]
        result['F_classes'].append({'blocks': sorted(sorted(B) for B in F), 'b4': b4, 'D': D, 'copies': c['copies'], 'degrees': degs,
                                    'need_l': need, 'labelled_line_sets': len(Ls), 'new_FL_classes': len(fl.classes) - before})
        if verbose:
            print(f"   F-class {k}: b4={b4} D={D} degrees={degs} copies={c['copies']}: need l>={need}: {len(Ls)} line sets, (F,L) classes so far {len(fl.classes)}  [{time.time()-T0:.0f}s]", flush=True)
    result['structures'] = fl.records()
    if verbose:
        print(f"TOTAL: stage1 {status}; {len(store.classes)} F-classes; {len(fl.classes)} (F,L) classes; time {time.time()-T0:.0f}s")
        for r in result['structures']:
            print("   count", r['count'], "l", len(r['lines']), "blocks", r['blocks'], "lines", r['lines'])
    if out:
        json.dump(result, open(out, 'w'), indent=0)
    return result

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('skel'); ap.add_argument('--thr', type=int, default=88); ap.add_argument('--out', default=None); ap.add_argument('--nolex', action='store_true'); ap.add_argument('--nocap10', action='store_true')
    a = ap.parse_args()
    if a.skel.startswith('['):
        skel = json.loads(a.skel)
    else:
        skel = json.load(open('/home/user/Solving-hard-math-problems-/experiments/506/n10-continue/skeletons_own.json'))[int(a.skel)]['blocks']
    run(skel, a.thr, a.out, lex=not a.nolex, cap10=not a.nocap10)
