"""Independent audit pipeline (own code): CP-SAT enumeration per big-block skeleton.

Model (Moebius reformulation): P = n points, blocks = circles-or-lines through >= 3 points; every triple in
exactly one block; rich blocks (size >= 4) pairwise share <= 2 points; lines pairwise share <= 1 point;
a line of size >= 4 is a rich block, a 3-line is a triple in no rich block.
circles = C(n,3) - D - ell,  D = sum_{rich} (C(k,3)-1),  ell = #lines.
Sylvester-Gallai caps (mode 'sg', o(m) = 1 only):
  capD: at each point p, rich blocks through p cover <= C(n-1,2) - 1 pairs of the other points
        (inversion about p: the images of the other n-1 points are not collinear, so some 2-point line exists)
  capL: lines cover <= C(n,2) - 1 pairs (P itself is not collinear)
  capQ: at each point p, the n-point set Q_p = images + p has >= 1 ordinary line: its >= 3-point lines are
        b \ {p} (b rich circle through p), b (b a rich line through p), {a', b', p} for a 3-line {p,a,b}.
"""
import itertools, sys, json, time, os
USE_Q = os.environ.get('NOQ') != '1'   # NOQ=1 disables the Q_p cap (to reproduce the other agents' weaker model)
from math import comb
from ortools.sat.python import cp_model
from canon import canon, aut_group_order

def caps_for(n, mode='sg'):
    if mode == 'sg':
        return dict(capD=comb(n-1, 2) - 1, capL=comb(n, 2) - 1)
    raise ValueError

# ---------------------------------------------------------------- skeletons (blocks of size >= 5)
def enumerate_skeletons(n, capD, verbose=True):
    """all families of blocks of size 5..n-1, pairwise sharing <= 2 points, per-point pair cap, up to iso"""
    subsets = [frozenset(c) for k in range(5, n) for c in itertools.combinations(range(n), k)]
    def ok(fam, B):
        if any(len(B & C) > 2 for C in fam): return False
        for p in B:
            used = sum(comb(len(C) - 1, 2) for C in fam if p in C) + comb(len(B) - 1, 2)
            if used > capD: return False
        return True
    level = {(): []}
    allfams = [[]]
    while level:
        nxt = {}
        for code, fam in level.items():
            for B in subsets:
                if B in fam or not ok(fam, B): continue
                fam2 = fam + [B]
                c = canon(n, fam2)[0]
                if c not in nxt: nxt[c] = fam2
        if verbose: print(f'  skeleton level {len(next(iter(nxt.values()))) if nxt else "-"}: {len(nxt)} classes', flush=True)
        allfams.extend(nxt.values())
        level = nxt
    return allfams

# ---------------------------------------------------------------- CP-SAT completion model
class Collector(cp_model.CpSolverSolutionCallback):
    def __init__(self, xvars, lvars):
        super().__init__(); self.x = xvars; self.l = lvars; self.sols = set(); self.n = 0
    def on_solution_callback(self):
        self.n += 1
        F4 = frozenset(S for S, v in self.x.items() if self.Value(v))
        L = frozenset(S for S, v in self.l.items() if self.Value(v))
        self.sols.add((F4, L))

def complete(n, fixed, target_Dl, caps, extra=None, time_limit=3600, verbose=False):
    """Enumerate all (F4, L): F4 = set of 4-blocks added to the fixed blocks, L = set of lines,
    with D + ell >= target_Dl and all constraints. extra: dict of options
      'no4through': points p such that no 4-block may pass through p (beyond the fixed ones)
      'maxdeg': max number of rich blocks through any point
    Returns (status_name, set of (F4, L))."""
    extra = extra or {}
    capD, capL = caps['capD'], caps['capL']
    fixed = [frozenset(B) for B in fixed]
    m = cp_model.CpModel()
    pts = range(n)
    fixed_triples = {T for B in fixed for T in itertools.combinations(sorted(B), 3)}
    cand4 = [frozenset(S) for S in itertools.combinations(pts, 4)
             if not any(T in fixed_triples for T in itertools.combinations(S, 3))]
    no4 = set(extra.get('no4through', ()))
    cand4 = [S for S in cand4 if not (S & no4)]
    x = {S: m.NewBoolVar(f'x{sorted(S)}') for S in cand4}
    free_triples = [T for T in itertools.combinations(pts, 3) if T not in fixed_triples]
    y3 = {frozenset(T): m.NewBoolVar(f'l{T}') for T in free_triples}
    y4 = {S: m.NewBoolVar(f'l{sorted(S)}') for S in cand4}
    yB = {B: m.NewBoolVar(f'l{sorted(B)}') for B in fixed}
    # every free triple in at most one rich 4-block, and 3-line only if uncovered
    sup = {frozenset(T): [] for T in free_triples}
    for S in cand4:
        for T in itertools.combinations(sorted(S), 3):
            sup[frozenset(T)].append(x[S])
    for T, lst in sup.items():
        m.Add(sum(lst) + y3[T] <= 1)
    for S in cand4:
        m.AddImplication(y4[S], x[S])
    lines = [(T, 3, v) for T, v in y3.items()] + [(S, 4, v) for S, v in y4.items()] + [(B, len(B), v) for B, v in yB.items()]
    # lines pairwise share <= 1 point
    for a, b in itertools.combinations(pts, 2):
        m.Add(sum(v for (Lset, k, v) in lines if a in Lset and b in Lset) <= 1)
    # capD at each point (rich blocks); capQ at each point (Q_p); optional maxdeg
    for p in pts:
        fx = sum(comb(len(B) - 1, 2) for B in fixed if p in B)
        richp = sum(3 * x[S] for S in cand4 if p in S)
        m.Add(richp + fx <= capD)
        qextra = sum((3 if k == 3 else k - 1) * v for (Lset, k, v) in lines if p in Lset)
        if USE_Q: m.Add(richp + fx + qextra <= capL)
        if 'maxdeg' in extra:
            m.Add(sum(x[S] for S in cand4 if p in S) + sum(1 for B in fixed if p in B) <= extra['maxdeg'])
        if 'max4' in extra and p in extra['max4'][0]:   # 4-block degree bound on selected points
            m.Add(sum(x[S] for S in cand4 if p in S) + sum(1 for B in fixed if p in B and len(B) == 4) <= extra['max4'][1])
    # capL
    m.Add(sum(comb(k, 2) * v for (Lset, k, v) in lines) <= capL)
    # target
    Dfix = sum(comb(len(B), 3) - 1 for B in fixed)
    m.Add(Dfix + 3 * sum(x.values()) + sum(v for (_, _, v) in lines) >= target_Dl)
    solver = cp_model.CpSolver()
    solver.parameters.enumerate_all_solutions = True
    solver.parameters.num_workers = 1
    solver.parameters.max_time_in_seconds = time_limit
    lv = {}
    for (Lset, k, v) in lines: lv[Lset] = v
    coll = Collector(x, lv)
    st = solver.Solve(m, coll)
    return solver.StatusName(st), coll.sols

# ---------------------------------------------------------------- exact line enumeration (own)
def line_sets(n, rich, capL, need):
    """all sets L of lines (rich blocks or uncovered triples), pairwise sharing <= 1 point, covering
    <= capL pairs, with Q_p cap at every point, |L| >= need. Returns list of frozensets."""
    rich = [frozenset(B) for B in rich]
    covered = {frozenset(T) for B in rich for T in itertools.combinations(sorted(B), 3)}
    cands = rich + [frozenset(T) for T in itertools.combinations(range(n), 3) if frozenset(T) not in covered]
    cands.sort(key=lambda S: (len(S), sorted(S)))
    richpairs = [sum(comb(len(B) - 1, 2) for B in rich if p in B) for p in range(n)]
    out = []
    def rec(i, chosen, pairs, qp):
        if len(chosen) + (len(cands) - i) < need: return
        if len(chosen) >= need: out.append(frozenset(chosen))
        for j in range(i, len(cands)):
            S = cands[j]
            if any(len(S & C) > 1 for C in chosen): continue
            np_ = pairs + comb(len(S), 2)
            if np_ > capL: continue
            ex = 3 if len(S) == 3 else len(S) - 1
            if USE_Q and any(qp[p] + ex > capL for p in S): continue
            qp2 = list(qp)
            for p in S: qp2[p] += ex
            chosen.append(S); rec(j + 1, chosen, np_, qp2); chosen.pop()
    rec(0, [], 0, list(richpairs))
    return out

# ---------------------------------------------------------------- derived real point sets & checks
O_TABLE = {3: 3, 4: 3, 5: 4, 6: 3, 7: 3, 8: 4, 9: 6, 10: 5, 11: 6, 12: 6, 13: 6, 14: 7}   # cited
T3_TABLE = {3: 1, 4: 1, 5: 2, 6: 4, 7: 6, 8: 7, 9: 10, 10: 12, 11: 16, 12: 19}          # Burr-Gruenbaum-Sloane 1974 (cited)

def derived_sets(n, rich, L):
    """the n+1 real point sets: Q_inf = P with lines L; Q_p = (P - p) + p' with the lines described above.
    Each returned as (m, list of >=3-point lines as frozensets of range(m))."""
    rich = [frozenset(B) for B in rich]; L = [frozenset(l) for l in L]
    out = {'inf': (n, [l for l in L])}
    for p in range(n):
        # relabel: other points keep labels via map, p' = index n-1 after compressing
        others = [q for q in range(n) if q != p]
        idx = {q: i for i, q in enumerate(others)}; pp = n - 1
        lines = []
        Lset = set(L)
        for B in rich:
            if p in B:
                img = frozenset(idx[q] for q in B if q != p)
                lines.append(img | {pp} if B in Lset else img)
        for l in L:
            if p in l and len(l) == 3:
                lines.append(frozenset(idx[q] for q in l if q != p) | {pp})
        out[p] = (n, lines)
    return out

def is_mk_config(m, lines8, S):
    """S (8 points) with restricted lines exactly 8 triples, each point on 3, pairwise <= 1 common -> (8_3)"""
    if len(lines8) != 8 or any(len(l) != 3 for l in lines8): return False
    deg = {v: 0 for v in S}
    for l in lines8:
        for v in l: deg[v] += 1
    if any(d != 3 for d in deg.values()): return False
    return all(len(a & b) <= 1 for a, b in itertools.combinations(lines8, 2))

def check_point_set(m, lines):
    """hereditary checks on a real m-point set with given >=3-point lines (as frozensets).
    returns dict tier -> None (pass) or a witness."""
    res = {'A': None, 'B': None, 'K': None, 'C': None}
    pts = list(range(m))
    for k in range(3, m + 1):
        for S in itertools.combinations(pts, k):
            Sset = frozenset(S)
            rl = [l & Sset for l in lines if len(l & Sset) >= 3]
            if any(len(l) == k for l in rl): continue          # S collinear: no SG statement
            covered = sum(comb(len(l), 2) for l in rl)
            ordinary = comb(k, 2) - covered
            if ordinary < 1 and res['A'] is None: res['A'] = (S, sorted(map(sorted, rl)))
            if ordinary < -(-3 * k // 7) and res['K'] is None: res['K'] = (S, ordinary)
            if k in O_TABLE and ordinary < O_TABLE[k] and res['C'] is None: res['C'] = ('o', S, ordinary)
            t3 = sum(1 for l in rl if len(l) == 3)
            if k in T3_TABLE and t3 > T3_TABLE[k] and res['C'] is None: res['C'] = ('t3', S, t3)
            if k == 8 and res['B'] is None and is_mk_config(m, rl, Sset): res['B'] = (S, sorted(map(sorted, rl)))
    return res

def full_check(n, rich, L):
    """returns dict tier -> list of (which derived set, witness) for tiers failing"""
    fails = {'A': [], 'B': [], 'K': [], 'C': []}
    for key, (m, lines) in derived_sets(n, rich, L).items():
        r = check_point_set(m, lines)
        for t, w in r.items():
            if w is not None: fails[t].append((key, w))
    return fails
