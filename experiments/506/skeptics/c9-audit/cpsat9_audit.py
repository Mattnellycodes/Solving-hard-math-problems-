"""Independent CP-SAT (OR-tools) model of the combinatorial relaxation for n = 9 (c(9) audit).

Variables: x_B (B a subset of size 4..8 is a block), y_S (S a subset of size 3..8 is a line).
Constraints:
  (I)   |B ∩ B'| >= 3  =>  x_B + x_B' <= 1
  (II)  |S ∩ S'| >= 2  =>  y_S + y_S' <= 1 ;  y_S <= x_S for |S| >= 4 ;  y_T + x_B <= 1 for a triple T ⊂ B
  (III) per point p: sum_{B ∋ p} C(|B|-1,2) x_B <= cap_point ;  sum_S C(|S|,2) y_S <= cap_lines
  optional MK:    per point p: sum_{|B|=4, B ∋ p} x_B <= 7           (Lemma 2.4, (8_3) not real)
  optional FANO:  derived: for every p, every 7-set S ⊂ P-{p}, every Fano plane on S with lines T_1..T_7:
                  sum_i sum_{B ⊇ T_i ∪ {p}} x_B <= 6 ;  at infinity: for every 7-set S ⊂ P and Fano plane:
                  sum_i sum_{S' ⊇ T_i} y_S' <= 6 ;  and AG(2,3): not all 12 lines of an affine plane on P.
                  (Hereditary Sylvester-Gallai; these are the only SG-closed configurations on <= 9 points
                  with all lines of size >= 3 -- see audit text.)
  optional forced block {0..m-1} with all other blocks of size <= m (symmetry breaking: relabel the
  largest block), or 'no block of size >= 5' (b5 = 0 case).
Objective: maximise D + ell = sum_B (C(|B|,3)-1) x_B + sum_S y_S ;  circles = 84 - (D + ell).

Usage: python3 cpsat9_audit.py MODE CAPS [time_limit]
  MODE ∈ {opt, opt_mk, opt_full, enum60, m4_mk, m6, m7, m8, m5_full}
  CAPS ∈ {strong (24,30), sg (27,35), none (28,36)}
"""
import sys, itertools, time
from math import comb
from ortools.sat.python import cp_model

N = 9
CAPS = {"strong": (24, 30), "sg": (27, 35), "none": (28, 36)}

def fano_planes(S):
    """all 30 Fano planes (sets of 7 triples) on the 7-set S."""
    S = list(S)
    tr = [frozenset(t) for t in itertools.combinations(S, 3)]
    out = []
    def rec(F, start):
        if len(F) == 7:
            out.append(list(F)); return
        for j in range(start, len(tr)):
            if all(len(tr[j] & u) <= 1 for u in F):
                rec(F + [tr[j]], j + 1)
    rec([], 0)
    return out

def build(caps, forced=None, maxsize=8, minsize=4, mk=False, fano=False, no5=False):
    cap_point, cap_lines = CAPS[caps]
    m = cp_model.CpModel()
    x = {frozenset(B): m.NewBoolVar("x" + "".join(map(str, B)))
         for k in range(4, maxsize + 1) for B in itertools.combinations(range(N), k)}
    y = {frozenset(S): m.NewBoolVar("y" + "".join(map(str, S)))
         for k in range(3, maxsize + 1) for S in itertools.combinations(range(N), k)}
    blocks = list(x); lines = list(y)
    for A, B in itertools.combinations(blocks, 2):
        if len(A & B) >= 3:
            m.Add(x[A] + x[B] <= 1)
    for A, B in itertools.combinations(lines, 2):
        if len(A & B) >= 2:
            m.Add(y[A] + y[B] <= 1)
    for S in lines:
        if len(S) >= 4:
            m.Add(y[S] <= x[S])
        else:
            for B in blocks:
                if S <= B:
                    m.Add(y[S] + x[B] <= 1)
    for p in range(N):
        m.Add(sum(comb(len(B) - 1, 2) * x[B] for B in blocks if p in B) <= cap_point)
    m.Add(sum(comb(len(S), 2) * y[S] for S in lines) <= cap_lines)
    if forced is not None:
        m.Add(x[frozenset(forced)] == 1)
    if no5:
        for B in blocks:
            if len(B) >= 5:
                m.Add(x[B] == 0)
    if mk:
        for p in range(N):
            m.Add(sum(x[B] for B in blocks if len(B) == 4 and p in B) <= 7)
    nf = 0
    if fano:
        for p in range(N):
            others = [q for q in range(N) if q != p]
            for S in itertools.combinations(others, 7):
                for plane in fano_planes(S):
                    m.Add(sum(x[B] for T in plane for B in blocks if (T | {p}) <= B) <= 6); nf += 1
        for S in itertools.combinations(range(N), 7):
            for plane in fano_planes(S):
                m.Add(sum(y[L] for T in plane for L in lines if T <= L) <= 6); nf += 1
        # AG(2,3) at infinity: 12 three-point lines covering all pairs -> at most 11 of them
        pts = list(range(N))
        # all affine planes AG(2,3) on 9 labelled points: generate from the standard one by permutations
        std = []
        for a in range(3):
            std.append({3 * a + b for b in range(3)})            # rows
            std.append({a + 3 * b for b in range(3)})            # columns
        for c in range(3):
            std.append({(3 * r + (r * 1 + c) % 3) for r in range(3)})   # diagonals slope 1
            std.append({(3 * r + (-r + c) % 3) for r in range(3)})      # diagonals slope -1
        std = [frozenset(L) for L in std]
        assert len(set(std)) == 12
        seen = set()
        for perm in itertools.permutations(pts):
            img = frozenset(frozenset(perm[i] for i in L) for L in std)
            if img in seen: continue
            seen.add(img)
            m.Add(sum(y[L] for L in img) <= 11); nf += 1
    obj = sum((comb(len(B), 3) - 1) * x[B] for B in blocks) + sum(y[S] for S in lines)
    return m, x, y, obj, nf

class Collector(cp_model.CpSolverSolutionCallback):
    def __init__(self, x, y):
        super().__init__(); self.x = x; self.y = y; self.sols = []
    def on_solution_callback(self):
        F = tuple(sorted(tuple(sorted(B)) for B, v in self.x.items() if self.Value(v)))
        L = tuple(sorted(tuple(sorted(S)) for S, v in self.y.items() if self.Value(v)))
        self.sols.append((F, L))

def main():
    mode, caps = sys.argv[1], sys.argv[2]
    tl = float(sys.argv[3]) if len(sys.argv) > 3 else 1800.0
    t0 = time.time()
    kw = {}
    if mode == "opt": pass
    elif mode == "opt_mk": kw = dict(mk=True)
    elif mode == "opt_full": kw = dict(mk=True, fano=True)
    elif mode == "enum60": kw = dict(forced=range(5))
    elif mode == "m5_full": kw = dict(forced=range(5), mk=True, fano=True)
    elif mode == "m4_mk": kw = dict(forced=range(4), no5=True, mk=True)
    elif mode == "m6": kw = dict(forced=range(6), maxsize=6)
    elif mode == "m7": kw = dict(forced=range(7), maxsize=7)
    elif mode == "m8": kw = dict(forced=range(8), maxsize=8)
    else: raise SystemExit("unknown mode")
    m, x, y, obj, nf = build(caps, **kw)
    print(f"mode={mode} caps={caps} {CAPS[caps]} build time {time.time()-t0:.1f}s, hereditary constraints {nf}", flush=True)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = tl
    if mode == "enum60":
        m.Add(obj >= 60)
        solver.parameters.enumerate_all_solutions = True
        solver.parameters.num_workers = 1
        col = Collector(x, y)
        st = solver.Solve(m, col)
        fams = {}
        for F, L in col.sols:
            fams.setdefault(F, []).append(L)
        print(f"status={solver.StatusName(st)} solutions(D+ell>=60, block {{0..4}} forced)={len(col.sols)} "
              f"distinct block families={len(fams)} time={time.time()-t0:.1f}s")
        for F, Ls in fams.items():
            D = sum(comb(len(B), 3) - 1 for B in F)
            print(f"  family sizes={sorted((len(B) for B in F), reverse=True)} D={D} line-sets={len(Ls)} "
                  f"ell values={sorted(set(len(L) for L in Ls))} blocks={[list(B) for B in F]}")
        return
    m.Maximize(obj)
    solver.parameters.num_workers = 2
    st = solver.Solve(m)
    print(f"status={solver.StatusName(st)} best D+ell={solver.ObjectiveValue():.0f} bound={solver.BestObjectiveBound():.0f} "
          f"=> min circles in this relaxation = {84 - solver.ObjectiveValue():.0f} (proved lower bound {84 - solver.BestObjectiveBound():.0f}) "
          f"time={time.time()-t0:.1f}s")
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        F = sorted(sorted(B) for B, v in x.items() if solver.Value(v))
        L = sorted(sorted(S) for S, v in y.items() if solver.Value(v))
        print("  blocks:", F); print("  lines:", L)
        degs = [sum(1 for B in F if p in B) for p in range(N)]
        print("  degrees:", degs)

if __name__ == "__main__":
    main()
