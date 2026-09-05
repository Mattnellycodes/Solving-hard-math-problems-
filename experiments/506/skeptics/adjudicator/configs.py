"""Adjudicator's own checks of the configuration facts used by the lower bounds.

(A) Every family of 8 triples on 8 points pairwise sharing <= 1 point is isomorphic to the
    Moebius-Kantor configuration MK = {i, i+1, i+3 mod 8}  (brute force + backtracking isomorphism).
    Also: no such family of 9 triples exists; every point of such a family has degree exactly 3.
(B) MK has no realisation by 8 distinct real points whose collinear triples are exactly its lines:
    frame {0,1,2,5} (no three on an MK line) -> the remaining incidences force c^2 - c + 1 = 0.
    We derive the coordinates from scratch (sympy), and ALSO run a frame-free numerical check.
(C) Fano plane: the unique 7-triple system on 7 points pairwise sharing <= 1 point covering all
    21 pairs; a real realisation would contradict Sylvester-Gallai (no ordinary line). We also show
    the classical algebraic obstruction det = 2.
(D) Hereditary-SG kill of the unique n=9 count-24 candidate (two 5-blocks sharing a point +
    12 four-blocks): at every degree-7 point p the derived lines restricted to S = P \\ {p, q}
    (q the degree-2 point) form a Fano plane, and S is inside no block through p.  Own ell_max
    recomputation for the candidate (expect 6, count 24).
(E) n=7 Fano-complement family: at every point the derived structure is a complete quadrilateral
    whose three diagonal quadruples are exactly the three blocks avoiding p (own check).
"""
from itertools import combinations, permutations
from math import comb
import sympy as sp

# ---------- (A) ----------
n = 8
triples = [frozenset(t) for t in combinations(range(n), 3)]
fams = []


def rec(F, start):
    if len(F) == 8:
        fams.append(list(F)); return
    for j in range(start, len(triples)):
        t = triples[j]
        if all(len(t & u) <= 1 for u in F):
            rec(F + [t], j + 1)


rec([], 0)
print("(A) labelled families of 8 triples on 8 points pairwise sharing <=1 point:", len(fams))
nine = 0


def rec9(F, start):
    global nine
    if len(F) == 9:
        nine += 1; return
    for j in range(start, len(triples)):
        t = triples[j]
        if all(len(t & u) <= 1 for u in F):
            rec9(F + [t], j + 1)


rec9([], 0)
print("(A) families of 9 such triples:", nine)
MK = [frozenset({i, (i + 1) % 8, (i + 3) % 8}) for i in range(8)]
MKset = set(MK)


def iso_to_MK(F):
    """backtracking search for a bijection phi: {0..7}->{0..7} with phi(F) = MK"""
    Fl = [tuple(sorted(t)) for t in F]
    lines_of = {p: [t for t in Fl if p in t] for p in range(8)}
    phi = {}
    used = set()

    def consistent():
        for t in Fl:
            if all(x in phi for x in t):
                if frozenset(phi[x] for x in t) not in MKset:
                    return False
        return True

    def bt(p):
        if p == 8:
            return True
        for q in range(8):
            if q in used:
                continue
            phi[p] = q; used.add(q)
            if consistent() and bt(p + 1):
                return True
            del phi[p]; used.discard(q)
        return False
    return bt(0)


alliso = all(iso_to_MK(F) for F in fams)
degs_ok = all(all(sum(1 for t in F if p in t) == 3 for p in range(8)) for F in fams)
print("(A) every family isomorphic to MK:", alliso, "; every point has degree exactly 3 in every family:", degs_ok)
print("(A) |Aut(MK)| = 8!/#labelled =", 40320 // len(fams))

# ---------- (B) ----------
a, b, c = sp.symbols('a b c')
frame = [0, 1, 2, 5]
print("(B) frame {0,1,2,5} contains an MK line:", any(frozenset(tr) in MKset for tr in combinations(frame, 3)))
p = {0: sp.Matrix([1, 0, 0]), 1: sp.Matrix([0, 1, 0]), 2: sp.Matrix([0, 0, 1]), 5: sp.Matrix([1, 1, 1])}
# derive p3 = line(0,1) ∩ line(2,5) from the lines {0,1,3} and {2,3,5}
l01 = p[0].cross(p[1]); l25 = p[2].cross(p[5])
p[3] = l01.cross(l25)
print("(B) p3 = line(0,1) x line(2,5) =", list(p[3]))
# p4 on line(1,2) (from {1,2,4}); p6 on line(5,0) (from {5,6,0}); p7 on line(0,2) (from {7,0,2}).
# A point on line(u,v) is s*u + r*v; it equals u iff r = 0, so (points distinct) we may scale r = 1.
p[4] = p[1] * c + p[2]   # (0, c, 1)   (r=1 on p2; p4 = p2 would need c... p4=(0,c,1) never equals p1; equals p2 iff c=0)
p[6] = p[5] + p[0] * (a - 1)   # (a, 1, 1)
p[7] = p[0] * b + p[2]   # (b, 0, 1)
eqs = []
for L in MK:
    i, j, k = sorted(L)
    e = sp.expand(sp.Matrix.hstack(p[i], p[j], p[k]).det())
    print(f"    line {sorted(L)}: incidence polynomial {e}")
    if e != 0:
        eqs.append(e)
G = sp.groebner(eqs, a, b, c, order='lex')
print("(B) Groebner basis (lex a>b>c):", list(G))
uni = [g for g in G if g.free_symbols <= {c}]
print("(B) univariate:", uni, " real roots:", [sp.real_roots(sp.Poly(g, c)) for g in uni])
sols = sp.solve(eqs, [a, b, c], dict=True)
print("(B) all complex solutions:", sols)
for s in sols:
    pts = {k: [sp.simplify(x.subs(s)) for x in v] for k, v in p.items()}
    distinct = len({tuple(v) for v in pts.values()}) == 8   # (projective distinctness up to scale checked loosely)
    print("     solution", s, "-> 8 distinct coordinate vectors:", distinct)
# Coincidence check: which parameter values would make a point coincide with a frame point?
print("(B) p4=(0,c,1) equals p2 iff c=0 ; p6=(a,1,1) never equals p5 unless a=1 ; p7=(b,0,1) equals p2 iff b=0 ; none of these are roots of c^2-c+1 (c=0? no).")
# frame-free numerical sanity: random least squares on 8 real points for the 8 collinearities
import random
random.seed(3)
import numpy as np
from scipy.optimize import least_squares


def resid(v):
    P = v.reshape(8, 2)
    r = []
    for L in MK:
        i, j, k = sorted(L)
        r.append((P[j, 0] - P[i, 0]) * (P[k, 1] - P[i, 1]) - (P[j, 1] - P[i, 1]) * (P[k, 0] - P[i, 0]))
    # normalisation: fix scale by penalising deviation of pairwise distances from being >= 0.3 and bounded
    return np.array(r)


bestsep = 0.0
for trial in range(300):
    v0 = np.random.RandomState(trial).uniform(-1, 1, 16)
    sol = least_squares(resid, v0, xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=2000)
    P = sol.x.reshape(8, 2)
    d = min(np.hypot(*(P[i] - P[j])) for i, j in combinations(range(8), 2))
    # collinearity of everything?
    M = np.c_[P, np.ones(8)]
    rank = np.linalg.matrix_rank(M, tol=1e-6)
    if np.max(np.abs(sol.fun)) < 1e-9 and rank == 3:
        bestsep = max(bestsep, d)
print("(B) numerical: best min-separation among zero-residual NON-collinear 8-point solutions over 300 restarts:", bestsep,
      "(a real realisation would give a value bounded away from 0)")

# ---------- (C) ----------
t7 = [frozenset(t) for t in combinations(range(7), 3)]
fano = []


def rec7(F, start):
    if len(F) == 7:
        fano.append(F); return
    for j in range(start, len(t7)):
        t = t7[j]
        if all(len(t & u) <= 1 for u in F):
            rec7(F + [t], j + 1)


rec7([], 0)
print("(C) labelled 7-triple systems on 7 points pairwise <=1:", len(fano),
      "; all cover all 21 pairs:", all(len({frozenset(pr) for t in F for pr in combinations(t, 2)}) == 21 for F in fano))
# classical algebraic obstruction: Fano lines {0,1,2},{0,3,4},{0,5,6},{1,3,5},{1,4,6},{2,3,6},{2,4,5}
FL = [(0, 1, 2), (0, 3, 4), (0, 5, 6), (1, 3, 5), (1, 4, 6), (2, 3, 6), (2, 4, 5)]
q = {3: sp.Matrix([1, 0, 0]), 4: sp.Matrix([0, 1, 0]), 5: sp.Matrix([0, 0, 1]), 6: sp.Matrix([1, 1, 1])}
q[0] = q[3].cross(q[4]).cross(q[5].cross(q[6]))
q[1] = q[3].cross(q[5]).cross(q[4].cross(q[6]))
q[2] = q[3].cross(q[6]).cross(q[4].cross(q[5]))
print("(C) forced Fano coordinates:", {k: list(v) for k, v in q.items()},
      "; det(p0,p1,p2) =", sp.Matrix.hstack(q[0], q[1], q[2]).det(), "(must be 0 for the 7th line: impossible over R)")

# ---------- (D) ----------
cand = [[0, 1, 2, 3, 4], [0, 1, 5, 6], [2, 3, 5, 6], [0, 2, 5, 7], [1, 3, 5, 7], [1, 2, 6, 7], [0, 3, 6, 7],
        [1, 2, 5, 8], [0, 3, 5, 8], [0, 2, 6, 8], [1, 3, 6, 8], [0, 1, 7, 8], [2, 3, 7, 8], [4, 5, 6, 7, 8]]
cand = [frozenset(B) for B in cand]
n9 = 9
assert all(len(A & B) <= 2 for A, B in combinations(cand, 2))
deg = {p: sum(1 for B in cand if p in B) for p in range(n9)}
print("(D) candidate degrees:", deg, " D =", sum(comb(len(B), 3) - 1 for B in cand))
for pnt in range(n9):
    if deg[pnt] < 7:
        continue
    derived = [B - {pnt} for B in cand if pnt in B]
    qdeg2 = [x for x in range(n9) if deg[x] == 2][0]
    S = frozenset(range(n9)) - {pnt, qdeg2}
    lines_S = [L & S for L in derived if len(L & S) >= 3]
    pairs = {frozenset(pr) for L in lines_S for pr in combinations(L, 2)}
    inside_block = any(S <= (B - {pnt}) for B in cand if pnt in B)
    print(f"    p={pnt}: S={sorted(S)} restricted lines={[sorted(L) for L in lines_S]} "
          f"#lines={len(lines_S)} sizes={[len(L) for L in lines_S]} pairs covered={len(pairs)}/21 "
          f"S inside one block through p: {inside_block}")


def ell_max9(F, cap_ln):
    tri = [frozenset(t) for t in combinations(range(n9), 3)]
    lc = list(F) + [t for t in tri if not any(t <= B for B in F)]
    best = [0]

    def rec_(i, chosen, pairs):
        best[0] = max(best[0], len(chosen))
        for j in range(i, len(lc)):
            if len(chosen) + len(lc) - j <= best[0]:
                return
            Sx = lc[j]
            if pairs + comb(len(Sx), 2) > cap_ln:
                continue
            if all(len(Sx & T) <= 1 for T in chosen):
                chosen.append(Sx); rec_(j + 1, chosen, pairs + comb(len(Sx), 2)); chosen.pop()
    rec_(0, [], 0)
    return best[0]


e = ell_max9(cand, 36)
print("(D) own ell_max (no cap) =", e, " count =", comb(9, 3) - 54 - e)

# ---------- (E) ----------
biplane = [frozenset(B) for B in [[0, 1, 2, 3], [0, 1, 4, 5], [0, 2, 4, 6], [0, 3, 5, 6], [1, 2, 5, 6], [1, 3, 4, 6], [2, 3, 4, 5]]]
ok_all = True
for pnt in range(7):
    lines = [B - {pnt} for B in biplane if pnt in B]
    avoid = [B for B in biplane if pnt not in B]
    others = frozenset(range(7)) - {pnt}
    quad = (len(lines) == 4 and all(len(L) == 3 for L in lines)
            and all(len(A & B) == 1 for A, B in combinations(lines, 2))
            and all(sum(1 for L in lines if x in L) == 2 for x in others))
    # diagonal quadruples: for each pair of disjoint... vertex P_ij = L_i ∩ L_j ; opposite pairs (P_ij, P_kl)
    idx = list(range(4))
    V = {}
    for i, j in combinations(idx, 2):
        V[(i, j)] = next(iter(lines[i] & lines[j]))
    opp = [({V[(0, 1)], V[(2, 3)]}), ({V[(0, 2)], V[(1, 3)]}), ({V[(0, 3)], V[(1, 2)]})]
    diag = {frozenset(opp[0] | opp[1]), frozenset(opp[0] | opp[2]), frozenset(opp[1] | opp[2])}
    same = diag == set(avoid)
    ok_all &= quad and same
    print(f"(E) p={pnt}: complete quadrilateral={quad}; blocks avoiding p are exactly the 3 diagonal quadruples={same}")
print("(E) reduction to the Angle Lemma holds at every point:", ok_all)
