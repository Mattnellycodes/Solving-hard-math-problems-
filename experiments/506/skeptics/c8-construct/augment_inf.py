"""Construction attack 3: Moebius sets containing infinity.
(a) 7 finite points + infinity = 8-point Moebius set Q; for every rich point O (pairwise
    intersections of blocks of Q, which include all 2-point lines since they pass through inf)
    count = |B(Q)| - deg_Q(O).  This reproduces the report's 17 (orthocentric + inf, O on
    altitude+circumcircle+nine-point circle) and looks for less.
(b) 8 finite points + infinity = 9-point Moebius set; drop one finite point q and invert about a
    rich point -> all 8-subsets containing infinity.
Rich points here: intersections of ALL circles/lines through >= 2 points (2-point lines are
blocks of Q; 2-point circles are not blocks but their intersections might still be rich)."""
import numpy as np, itertools, math
import polyhedra as ph
from robust import normalise
def blocks_with_inf(P, tol=1e-9):
    """blocks of P ∪ {inf}: circles through >=3 pts, lines through >=2 pts (lines get member 'inf')."""
    P = np.array(P, float); n = len(P); B = ph.blocks_planar(P, tol=tol)
    # add 2-point lines (not already a >=3-point line)
    lines3 = [b[1] for b in B if b[0] == "L"]
    for i, j in itertools.combinations(range(n), 2):
        if any({i, j} <= l for l in lines3): continue
        B.append(("L", frozenset([i, j]), (P[i], P[j])))
    return B
def count_inf(P, O, B, tol=1e-7):
    deg = sum(1 for b in B if ph.on_block(b, O[0], O[1], tol=tol))
    return len(B) - deg, deg
def best_inv(P, B):
    P = np.array(P, float); n = len(P)
    best = (10**9, None, 0); cands = []
    for p, q in itertools.combinations(B, 2):
        cands += ph.intersections(p, q)
    for O in cands:
        if np.linalg.norm(O) > 1e4: continue
        if any(np.linalg.norm(O - P[t]) < 1e-6 for t in range(n)): continue
        c, d = count_inf(P, O, B)
        c2, d2 = count_inf(P, O, B, tol=1e-10)
        if c != c2: continue
        if c < best[0]: best = (c, O, d)
    return best
def attack7(name, P7):
    P7 = normalise(P7); B = blocks_with_inf(P7)
    c, O, d = best_inv(P7, B)
    print(f"{name} + inf: |B|={len(B)}, best inversion count {c} (centre {np.round(O,5) if O is not None else None}, degree {d})", flush=True)
    return c, O
def attack8(name, P8):
    P8 = normalise(P8); best = 99
    for q in range(8):
        P7 = np.delete(P8, q, axis=0); B = blocks_with_inf(P7); c, O, d = best_inv(P7, B)
        best = min(best, c)
    print(f"{name}: 8 finite + inf, best over dropped point & inversions: {best}", flush=True)
    return best
if __name__ == "__main__":
    A = np.array([0, 3.]); Bp = np.array([1, 0.]); C = np.array([3, 0.]); H = np.array([0, -1.])
    D = np.array([6/5, -3/5]); E = np.array([2, 1.]); F = np.array([0, 0.])
    orth7 = np.array([A, Bp, C, H, D, E, F])
    eq7 = np.array([[0, 0], [1, 0], [0.5, math.sqrt(3)/2], [0.5, 0], [0.25, math.sqrt(3)/4], [0.75, math.sqrt(3)/4], [0.5, math.sqrt(3)/6]])
    attack7("orthocentric7 (report)", orth7)
    attack7("equilateral+midpoints+centroid", eq7)
    rng = np.random.default_rng(11)
    def foot(p, a, b):
        d = b - a; return a + ((p - a) @ d) / (d @ d) * d
    res = []
    for t in range(40):
        A = rng.normal(size=2); Bp = rng.normal(size=2); C = rng.normal(size=2)
        M = np.array([[2*(Bp[0]-A[0]), 2*(Bp[1]-A[1])], [2*(C[0]-A[0]), 2*(C[1]-A[1])]])
        O = np.linalg.solve(M, np.array([Bp@Bp - A@A, C@C - A@A])); H = A + Bp + C - 2*O
        P7 = np.array([A, Bp, C, H, foot(A, Bp, C), foot(Bp, A, C), foot(C, A, Bp)])
        P7 = normalise(P7); B = blocks_with_inf(P7); c, O, d = best_inv(P7, B); res.append(c)
    print("random orthocentric 7 + inf: counts", sorted(set(res)))
    # 7-point sets from the 5x5 grid with 11 circles (Moebius-optimal): 3x3 grid minus 2 corners... use known: (0,0),(0,3),(1,1),(1,2),(2,1),(2,2),(3,0)
    g7 = np.array([(0,0),(0,3),(1,1),(1,2),(2,1),(2,2),(3,0)], float)
    attack7("grid7 (0,0),(0,3),(1,1),(1,2),(2,1),(2,2),(3,0)", g7)
    # random 7-subsets of small grids + inf
    best = 99; import random
    random.seed(5)
    pts = [(x, y) for x in range(5) for y in range(5)]
    for t in range(3000):
        S = random.sample(pts, 7); P7 = np.array(S, float)
        if len({p[0] for p in S}) == 1 or len({p[1] for p in S}) == 1: continue
        B = blocks_with_inf(normalise(P7)); c, O, d = best_inv(normalise(P7), B)
        if c < best: best = c; print("  grid7+inf new best", c, S, flush=True)
    print("random 7-subsets of 5x5 + inf: best", best)
    rec17 = np.array([(0,0),(0,5),(0,10),(0,15),(3,6),(3,9),(5,5),(5,10)], float)
    sq18 = np.array([(1,0),(0,1),(-1,0),(0,-1),(2,0),(0,2),(-2,0),(0,-2)], float)
    attack8("record17", rec17); attack8("two squares", sq18)
    # 3x3 grid minus one point, and 3x3 grid subsets
    g9 = np.array([(x, y) for x in range(3) for y in range(3)], float)
    for q in range(9):
        attack8(f"3x3 grid minus {q}", np.delete(g9, q, axis=0))
