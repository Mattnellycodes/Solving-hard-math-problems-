"""Construction attack 2: rich-point augmentation.
Start from good 7-point sets (11 circles) and the 8-point 17/18-sets; candidate new points =
all pairwise intersections of blocks, all circle centres, midpoints, reflections of points in
lines, orthocentres/circumcentres of triples.  For 7+1: evaluate the 8-set (Moebius-optimised).
For 8+1 = 9-point Moebius sets Q: evaluate every 8-subset (Moebius-optimised) — this covers
'subsets of rich 9-point configurations'."""
import numpy as np, itertools, math
import polyhedra as ph
from robust import evaluate, normalise
def rich_points(P):
    P = normalise(P); B = ph.blocks_planar(P); pts = []
    for p, q in itertools.combinations(B, 2):
        pts += ph.intersections(p, q)
    for b in B:
        if b[0] == "C": pts.append(np.array(b[2][:2]))
    for i, j in itertools.combinations(range(len(P)), 2):
        pts.append((P[i] + P[j]) / 2)
    for i, j, k in itertools.permutations(range(len(P)), 3):
        a, b, c = P[i], P[j], P[k]; d = b - a
        if np.linalg.norm(d) < 1e-9: continue
        t = ((c - a) @ d) / (d @ d); foot = a + t * d
        pts.append(foot); pts.append(2 * foot - c)
        # circumcentre / orthocentre of triangle
        M = np.array([[2*(b[0]-a[0]), 2*(b[1]-a[1])], [2*(c[0]-a[0]), 2*(c[1]-a[1])]])
        rhs = np.array([b@b - a@a, c@c - a@a])
        if abs(np.linalg.det(M)) > 1e-9:
            O = np.linalg.solve(M, rhs); pts.append(O); pts.append(a + b + c - 2 * O)
    # dedupe
    out = []
    for p in pts:
        if np.linalg.norm(p) > 50: continue
        if any(np.linalg.norm(p - q) < 1e-7 for q in out): continue
        if any(np.linalg.norm(p - q) < 1e-7 for q in P): continue
        out.append(p)
    return P, out
def best_add(name, P7, target):
    P, cand = rich_points(P7)
    print(f"{name}: {len(cand)} candidate points", flush=True)
    best = 99; bestP = None
    for c in cand:
        Q = np.vstack([P, c])
        r = evaluate(Q)
        if r is None or r[0] == "UNSTABLE": continue
        if r[0] < best: best = r[0]; bestP = Q
        if r[0] <= target: print("  !!", r[0], "euclid", r[1], "with new point", np.round(c, 5), flush=True)
    print(f"  best 8-point count after augmentation: {best}")
    return bestP
def best_sub(name, P9):
    best = 99
    for S in itertools.combinations(range(len(P9)), 8):
        r = evaluate(P9[list(S)])
        if r is None or r[0] == "UNSTABLE": continue
        if r[0] < best: best = r[0]; bs = S
    print(f"{name}: best 8-subset count {best}")
    return best
if __name__ == "__main__":
    # 7-point record sets
    A = np.array([0, 3.]); Bp = np.array([1, 0.]); C = np.array([3, 0.]); H = np.array([0, -1.])
    D = np.array([6/5, -3/5]); E = np.array([2, 1.]); F = np.array([0, 0.])
    orth7 = np.array([A, Bp, C, H, D, E, F])
    eq7 = np.array([[0, 0], [1, 0], [0.5, math.sqrt(3)/2], [0.5, 0], [0.25, math.sqrt(3)/4], [0.75, math.sqrt(3)/4], [0.5, math.sqrt(3)/6]])
    print("7-point records:", evaluate(orth7)[:2], evaluate(eq7)[:2])
    best_add("orthocentric 7 (A,B,C,H,D,E,F)", orth7, 16)
    best_add("equilateral+midpoints+centroid", eq7, 16)
    # generic orthocentric systems (random triangle): A,B,C,H,D,E,F
    rng = np.random.default_rng(3)
    for t in range(6):
        A = rng.normal(size=2); Bp = rng.normal(size=2); C = rng.normal(size=2)
        def foot(p, a, b):
            d = b - a; return a + ((p - a) @ d) / (d @ d) * d
        Dp = foot(A, Bp, C); Ep = foot(Bp, A, C); Fp = foot(C, A, Bp)
        M = np.array([[2*(Bp[0]-A[0]), 2*(Bp[1]-A[1])], [2*(C[0]-A[0]), 2*(C[1]-A[1])]])
        O = np.linalg.solve(M, np.array([Bp@Bp - A@A, C@C - A@A])); H = A + Bp + C - 2*O
        P7 = np.array([A, Bp, C, H, Dp, Ep, Fp])
        print("random orthocentric 7:", evaluate(P7)[:2])
        best_add(f"random orthocentric {t}", P7, 16)
    # 8-point 17-sets and 18-sets -> 9 -> subsets
    rec17 = np.array([(0,0),(0,5),(0,10),(0,15),(3,6),(3,9),(5,5),(5,10)], float)
    sq18 = np.array([(1,0),(0,1),(-1,0),(0,-1),(2,0),(0,2),(-2,0),(0,-2)], float)
    for name, P8 in [("record17", rec17), ("two squares 18", sq18)]:
        P, cand = rich_points(P8)
        print(f"{name}: {len(cand)} candidate 9th points", flush=True)
        best = 99
        for c in cand[:400]:
            Q = np.vstack([P, c])
            for S in itertools.combinations(range(9), 8):
                if 8 not in S: continue
                r = evaluate(Q[list(S)])
                if r is None or r[0] == "UNSTABLE": continue
                if r[0] < best: best = r[0]
                if r[0] <= 16: print("  !!", r, S, np.round(c, 5), flush=True)
        print(f"  best 8-subset containing the new point: {best}")
