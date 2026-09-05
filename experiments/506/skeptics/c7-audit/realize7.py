"""Skeptic c7-audit: direct numerical attempt to realise the Fano-complement circle structure by
7 points of R^2 (independent of the Angle Lemma).  Unknowns: 7 points, gauge p0=(0,0), p1=(1,0)
(similarity), 10 unknowns, 7 quadruple-concyclicity equations (4x4 determinants).  If a real
non-degenerate realisation existed, generic Levenberg-Marquardt runs would converge to it (the
solution set would be a >= 3-dimensional real variety in the gauge-fixed space).  We record every
converged solution and test non-degeneracy: distinct points, no block containing >= 5 points,
not all concyclic/collinear.
"""
import numpy as np, itertools, sys
F = [[0,1,2,3],[0,1,4,5],[0,2,4,6],[0,3,5,6],[1,2,5,6],[1,3,4,6],[2,3,4,5]]
def unpack(w):
    P = np.zeros((7, 2)); P[1, 0] = 1.0; P[2:] = w.reshape(5, 2); return P
def concyc_det(P, idx):
    M = np.array([[P[i,0], P[i,1], P[i,0]**2 + P[i,1]**2, 1.0] for i in idx]); return np.linalg.det(M)
def resid(w):
    P = unpack(w); return np.array([concyc_det(P, B) for B in F])
def jac(w, h=1e-7):
    r0 = resid(w); J = np.zeros((7, w.size))
    for k in range(w.size):
        w2 = w.copy(); w2[k] += h; J[:, k] = (resid(w2) - r0) / h
    return J
def lm(w, iters=200):
    lam = 1e-3
    for _ in range(iters):
        r = resid(w); J = jac(w)
        A = J.T @ J + lam * np.eye(w.size); g = J.T @ r
        try: d = np.linalg.solve(A, -g)
        except np.linalg.LinAlgError: lam *= 10; continue
        r2 = resid(w + d)
        if np.sum(r2**2) < np.sum(r**2): w = w + d; lam = max(lam / 3, 1e-12)
        else: lam *= 5
        if np.sum(resid(w)**2) < 1e-28: break
    return w
def classify(P, tol=1e-6):
    # scale-normalise
    sc = max(np.ptp(P[:,0]), np.ptp(P[:,1]), 1e-300)
    dmin = min(np.linalg.norm(P[i]-P[j]) for i, j in itertools.combinations(range(7), 2)) / sc
    if dmin < 1e-4: return "degenerate: coincident points"
    # points per circle: for each block, count all points on that circle
    def circle(P3):
        (x1,y1),(x2,y2),(x3,y3) = P3
        d = 2*(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))
        if abs(d) < 1e-12*sc**2: return None
        ux = ((x1**2+y1**2)*(y2-y3)+(x2**2+y2**2)*(y3-y1)+(x3**2+y3**2)*(y1-y2))/d
        uy = ((x1**2+y1**2)*(x3-x2)+(x2**2+y2**2)*(x1-x3)+(x3**2+y3**2)*(x2-x1))/d
        return ux, uy, np.hypot(x1-ux, y1-uy)
    for B in F:
        c = circle(P[B[:3]])
        if c is None: return "degenerate: collinear block"
        ux, uy, r = c
        on = sum(1 for i in range(7) if abs(np.hypot(P[i,0]-ux, P[i,1]-uy) - r) < tol * sc)
        if on >= 5: return "degenerate: %d points on a block circle" % on
    return "NON-DEGENERATE REALISATION"
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
N = int(sys.argv[2]) if len(sys.argv) > 2 else 400
summary = {}
best_nondeg = None
for run in range(N):
    w = rng.normal(size=10) * rng.choice([0.5, 1, 3])
    w = lm(w)
    P = unpack(w); res = np.sqrt(np.sum(resid(w)**2))
    if res > 1e-10:
        summary["not converged"] = summary.get("not converged", 0) + 1; continue
    tag = classify(P)
    summary[tag] = summary.get(tag, 0) + 1
    if tag.startswith("NON"):
        best_nondeg = P
print("LM runs:", N, "summary:", summary)
if best_nondeg is not None:
    print("!!! candidate realisation found:", best_nondeg)
