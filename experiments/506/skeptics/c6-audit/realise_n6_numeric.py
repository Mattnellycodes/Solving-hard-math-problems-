"""Numerical attempt (NOT a proof; consistency check only) to realise the unique n = 6 combinatorial
candidate with 7 circles: points 0..5 with
   lines   {0,2,4}, {0,3,5}, {1,2,5}, {1,3,4}      (collinear triples)
   circles {0,1,2,3}, {0,1,4,5}, {2,3,4,5}          (concyclic quadruples)
Unknowns: coordinates of the 6 points with the similarity normalisation p0 = (0,0), p2 = (1,0)
(8 unknowns), 7 polynomial equations (4 collinearity determinants, 3 concyclicity determinants).
Random-restart Levenberg-Marquardt / Gauss-Newton; a genuine solution must have residual ~ 0 AND be
non-degenerate (distinct points, four distinct lines).  We report the best non-degenerate residual.
"""
import numpy as np
import random

LINES = [(0, 2, 4), (0, 3, 5), (1, 2, 5), (1, 3, 4)]
CIRCS = [(0, 1, 2, 3), (0, 1, 4, 5), (2, 3, 4, 5)]


def unpack(z):
    P = np.zeros((6, 2))
    P[2] = (1.0, 0.0)
    P[1] = z[0:2]; P[3] = z[2:4]; P[4] = z[4:6]; P[5] = z[6:8]
    return P


def residuals(z):
    P = unpack(z)
    r = []
    for i, j, k in LINES:
        r.append((P[j, 0] - P[i, 0]) * (P[k, 1] - P[i, 1]) - (P[j, 1] - P[i, 1]) * (P[k, 0] - P[i, 0]))
    for q in CIRCS:
        M = np.array([[P[i, 0], P[i, 1], P[i, 0] ** 2 + P[i, 1] ** 2, 1.0] for i in q])
        r.append(np.linalg.det(M))
    return np.array(r)


def jac(z, h=1e-7):
    r0 = residuals(z)
    J = np.zeros((len(r0), len(z)))
    for i in range(len(z)):
        dz = np.zeros_like(z); dz[i] = h
        J[:, i] = (residuals(z + dz) - residuals(z - dz)) / (2 * h)
    return J


def min_pair_dist(P):
    d = [np.linalg.norm(P[i] - P[j]) for i in range(6) for j in range(i + 1, 6)]
    return min(d)


def lm(z, iters=200):
    lam = 1e-3
    for _ in range(iters):
        r = residuals(z); J = jac(z)
        A = J.T @ J + lam * np.eye(len(z)); g = J.T @ r
        try:
            step = -np.linalg.solve(A, g)
        except np.linalg.LinAlgError:
            break
        z2 = z + step
        if np.sum(residuals(z2) ** 2) < np.sum(r ** 2):
            z = z2; lam = max(lam / 3, 1e-12)
        else:
            lam = min(lam * 5, 1e6)
    return z


if __name__ == "__main__":
    random.seed(7); np.random.seed(7)
    best = []
    for trial in range(3000):
        z = np.random.uniform(-3, 3, 8)
        z = lm(z)
        P = unpack(z)
        res = np.sqrt(np.sum(residuals(z) ** 2))
        md = min_pair_dist(P)
        best.append((res, md, trial))
    best.sort()
    print("smallest residuals (residual, min pairwise distance, trial) -- degenerate = collapsed points:")
    for b in best[:10]:
        print("   ", b)
    nondeg = [b for b in best if b[1] > 1e-2]
    print("best residual among solutions with min pairwise distance > 1e-2:", nondeg[0] if nondeg else None)
