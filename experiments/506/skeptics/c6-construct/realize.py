"""Try to REALISE the unique abstract 6-point structure with 7 circles:
circles {0,1,2,3},{0,1,4,5},{2,3,4,5}; lines {0,2,4},{0,3,5},{1,2,5},{1,3,4}.
(a) random-restart least squares (scipy) in a point parametrisation, with a barrier against
    coincident points;  (b) the same in a 4-line parametrisation (6 vertices = pairwise
    intersections, only the 3 concyclicity residuals remain);  (c) z3 nlsat decision procedure.
"""
import numpy as np, sys, time
from scipy.optimize import least_squares

CIRC = [(0, 1, 2, 3), (0, 1, 4, 5), (2, 3, 4, 5)]
LINES = [(0, 2, 4), (0, 3, 5), (1, 2, 5), (1, 3, 4)]
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)


def col(P, a, b, c):
    return (P[b, 0] - P[a, 0]) * (P[c, 1] - P[a, 1]) - (P[b, 1] - P[a, 1]) * (P[c, 0] - P[a, 0])


def cyc(P, a, b, c, d):
    M = np.array([[P[i, 0], P[i, 1], P[i, 0] ** 2 + P[i, 1] ** 2, 1.0] for i in (a, b, c, d)])
    return np.linalg.det(M)


def unpack(v):
    P = np.zeros((6, 2)); P[1] = (1, 0); P[2:] = v.reshape(4, 2); return P


def resid_pts(v):
    P = unpack(v)
    r = [col(P, *t) for t in LINES] + [cyc(P, *q) for q in CIRC]
    # barrier: keep points apart and configuration bounded (scale-free thanks to p0,p1 fixed)
    for i in range(6):
        for j in range(i + 1, 6):
            d2 = np.sum((P[i] - P[j]) ** 2)
            r.append(1e-3 / (d2 + 1e-12))
    return np.array(r)


def minsep(P):
    return min(np.linalg.norm(P[i] - P[j]) for i in range(6) for j in range(i + 1, 6))


print("=== (a) point parametrisation, 2000 random restarts")
best = []
for trial in range(2000):
    v0 = rng.normal(size=8) * 2
    sol = least_squares(resid_pts, v0, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
    P = unpack(sol.x)
    geo = np.array([col(P, *t) for t in LINES] + [cyc(P, *q) for q in CIRC])
    best.append((np.max(np.abs(geo)), minsep(P), P.copy()))
best.sort(key=lambda t: t[0])
for g, s, P in best[:5]:
    print(f"  max|incidence residual|={g:.3e}  min point separation={s:.3e}  points={np.round(P,4).tolist()}")

print("=== (b) line parametrisation: lines y = m_i x + b_i, 5000 random restarts")


def verts(w):
    m = w[:4]; b = w[4:]
    P = {}
    for i in range(4):
        for j in range(i + 1, 4):
            x = (b[j] - b[i]) / (m[i] - m[j]); P[(i, j)] = np.array([x, m[i] * x + b[i]])
    return P


def resid_lines(w):
    P = verts(w)
    # diagonal quadruples: {P12,P34,P13,P24}, {P12,P34,P14,P23}, {P13,P24,P14,P23} (0-based pairs)
    quads = [((0, 1), (2, 3), (0, 2), (1, 3)), ((0, 1), (2, 3), (0, 3), (1, 2)), ((0, 2), (1, 3), (0, 3), (1, 2))]
    r = []
    for q in quads:
        M = np.array([[P[k][0], P[k][1], P[k] @ P[k], 1.0] for k in q]); r.append(np.linalg.det(M))
    # keep slopes apart (non-parallel) and offsets bounded
    m = w[:4]
    for i in range(4):
        for j in range(i + 1, 4):
            r.append(1e-4 / ((m[i] - m[j]) ** 2 + 1e-12))
    r.append(1e-3 * np.sum(w[4:] ** 2))
    return np.array(r)


best = []
for trial in range(5000):
    w0 = np.concatenate([np.tan(rng.uniform(-1.5, 1.5, 4)), rng.normal(size=4)])
    try:
        sol = least_squares(resid_lines, w0, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
    except Exception:
        continue
    P = verts(sol.x)
    pts = np.array(list(P.values()))
    quads = [((0, 1), (2, 3), (0, 2), (1, 3)), ((0, 1), (2, 3), (0, 3), (1, 2)), ((0, 2), (1, 3), (0, 3), (1, 2))]
    geo = []
    for q in quads:
        M = np.array([[P[k][0], P[k][1], P[k] @ P[k], 1.0] for k in q]); geo.append(np.linalg.det(M))
    m = sol.x[:4]
    ang = np.arctan(m)
    minang = min(abs(np.sin(ang[i] - ang[j])) for i in range(4) for j in range(i + 1, 4))
    sep = min(np.linalg.norm(pts[i] - pts[j]) for i in range(6) for j in range(i + 1, 6))
    scale = np.max(np.abs(pts)) + 1
    best.append((np.max(np.abs(geo)) / scale ** 4, minang, sep / scale, np.round(np.degrees(ang), 3).tolist()))
best.sort(key=lambda t: t[0])
for g, a, s, ang in best[:8]:
    print(f"  scaled max|concyclicity residual|={g:.3e}  min|sin(angle between lines)|={a:.3e}  min sep/scale={s:.3e}  angles(deg)={ang}")

print("=== (c) z3 nlsat over the reals")
import z3
X = [z3.Real(f"x{i}") for i in range(6)]; Y = [z3.Real(f"y{i}") for i in range(6)]


def zcol(a, b, c):
    return (X[b] - X[a]) * (Y[c] - Y[a]) - (Y[b] - Y[a]) * (X[c] - X[a])


def zcyc(a, b, c, d):
    # 4x4 determinant expanded via translation to a: rows (x,y,x^2+y^2) relative to a
    rows = []
    for i in (b, c, d):
        dx, dy = X[i] - X[a], Y[i] - Y[a]
        rows.append((dx, dy, dx * dx + dy * dy))
    (a1, a2, a3), (b1, b2, b3), (c1, c2, c3) = rows
    return a1 * (b2 * c3 - b3 * c2) - a2 * (b1 * c3 - b3 * c1) + a3 * (b1 * c2 - b2 * c1)


s = z3.SolverFor("QF_NRA")
s.add(X[0] == 0, Y[0] == 0, X[1] == 1, Y[1] == 0)   # similarity normalisation
for t in LINES:
    s.add(zcol(*t) == 0)
for q in CIRC:
    s.add(zcyc(*q) == 0)
    s.add(zcol(q[0], q[1], q[2]) != 0)   # genuinely a circle, not a line
for i in range(6):
    for j in range(i + 1, 6):
        s.add(z3.Or(X[i] != X[j], Y[i] != Y[j]))
t0 = time.time(); res = s.check(); print("  z3 result:", res, f"({time.time()-t0:.1f}s)")
if res == z3.sat:
    print("  MODEL:", s.model())

print("=== (c') z3 without the 'genuinely a circle' constraints (all 6 collinear allowed?)")
s2 = z3.SolverFor("QF_NRA")
s2.add(X[0] == 0, Y[0] == 0, X[1] == 1, Y[1] == 0)
for t in LINES:
    s2.add(zcol(*t) == 0)
for q in CIRC:
    s2.add(zcyc(*q) == 0)
for i in range(6):
    for j in range(i + 1, 6):
        s2.add(z3.Or(X[i] != X[j], Y[i] != Y[j]))
t0 = time.time(); res = s2.check(); print("  z3 result:", res, f"({time.time()-t0:.1f}s)")
if res == z3.sat:
    print("  MODEL:", s2.model())
