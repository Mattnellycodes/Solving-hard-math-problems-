"""Construction attack 1: stereographic projections / inversions of point sets with rich
Moebius structure.  For a Moebius set Q of 8 points (on the sphere), the best planar count is
|B(Q)| - max_O deg(O) over all sphere points O (O in Q gives 7 points; we want O not in Q).
Candidates: cube-family (boxes), square antiprism, triangular prism + 2 apexes, two squares
rotated, 8 points of a regular n-gon (n>8), orthocentric systems + feet + O, vertices+centres.
Also planar sets: for each we evaluate all inversion centres (pairwise block intersections)."""
import numpy as np, itertools, math, sys
from realise import count_float

def blocks_planar(P, tol=1e-8):
    n = len(P); res = []; done = []
    for i, j, k in itertools.combinations(range(n), 3):
        mk = {i, j, k}
        if any(mk <= d for d in done): continue
        a, b, c = P[i], P[j], P[k]
        d = 2*(a[0]*(b[1]-c[1]) + b[0]*(c[1]-a[1]) + c[0]*(a[1]-b[1]))
        scale = max(np.linalg.norm(b-a), np.linalg.norm(c-a), np.linalg.norm(c-b))
        if abs(d) < tol*scale*scale*10:
            dx, dy = b-a
            memb = {t for t in range(n) if abs((P[t][0]-a[0])*dy - (P[t][1]-a[1])*dx) < tol*scale*(1+np.linalg.norm(P[t]-a))}
            res.append(("L", frozenset(memb), (a, b)))
        else:
            a2 = a@a; b2 = b@b; c2 = c@c
            ux = (a2*(b[1]-c[1]) + b2*(c[1]-a[1]) + c2*(a[1]-b[1]))/d
            uy = (a2*(c[0]-b[0]) + b2*(a[0]-c[0]) + c2*(b[0]-a[0]))/d
            r = math.hypot(a[0]-ux, a[1]-uy)
            memb = {t for t in range(n) if abs(math.hypot(P[t][0]-ux, P[t][1]-uy) - r) < tol*(1+r)}
            res.append(("C", frozenset(memb), (ux, uy, r)))
        done.append(memb)
    return res

def on_block(b, x, y, tol=1e-7):
    if b[0] == "L":
        a, c = b[2]; dx, dy = c-a
        return abs((x-a[0])*dy - (y-a[1])*dx) < tol*(1+abs(dx)+abs(dy))*(1+abs(x)+abs(y))
    ux, uy, r = b[2]
    return abs(math.hypot(x-ux, y-uy) - r) < tol*(1+r)

def intersections(p, q):
    if p[0] == "L" and q[0] == "L":
        a, b = p[2]; c, d = q[2]; d1 = b-a; d2 = d-c; den = d1[0]*d2[1]-d1[1]*d2[0]
        if abs(den) < 1e-12: return []
        t = ((c[0]-a[0])*d2[1] - (c[1]-a[1])*d2[0])/den
        return [a + t*d1]
    if p[0] == "L" or q[0] == "L":
        L, C = (p, q) if p[0] == "L" else (q, p)
        a, b = L[2]; d = (b-a)/np.linalg.norm(b-a); ux, uy, r = C[2]
        f = a - np.array([ux, uy]); t0 = -(f@d); h = f + t0*d; disc = r*r - h@h
        if disc < -1e-9: return []
        s = math.sqrt(max(disc, 0)); return [a + (t0+s)*d, a + (t0-s)*d]
    (ux, uy, r1), (vx, vy, r2) = p[2], q[2]
    dx, dy = vx-ux, vy-uy; d = math.hypot(dx, dy)
    if d < 1e-12: return []
    a = (r1*r1 - r2*r2 + d*d)/(2*d); h2 = r1*r1 - a*a
    if h2 < -1e-9: return []
    h = math.sqrt(max(h2, 0)); mx, my = ux + a*dx/d, uy + a*dy/d
    return [np.array([mx + h*dy/d, my - h*dx/d]), np.array([mx - h*dy/d, my + h*dx/d])]

def moebius_best(P):
    """Return (best count, |B|, best centre, degree) for the planar set P: |B| - max external degree.
    Also returns the Euclidean count = |B| - #lines."""
    B = blocks_planar(P); nb = len(B); nl = sum(1 for b in B if b[0] == "L")
    best = (nb - nl, None, 0)
    for p, q in itertools.combinations(B, 2):
        for pt in intersections(p, q):
            if any(np.linalg.norm(pt - P[t]) < 1e-6 for t in range(len(P))): continue
            deg = sum(1 for b in B if on_block(b, pt[0], pt[1]))
            if nb - deg < best[0]: best = (nb - deg, pt, deg)
    return best, nb, nb - nl

def invert(P, O, R2=1.0):
    Q = []
    for p in P:
        v = p - O; Q.append(O + R2*v/(v@v))
    return np.array(Q)

def report(name, P):
    P = np.array(P, float)
    (bc, O, deg), nb, euc = moebius_best(P)
    print(f"{name}: |B|={nb} euclid={euc} best-inversion count={bc} (centre {None if O is None else np.round(O,4)}, degree {deg})")
    return bc

if __name__ == "__main__":
    results = []
    # 1. Box family (cube structure): points = stereographic projections of box vertices for many boxes
    def box_proj(a, b, c, view):
        # vertices of box (+-a,+-b,+-c) on sphere radius R; stereographic projection from a generic
        # sphere point = inversion; instead: project vertices to plane via a Moebius map: we realise the
        # cube structure directly by picking a generic sphere point N and projecting.
        V = np.array(list(itertools.product([-a, a], [-b, b], [-c, c])), float)
        R = math.sqrt(a*a + b*b + c*c); V /= R
        N = np.array(view, float); N /= np.linalg.norm(N)
        # rotate N to (0,0,1)
        z = N; xx = np.cross(z, [0.3, 0.5, 0.8]); xx /= np.linalg.norm(xx); yy = np.cross(z, xx)
        M = np.array([xx, yy, z])
        W = V @ M.T
        return np.array([[w[0]/(1-w[2]), w[1]/(1-w[2])] for w in W])
    rng = np.random.default_rng(1)
    best = 99
    for trial in range(300):
        a, b, c = rng.uniform(0.5, 2, 3); view = rng.normal(size=3)
        P = box_proj(a, b, c, view)
        (bc, O, deg), nb, euc = moebius_best(P)
        best = min(best, bc)
    print("random boxes (cube structure), best over 300 random boxes and all inversion centres:", best)
    # 2. Regular polygon subsets: 8 vertices of a regular m-gon (m=8..16), plus inversions
    for m in range(8, 17):
        for S in itertools.combinations(range(m), 8):
            if S[0] != 0: continue
            P = np.array([[math.cos(2*math.pi*k/m), math.sin(2*math.pi*k/m)] for k in S])
            (bc, O, deg), nb, euc = moebius_best(P)
            if bc <= 18: print("regular", m, S, "best", bc, "euclid", euc)
    # 3. Square antiprism / other polyhedra with 8 vertices on a sphere, projected
    def sphere_proj(V, view):
        V = np.array(V, float); V /= np.linalg.norm(V, axis=1)[:, None]
        N = np.array(view, float); N /= np.linalg.norm(N)
        z = N; xx = np.cross(z, [0.3, 0.5, 0.8]); xx /= np.linalg.norm(xx); yy = np.cross(z, xx)
        W = V @ np.array([xx, yy, z]).T
        return np.array([[w[0]/(1-w[2]), w[1]/(1-w[2])] for w in W])
    polys = {}
    for h in [0.3, 0.5, 0.7071, 1.0, 1.5]:
        polys[f"square antiprism h={h}"] = [(math.cos(k*math.pi/2), math.sin(k*math.pi/2), h) for k in range(4)] + [(math.cos(k*math.pi/2+math.pi/4), math.sin(k*math.pi/2+math.pi/4), -h) for k in range(4)]
        polys[f"square prism h={h}"] = [(math.cos(k*math.pi/2), math.sin(k*math.pi/2), h) for k in range(4)] + [(math.cos(k*math.pi/2), math.sin(k*math.pi/2), -h) for k in range(4)]
        polys[f"triangular prism+2 apex h={h}"] = [(math.cos(2*k*math.pi/3), math.sin(2*k*math.pi/3), h) for k in range(3)] + [(math.cos(2*k*math.pi/3), math.sin(2*k*math.pi/3), -h) for k in range(3)] + [(0,0,1),(0,0,-1)]
        polys[f"hexagon+2 poles h={h}"] = [(math.cos(k*math.pi/3), math.sin(k*math.pi/3), h) for k in range(6)] + [(0,0,1),(0,0,-1)]
        polys[f"twisted trig prism+2 h={h}"] = [(math.cos(2*k*math.pi/3), math.sin(2*k*math.pi/3), h) for k in range(3)] + [(math.cos(2*k*math.pi/3+math.pi/3), math.sin(2*k*math.pi/3+math.pi/3), -h) for k in range(3)] + [(0,0,1),(0,0,-1)]
    for name, V in polys.items():
        b = 99
        for trial in range(40):
            P = sphere_proj(V, rng.normal(size=3))
            (bc, O, deg), nb, euc = moebius_best(P)
            b = min(b, bc)
        print(f"{name}: |B|={nb} best over views+inversions = {b}")
