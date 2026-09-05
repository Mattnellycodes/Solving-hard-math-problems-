"""Inversion / augmentation attack (floating point, tolerance-based; hits are re-verified exactly when
rational).  For a 9-point set P (as a Mobius set on the sphere) circles(P inverted about O) =
|B(P)| - deg(O) for O not in P, where deg(O) = number of blocks of P (circles or lines with >= 3 points)
through O.  So for every candidate P we compute all blocks, all pairwise intersection points of
blocks, the maximal number of blocks concurrent at such a point, and hence the best achievable count
from P.  Candidates: structured 9-point sets (antipodal, concentric polygons + centre, grids with
few blocks, orthocentric/cube record sets + a 9th point, Pappus, etc.)."""
import numpy as np, itertools, math, sys
from mycount import analyse
TOL = 1e-7

def block_of(p, q, r):
    """circle (cx,cy,R) or line (nx,ny,d) through p,q,r; returns ('C',cx,cy,R) or ('L',nx,ny,d) or None."""
    ax, ay = p; bx, by = q; cx, cy = r
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    if abs(d) < 1e-9 * (1 + abs(ax) + abs(bx) + abs(cx)) ** 2:
        nx, ny = by - ay, ax - bx
        nn = math.hypot(nx, ny)
        if nn < 1e-12: return None
        nx, ny = nx / nn, ny / nn
        dd = nx * ax + ny * ay
        if nx < 0 or (abs(nx) < 1e-12 and ny < 0): nx, ny, dd = -nx, -ny, -dd
        return ('L', nx, ny, dd)
    a2 = ax * ax + ay * ay; b2 = bx * bx + by * by; c2 = cx * cx + cy * cy
    ux = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d
    uy = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d
    return ('C', ux, uy, math.hypot(ax - ux, ay - uy))

def same(b1, b2):
    return b1[0] == b2[0] and all(abs(x - y) < TOL * (1 + abs(x)) for x, y in zip(b1[1:], b2[1:]))

def blocks_of(P):
    bl = []   # list of (block, set of point indices)
    for i, j, k in itertools.combinations(range(len(P)), 3):
        b = block_of(P[i], P[j], P[k])
        if b is None: continue
        for t, (b2, s) in enumerate(bl):
            if same(b, b2): s.update((i, j, k)); break
        else:
            bl.append((b, {i, j, k}))
    return bl

def on_block(pt, b):
    x, y = pt
    if b[0] == 'C':
        return abs(math.hypot(x - b[1], y - b[2]) - b[3]) < 1e-6 * (1 + b[3])
    return abs(b[1] * x + b[2] * y - b[3]) < 1e-6 * (1 + abs(b[3]))

def intersections(b1, b2):
    pts = []
    if b1[0] == 'L' and b2[0] == 'L':
        A = np.array([[b1[1], b1[2]], [b2[1], b2[2]]]); rhs = np.array([b1[3], b2[3]])
        if abs(np.linalg.det(A)) > 1e-9: pts.append(tuple(np.linalg.solve(A, rhs)))
        return pts
    if b1[0] == 'C' and b2[0] == 'L': b1, b2 = b2, b1
    if b1[0] == 'L':   # line-circle
        nx, ny, d = b1[1:]; cx, cy, R = b2[1:]
        # foot of perpendicular from centre
        t = d - (nx * cx + ny * cy); fx, fy = cx + t * nx, cy + t * ny
        h2 = R * R - t * t
        if h2 < -1e-9: return pts
        h = math.sqrt(max(h2, 0.0))
        for s in (+1, -1): pts.append((fx + s * h * (-ny), fy + s * h * nx))
        return pts
    # circle-circle
    x1, y1, r1 = b1[1:]; x2, y2, r2 = b2[1:]
    dd = math.hypot(x2 - x1, y2 - y1)
    if dd < 1e-12 or dd > r1 + r2 + 1e-9 or dd < abs(r1 - r2) - 1e-9: return pts
    a = (r1 * r1 - r2 * r2 + dd * dd) / (2 * dd); h2 = r1 * r1 - a * a
    h = math.sqrt(max(h2, 0.0))
    mx, my = x1 + a * (x2 - x1) / dd, y1 + a * (y2 - y1) / dd
    for s in (+1, -1): pts.append((mx + s * h * (-(y2 - y1) / dd), my + s * h * ((x2 - x1) / dd)))
    return pts

def best_inversion(P, verbose=False):
    P = [(float(x), float(y)) for x, y in P]
    bl = blocks_of(P)
    nb = len(bl)
    if any(len(s) == len(P) for _, s in bl): return None
    best = 0; bestO = None
    cands = []
    for (b1, _), (b2, _) in itertools.combinations(bl, 2):
        cands.extend(intersections(b1, b2))
    for O in cands:
        if any(math.hypot(O[0] - x, O[1] - y) < 1e-6 for x, y in P): continue
        deg = sum(1 for b, _ in bl if on_block(O, b))
        if deg > best: best, bestO = deg, O
    nlines = sum(1 for b, _ in bl if b[0] == 'L')
    return nb, nlines, best, bestO, nb - best

if __name__ == "__main__":
    from fractions import Fraction as Fr
    pools = {}
    def circ(k, r=1.0, phase=0.0, cx=0.0, cy=0.0):
        return [(cx + r * math.cos(2 * math.pi * i / k + phase), cy + r * math.sin(2 * math.pi * i / k + phase)) for i in range(k)]
    pools["antipodal (8 on circle + centre)"] = circ(8) + [(0, 0)]
    pools["regular 9-gon minus? no: 9-gon+..."] = circ(9)[:8] + [(0, 0)]
    pools["two squares + centre"] = circ(4) + circ(4, 2) + [(0, 0)]
    pools["two squares rotated + centre"] = circ(4) + circ(4, 2, math.pi / 4) + [(0, 0)]
    pools["two squares (r=sqrt2 rot) + centre"] = circ(4) + circ(4, math.sqrt(2), math.pi / 4) + [(0, 0)]
    pools["3x3 grid"] = [(x, y) for x in range(3) for y in range(3)]
    pools["hexagon + centre + 2"] = circ(6) + [(0, 0), (2, 0), (-2, 0)]
    pools["two triangles + centre + 2"] = circ(3) + circ(3, 2, math.pi / 3) + [(0, 0), (0.5, 0), (-0.5, 0)]
    pools["cube record + centre-ish"] = [(0,0),(0,5),(0,10),(0,15),(3,6),(3,9),(5,5),(5,10),(0,7.5)]
    pools["orthocentric 7 + 2"] = [(0,0),(20,0),(5,15),(5,0),(2,6),(10,10),(5,5),(10,7.5),(20/3,5)]
    pools["pappus"] = [(0,0),(1,0),(2,0),(0,1),(1,1),(2,1),(0.5,0.5),(1.5,0.5),(1.0,2/3)]
    pools["3 concentric triangles"] = circ(3) + circ(3, 2) + circ(3, 3)
    pools["3 conc. triangles alt"] = circ(3) + circ(3, 2, math.pi/3) + circ(3, 3)
    pools["square+centre+4 midpoints"] = [(1,1),(-1,1),(-1,-1),(1,-1),(0,0),(1,0),(0,1),(-1,0),(0,-1)]
    pools["nine-point-ish: triangle, feet, midpoints"] = None
    # triangle A,B,C, feet D,E,F, midpoints of sides
    A, B, C = (0, 3), (1, 0), (3, 0)
    def foot(Pp, Q, R):
        qx, qy = Q; rx, ry = R; px, py = Pp
        t = ((px - qx) * (rx - qx) + (py - qy) * (ry - qy)) / ((rx - qx) ** 2 + (ry - qy) ** 2)
        return (qx + t * (rx - qx), qy + t * (ry - qy))
    pools["nine-point-ish: triangle, feet, midpoints"] = [A, B, C, foot(A, B, C), foot(B, A, C), foot(C, A, B),
        ((A[0]+B[0])/2,(A[1]+B[1])/2), ((B[0]+C[0])/2,(B[1]+C[1])/2), ((A[0]+C[0])/2,(A[1]+C[1])/2)]
    H = (0, -1)
    pools["orthocentric ABCH + feet + midpoints(2)"] = [A, B, C, H, foot(A, B, C), foot(B, A, C), foot(C, A, B), ((A[0]+H[0])/2,(A[1]+H[1])/2), ((B[0]+C[0])/2,(B[1]+C[1])/2)]
    for name, P in pools.items():
        r = best_inversion(P)
        if r is None: print(f"{name}: degenerate"); continue
        nb, nl, deg, O, cnt = r
        print(f"{name}: |B|={nb} lines={nl} circles(as is)={nb-nl}  max concurrency at a non-point={deg}  -> best count {cnt}")
