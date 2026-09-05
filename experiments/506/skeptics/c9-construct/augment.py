"""Augmentation attack: base 8-point sets with rich block structure; 9th point placed at every
intersection point of two blocks of the base (points where many blocks concur), then the best
inversion is applied.  Also feeds the low-count 9-sets found on grids into the inversion optimiser."""
import math, itertools, sys
from invert_attack import blocks_of, intersections, best_inversion, on_block

def record8(b):
    A = (0.0, 3.0); B = (b, 0.0); C = (3.0 / b, 0.0)
    H = (0.0, -1.0)   # orthocentre of A,B,C with B*C x-coords product 3: H=(0, -bc/3)= (0,-1)
    def foot(P, Q, R):
        qx, qy = Q; rx, ry = R; px, py = P
        t = ((px - qx) * (rx - qx) + (py - qy) * (ry - qy)) / ((rx - qx) ** 2 + (ry - qy) ** 2)
        return (qx + t * (rx - qx), qy + t * (ry - qy))
    D, E, F = foot(A, B, C), foot(B, A, C), foot(C, A, B)
    O = (0.0, 1.0)
    pts = [A, B, C, H, D, E, F, O]
    # invert about O with radius 1
    out = []
    for x, y in pts:
        dx, dy = x - O[0], y - O[1]; r2 = dx * dx + dy * dy
        if r2 < 1e-12: continue
        out.append((O[0] + dx / r2, O[1] + dy / r2))
    return pts, out   # Mobius 8-set with O in it (7 finite pts + O), and the 8-point inverted record set

def circ(k, r=1.0, phase=0.0):
    return [(r * math.cos(2 * math.pi * i / k + phase), r * math.sin(2 * math.pi * i / k + phase)) for i in range(k)]

bases = {}
for b in (1.0, 0.5, 1.5, 2/3, 1.2):
    pts, inv = record8(b)
    bases[f"record8 b={b}"] = inv
    bases[f"orthocentric+feet+O b={b}"] = pts
bases["two squares"] = circ(4) + circ(4, 2)
bases["two squares rot"] = circ(4) + circ(4, 2, math.pi / 4)
bases["antipodal 8"] = circ(8)
bases["cube projection (generic)"] = [(x + 0.3 * z, y + 0.55 * z) for x in (0, 1) for y in (0, 1) for z in (0, 1)]
bases["3x3 grid minus corner"] = [(x, y) for x in range(3) for y in range(3)][:8]
bases["hexagon+centre+1"] = circ(6) + [(0, 0), (2, 0)]
bases["square+centre+3 mid"] = [(1,1),(-1,1),(-1,-1),(1,-1),(0,0),(1,0),(0,1),(-1,0)]

overall = []
for name, P in bases.items():
    bl = blocks_of(P)
    cands = []
    for (b1, _), (b2, _) in itertools.combinations(bl, 2):
        cands.extend(intersections(b1, b2))
    # dedupe
    uniq = []
    for c in cands:
        if any(math.hypot(c[0] - x, c[1] - y) < 1e-6 for x, y in P): continue
        if any(math.hypot(c[0] - u[0], c[1] - u[1]) < 1e-6 for u in uniq): continue
        uniq.append(c)
    best = (10 ** 9, None)
    for X in uniq:
        r = best_inversion(P + [X])
        if r is None: continue
        nb, nl, deg, O, cnt = r
        if cnt < best[0]: best = (cnt, (X, nb, nl, deg))
    print(f"{name}: base |B|={len(bl)}, {len(uniq)} candidate 9th points; best 9-point count after optimal inversion = {best[0]}  detail={best[1]}")
    overall.append(best[0])
print("overall minimum:", min(overall))
