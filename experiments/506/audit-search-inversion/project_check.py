"""Cross-check of sphere_search hits via an independent path: numerical centre search over pairwise
circle intersections (mpmath, 40 digits), stereographic projection to the plane, planar block count."""
import mpmath as mp
from itertools import combinations
from pent_check import blocks
from sphere_audit import rich_planes
mp.mp.dps = 40
def check(S, N, expect):
    S = [tuple(map(int, p)) for p in S]
    planes = list(rich_planes(S).keys())
    R = len(planes)
    best = (0, None)
    for k1, k2 in combinations(planes, 2):
        n1 = mp.matrix(k1[:3]); n2 = mp.matrix(k2[:3])
        u = mp.matrix([n1[1]*n2[2]-n1[2]*n2[1], n1[2]*n2[0]-n1[0]*n2[2], n1[0]*n2[1]-n1[1]*n2[0]])
        if mp.norm(u) == 0: continue
        # point on the line: least squares
        M = mp.matrix([[k1[0],k1[1],k1[2]],[k2[0],k2[1],k2[2]],[u[0],u[1],u[2]]])
        p0 = mp.lu_solve(M, mp.matrix([k1[3], k2[3], 0]))
        a = (u.T*u)[0]; b = 2*(p0.T*u)[0]; c = (p0.T*p0)[0] - N
        disc = b*b - 4*a*c
        if disc < 0: continue
        for sg in (1, -1):
            t = (-b + sg*mp.sqrt(disc))/(2*a)
            Q = p0 + t*u
            if any(mp.norm(Q - mp.matrix(p)) < mp.mpf('1e-20') for p in S): continue
            L = sum(1 for pl in planes if abs(pl[0]*Q[0]+pl[1]*Q[1]+pl[2]*Q[2]-pl[3]) < mp.mpf('1e-20'))
            if L > best[0]: best = (L, Q)
    L, Q = best
    # stereographic projection from Q to the plane Q.x = 0
    e1 = mp.matrix([Q[1], -Q[0], 0]); 
    if mp.norm(e1) < mp.mpf('1e-10'): e1 = mp.matrix([0, Q[2], -Q[1]])
    e1 = e1/mp.norm(e1)
    e2 = mp.matrix([Q[1]*e1[2]-Q[2]*e1[1], Q[2]*e1[0]-Q[0]*e1[2], Q[0]*e1[1]-Q[1]*e1[0]]); e2 = e2/mp.norm(e2)
    pts = []
    for p in S:
        p = mp.matrix(p); d = (Q.T*(p - Q))[0]
        t = -N/d
        pp = Q + t*(p - Q)
        pts.append(((e1.T*pp)[0], (e2.T*pp)[0]))
    B = blocks(pts)
    circ = sum(1 for b in B if b[0] == 'C'); lin = sum(1 for b in B if b[0] == 'L')
    print(f"N={N} n={len(S)}: 3D: R={R}, best L={L} -> {R-L};  planar recount: circles={circ}, lines={lin}; expected {expect}; OK={circ==expect and R-L==expect}")
check([(-2, -1, 0), (-2, 0, -1), (-2, 0, 1), (-1, 0, -2), (-1, 0, 2), (1, 0, -2), (1, 0, 2), (2, 0, -1), (2, 0, 1)], 5, 25)
check([(-2, -1, 0), (-2, 0, -1), (-1, -2, 0), (-1, 0, -2), (0, -2, -1), (0, -1, -2), (0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0)], 5, 39)
check([(-2, -1, -1), (-2, -1, 1), (-1, -2, -1), (-1, -2, 1), (-1, -1, -2), (-1, -1, 2), (-1, 1, -2), (-1, 1, 2), (-1, 2, -1), (-1, 2, 1)], 6, 45)
