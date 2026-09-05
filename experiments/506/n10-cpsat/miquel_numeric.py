"""Numerical verification of Miquel's theorem in the cube form used by filters.py, with random data:
A1..A4 concyclic; B1 random; B2 on circle(A1,A2,B1); B3 on circle(A2,A3,B2); B4 = circle(A3,A4,B3) ∩
circle(A4,A1,B1) (second point); then B1..B4 must be concyclic.  Also a version where the 'A-circle' is
a line (A1..A4 collinear, i.e. the face passes through infinity)."""
import numpy as np, math, itertools, sys
sys.path.insert(0, '.')
from realise import structure_of
rng = np.random.default_rng(3)
def circ(p, q, r):
    ax, ay = p; bx, by = q; cx, cy = r
    d = 2*(ax*(by-cy)+bx*(cy-ay)+cx*(ay-by))
    a2=ax*ax+ay*ay; b2=bx*bx+by*by; c2=cx*cx+cy*cy
    ux=(a2*(by-cy)+b2*(cy-ay)+c2*(ay-by))/d; uy=(a2*(cx-bx)+b2*(ax-cx)+c2*(bx-ax))/d
    return np.array([ux,uy]), math.hypot(ax-ux, ay-uy)
def on_circle(c, r, t): return c + r*np.array([math.cos(t), math.sin(t)])
def second_intersection(c1, r1, c2, r2, known):
    d = np.linalg.norm(c2-c1); a = (r1*r1 - r2*r2 + d*d)/(2*d); h = math.sqrt(max(r1*r1-a*a, 0))
    m = c1 + a*(c2-c1)/d; perp = np.array([-(c2-c1)[1], (c2-c1)[0]])/d
    p1, p2 = m + h*perp, m - h*perp
    return p2 if np.linalg.norm(p1-known) < np.linalg.norm(p2-known) else p1
def concyclic(pts):
    M = np.array([[x*x+y*y, x, y, 1.0] for x, y in pts]); return abs(np.linalg.det(M))
for trial in range(5):
    if trial < 3:
        cA, rA = rng.normal(size=2), 1 + rng.random()
        A = [on_circle(cA, rA, t) for t in sorted(rng.uniform(0, 2*math.pi, 4))]
    else:  # A-face collinear (through infinity)
        A = [np.array([t, 0.3*t + 1.0]) for t in sorted(rng.uniform(-2, 2, 4))]
    B1 = rng.normal(size=2)*2
    c12, r12 = circ(A[0], A[1], B1); B2 = on_circle(c12, r12, rng.uniform(0, 2*math.pi))
    c23, r23 = circ(A[1], A[2], B2); B3 = on_circle(c23, r23, rng.uniform(0, 2*math.pi))
    c34, r34 = circ(A[2], A[3], B3); c41, r41 = circ(A[3], A[0], B1)
    B4 = second_intersection(c34, r34, c41, r41, known=np.array([1e9,1e9]))
    # pick the intersection that is not (numerically) A[3] or A[0]... both circles pass through A[3]? no: c34 has A3,A4; c41 has A4,A1: common point A4=A[3]
    B4 = second_intersection(c34, r34, c41, r41, known=A[3])
    print(f"trial {trial}: faces A:{concyclic(A):.1e}  B-face determinant: {concyclic([B1,B2,B3,B4]):.2e}  (scale {np.max(np.abs([B1,B2,B3,B4])):.1f})")
