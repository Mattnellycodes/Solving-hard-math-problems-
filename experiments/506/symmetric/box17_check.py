"""Exact check: rectangular box a x b x c on its circumsphere with c^2 = 3 b^2, stereographic
projection from Q=(a,2b,0) (a sphere point on 3 of the 20 vertex planes) gives 17 circles.
Uses sympy exact arithmetic (sqrt(3)), grouping circles by exact centre/radius^2."""
import sympy as sp, itertools, sys
sys.path.insert(0, '/home/user/Solving-hard-math-problems-/experiments/506')
from collections import defaultdict

def stereo(P, Q, R2):
    """stereographic projection from pole Q (on sphere x^2+y^2+z^2=R2) onto the plane through the
    centre orthogonal to Q: maps sphere -> plane; circles through Q -> lines."""
    Q = sp.Matrix(Q); P = sp.Matrix(P)
    # orthonormal basis of Q^perp
    u = Q.cross(sp.Matrix([0, 0, 1]))
    if u.norm() == 0: u = Q.cross(sp.Matrix([1, 0, 0]))
    u = u / u.norm(); v = Q.cross(u) / Q.norm()
    n = Q / Q.norm()
    t = (P - Q); h = n.dot(t)   # signed height relative to Q (negative)
    # project along line Q->P onto plane n.x = 0 : Q + s (P-Q) with n.(Q + s(P-Q)) = 0
    s = -n.dot(Q) / h
    X = Q + s * (P - Q)
    return (sp.nsimplify(sp.simplify(u.dot(X))), sp.nsimplify(sp.simplify(v.dot(X))))

def circ(a, b, c):
    (ax, ay), (bx, by), (cx, cy) = a, b, c
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    d = sp.simplify(d)
    if d == 0: return None
    a2 = ax**2 + ay**2; b2 = bx**2 + by**2; c2 = cx**2 + cy**2
    ux = sp.simplify((a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d)
    uy = sp.simplify((a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d)
    r2 = sp.simplify((ax - ux)**2 + (ay - uy)**2)
    return (sp.radsimp(ux), sp.radsimp(uy), sp.radsimp(r2))

def analyze_exact(pts):
    circles = defaultdict(set); coll = 0
    for i, j, k in itertools.combinations(range(len(pts)), 3):
        c = circ(pts[i], pts[j], pts[k])
        if c is None: coll += 1
        else: circles[c].update((i, j, k))
    return circles, coll

a, b = sp.Integer(1), sp.Integer(1); c = sp.sqrt(3)
R2 = a**2 + b**2 + c**2
V = [(sa * a, sb * b, sc * c) for sa in (1, -1) for sb in (1, -1) for sc in (1, -1)]
Q = (a, 2 * b, 0)
assert sp.simplify(Q[0]**2 + Q[1]**2 + Q[2]**2 - R2) == 0
pts = [stereo(p, Q, R2) for p in V]
print("projected points:")
for p in pts: print("  ", p, [float(x) for x in p])
circles, coll = analyze_exact(pts)
sizes = sorted(len(s) for s in circles.values())
print("circles:", len(circles), "collinear triples:", coll, "sizes:", sizes)
# lines
lines = defaultdict(set)
for i, j, k in itertools.combinations(range(len(pts)), 3):
    if circ(pts[i], pts[j], pts[k]) is None:
        lines[frozenset([i, j, k])] |= {i, j, k}
print("collinear triples list:", [sorted(s) for s in lines])
# non-degeneracy
print("all concyclic?", any(len(s) == len(pts) for s in circles.values()))
