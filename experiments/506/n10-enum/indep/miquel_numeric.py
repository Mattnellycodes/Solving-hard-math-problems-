"""Numeric sanity check of Miquel's theorem as used in mo.miquel_violations: 8 distinct points, faces
(0123),(4567),(0145),(2367),(0246),(1357); if five faces are concyclic so is the sixth.  Random
instances: points 0,1,2,3 on a circle, 4 random, 5 on circle(0,1,4), 6 on circle(0,2,4),
7 = second intersection of circle(2,3,6) and circle(1,3,5); test face (4,5,6,7) via the cross-ratio."""
import random, cmath
def circle3(a, b, c):
    # returns (centre, radius) of circle through a,b,c (complex)
    w = (c - a) / (b - a)
    if abs(w.imag) < 1e-12: return None
    z = (a - b) * (w - abs(w) ** 2) / (2j * w.imag) - a  # standard formula; centre = -z
    centre = -z
    return centre, abs(a - centre)
def on_circle(cen, r, t):
    return cen + r * cmath.exp(1j * t)
def second_intersection(c1, c2, shared):
    (o1, r1), (o2, r2) = c1, c2
    # radical line reflection of the shared point
    d = o2 - o1
    if abs(d) < 1e-12: return None
    # reflect 'shared' across the line of centres
    u = d / abs(d)
    p = (shared - o1) / u
    q = o1 + u * p.conjugate()
    return q
def concyclic(a, b, c, d):
    cr = ((a - c) * (b - d)) / ((a - d) * (b - c))
    return abs(cr.imag)
random.seed(1)
worst = 0
for trial in range(2000):
    cen = complex(random.uniform(-1, 1), random.uniform(-1, 1)); r = random.uniform(0.5, 2)
    P = [on_circle(cen, r, random.uniform(0, 6.283)) for _ in range(4)]
    P.append(complex(random.uniform(-2, 2), random.uniform(-2, 2)))
    c = circle3(P[0], P[1], P[4]); P.append(on_circle(*c, random.uniform(0, 6.283)))
    c = circle3(P[0], P[2], P[4]); P.append(on_circle(*c, random.uniform(0, 6.283)))
    c1 = circle3(P[2], P[3], P[6]); c2 = circle3(P[1], P[3], P[5])
    if c1 is None or c2 is None: continue
    p7 = second_intersection(c1, c2, P[3])
    if p7 is None: continue
    P.append(p7)
    if min(abs(P[i] - P[j]) for i in range(8) for j in range(i + 1, 8)) < 1e-3: continue
    # verify the five faces numerically, then test the sixth
    faces = [(0,1,2,3),(0,1,4,5),(2,3,6,7),(0,2,4,6),(1,3,5,7)]
    if max(concyclic(*[P[i] for i in f]) for f in faces) > 1e-8: continue
    worst = max(worst, concyclic(P[4], P[5], P[6], P[7]))
print("max |Im cross-ratio| of the sixth face over random Miquel cubes:", worst)
