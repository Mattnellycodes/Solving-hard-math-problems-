"""Verify the cross-ratio identity behind Miquel's theorem (auditor's own check).

Cube vertices z0..z7: 0:(000) 1:(100) 2:(110) 3:(010) 4:(001) 5:(101) 6:(111) 7:(011).
cr(a,b,c,d) = ((a-b)(c-d)) / ((b-c)(d-a)) is real iff the four distinct points are concyclic or
collinear (it is a cross-ratio of the four points).  Numerator edges: ab, cd; denominator: bc, da.
Choose for each face the traversal so that x-edges are numerators on the xy-faces, y-edges on the
yz-faces, z-edges on the xz-faces; then every edge occurs once in a numerator and once in a
denominator, so the product over the six faces is +-1: if five faces have real cross-ratio, so has
the sixth."""
import sympy as sp
import cmath, random

z = sp.symbols('z0:8')


def cr(a, b, c, d):
    return ((a - b) * (c - d)) / ((b - c) * (d - a))


faces = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 3, 7, 4), (1, 2, 6, 5), (4, 0, 1, 5), (7, 3, 2, 6)]
prod = sp.Integer(1)
for f in faces:
    prod *= cr(*[z[i] for i in f])
print('product simplifies to:', sp.cancel(prod))
pts = [cmath.exp(1j * t) for t in (0.3, 1.1, 2.0, 4.5)]
print('cr of 4 concyclic points:', cr(*pts))
pts = [1.0, 2.0, 3.5, -4.0]
print('cr of 4 collinear points:', cr(*pts))
pts = [complex(random.random(), random.random()) for _ in range(4)]
print('cr of 4 random points:', cr(*pts))
# each face quadruple is a face of common.CUBE_FACES cube
import common
cf = set()
for a, b in common.CUBE_FACES:
    cf.add(frozenset(a)); cf.add(frozenset(b))
assert cf == set(frozenset(f) for f in faces)
print('faces agree with common.CUBE_FACES')
