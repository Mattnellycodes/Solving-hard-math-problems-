#!/usr/bin/env python3
"""miquel_identity.py -- exact proof of Miquel's theorem used in miquel.py, via a cross-ratio identity.
Cube vertices p[abc], a,b,c in {0,1}.  For a face we use the cross-ratio of its four vertices in the cyclic
order around the face.  We search (numerically) for a choice of exponents e_f in {+1,-1} and of the cross-ratio
'type' per face such that the product over the six faces is identically 1, then verify the identity EXACTLY
with sympy as an identity of rational functions in the 8 complex coordinates.  Since a quadruple of distinct
points is concyclic (or collinear) iff its cross-ratio is real, five real factors force the sixth to be real:
Miquel's theorem for eight distinct points of the inversive plane, with no further genericity assumption."""
import itertools, random, cmath, sympy as sp

verts = list(itertools.product((0, 1), repeat=3))
def face_cycle(i, c):
    f = [v for v in verts if v[i] == c]
    # cyclic order around the face: (0,0),(0,1),(1,1),(1,0) in the two free coordinates
    others = [k for k in range(3) if k != i]
    order = [(0, 0), (0, 1), (1, 1), (1, 0)]
    cyc = []
    for o in order:
        v = [None] * 3; v[i] = c; v[others[0]] = o[0]; v[others[1]] = o[1]
        cyc.append(verts.index(tuple(v)))
    return cyc
faces = [face_cycle(i, c) for i in range(3) for c in (0, 1)]
def cr(a, b, c, d): return (a - c) * (b - d) / ((a - d) * (b - c))
types = [lambda a, b, c, d: cr(a, b, c, d), lambda a, b, c, d: cr(a, c, b, d), lambda a, b, c, d: cr(a, b, d, c)]
random.seed(0)
samples = [[complex(random.uniform(-2, 2), random.uniform(-2, 2)) for _ in range(8)] for _ in range(4)]
found = None
for tsel in itertools.product(range(3), repeat=6):
    for esel in itertools.product((1, -1), repeat=6):
        ok = True
        for p in samples:
            prod = 1
            for f, t, e in zip(faces, tsel, esel):
                prod *= types[t](*[p[k] for k in f]) ** e
            if abs(prod - 1) > 1e-9 and abs(prod + 1) > 1e-9:
                ok = False; break
        if ok:
            found = (tsel, esel); break
    if found: break
print("numerically found identity:", found)
tsel, esel = found
z = sp.symbols('z0:8')
prod = sp.Integer(1)
for f, t, e in zip(faces, tsel, esel):
    a, b, c, d = [z[k] for k in f]
    if t == 0: val = (a - c) * (b - d) / ((a - d) * (b - c))
    elif t == 1: val = (a - b) * (c - d) / ((a - d) * (c - b))
    else: val = (a - d) * (b - c) / ((a - c) * (b - d))
    prod *= val ** e
print("exact simplification of the product of the six face cross-ratios:", sp.simplify(prod))
