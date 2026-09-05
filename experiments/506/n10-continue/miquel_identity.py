"""Own verification of the Miquel cross-ratio identity used by the Miquel filter.
Cube vertices (x,y,z) ∈ {0,1}^3 carry complex symbols.  For a face with vertices in cyclic order (a,b,c,d),
E = (a-b)(c-d)/((b-c)(d-a)) = cross-ratio (a,c;b,d): real iff a,b,c,d concyclic/collinear (distinct points).
Orientation: faces x=const use (y-edges)/(z-edges), y=const use (z-edges)/(x-edges), z=const use (x-edges)/(y-edges),
so every edge of the cube occurs once in a numerator and once in a denominator: the product of the six E is
identically ±1.  Hence if five faces are concyclic (five E real and nonzero) the sixth E is real: Miquel's theorem
for 8 distinct points of the Möbius plane (cross-ratios are Möbius invariant, so ∞ may be one of the points)."""
import itertools, sympy as sp
verts = list(itertools.product((0, 1), repeat=3))
z = {v: sp.Symbol('z%d%d%d' % v) for v in verts}
def face_E(axis, val):
    f = [v for v in verts if v[axis] == val]
    j, k = [c for c in range(3) if c != axis]          # the face's two edge directions, j < k
    # cyclic order: (j,k) = (0,0) -> (1,0) -> (1,1) -> (0,1)
    order = [(0, 0), (1, 0), (1, 1), (0, 1)]
    a, b, c, d = [next(v for v in f if (v[j], v[k]) == o) for o in order]
    E = (z[a] - z[b]) * (z[c] - z[d]) / ((z[b] - z[c]) * (z[d] - z[a]))
    # edges ab, cd are j-edges; bc, da are k-edges -> E = (j-edges)/(k-edges)
    want = {0: (1, 2), 1: (2, 0), 2: (0, 1)}[axis]    # (numerator direction, denominator direction)
    if (j, k) != want:
        E = 1 / E
    return E, (a, b, c, d)
prod = sp.Integer(1)
for axis in range(3):
    for val in (0, 1):
        E, (a, b, c, d) = face_E(axis, val)
        cr = (z[a] - z[b]) * (z[c] - z[d]) / ((z[a] - z[d]) * (z[c] - z[b]))   # (a,c;b,d)
        assert sp.simplify(E - cr) == 0 or sp.simplify(E - 1 / cr) == 0
        prod *= E
print("product of the six face cross-ratios =", sp.simplify(prod))
