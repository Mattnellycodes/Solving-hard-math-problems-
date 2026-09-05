"""Independent check of the (8_3) ingredient.
(1) every family of 8 triples on 8 points pairwise sharing <= 1 point is one isomorphism class (Moebius-Kantor).
(2) no 8 distinct points of the real projective plane have exactly these 8 collinear triples among their
    >=3-point lines: fix a projective frame on 4 points containing no MK triple, determine the other points
    successively as intersections of two MK lines through known points (homogeneous coordinates: no chart issue),
    with one free parameter t where needed; the remaining incidences are polynomial equations in the free
    parameters; check that they have no real solutions with distinct points."""
import itertools, sympy as sp
from canon import canon
pts = range(8)
triples = [frozenset(T) for T in itertools.combinations(pts, 3)]
fams = []
def rec(i, fam):
    if len(fam) == 8: fams.append(list(fam)); return
    if len(triples) - i < 8 - len(fam): return
    for j in range(i, len(triples)):
        T = triples[j]
        if all(len(T & U) <= 1 for U in fam):
            fam.append(T); rec(j + 1, fam); fam.pop()
rec(0, [])
codes = {canon(8, f)[0] for f in fams}
print(f'(1) {len(fams)} labelled families of 8 pairwise <=1-sharing triples on 8 points; isomorphism classes: {len(codes)}')
MK = fams[0]
print('   representative:', sorted(sorted(T) for T in MK))
# (2) realisability
frame = next(S for S in itertools.combinations(pts, 4) if not any(T <= set(S) for T in MK))
print('   frame (no MK triple inside):', frame)
coords = {frame[0]: sp.Matrix([1, 0, 0]), frame[1]: sp.Matrix([0, 1, 0]), frame[2]: sp.Matrix([0, 0, 1]), frame[3]: sp.Matrix([1, 1, 1])}
params = []
used_lines = set()
def line_through(a, b): return coords[a].cross(coords[b])
progress = True
while len(coords) < 8 and progress:
    progress = False
    for p in pts:
        if p in coords: continue
        Ls = [T for T in MK if p in T and all(q in coords for q in T if q != p)]
        if len(Ls) >= 2:
            (a1, b1), (a2, b2) = [tuple(q for q in T if q != p) for T in Ls[:2]]
            P = line_through(a1, b1).cross(line_through(a2, b2))
            coords[p] = P; used_lines.update(Ls[:2]); progress = True
            print(f'   point {p} = intersection of lines {sorted(Ls[0])} and {sorted(Ls[1])}')
            break
        if len(Ls) == 1:
            (a, b) = tuple(q for q in Ls[0] if q != p)
            t = sp.Symbol(f't{p}', real=True); params.append(t)
            coords[p] = coords[a] + t * coords[b]; used_lines.add(Ls[0]); progress = True
            print(f'   point {p} = {a} + t{p} * {b} on line {sorted(Ls[0])} (free parameter)')
            break
assert len(coords) == 8
eqs = []
for T in MK:
    if T in used_lines: continue
    a, b, c = sorted(T)
    eqs.append(sp.expand(sp.Matrix.hstack(coords[a], coords[b], coords[c]).det()))
print('   parameters:', params, '; closing equations:', [sp.factor(e) for e in eqs])
sols = sp.solve(eqs, params, dict=True) if params else []
print('   solutions over C:', sols)
# distinctness / degeneracy for each solution: points must be distinct; every solution must be non-real
real_ok = []
for s in sols:
    vals = {p: coords[p].subs(s) for p in pts}
    distinct = all(sp.Matrix.hstack(vals[a], vals[b]).rank() == 2 for a, b in itertools.combinations(pts, 2))
    isreal = all(sp.simplify(sp.im(v)) == 0 for v in s.values())
    print(f'   solution {s}: real={isreal}, points distinct={distinct}')
    if isreal and distinct: real_ok.append(s)
if not params:
    print('   no free parameter: consistency of closing equations:', [sp.simplify(e) == 0 for e in eqs])
print('(2) real realisations with distinct points:', len(real_ok), '-> lemma', 'VERIFIED' if not real_ok else 'FAILS')
