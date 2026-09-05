"""Necessary realisability conditions for an abstract structure (F, L) on P = {0..9} (own implementation).
F = rich blocks (size >= 4), L = lines (size >= 3).  The 11-point Möbius set is P ∪ {INF}, INF = 10, with blocks
F ∪ {S ∪ {INF} : S ∈ L} (plus the 3-point blocks = uncovered triples, which are implicit).

 (1) hereditary Sylvester–Gallai.  Real linear spaces considered:
       - derived structure at p ∈ P (inversion about p): points P - {p}, lines {B - p : p ∈ B ∈ F};
       - the same for the 11-point set: points (P - {p}) ∪ {INF}, lines as above plus {(S - p) ∪ {INF} : p ∈ S ∈ L};
       - the line structure of P: points P, lines L.
     For every subset S (|S| >= 3) not contained in one line, some pair of S must lie on no line meeting S in
     >= 3 points ('ordinary line of S').  Violating S = 'SG-closed'.
 (2) hereditary orchard: #{lines meeting S in exactly 3 points} <= t3(|S|); t3(7..9) = 6, 7, 10 are proved
     (facts.py: Fano, (8_3), AG(2,3)-minus-line); t3(10) = 12 [BGS 1974, CITED] is applied only when cited=True.
 (3) Miquel closure on the 11-point set: for 8 points labelled by the cube vertices, if 5 faces are 'coblock'
     (contained in a block) the 6th must be (facts.py F5).  A structure with exactly 5 coblock faces of some cube is
     not exactly realisable.
 NOTE: the naive 'bundle theorem' (five of the six quadruples Pi ∪ Pj of four disjoint pairs coblock => sixth)
     is FALSE in degenerate position (three pairs on one circle; the antipodal 10-point set violates it), so no
     bundle filter is used.
"""
import itertools
from math import comb

INF = 10
T3_PROVED = {3: 1, 4: 1, 5: 2, 6: 4, 7: 6, 8: 7, 9: 10}
T3_CITED = {10: 12, 11: 16, 12: 19}


def derived(F, p):
    return [frozenset(B) - {p} for B in F if p in B]


def derived11(F, L, p):
    """>= 3-point lines of the real 10-point set obtained by inverting P ∪ {∞} about p: a rich block B ∋ p gives
    B - p, extended by ∞' (= INF) iff B is also a line; a 3-line S ∋ p gives (S - p) ∪ {∞'}.  (A rich line is ONE
    line of the derived set.)"""
    Fs = {frozenset(B) for B in F}; Ls = {frozenset(S) for S in L}
    out = []
    for B in Fs:
        if p in B:
            out.append((B - {p}) | ({INF} if B in Ls else frozenset()))
    for S in Ls:
        if p in S and S not in Fs:
            out.append((S - {p}) | {INF})
    return out


def sg_closed(points, lines):
    """minimal SG-closed subsets of the linear space (points, lines); lines pairwise share <= 1 point."""
    pts = sorted(points)
    lines = [frozenset(l) for l in lines]
    pair_line = {}
    for l in lines:
        for a, b in itertools.combinations(sorted(l), 2):
            pair_line[(a, b)] = l
    viol = []
    for r in range(3, len(pts) + 1):
        for S in itertools.combinations(pts, r):
            Sset = frozenset(S)
            if any(v <= Sset for v in viol):
                continue
            if any(Sset <= l for l in lines):
                continue
            ordinary = False
            for a, b in itertools.combinations(S, 2):
                l = pair_line.get((a, b))
                if l is None or len(l & Sset) == 2:
                    ordinary = True; break
            if not ordinary:
                viol.append(Sset)
    return viol


def orchard_excess(points, lines, cited=False):
    pts = sorted(points); lines = [frozenset(l) for l in lines]
    t3 = dict(T3_PROVED); 
    if cited: t3.update(T3_CITED)
    out = []
    for r in range(7, len(pts) + 1):
        if r not in t3: continue
        for S in itertools.combinations(pts, r):
            Sset = frozenset(S)
            c = sum(1 for l in lines if len(l & Sset) == 3)
            if c > t3[r]:
                out.append((S, c, t3[r]))
    return out


def all_linear_spaces(F, L):
    """(tag, points, lines) for every real linear space that must exist."""
    P = frozenset(range(10))
    out = []
    for p in range(10):
        out.append((('derived', p), P - {p}, derived(F, p)))
        out.append((('derived11', p), (P - {p}) | {INF}, derived11(F, L, p)))
    out.append((('lines',), P, [frozenset(S) for S in L]))
    return out


def hereditary_sg_report(F, L):
    rep = {}
    for tag, pts, lines in all_linear_spaces(F, L):
        v = sg_closed(pts, lines)
        if v:
            rep[tag] = [sorted(s) for s in v]
    return rep


def orchard_report(F, L, cited=False):
    rep = {}
    for tag, pts, lines in all_linear_spaces(F, L):
        v = orchard_excess(pts, lines, cited)
        if v:
            rep[tag] = v[:3]
    return rep


_CUBES = None
def cube_face_systems():
    """the 840 distinct face systems of a cube on labels 0..7 (each = 6 faces as frozensets of 4 labels)."""
    global _CUBES
    if _CUBES is None:
        verts = list(itertools.product((0, 1), repeat=3))
        faces = [[verts.index(v) for v in verts if v[i] == b] for i in range(3) for b in (0, 1)]
        seen = set()
        for perm in itertools.permutations(range(8)):
            fs = frozenset(frozenset(perm[j] for j in f) for f in faces)
            seen.add(fs)
        _CUBES = [tuple(fs) for fs in seen]
    return _CUBES


def coblock_table(F, L):
    blocks = [frozenset(B) for B in F] + [frozenset(S) | {INF} for S in L]
    pts = list(range(11))
    tab = {}
    for Q in itertools.combinations(pts, 4):
        Q = frozenset(Q)
        tab[Q] = any(Q <= B for B in blocks)
    return tab


def miquel_violations(F, L, tab=None):
    tab = tab or coblock_table(F, L)
    viol = []
    for S in itertools.combinations(range(11), 8):
        for fs in cube_face_systems():
            faces = [frozenset(S[j] for j in f) for f in fs]
            c = sum(tab[f] for f in faces)
            if c == 5:
                viol.append(tuple(sorted(f) for f in faces))
    return viol


def full_report(F, L, cited=False):
    F = [frozenset(B) for B in F]; L = [frozenset(S) for S in L]
    hs = hereditary_sg_report(F, L)
    orc = orchard_report(F, L, cited)
    tab = coblock_table(F, L)
    mq = miquel_violations(F, L, tab)
    kills = []
    if hs: kills.append('hereditarySG')
    if orc: kills.append('orchard' + ('(cited t3(10))' if cited else ''))
    if mq: kills.append('Miquel')
    return {'hereditary_sg': {str(k): v for k, v in hs.items()}, 'orchard': {str(k): v for k, v in orc.items()},
            'n_miquel': len(mq), 'miquel_example': mq[:1], 'kills': kills}


if __name__ == '__main__':
    # self-tests
    FANO = [{0,1,2},{0,3,4},{0,5,6},{1,3,5},{1,4,6},{2,3,6},{2,4,5}]
    assert sg_closed(range(7), FANO) == [frozenset(range(7))]
    MK = [{i % 8, (i+1) % 8, (i+3) % 8} for i in range(8)]
    assert sg_closed(range(8), MK) == [] and orchard_excess(range(8), MK) != []
    assert len(cube_face_systems()) == 840
    # the antipodal 10-point structure: 9-block + 4 lines through point 0 -> must pass everything
    F = [frozenset(range(1, 10))]; L = [frozenset({0, 1, 2}), frozenset({0, 3, 4}), frozenset({0, 5, 6}), frozenset({0, 7, 8})]
    print("antipodal structure:", full_report(F, L, cited=True))
    # 8-point cube structure (SQS(8) minus a parallel class) padded with 2 isolated points: Miquel-closed
    verts = list(itertools.product((0, 1), repeat=3)); idx = {v: i for i, v in enumerate(verts)}
    faces = [frozenset(idx[v] for v in verts if v[i] == b) for i in range(3) for b in (0, 1)]
    diag = [frozenset(idx[v] for v in verts if (v[i] + v[j]) % 2 == s) for i, j in itertools.combinations(range(3), 2) for s in (0, 1)]
    cube = faces + diag
    r = full_report(cube, [], cited=True)
    print("cube structure (Miquel-closed expected):", r['n_miquel'], r['kills'])
    r = full_report(cube[:-1], [], cited=True)
    print("cube minus one block (Miquel violation expected):", r['n_miquel'], r['kills'])
