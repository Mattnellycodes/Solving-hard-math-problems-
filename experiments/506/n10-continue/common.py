"""Erdős #506, n = 10 — shared utilities for the n10-continue agent (2026-09-05).  Own code; nothing imported
from the other experiment directories.

Möbius setting.  P = {0..9}.  A rich block is a subset of size >= 4 (points of P on one circle-or-line), a line
is a block through infinity (collinear subset of size >= 3).  circles(P) = C(10,3) - D - l where
D = sum_{rich B} (C(|B|,3) - 1) and l = number of lines.  Two blocks share <= 2 points, two lines <= 1 point,
a line of size >= 4 is a rich block, a 3-line is a triple in no rich block.

Derived structures.  Inverting P ∪ {∞} about p ∈ P gives a real 10-point set Q_p = (P - p) ∪ {∞'} whose
>= 3-point lines are: for a rich block B ∋ p the set B - p, extended by ∞' iff B is a line; for a 3-line
S ∋ p the set (S - p) ∪ {∞'}.  Dropping ∞' gives the 9-point derived set of P alone.  The line structure of P
itself (points P, lines L) is a third real linear space.  Every real linear space must satisfy:
 (HSG) hereditary Sylvester–Gallai: every subset S with |S| >= 3 that is NOT contained in one line has an
       ordinary line, i.e. a pair of S on no line meeting S in >= 3 points;
 (ORCH) at most t3(r) lines meeting an r-subset in exactly 3 points, t3(7..9) = 6, 7, 10 proved in this session
       (Fano by SG, (8_3) unrealisable, 11-packings on 9 points contain an (8_3)), t3(10) <= 13 by parity
       (proved below), t3(10) = 12 CITED (Burr–Grünbaum–Sloane 1974);
 (MIQ) Miquel: 8 distinct points of the Möbius plane labelled by a cube never have exactly 5 coblock faces.
"""
import itertools, math
from collections import Counter
import networkx as nx

N = 10
INF = 10
T3_PROVED = {3: 1, 4: 1, 5: 2, 6: 4, 7: 6, 8: 7, 9: 10, 10: 13}   # 10: parity/packing bound (see t3_10_parity)
T3_CITED = {10: 12, 11: 16, 12: 19}                              # Burr–Grünbaum–Sloane 1974


def deficit(k):
    return math.comb(k, 3) - 1


def count(F, L, n=N):
    return math.comb(n, 3) - sum(deficit(len(B)) for B in F) - len(L)


def fs(x):
    return frozenset(x)


# ------------------------------------------------------------------------------------------ isomorphism
def incidence_graph(F, L, n=N):
    G = nx.Graph()
    for p in range(n):
        G.add_node(('p', p), c=0)
    for i, B in enumerate(F):
        G.add_node(('b', i), c=10 + len(B))
        for p in B:
            G.add_edge(('b', i), ('p', p))
    for i, S in enumerate(L):
        G.add_node(('l', i), c=30 + len(S))
        for p in S:
            G.add_edge(('l', i), ('p', p))
    return G


def invariant(F, L, n=N):
    prof = []
    for p in range(n):
        c = Counter()
        for B in F:
            if p in B:
                c[('b', len(B))] += 1
        for S in L:
            if p in S:
                c[('l', len(S))] += 1
        prof.append(tuple(sorted(c.items())))
    pairs = Counter()
    for a, b in itertools.combinations(range(n), 2):
        bs = tuple(sorted(len(B) for B in F if a in B and b in B))
        ls = tuple(sorted(len(S) for S in L if a in S and b in S))
        pairs[(bs, ls)] += 1
    G = incidence_graph(F, L, n)
    return (tuple(sorted(len(B) for B in F)), tuple(sorted(len(S) for S in L)), tuple(sorted(prof)),
            tuple(sorted(pairs.items())), nx.weisfeiler_lehman_graph_hash(G, node_attr='c', iterations=3))


NM = nx.algorithms.isomorphism.categorical_node_match('c', 0)


class ClassStore:
    """isomorphism classes of structures (F, L) on n points (invariant buckets + VF2)."""
    def __init__(self, n=N):
        self.n = n; self.buckets = {}; self.classes = []; self.total = 0

    def add(self, F, L=()):
        F = [fs(B) for B in F]; L = [fs(S) for S in L]
        self.total += 1
        inv = invariant(F, L, self.n)
        G = incidence_graph(F, L, self.n)
        for idx in self.buckets.get(inv, []):
            if nx.is_isomorphic(self.classes[idx]['G'], G, node_match=NM):
                self.classes[idx]['copies'] += 1
                return idx, False
        self.buckets.setdefault(inv, []).append(len(self.classes))
        self.classes.append({'F': F, 'L': L, 'copies': 1, 'G': G})
        return len(self.classes) - 1, True

    def records(self):
        return [{'blocks': sorted(sorted(B) for B in c['F']), 'lines': sorted(sorted(S) for S in c['L']),
                 'copies': c['copies'], 'count': count(c['F'], c['L'], self.n)} for c in self.classes]


def automorphisms(F, L=(), n=N):
    G = incidence_graph([fs(B) for B in F], [fs(S) for S in L], n)
    GM = nx.algorithms.isomorphism.GraphMatcher(G, G, node_match=NM)
    return [tuple(iso[('p', p)][1] for p in range(n)) for iso in GM.isomorphisms_iter()]


def apply_perm(g, fam):
    return frozenset(fs(g[p] for p in B) for B in fam)


# ------------------------------------------------------------------------------------------ linear spaces
def derived9(F, p):
    return [fs(B) - {p} for B in F if p in B]


def derived10(F, L, p):
    Fs = {fs(B) for B in F}; Ls = {fs(S) for S in L}
    out = []
    for B in Fs:
        if p in B:
            out.append((B - {p}) | ({INF} if B in Ls else frozenset()))
    for S in Ls:
        if p in S and S not in Fs:
            out.append((S - {p}) | {INF})
    return out


def linear_spaces(F, L, n=N):
    """all real linear spaces forced by (F, L): (tag, points, lines)."""
    P = frozenset(range(n))
    out = []
    for p in range(n):
        out.append((('derived9', p), P - {p}, derived9(F, p)))
        out.append((('derived10', p), (P - {p}) | {INF}, derived10(F, L, p)))
    out.append((('lines',), P, [fs(S) for S in L]))
    return out


def check_linear_space(points, lines):
    lines = [fs(l) for l in lines]
    for a, b in itertools.combinations(lines, 2):
        assert len(a & b) <= 1, (a, b)
    assert all(len(l) >= 3 and l <= set(points) for l in lines)


def sg_closed_subsets(points, lines, minimal_only=True):
    """subsets S (|S| >= 3) of the linear space that are not contained in one line and have no ordinary line
    (pair of S lying on no line that meets S in >= 3 points).  Any such S contradicts Sylvester–Gallai."""
    pts = sorted(points)
    lines = [fs(l) for l in lines]
    pair_line = {}
    for l in lines:
        for a, b in itertools.combinations(sorted(l), 2):
            pair_line[(a, b)] = l
    viol = []
    for r in range(3, len(pts) + 1):
        for S in itertools.combinations(pts, r):
            Sset = frozenset(S)
            if minimal_only and any(v <= Sset for v in viol):
                continue
            if any(Sset <= l for l in lines):
                continue
            ordinary = False
            for a, b in itertools.combinations(S, 2):
                l = pair_line.get((a, b))
                if l is None or len(l & Sset) == 2:
                    ordinary = True
                    break
            if not ordinary:
                viol.append(Sset)
    return viol


def ordinary_lines(points, lines):
    """number of ordinary lines (pairs on no >= 3-point line) of the whole space."""
    lines = [fs(l) for l in lines]
    cov = sum(math.comb(len(l), 2) for l in lines)
    return math.comb(len(points), 2) - cov


def orchard_excess(points, lines, cited=False):
    pts = sorted(points); lines = [fs(l) for l in lines]
    t3 = dict(T3_PROVED)
    if cited:
        t3.update(T3_CITED)
    out = []
    for r in range(7, len(pts) + 1):
        if r not in t3:
            continue
        for S in itertools.combinations(pts, r):
            Sset = frozenset(S)
            c = sum(1 for l in lines if len(l & Sset) == 3)
            if c > t3[r]:
                out.append((S, c, t3[r]))
    return out


def t3_10_parity():
    """13 triples on 10 points pairwise sharing <= 1 point leave 6 uncovered pairs; every point has odd
    leave-degree (9 - 2k); 14 triples would leave 3 pairs, but 10 odd degrees need >= 10 endpoint incidences > 6.
    Hence t3(10) <= 13 by pure packing (used as the proved value)."""
    return 45 - 3 * 14 < 5


_CUBES = None
def cube_face_systems():
    global _CUBES
    if _CUBES is None:
        verts = list(itertools.product((0, 1), repeat=3))
        faces = [[verts.index(v) for v in verts if v[i] == b] for i in range(3) for b in (0, 1)]
        seen = set()
        for perm in itertools.permutations(range(8)):
            seen.add(frozenset(frozenset(perm[j] for j in f) for f in faces))
        _CUBES = [tuple(x) for x in seen]
        assert len(_CUBES) == 840
    return _CUBES


def coblock_table(F, L, n=N):
    blocks = [fs(B) for B in F] + [fs(S) | {INF} for S in L]
    tab = {}
    for Q in itertools.combinations(range(n + 1), 4):
        Q = frozenset(Q)
        tab[Q] = any(Q <= B for B in blocks)
    return tab


def miquel_violations(F, L, n=N, first_only=False):
    tab = coblock_table(F, L, n)
    viol = []
    cubes = cube_face_systems()
    for S in itertools.combinations(range(n + 1), 8):
        for faces in cubes:
            c = 0
            for f in faces:
                if tab[frozenset(S[j] for j in f)]:
                    c += 1
            if c == 5:
                viol.append(tuple(sorted(S[j] for j in f) for f in faces))
                if first_only:
                    return viol
    return viol


def full_report(F, L, cited=False, miquel=True, n=N):
    F = [fs(B) for B in F]; L = [fs(S) for S in L]
    rep = {'hsg': {}, 'orchard': {}, 'ordinary': {}, 'miquel': None}
    for tag, pts, lines in linear_spaces(F, L, n):
        check_linear_space(pts, lines)
        v = sg_closed_subsets(pts, lines)
        if v:
            rep['hsg'][str(tag)] = [sorted(s) for s in v[:2]]
        o = orchard_excess(pts, lines, cited)
        if o:
            rep['orchard'][str(tag)] = o[:2]
        if cited:
            # cited minimum numbers of ordinary lines o(9) = 6, o(10) = 5 (Crowe–McKee 1968)
            m = len(pts)
            need = {9: 6, 10: 5}[m]
            if ordinary_lines(pts, lines) < need:
                rep['ordinary'][str(tag)] = (ordinary_lines(pts, lines), need)
    if miquel:
        mv = miquel_violations(F, L, n, first_only=True)
        rep['miquel'] = mv[0] if mv else None
    kills = []
    if rep['hsg']: kills.append('HSG')
    if rep['orchard']: kills.append('ORCH' + ('(cited)' if cited else ''))
    if rep['ordinary']: kills.append('ORD(cited)')
    if rep['miquel']: kills.append('MIQ')
    rep['kills'] = kills
    return rep


if __name__ == '__main__':
    # self-tests
    FANO = [{0, 1, 2}, {0, 3, 4}, {0, 5, 6}, {1, 3, 5}, {1, 4, 6}, {2, 3, 6}, {2, 4, 5}]
    assert sg_closed_subsets(range(7), FANO) == [frozenset(range(7))]
    MK = [{i % 8, (i + 1) % 8, (i + 3) % 8} for i in range(8)]
    assert sg_closed_subsets(range(8), MK) == [] and orchard_excess(range(8), MK)
    # a 7-point line with an extra point: S = the line is collinear -> not a violation
    assert sg_closed_subsets(range(8), [set(range(7))]) == []
    assert t3_10_parity()
    # antipodal 10-point structure (9-block + 4 chords through the centre) passes everything
    F = [frozenset(range(1, 10))]; L = [{0, 1, 2}, {0, 3, 4}, {0, 5, 6}, {0, 7, 8}]
    r = full_report(F, L, cited=True)
    print("antipodal:", r['kills'], "count", count(F, L))
    assert r['kills'] == []
    # cube structure on 8 points (+2 isolated) is Miquel-closed; minus a block it is not
    verts = list(itertools.product((0, 1), repeat=3)); idx = {v: i for i, v in enumerate(verts)}
    faces = [frozenset(idx[v] for v in verts if v[i] == b) for i in range(3) for b in (0, 1)]
    diag = [frozenset(idx[v] for v in verts if (v[i] + v[j]) % 2 == s) for i, j in itertools.combinations(range(3), 2) for s in (0, 1)]
    cube = faces + diag
    assert miquel_violations(cube, [], first_only=True) == []
    assert miquel_violations(cube[:-1], [], first_only=True) != []
    print("common.py self-tests OK")
